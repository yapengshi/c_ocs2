/*
 * Implementation of SLQP_BASE.h
 *
 *  Created on: October 7, 2016
 *      Author: mgiftthaler@ethz.ch
 */


namespace ocs2{

template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP_BASE<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::solveSequentialRiccatiEquations(const scalar_t& learningRate)
{

	LinearInterpolation<state_matrix_t, Eigen::aligned_allocator<state_matrix_t> > SmFunc;

	// final value for the last Riccati equations
	typename RiccatiEquations_t::s_vector_t allSsFinal;
	allSsFinal.setZero();
	// final value for the last error equation
	output_vector_t SveFinal = output_vector_t::Zero();

	for (int i=NUM_SUBSYSTEMS-1; i>=0; i--) {

		// final cost of the subsystem is added to the following subsystem solution
		typename RiccatiEquations_t::s_vector_t allCostFinal;
		RiccatiEquations_t::convert2Vector(QmFinalStock_[i], QvFinalStock_[i], qFinalStock_[i], allCostFinal);
		allSsFinal += allCostFinal;

		// set data for Riccati equations
		std::shared_ptr<RiccatiEquations_t> riccatiEquationsPtr( new RiccatiEquations_t() );
		riccatiEquationsPtr->setData(learningRate, i, switchingTimes_[i], switchingTimes_[i+1],
				&nominalTimeTrajectoriesStock_[i],
				&AmConstrainedTrajectoryStock_[i], &BmTrajectoryStock_[i],
				&qTrajectoryStock_[i], &QvConstrainedTrajectoryStock_[i], &QmConstrainedTrajectoryStock_[i],
				&RvTrajectoryStock_[i], &RmInverseTrajectoryStock_[i], &RmConstrainedTrajectoryStock_[i], &PmTrajectoryStock_[i]);

		// max number of steps of integration
		size_t maxNumSteps = options_.maxNumStepsPerSecond_ * std::max( 1.0, switchingTimes_[i+1]-switchingTimes_[i] );

		std::vector<double> normalizedTimeTrajectory;
		std::vector<typename RiccatiEquations_t::s_vector_t, Eigen::aligned_allocator<typename RiccatiEquations_t::s_vector_t> > allSsTrajectory;

		switch(options_.RiccatiIntegratorType_){

		case DIMENSIONS::RICCATI_INTEGRATOR_TYPE::ODE45 : {
			ODE45<RiccatiEquations_t::S_DIM_> riccati_integrator (riccatiEquationsPtr);
			riccati_integrator.integrate(allSsFinal, i, i+1, allSsTrajectory, normalizedTimeTrajectory, 1e-5, options_.AbsTolODE_, options_.RelTolODE_, maxNumSteps);
			break;
		}
		/*note: this case is not yet working. It would most likely work if we had an adaptive time adams-bashforth integrator */
		case DIMENSIONS::RICCATI_INTEGRATOR_TYPE::ADAMS_BASHFORTH : {
			const size_t order = 4;
			IntegratorAdamsBashforth<RiccatiEquations_t::S_DIM_,order> riccati_integrator (riccatiEquationsPtr);
			riccati_integrator.integrate(allSsFinal, i,	i+1, options_.adams_integrator_dt_, allSsTrajectory, normalizedTimeTrajectory); // fixed time step
			break;
		}
		case DIMENSIONS::RICCATI_INTEGRATOR_TYPE::BULIRSCH_STOER : {
			IntegratorBulirschStoer<RiccatiEquations_t::S_DIM_> riccati_integrator (riccatiEquationsPtr);
			riccati_integrator.integrate(allSsFinal, i, i+1, allSsTrajectory, normalizedTimeTrajectory, 1e-5, options_.AbsTolODE_, options_.RelTolODE_, maxNumSteps);
			break;
		}
		default:
			throw (std::runtime_error("Riccati equation integrator type specified wrongly in solveSequentialRiccatiEquations()"));
		}

		// denormalizing time and constructing 'Sm', 'Sv', and 's'
		int N = normalizedTimeTrajectory.size();
		SsTimeTrajectoryStock_[i].resize(N);
		SmTrajectoryStock_[i].resize(N);
		SvTrajectoryStock_[i].resize(N);
		sTrajectoryStock_[i].resize(N);
		for (int k=0; k<N; k++) {
			RiccatiEquations_t::convert2Matrix(allSsTrajectory[N-1-k], SmTrajectoryStock_[i][k], SvTrajectoryStock_[i][k], sTrajectoryStock_[i][k]);
			SsTimeTrajectoryStock_[i][k] = (switchingTimes_[i]-switchingTimes_[i+1])*(normalizedTimeTrajectory[N-1-k]-i) + switchingTimes_[i+1];
		}  // end of k loop


		// testing the numerical stability of the Riccati equations
		for (int k=N-1; k>=0; k--) {
			try {
				if (SmTrajectoryStock_[i][k].hasNaN())  throw std::runtime_error("Sm is unstable.");
				if (SvTrajectoryStock_[i][k].hasNaN())  throw std::runtime_error("Sv is unstable.");
				if (sTrajectoryStock_[i][k].hasNaN())   throw std::runtime_error("s is unstable.");
			}
			catch(const std::exception& error)
			{
				std::cerr << "what(): " << error.what() << " at time " << SsTimeTrajectoryStock_[i][k] << " [sec]." << std::endl;
				for (int kp=k; kp<k+10; kp++)  {
					if (kp >= N) continue;
					std::cerr << "Sm[" << SsTimeTrajectoryStock_[i][kp] << "]:\n"<< SmTrajectoryStock_[i][kp].norm() << std::endl;
					std::cerr << "Sv[" << SsTimeTrajectoryStock_[i][kp] << "]:\t"<< SvTrajectoryStock_[i][kp].transpose().norm() << std::endl;
					std::cerr << "s["  << SsTimeTrajectoryStock_[i][kp] << "]: \t"<< sTrajectoryStock_[i][kp].transpose().norm() << std::endl;
				}
				exit(0);
			}
		}  // end of k loop

		// set the final value for next Riccati equation
		allSsFinal = allSsTrajectory.back();

		/*
		 * Type_1 constraints error correction compensation
		 */

		// Skip calculation of the error correction term Sve if the constrained simulation is used for forward simulation
		if (options_.simulationIsConstrained_) {
			SveTrajectoryStock_[i].resize(N);
			for (int k=0; k<N; k++)
				SveTrajectoryStock_[i][k].setZero();
			continue;
		}

		// Calculating the coefficients of the error equation
		SmFunc.setTimeStamp( &(SsTimeTrajectoryStock_[i]) );
		SmFunc.setData( &(SmTrajectoryStock_[i]) );

		output_vector_array_t GvTrajectory(nominalTimeTrajectoriesStock_[i].size());
		state_matrix_array_t  GmTrajectory(nominalTimeTrajectoriesStock_[i].size());

		for (int k=nominalTimeTrajectoriesStock_[i].size()-1; k>=0; k--) {

			state_matrix_t Sm;
			SmFunc.interpolate(nominalTimeTrajectoriesStock_[i][k], Sm);

			control_feedback_t Lm = RmInverseTrajectoryStock_[i][k]*(PmTrajectoryStock_[i][k]+BmTrajectoryStock_[i][k].transpose()*Sm);

			GmTrajectory[k] = AmConstrainedTrajectoryStock_[i][k] -
					BmTrajectoryStock_[i][k]*RmInverseTrajectoryStock_[i][k]*RmConstrainedTrajectoryStock_[i][k]*Lm;

			GvTrajectory[k] = (CmProjectedTrajectoryStock_[i][k]-Lm).transpose()*
					RmTrajectoryStock_[i][k]*EvProjectedTrajectoryStock_[i][k];

		}  // end of k loop

		// set data for error equations
		std::shared_ptr<ErrorEquation_t> errorEquationPtr( new ErrorEquation_t() );
		errorEquationPtr->setData(i, switchingTimes_[i], switchingTimes_[i+1],&nominalTimeTrajectoriesStock_[i], &GvTrajectory, &GmTrajectory);

		// integrating the Riccati equations
		ODE45<OUTPUT_DIM> errorOde45(errorEquationPtr);
		output_vector_array_t SveTrajectory;
		errorOde45.integrate(SveFinal, normalizedTimeTrajectory, SveTrajectory, 1e-3, options_.AbsTolODE_, options_.RelTolODE_);

		// reset the final value for next Riccati equation
		SveFinal = SveTrajectory.back();

		if(SveTrajectory.size() != N)
			throw std::runtime_error("sve traj size not equal to N");

		SveTrajectoryStock_[i].resize(N);
		for (int k=0; k<N; k++) {
			SveTrajectoryStock_[i][k] = SveTrajectory[N-1-k];

			// testing the numerical stability of the Riccati error equation
			try {
				if (SveTrajectoryStock_[i][k].hasNaN())  throw std::runtime_error("Sve is unstable");
			}
			catch(const std::exception& error) 	{
				std::cerr << "what(): " << error.what() << " at time " << SsTimeTrajectoryStock_[i][k] << " [sec]." << std::endl;
				for (int kp=k; kp<N; kp++){
					std::cerr << "Sve[" << SsTimeTrajectoryStock_[i][kp] << "]:\t"<< SveTrajectoryStock_[i][kp].transpose().norm() << std::endl;
				}
				for(size_t kp = 0; kp<nominalTimeTrajectoriesStock_[i].size()-1; kp++){
					std::cerr << "Gm[" << SsTimeTrajectoryStock_[i][kp] << "]:\t"<< GmTrajectory[kp].transpose().norm() << std::endl;
					std::cerr << "Gv[" << SsTimeTrajectoryStock_[i][kp] << "]:\t"<< GvTrajectory[kp].transpose().norm() << std::endl;
				}

				exit(0);
			}
		}

	}  // end of i loop
}


}

