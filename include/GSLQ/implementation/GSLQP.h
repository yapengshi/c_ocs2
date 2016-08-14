/*
 * Implementation of GSLQP.h
 *
 *  Created on: Jan 5, 2016
 *      Author: farbod
 */


namespace ocs2{

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * Forward integrate the system dynamics with given controller:
 * 		inputs:
 * 			+ initState: initial state at time switchingTimes_[0]
 * 			+ controllersStock: controller for each subsystem
 * 		outputs:
 * 			+ timeTrajectoriesStock:  rollout simulated time steps
 * 			+ stateTrajectoriesStock: rollout states
 * 			+ inputTrajectoriesStock: rollout control inputs
 * 			+ (optional) outputTrajectoriesStock: rollout outputs
 * 			+ (optional) nc1TrajectoriesStock: number of active constraints at each time step
 * 			+ (optional) EvTrajectoryStock: value of the constraint (if the rollout is constrained the value is
 * 											always zero otherwise it is nonzero)
 */

/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::rollout(const state_vector_t& initState,
		const std::vector<controller_t>& controllersStock,
		std::vector<scalar_array_t>& timeTrajectoriesStock,
		std::vector<state_vector_array_t>& stateTrajectoriesStock,
		std::vector<control_vector_array_t>& inputTrajectoriesStock,
		std::vector<output_vector_array_t>& outputTrajectoriesStock)  {

	slqp_->rollout(initState,
			controllersStock, timeTrajectoriesStock, stateTrajectoriesStock, inputTrajectoriesStock, outputTrajectoriesStock);

}


/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::rollout(const state_vector_t& initState,
		const std::vector<controller_t>& controllersStock,
		std::vector<scalar_array_t>& timeTrajectoriesStock,
		std::vector<state_vector_array_t>& stateTrajectoriesStock,
		std::vector<control_vector_array_t>& inputTrajectoriesStock)  {

	slqp_->rollout(initState,
			controllersStock, timeTrajectoriesStock, stateTrajectoriesStock, inputTrajectoriesStock);
}

/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::rollout(const state_vector_t& initState,
		const std::vector<controller_t>& controllersStock,
		std::vector<scalar_array_t>& timeTrajectoriesStock,
		std::vector<state_vector_array_t>& stateTrajectoriesStock,
		std::vector<control_vector_array_t>& inputTrajectoriesStock,
		std::vector<output_vector_array_t>& outputTrajectoriesStock,
		std::vector<std::vector<size_t> >& nc1TrajectoriesStock,
		std::vector<constraint1_vector_array_t>& EvTrajectoryStock)  {

	slqp_->rollout(initState,
			controllersStock, timeTrajectoriesStock, stateTrajectoriesStock, inputTrajectoriesStock,
			outputTrajectoriesStock, nc1TrajectoriesStock, EvTrajectoryStock);
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * compute the cost for a given rollout
 * 		inputs:
 * 			+ timeTrajectoriesStock:  rollout simulated time steps
 * 			+ outputTrajectoriesStock: rollout outputs
 * 			+ inputTrajectoriesStock: rollout control inputs
 *
 * 		outputs:
 * 			+ totalCost: the total cost of the trajectory
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateCostFunction(
		const std::vector<scalar_array_t>& timeTrajectoriesStock,
		const std::vector<output_vector_array_t>& outputTrajectoriesStock,
		const std::vector<control_vector_array_t>& inputTrajectoriesStock,
		scalar_t& totalCost)  {

	slqp_->calculateCostFunction(timeTrajectoriesStock, outputTrajectoriesStock, inputTrajectoriesStock,
			totalCost);
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * compute the merit function for given rollout
 * 		inputs:
 * 			+ timeTrajectoriesStock: simulation time trajectory
 * 			+ nc1TrajectoriesStock: rollout's number of active constraints in each time step
 * 			+ EvTrajectoryStock: rollout's constraints value
 * 			+ lagrangeTrajectoriesStock: constraint Lagrange multiplier for the given rollout
 * 			+ totalCost: the total cost of the trajectory
 *
 * 		outputs:
 * 			+ meritFuntionValue: the merit function value
 * 			+ constraintISE: integral of Square Error (ISE)
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateMeritFunction(
		const std::vector<scalar_array_t>& timeTrajectoriesStock,
		const std::vector<std::vector<size_t> >& nc1TrajectoriesStock,
		const std::vector<constraint1_vector_array_t>& EvTrajectoryStock,
		const std::vector<std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd> > >&  lagrangeTrajectoriesStock,
		const scalar_t& totalCost,
		scalar_t& meritFuntionValue,
		scalar_t& constraintISE)  {

	slqp_->calculateMeritFunction(timeTrajectoriesStock, nc1TrajectoriesStock, EvTrajectoryStock, lagrangeTrajectoriesStock, totalCost,
			meritFuntionValue, constraintISE);
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * Constraint's Integral of Squared Error (ISE)
 * 		Inputs:
 * 			+ timeTrajectoriesStock: simulation time trajectory
 * 			+ nc1TrajectoriesStock: rollout's number of active constraints in each time step
 * 			+ EvTrajectoriesStock: rollout's constraints value
 * 		Output:
 * 			+ constraintISE: integral of Square Error (ISE)
 * 		Return:
 * 			+ maximum constraint norm
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
double GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateConstraintISE(
		const std::vector<scalar_array_t>& timeTrajectoriesStock,
		const std::vector<std::vector<size_t>>& nc1TrajectoriesStock,
		const std::vector<constraint1_vector_array_t>& EvTrajectoriesStock,
		scalar_t& constraintISE)  {

	return slqp_->calculateConstraintISE(timeTrajectoriesStock, nc1TrajectoriesStock, EvTrajectoriesStock,
			constraintISE);
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * transform the local value function derivatives to the global one.
 * 		it manipulates the following member variables:
 * 			+ nablasTrajectoryStock_
 * 			+ nablaSvTrajectoryStock_
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::transformLocalValueFuntionDerivative2Global()  {

	LinearInterpolation<output_vector_t,Eigen::aligned_allocator<output_vector_t> > nominalOutputFunc;

	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		nominalOutputFunc.setTimeStamp( &(slqp_->nominalTimeTrajectoriesStock_[i]) );
		nominalOutputFunc.setData( &(slqp_->nominalOutputTrajectoriesStock_[i]) );

		for (int k=0; k<slqp_->SsTimeTrajectoryStock_[i].size(); k++) {

			output_vector_t nominalOutput;
			nominalOutputFunc.interpolate(slqp_->SsTimeTrajectoryStock_[i][k], nominalOutput);

			for (int j=0; j<NUM_SUBSYSTEMS-1; j++)  {

				nablasTrajectoryStock_[i][k][j] += - nominalOutput.transpose()*nablaSvTrajectoryStock_[i][k][j] +
						0.5*nominalOutput.transpose()*nablaSmTrajectoryStock_[i][k][j]*nominalOutput;
				nablaSvTrajectoryStock_[i][k][j]+= - nablaSmTrajectoryStock_[i][k][j]*nominalOutput;
			}  // end of j loop
		}  // end of k loop
	}  // end of i loop
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * get the calculated rollout's sensitivity to switchingTimes
 * 		outputs:
 * 			+ sensitivityTimeTrajectoriesStock: time stamps of the sensitivity values
 * 			+ sensitivityOutputTrajectoriesStock: output trajectory sensitivity to the switching times
 * 			+ sensitivityInputTrajectoriesStock: control input trajectory sensitivity to the switching times
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getRolloutSensitivity2SwitchingTime(
		std::vector<scalar_array_t>& sensitivityTimeTrajectoriesStock,
		std::vector<nabla_output_matrix_array_t>& sensitivityOutputTrajectoriesStock,
		std::vector<nabla_input_matrix_array_t>& sensitivityInputTrajectoriesStock)  {

	sensitivityTimeTrajectoriesStock   = sensitivityTimeTrajectoryStock_;
	sensitivityOutputTrajectoriesStock = nablaOutputTrajectoryStock_;
	sensitivityInputTrajectoriesStock  = nablaInputTrajectoryStock_;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * get the calculated optimal controller structure
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getController(std::vector<controller_t>& controllersStock) {

	slqp_->getController(controllersStock);
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * calculate the value function for the given time and output vector
 * 		inputs
 * 			+ time: inquiry time
 * 			+ output: inquiry output
 *
 * 		output:
 * 			+ valueFuntion: value function at the inquiry time and output
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getValueFuntion(const scalar_t& time, const output_vector_t& output, scalar_t& valueFuntion)  {

	slqp_->getValueFuntion(time, output, valueFuntion);
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * calculate the cost function at the initial time
 * 		inputs
 * 			+ initOutput: initial output
 *
 * 		output:
 * 			+ cost function value
 * 			+ cost function value plus the constraint ISE multiplied by pho
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getCostFuntion(const output_vector_t& initOutput, scalar_t& costFunction, scalar_t& constriantCostFunction)  {

	slqp_->getCostFuntion(initOutput, costFunction, constriantCostFunction);
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * calculate the value function's derivatives w.r.t. switchingTimes at the initial time
 * 		inputs
 * 			+ initOutput: initial output
 *
 * 		output:
 * 			+ valueFuntionDerivative: cost function' derivatives w.r.t. switchingTimes for given initial output vector
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getValueFuntionDerivative(const output_vector_t& initOutput,
		Eigen::Matrix<double,NUM_SUBSYSTEMS-1,1>& valueFuntionDerivative)  {

	for (int j=0; j<NUM_SUBSYSTEMS-1; j++)  {

		state_matrix_t dSm  = nablaSmTrajectoryStock_[0][0][j];
		output_vector_t dSv = nablaSvTrajectoryStock_[0][0][j];
		eigen_scalar_t ds   = nablasTrajectoryStock_[0][0][j];

		valueFuntionDerivative(j) = (ds + initOutput.transpose()*dSv + 0.5*initOutput.transpose()*dSm*initOutput).eval()(0);
	}
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * calculate the cost function's derivatives w.r.t. switchingTimes
 *
 * 		output:
 * 			+ costFunctionDerivative: cost function' derivatives w.r.t. switchingTimes for given initial output vector
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getCostFuntionDerivative(Eigen::Matrix<double,NUM_SUBSYSTEMS-1,1>& costFunctionDerivative)  {

	costFunctionDerivative = nominalCostFuntionDerivative_;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * get the optimal state, output, and input trajectories
 * 		output
 * 			+ nominalTimeTrajectoriesStock_: optimal time trajectory
 * 			+ nominalStateTrajectoriesStock_: optimal state trajectory
 * 			+ nominalInputTrajectoriesStock_: optimal control input trajectory
 * 			+ nominalOutputTrajectoriesStock_: optimal output trajectory
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getNominalTrajectories(std::vector<scalar_array_t>& nominalTimeTrajectoriesStock,
		std::vector<state_vector_array_t>& nominalStateTrajectoriesStock,
		std::vector<control_vector_array_t>& nominalInputTrajectoriesStock,
		std::vector<output_vector_array_t>& nominalOutputTrajectoriesStock)   {

	slqp_->getNominalTrajectories(nominalTimeTrajectoriesStock, nominalStateTrajectoriesStock,
			nominalInputTrajectoriesStock, nominalOutputTrajectoriesStock);
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * solve the SLQ Riccati differential equations plus its derivatives differential equations :
 * 		input:
 * 			+ learningRate: the feeadforward learningRate
 *
 * 		uses:
 * 			+ linearized dynamics
 * 			+ quadratized cost
 * 			+ nominal system sensitivity analysis
 * 			+ SsTimeTrajectoryStock_: time stamp
 * 			V(t,y) = y^T*Sm*y + y^T*(Sv) + s
 * 			+ SmTrajectoryStock_: Sm matrix
 * 			+ SvTrajectoryStock_: Sv vector
 *
 * 		modifies:

 * 			dV(t,y) = y^T*dSm*y + y^T*(dSv) + ds
 * 			+ nablaSmTrajectoryStock_: dSm
 * 			+ nablaSvTrajectoryStock_: dSv
 * 			+ nablasTrajectoryStock_: ds
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::solveSensitivityRiccatiEquations(const scalar_t& learningRate)  {

	// final value for the last Riccati equations
	typename SensitivityRiccatiEquations_t::all_s_vector_t allSsFinal;
	nabla_Sm_t nablaQmFinal;  nablaQmFinal.fill(state_matrix_t::Zero());
	nabla_Sv_t nablaQvFinal;  nablaQvFinal.fill(output_vector_t::Zero());
	nabla_s_t  nablaqFinal;   nablaqFinal.fill(eigen_scalar_t::Zero());
	for (int i=0; i<NUM_SUBSYSTEMS-1; i++) {
		nablaQvFinal[i] = nablaQvFinal_.template block<OUTPUT_DIM,1>(0,i);
		nablaqFinal[i](0) = nablaqFinal_(i);
	}
	SensitivityRiccatiEquations_t::convert2Vector(nablaQmFinal, nablaQvFinal, nablaqFinal, allSsFinal);

	for (int i=NUM_SUBSYSTEMS-1; i>=0; i--) {

		// set data for Riccati equations
		std::shared_ptr<SensitivityRiccatiEquations_t> riccatiEquationsPtr( new SensitivityRiccatiEquations_t());
		riccatiEquationsPtr->setData(learningRate,
				i, switchingTimes_[i], switchingTimes_[i+1],
				&slqp_->SsTimeTrajectoryStock_[i], &slqp_->SmTrajectoryStock_[i], &slqp_->SvTrajectoryStock_[i],
				&slqp_->nominalTimeTrajectoriesStock_[i],
				&slqp_->AmConstrainedTrajectoryStock_[i], &slqp_->BmTrajectoryStock_[i],
				&slqp_->qTrajectoryStock_[i], &slqp_->QvConstrainedTrajectoryStock_[i], &slqp_->QmConstrainedTrajectoryStock_[i],
				&slqp_->RvTrajectoryStock_[i], &slqp_->RmInverseTrajectoryStock_[i], &slqp_->RmConstrainedTrajectoryStock_[i], &slqp_->PmTrajectoryStock_[i],
				&sensitivityTimeTrajectoryStock_[i], &nablaqTrajectoryStock_[i],
				&nablaQvTrajectoryStock_[i], &nablaRvTrajectoryStock_[i]);

		// normalized integration time based on SsTimeTrajectoryStock_
		int N = slqp_->SsTimeTrajectoryStock_[i].size();
		scalar_array_t normalizedTimeTrajectory(N);
		for (int k=0; k<N; k++)
			normalizedTimeTrajectory[N-1-k] = (slqp_->SsTimeTrajectoryStock_[i][k]-switchingTimes_[i+1])/(switchingTimes_[i]-switchingTimes_[i+1]);

		// integrating the Riccati equations
		ODE45<SensitivityRiccatiEquations_t::S_DIM_*(NUM_SUBSYSTEMS-1)> ode45(riccatiEquationsPtr);

		std::vector<typename SensitivityRiccatiEquations_t::all_s_vector_t, Eigen::aligned_allocator<typename SensitivityRiccatiEquations_t::all_s_vector_t> > allSsTrajectory;
		ode45.integrate(allSsFinal, normalizedTimeTrajectory, allSsTrajectory,
				1e-3, options_.AbsTolODE_, options_.RelTolODE_);

		// construct 'Sm', 'Sv', and 's'
		nablaSmTrajectoryStock_[i].resize(N);
		nablaSvTrajectoryStock_[i].resize(N);
		nablasTrajectoryStock_[i].resize(N);
		for (int k=0; k<N; k++)
			SensitivityRiccatiEquations_t::convert2Matrix(allSsTrajectory[N-1-k], nablaSmTrajectoryStock_[i][k], nablaSvTrajectoryStock_[i][k], nablasTrajectoryStock_[i][k]);

		// reset the final value for next Riccati equation
		allSsFinal = allSsTrajectory.back();

	}  // end of i loop

}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * calculates sensitivity controller feedback part (constrained feedback):
 * 		This method uses the following variables:
 * 			+ constrained, linearized model
 * 			+ constrained, quadratized cost
 *
 * 		output:
 * 			+ sensitivityControllersStock: the sensitivity controller
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateSensitivityControllerFeedback(
		std::vector<sensitivity_controller_t>& sensitivityControllersStock) {


	LinearInterpolation<control_gain_matrix_t,Eigen::aligned_allocator<control_gain_matrix_t> > BmFunc;
	LinearInterpolation<control_feedback_t,Eigen::aligned_allocator<control_feedback_t> > PmFunc;
	LinearInterpolation<control_matrix_t,Eigen::aligned_allocator<control_matrix_t> >     RmInverseFunc;
	LinearInterpolation<control_feedback_t,Eigen::aligned_allocator<control_feedback_t> > CmProjectedFunc;
	LinearInterpolation<control_matrix_t,Eigen::aligned_allocator<control_matrix_t> >     DmProjectedFunc;

	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		BmFunc.setTimeStamp( &(slqp_->nominalTimeTrajectoriesStock_[i]) );
		BmFunc.setData( &(slqp_->BmTrajectoryStock_[i]) );

		PmFunc.setTimeStamp( &(slqp_->nominalTimeTrajectoriesStock_[i]) );
		PmFunc.setData( &(slqp_->PmTrajectoryStock_[i]) );

		RmInverseFunc.setTimeStamp( &(slqp_->nominalTimeTrajectoriesStock_[i]) );
		RmInverseFunc.setData( &(slqp_->RmInverseTrajectoryStock_[i]) );

		CmProjectedFunc.setTimeStamp( &(slqp_->nominalTimeTrajectoriesStock_[i]) );
		CmProjectedFunc.setData( &(slqp_->CmProjectedTrajectoryStock_[i]) );

		DmProjectedFunc.setTimeStamp( &(slqp_->nominalTimeTrajectoriesStock_[i]) );
		DmProjectedFunc.setData( &(slqp_->DmProjectedTrajectoryStock_[i]) );

		sensitivityControllersStock[i].time_ = slqp_->SsTimeTrajectoryStock_[i];

		size_t N = slqp_->SsTimeTrajectoryStock_[i].size();
		sensitivityControllersStock[i].k_.resize(N);

		for (int k=0; k<N; k++) {

			const double& time = slqp_->SsTimeTrajectoryStock_[i][k];

			control_gain_matrix_t Bm;
			BmFunc.interpolate(time, Bm);
			size_t greatestLessTimeStampIndex = BmFunc.getGreatestLessTimeStampIndex();
			control_feedback_t Pm;
			PmFunc.interpolate(time, Pm, greatestLessTimeStampIndex);
			control_matrix_t RmInverse;
			RmInverseFunc.interpolate(time, RmInverse, greatestLessTimeStampIndex);
			control_feedback_t CmProjected;
			CmProjectedFunc.interpolate(time, CmProjected, greatestLessTimeStampIndex);
			control_matrix_t DmProjected;
			DmProjectedFunc.interpolate(time, DmProjected, greatestLessTimeStampIndex);

			control_feedback_t Lm  = RmInverse * (Pm + Bm.transpose()*slqp_->SmTrajectoryStock_[i][k]);
			control_matrix_t DmNullProjection = control_matrix_t::Identity()-DmProjected;
			sensitivityControllersStock[i].k_[k] = -DmNullProjection*Lm - CmProjected;

			// checking the numerical stability of the controller parameters
			try {
				if (sensitivityControllersStock[i].k_[k] != sensitivityControllersStock[i].k_[k])
					throw std::runtime_error("sensitivityController feedback gains are unstable.");
			}
			catch(const std::exception& error)  {
			    std::cerr << "what(): " << error.what() << " at time " << sensitivityControllersStock[i].time_[k] << " [sec]." << std::endl;
			}

		}  // end of k loop
	}  // end of i loop

}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * calculate the sensitivity of the control input increment to switchingTimes based on the LQ method
 * 		input & output:
 * 			+ sensitivityControllersStock
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateLQSensitivityControllerForward(
		std::vector<sensitivity_controller_t>& sensitivityControllersStock)  {

	LinearInterpolation<control_matrix_t,Eigen::aligned_allocator<control_matrix_t> > RmInverseFunc;
	LinearInterpolation<control_gain_matrix_t,Eigen::aligned_allocator<control_gain_matrix_t> > BmFunc;
	LinearInterpolation<nabla_input_matrix_t,Eigen::aligned_allocator<nabla_input_matrix_t> > nabla_RvFunc;

	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		// set data
		BmFunc.setTimeStamp(&slqp_->nominalTimeTrajectoriesStock_[i]);
		BmFunc.setData(&slqp_->BmTrajectoryStock_[i]);
		RmInverseFunc.setTimeStamp(&slqp_->nominalTimeTrajectoriesStock_[i]);
		RmInverseFunc.setData(&slqp_->RmInverseTrajectoryStock_[i]);
		nabla_RvFunc.setTimeStamp(&sensitivityTimeTrajectoryStock_[i]);
		nabla_RvFunc.setData(&nablaRvTrajectoryStock_[i]);

		// resizing the
		size_t N = sensitivityControllersStock[i].time_.size();
		sensitivityControllersStock[i].uff_.resize(N);

		for (size_t k=0; k<N; k++) {

			// time
			const double& t = sensitivityControllersStock[i].time_[k];

			// Bm
			control_gain_matrix_t Bm;
			BmFunc.interpolate(t, Bm);
			size_t greatestLessTimeStampIndex = BmFunc.getGreatestLessTimeStampIndex();
			// RmInverse
			control_matrix_t RmInverse;
			RmInverseFunc.interpolate(t, RmInverse, greatestLessTimeStampIndex);

			// nabla_Rv
			nabla_input_matrix_t nabla_Rv;
			nabla_RvFunc.interpolate(t, nabla_Rv);

			// nabla_sv
			nabla_output_matrix_t nabla_Sv;
			for (size_t j=0; j<NUM_SUBSYSTEMS-1; j++)
				nabla_Sv.col(j) = nablaSvTrajectoryStock_[i][k][j];

			sensitivityControllersStock[i].uff_[k] = -RmInverse*(nabla_Rv+Bm.transpose()*nabla_Sv);
		}
	}
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * calculate the sensitivity of the control input increment to switchingTimes based on the BVP method
 * 		inputs
 * 			+ switchingTimeIndex: the index of the switching time which the cost derivative will be calculated
 * 			+ SvTrajectoriesStock: sweeping method S vector
 *
 * 		output:
 * 			+ sensitivityControllersStock
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateBVPSensitivityControllerForward(
		const size_t& switchingTimeIndex,
		const std::vector<output_vector_array_t>& SvTrajectoriesStock,
		std::vector<sensitivity_controller_t>& sensitivityControllersStock)  {

	if (switchingTimeIndex < 1)  throw std::runtime_error("The initial switching time (startTime) is fixed and cost function derivative is not defined.");

	LinearInterpolation<control_matrix_t,Eigen::aligned_allocator<control_matrix_t> > RmInverseFunc;
	LinearInterpolation<control_gain_matrix_t,Eigen::aligned_allocator<control_gain_matrix_t> > BmFunc;

	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		// set data
		BmFunc.setTimeStamp(&slqp_->nominalTimeTrajectoriesStock_[i]);
		BmFunc.setData(&slqp_->BmTrajectoryStock_[i]);
		RmInverseFunc.setTimeStamp(&slqp_->nominalTimeTrajectoriesStock_[i]);
		RmInverseFunc.setData(&slqp_->RmInverseTrajectoryStock_[i]);

		// resizing the
		size_t N = slqp_->SsTimeTrajectoryStock_[i].size();
		sensitivityControllersStock[i].uff_.resize(N);

		for (size_t k=0; k<N; k++) {

			// time
			const double& t = slqp_->SsTimeTrajectoryStock_[i][k];

			// Bm
			control_gain_matrix_t Bm;
			BmFunc.interpolate(t, Bm);
			size_t greatestLessTimeStampIndex = BmFunc.getGreatestLessTimeStampIndex();
			// RmInverse
			control_matrix_t RmInverse;
			RmInverseFunc.interpolate(t, RmInverse, greatestLessTimeStampIndex);

			sensitivityControllersStock[i].uff_[k].col(switchingTimeIndex-1) = -RmInverse*Bm.transpose()*SvTrajectoriesStock[i][k]; //FIXME

		}  // end of k loop
	}  // end of i loop

}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * calculates the sensitivity of the rollout and LQ model to the switchingTimes
 * 		inputs:
 * 			+ sensitivityControllersStock
 * 		outputs:
 * 			+ sensitivityTimeTrajectoryStock: time stamp
 * 			+ nablaOutputTrajectoryStock: dy
 * 			+ nablaInputTrajectoryStock: du
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::rolloutSensitivity2SwitchingTime(
		const std::vector<sensitivity_controller_t>& sensitivityControllersStock,
		std::vector<scalar_array_t>& sensitivityTimeTrajectoryStock,
		std::vector<nabla_output_matrix_array_t>& nablaOutputTrajectoryStock,
		std::vector<nabla_input_matrix_array_t>& nablaInputTrajectoryStock)  {

	std::shared_ptr<RolloutSensitivityEquations_t> rolloutSensitivityEquationsPtr( new RolloutSensitivityEquations_t() );
	nabla_output_vector_t nabla_YmInit;
	RolloutSensitivityEquations_t::convert2Vector(nabla_output_matrix_t::Zero(), nabla_YmInit);

	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		// initialize subsystem i
		slqp_->getSubsystemDynamicsPtrStock()[i]->initializeModel(switchingTimes_, slqp_->nominalStateTrajectoriesStock_[i].front(), i, "GSLQP");

		rolloutSensitivityEquationsPtr->setData(i, switchingTimes_, &sensitivityControllersStock[i],
				&slqp_->nominalTimeTrajectoriesStock_[i], &nominalOutputTimeDerivativeTrajectoriesStock_[i],
				&slqp_->AmTrajectoryStock_[i], &slqp_->BmTrajectoryStock_[i]);

		// integrating
		scalar_array_t normalizedSensitivityTimeTrajectory;
		nabla_output_vector_array_t sensitivityOutputTrajectory;
		ODE45<(NUM_SUBSYSTEMS-1)*OUTPUT_DIM> ode45(rolloutSensitivityEquationsPtr);
		ode45.integrate(nabla_YmInit, 0.0, 1.0, sensitivityOutputTrajectory, normalizedSensitivityTimeTrajectory,
				1e-3, options_.AbsTolODE_, options_.RelTolODE_);

		// denormalizing time and constructing SensitivityStateTrajectory and computing control trajectory sensitivity for subsystem i
		int N = normalizedSensitivityTimeTrajectory.size();
		sensitivityTimeTrajectoryStock[i].resize(N);
		nablaOutputTrajectoryStock[i].resize(N);
		nablaInputTrajectoryStock[i].resize(N);
		for (int k=0; k<N; k++) {

			sensitivityTimeTrajectoryStock[i][k] = switchingTimes_[i] + (switchingTimes_[i+1]-switchingTimes_[i])*normalizedSensitivityTimeTrajectory[k];
			RolloutSensitivityEquations_t::convert2Matrix(sensitivityOutputTrajectory[k], nablaOutputTrajectoryStock[i][k]);
			rolloutSensitivityEquationsPtr->computeInputSensitivity(sensitivityTimeTrajectoryStock[i][k], nablaOutputTrajectoryStock[i][k],
					nablaInputTrajectoryStock[i][k]);
		}

		// reset the initial state
		nabla_YmInit = sensitivityOutputTrajectory.back();

	}  // end of i loop
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * approximate nominal LQ problem sensitivity to switching times
 * 		modifies:
 * 			+ nablaqTrajectoryStock_: dq
 * 			+ nablaQvTrajectoryStock_: dQv
 * 			+ nablaRvTrajectoryStock_: dRv
 * 			+ nablaqFinal_: dq_f
 * 			+ nablaQvFinal_: dQv_f
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::approximateNominalLQPSensitivity2SwitchingTime() {

	// calculate nabla_q, nabla_Qv, nabla_Rv
	LinearInterpolation<output_vector_t,Eigen::aligned_allocator<output_vector_t> > QvFunc;
	LinearInterpolation<state_matrix_t,Eigen::aligned_allocator<state_matrix_t> > QmFunc;
	LinearInterpolation<control_vector_t,Eigen::aligned_allocator<control_vector_t> > RvFunc;
	LinearInterpolation<control_matrix_t,Eigen::aligned_allocator<control_matrix_t> > RmFunc;
	LinearInterpolation<control_feedback_t,Eigen::aligned_allocator<control_feedback_t> > PmFunc;

	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		QvFunc.setTimeStamp( &(slqp_->nominalTimeTrajectoriesStock_[i]) );
		QvFunc.setData( &(slqp_->QvTrajectoryStock_[i]) );
		QmFunc.setTimeStamp( &(slqp_->nominalTimeTrajectoriesStock_[i]) );
		QmFunc.setData( &(slqp_->QmTrajectoryStock_[i]) );
		RvFunc.setTimeStamp( &(slqp_->nominalTimeTrajectoriesStock_[i]) );
		RvFunc.setData( &(slqp_->RvTrajectoryStock_[i]) );
		RmFunc.setTimeStamp( &(slqp_->nominalTimeTrajectoriesStock_[i]) );
		RmFunc.setData( &(slqp_->RmTrajectoryStock_[i]) );
		PmFunc.setTimeStamp( &(slqp_->nominalTimeTrajectoriesStock_[i]) );
		PmFunc.setData( &(slqp_->PmTrajectoryStock_[i]) );

		int N = sensitivityTimeTrajectoryStock_[i].size();
		nablaqTrajectoryStock_[i].resize(N);
		nablaQvTrajectoryStock_[i].resize(N);
		nablaRvTrajectoryStock_[i].resize(N);

		for (int k=0; k<N; k++) {

			control_matrix_t Rm;
			RmFunc.interpolate(sensitivityTimeTrajectoryStock_[i][k], Rm);
			size_t greatestLessTimeStampIndex = RmFunc.getGreatestLessTimeStampIndex();
			output_vector_t Qv;
			QvFunc.interpolate(sensitivityTimeTrajectoryStock_[i][k], Qv, greatestLessTimeStampIndex);
			state_matrix_t Qm;
			QmFunc.interpolate(sensitivityTimeTrajectoryStock_[i][k], Qm, greatestLessTimeStampIndex);
			control_vector_t Rv;
			RvFunc.interpolate(sensitivityTimeTrajectoryStock_[i][k], Rv, greatestLessTimeStampIndex);
			control_feedback_t Pm;
			PmFunc.interpolate(sensitivityTimeTrajectoryStock_[i][k], Pm, greatestLessTimeStampIndex);

			nablaqTrajectoryStock_[i][k]  = Qv.transpose()*nablaOutputTrajectoryStock_[i][k] + Rv.transpose()*nablaInputTrajectoryStock_[i][k];
			nablaQvTrajectoryStock_[i][k] = Qm*nablaOutputTrajectoryStock_[i][k] + Pm.transpose()*nablaInputTrajectoryStock_[i][k];
			nablaRvTrajectoryStock_[i][k] = Pm*nablaOutputTrajectoryStock_[i][k] + Rm*nablaInputTrajectoryStock_[i][k];
		}

		if (i==NUM_SUBSYSTEMS-1)  {
			output_vector_t Qv = slqp_->QvTrajectoryStock_[NUM_SUBSYSTEMS-1].back();
			state_matrix_t Qm  = slqp_->QmTrajectoryStock_[NUM_SUBSYSTEMS-1].back();

			nablaqFinal_  = Qv.transpose()*nablaOutputTrajectoryStock_[NUM_SUBSYSTEMS-1].back();
			nablaQvFinal_ = Qm*nablaOutputTrajectoryStock_[NUM_SUBSYSTEMS-1].back();
		}
	}  // end of i loop
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * calculate the nominal output time derivative
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateOutputTimeDerivative()  {

	for (size_t i=0; i<NUM_SUBSYSTEMS; i++) {

		slqp_->getSubsystemDynamicsPtrStock()[i]->initializeModel(switchingTimes_, slqp_->nominalStateTrajectoriesStock_[i].front(), i, "GSLQP");

		size_t N = slqp_->nominalTimeTrajectoriesStock_[i].size();
		nominalOutputTimeDerivativeTrajectoriesStock_[i].resize(N);

		for (size_t k=0; k<N; k++) {
			const scalar_t& 	    t = slqp_->nominalTimeTrajectoriesStock_[i][k];
			const state_vector_t&   x = slqp_->nominalStateTrajectoriesStock_[i][k];
			const control_vector_t& u = slqp_->nominalInputTrajectoriesStock_[i][k];
			state_vector_t dxdt;
			slqp_->getSubsystemDynamicsPtrStock()[i]->computeDerivative(t, x, u, dxdt);
			nominalOutputTimeDerivativeTrajectoriesStock_[i][k] = slqp_->getSubsystemDynamicsPtrStock()[i]->computeOutputStateDerivative(t, x, u) * dxdt;

		}  // end of k loop
	}  // end of i loop

}



/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * solve sensitivity BVP (the boundary value problem of the sensitivity equations for a given switching time)
 * 		inputs
 * 			+ switchingTimeIndex: the index of the switching time which the cost derivative will be calculated
 * 			+ timeTrajectoriesStock: time stamp
 *
 * 		outputs
 * 			+ MmTrajectoriesStock: sweeping method M matrix
 * 			+ SvTrajectoriesStock: sweeping method S vector
 *
 * 		uses
 * 			+ linearized dynamics
 * 			+ quadratized cost
 *
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::solveSensitivityBVP(
		const size_t& switchingTimeIndex,
		const std::vector<scalar_array_t>& timeTrajectoriesStock,
		std::vector<state_matrix_array_t>& MmTrajectoriesStock,
		std::vector<output_vector_array_t>& SvTrajectoriesStock)  {

	if (switchingTimeIndex < 1)  throw std::runtime_error("The initial switching time (startTime) is fixed and cost function derivative is not defined.");

	// calculate the BVP coefficients
	output_vector_array_t bvpGvPositivetraTrajectory(slqp_->nominalTimeTrajectoriesStock_[switchingTimeIndex-1].size());
	output_vector_array_t bvpQvPositivetraTrajectory(slqp_->nominalTimeTrajectoriesStock_[switchingTimeIndex-1].size());
	output_vector_array_t bvpGvNegativetraTrajectory(slqp_->nominalTimeTrajectoriesStock_[switchingTimeIndex].size());
	output_vector_array_t bvpQvNegativetraTrajectory(slqp_->nominalTimeTrajectoriesStock_[switchingTimeIndex].size());

	const double scalingFactor = 1/(switchingTimes_[switchingTimeIndex]-switchingTimes_[switchingTimeIndex-1]);

	for (size_t k=0; k<slqp_->nominalTimeTrajectoriesStock_[switchingTimeIndex-1].size(); k++) {
		bvpGvPositivetraTrajectory[k] = scalingFactor * nominalOutputTimeDerivativeTrajectoriesStock_[switchingTimeIndex-1][k];
		bvpQvPositivetraTrajectory[k] = scalingFactor * (slqp_->QvTrajectoryStock_[switchingTimeIndex-1][k]+
				slqp_->AmTrajectoryStock_[switchingTimeIndex-1][k].transpose()*slqp_->nominalcostateTrajectoriesStock_[switchingTimeIndex-1][k]);
	}

	for (size_t k=0; k<slqp_->nominalTimeTrajectoriesStock_[switchingTimeIndex].size(); k++) {
		bvpGvNegativetraTrajectory[k] = -scalingFactor * nominalOutputTimeDerivativeTrajectoriesStock_[switchingTimeIndex][k];
		bvpQvNegativetraTrajectory[k] = -scalingFactor * (slqp_->QvTrajectoryStock_[switchingTimeIndex][k]+
				slqp_->AmTrajectoryStock_[switchingTimeIndex][k].transpose()*slqp_->nominalcostateTrajectoriesStock_[switchingTimeIndex][k]);
	}

	SolveBVP<OUTPUT_DIM, INPUT_DIM> bvpSolver;
	output_vector_t SvFinal = output_vector_t::Zero();
	state_matrix_t  MmFinal = slqp_->QmFinal_;

	for (int i=NUM_SUBSYSTEMS-1; i>=0; i--) {

		const output_vector_array_t* GvPtr;
		const output_vector_array_t* QvPtr;
		if (i==switchingTimeIndex-1) {
			GvPtr = &bvpGvPositivetraTrajectory;
			QvPtr = &bvpQvPositivetraTrajectory;
		} else if (i==switchingTimeIndex) {
			GvPtr = &bvpGvNegativetraTrajectory;
			QvPtr = &bvpQvNegativetraTrajectory;
		} else {
			GvPtr = NULL;
			QvPtr = NULL;
		}

		// set the general BVP solver coefficient
		bvpSolver.setData(&slqp_->nominalTimeTrajectoriesStock_[i],
				&slqp_->AmConstrainedTrajectoryStock_[i], NULL,  &slqp_->BmTrajectoryStock_[i], GvPtr,
				QvPtr, &slqp_->QmConstrainedTrajectoryStock_[i], &slqp_->PmTrajectoryStock_[i],
				NULL, &slqp_->RmConstrainedTrajectoryStock_[i], &slqp_->RmInverseTrajectoryStock_[i]);

		// solve BVP for the given time trajectory
		bvpSolver.solve(timeTrajectoriesStock[i], SvFinal, MmFinal,
				MmTrajectoriesStock[i], SvTrajectoriesStock[i],
				options_.AbsTolODE_, options_.RelTolODE_);

		//set the final value of the previous subsystem solver to the starting time value of the current subsystem's solution
		SvFinal = SvTrajectoriesStock[i].front();
		MmFinal = MmTrajectoriesStock[i].front();

		for (size_t k=0; k<timeTrajectoriesStock[i].size(); k++)
			if (!MmTrajectoriesStock[i][k].isApprox(slqp_->SmTrajectoryStock_[i][k], 1e-3)) {
				std::cerr << "In solveSensitivityBVP, Mm and Sm do not match." << std::endl;
				std::cerr << "Mm[" << i << "][" << k << "]\n" <<  MmTrajectoriesStock[i][k] << std::endl;
				std::cerr << "Sm[" << i << "][" << k << "]\n" <<  slqp_->SmTrajectoryStock_[i][k] << std::endl;
			}

	}  // end of i loop

}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * calculates cost function derivative based on BVP solution
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateBVPCostFunctionDerivative(
		Eigen::Matrix<double,NUM_SUBSYSTEMS-1,1>& costFunctionDerivative)  {

	// final time
	costFunctionDerivative = nablaOutputTrajectoryStock_.back().back().transpose() * slqp_->QvFinal_;

	for (size_t i=0; i<NUM_SUBSYSTEMS; i++) {

		LinearInterpolation<eigen_scalar_t, Eigen::aligned_allocator<eigen_scalar_t> > qFunc(
				&slqp_->nominalTimeTrajectoriesStock_[i], &slqp_->qTrajectoryStock_[i]);
		LinearInterpolation<output_vector_t, Eigen::aligned_allocator<output_vector_t> > QvFunc(
				&slqp_->nominalTimeTrajectoriesStock_[i], &slqp_->QvTrajectoryStock_[i]);
		LinearInterpolation<control_vector_t, Eigen::aligned_allocator<control_vector_t> > RvFunc(
				&slqp_->nominalTimeTrajectoriesStock_[i], &slqp_->RvTrajectoryStock_[i]);

		Eigen::Matrix<double,NUM_SUBSYSTEMS-1,1> previousIntermediatecostFunctionDev;
		Eigen::Matrix<double,NUM_SUBSYSTEMS-1,1> currentIntermediatecostFunctionDev;

		for (int k=0; k<sensitivityTimeTrajectoryStock_[i].size(); k++) {

			const double& t = sensitivityTimeTrajectoryStock_[i][k];

			output_vector_t Qv;
			QvFunc.interpolate(t, Qv);
			size_t greatestLessTimeStampIndex = QvFunc.getGreatestLessTimeStampIndex();
			control_vector_t Rv;
			RvFunc.interpolate(t, Rv, greatestLessTimeStampIndex);
			eigen_scalar_t q;
			qFunc.interpolate(t, q, greatestLessTimeStampIndex);

			Eigen::Matrix<double,NUM_SUBSYSTEMS-1,1> coeff;
			for (int j=1; j<NUM_SUBSYSTEMS; j++)
				if (i==j-1)
					coeff(j-1) = +1.0;
				else if (i==j)
					coeff(j-1) = -1.0;
				else
					coeff(j-1) = 0.0;

			if (k>0)
				previousIntermediatecostFunctionDev = currentIntermediatecostFunctionDev;

			currentIntermediatecostFunctionDev = coeff*q/(switchingTimes_[i+1]-switchingTimes_[i]) +
					(nablaOutputTrajectoryStock_[i][k].transpose()*Qv + nablaInputTrajectoryStock_[i][k].transpose()*Rv);

			if (k>0)
				costFunctionDerivative += 0.5*(sensitivityTimeTrajectoryStock_[i][k]-sensitivityTimeTrajectoryStock_[i][k-1]) *
					(currentIntermediatecostFunctionDev+previousIntermediatecostFunctionDev);

		}  // end of k loop
	}  // end of i loop

}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * run the SLQ algorithm for a given state and switching times
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::runLQBasedMethod(const state_vector_t& initState, const std::vector<scalar_t>& switchingTimes)  {

	switchingTimes_ = switchingTimes;
	initState_ = initState;

	// make sure that the minimum difference between to successive switching times is at least options_.minSimulationTimeDuration_
	for (size_t i=0; i<NUM_SUBSYSTEMS; i++)
		if (switchingTimes_[i+1]-switchingTimes_[i] < options_.minSimulationTimeDuration_) {
			if (i+1 == NUM_SUBSYSTEMS)
				std::cerr << "WARNING: The minimum simulation time between the last subsystem's stratTime and finalTime should be at least "
						<< options_.minSimulationTimeDuration_ << "." << std::endl;
			switchingTimes_[i+1] = switchingTimes_[i]+options_.minSimulationTimeDuration_;
		}

	// run the SLQ algorithm
	slqp_->run(initState, switchingTimes_);

	if (options_.dispayGSLQP_)  std::cerr << "\n#### Calculating cost function sensitivity ..." << std::endl;

	// calculate output time derivative
	calculateOutputTimeDerivative();

	// calculate sensitivity controller feedback part
	calculateSensitivityControllerFeedback(nominalSensitivityControllersStock_);
	// set sensitivity controller feedforward part to zero
	for (size_t i=0; i<NUM_SUBSYSTEMS; i++) {
		size_t N = nominalSensitivityControllersStock_[i].time_.size();
		nominalSensitivityControllersStock_[i].uff_.resize(N);
		for (size_t k=0; k<N; k++)  nominalSensitivityControllersStock_[i].uff_[k].setZero();
	}


	for (size_t j=0; j<3; j++) {

		// calculate nominal rollout sensitivity to switching times
		rolloutSensitivity2SwitchingTime(nominalSensitivityControllersStock_,
				sensitivityTimeTrajectoryStock_, nablaOutputTrajectoryStock_, nablaInputTrajectoryStock_);

		// approximate the nominal LQP sensitivity to switching times
		approximateNominalLQPSensitivity2SwitchingTime();

		// solve Riccati equations
		solveSensitivityRiccatiEquations(0.0 /*learningRateStar*/); // prevents the changes in the nominal trajectories and just update the gains

		// calculate sensitivity controller feedforward part
		calculateLQSensitivityControllerForward(nominalSensitivityControllersStock_);
	}

	// transform from local value function derivatives to global representation
	transformLocalValueFuntionDerivative2Global();

	// calculate the cost function derivatives w.r.t. switchingTimes
	getValueFuntionDerivative(slqp_->nominalOutputTrajectoriesStock_[0][0], nominalCostFuntionDerivative_);

}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * run the SLQ algorithm for a given state and switching times based on the BVP method
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::run(const state_vector_t& initState, const std::vector<scalar_t>& switchingTimes)  {

	if (options_.useLQForDerivatives_==true) {
		runLQBasedMethod(initState, switchingTimes);
		return;
	}

	switchingTimes_ = switchingTimes;
	initState_ = initState;

	// make sure that the minimum difference between to successive switching times is at least options_.minSimulationTimeDuration_
	for (size_t i=0; i<NUM_SUBSYSTEMS; i++)
		if (switchingTimes_[i+1]-switchingTimes_[i] < options_.minSimulationTimeDuration_) {
			if (i+1 == NUM_SUBSYSTEMS)
				std::cerr << "WARNING: The minimum simulation time between the last subsystem's stratTime and finalTime should be at least "
					<< options_.minSimulationTimeDuration_ << "." << std::endl;
			switchingTimes_[i+1] = switchingTimes_[i]+options_.minSimulationTimeDuration_;
		}

	// run the SLQ algorithm
	slqp_->run(initState, switchingTimes_);

	// calculate output time derivative
	calculateOutputTimeDerivative();

//	{
//		SolveBVP<OUTPUT_DIM, INPUT_DIM> bvpSolver;
//		std::vector<state_matrix_array_t> MmTrajectoriesStock(NUM_SUBSYSTEMS);
//		std::vector<state_vector_array_t> SvTrajectoriesStock(NUM_SUBSYSTEMS);
//		state_vector_t SvFinal = slqp_->QvFinal_;
//		state_matrix_t MmFinal = slqp_->QmFinal_;
//		for (int i=NUM_SUBSYSTEMS-1; i>=0; i--) {
//
//			bvpSolver.setData(&slqp_->nominalTimeTrajectoriesStock_[i],
//					&slqp_->AmConstrainedTrajectoryStock_[i], NULL,  &slqp_->BmTrajectoryStock_[i], NULL,
//					&slqp_->QvConstrainedTrajectoryStock_[i], &slqp_->QmConstrainedTrajectoryStock_[i], &slqp_->PmTrajectoryStock_[i],
//					&slqp_->RvTrajectoryStock_[i], &slqp_->RmConstrainedTrajectoryStock_[i], &slqp_->RmInverseTrajectoryStock_[i]);
//
//			bvpSolver.solve(slqp_->SsTimeTrajectoryStock_[i], SvFinal, MmFinal,
//					MmTrajectoriesStock[i], SvTrajectoriesStock[i],
//					options_.AbsTolODE_, options_.RelTolODE_);
//
//			SvFinal = SvTrajectoriesStock[i].front();
//			MmFinal = MmTrajectoriesStock[i].front();
//
//			for (size_t k=0; k<slqp_->SsTimeTrajectoryStock_[i].size(); k++) {
//				if (!MmTrajectoriesStock[i][k].isApprox(slqp_->SmTrajectoryStock_[i][k], 1e-3)) {
//					std::cout << "Mm[" << i << "][" << k << "]\n" <<  MmTrajectoriesStock[i][k] << std::endl;
//					std::cout << "Sm[" << i << "][" << k << "]\n" <<  slqp_->SmTrajectoryStock_[i][k] << std::endl;
//				}
//				if (!SvTrajectoriesStock[i][k].isApprox(slqp_->SvTrajectoryStock_[i][k], 1e-3)) {
//					std::cout << "Sv[" << i << "][" << k << "]\n" <<  SvTrajectoriesStock[i][k] << std::endl;
//					std::cout << "Sv[" << i << "][" << k << "]\n" <<  slqp_->SvTrajectoryStock_[i][k] << std::endl;
//				}
//			}
//		}
//	}

	// calculate sensitivity controller feedback part
	calculateSensitivityControllerFeedback(nominalSensitivityControllersStock_);

	// for each switching time solve BVP
	for (int j=1; j<NUM_SUBSYSTEMS; j++)  {

		std::vector<state_matrix_array_t>  MmTrajectoriesStock(NUM_SUBSYSTEMS);
		std::vector<output_vector_array_t> SvTrajectoriesStock(NUM_SUBSYSTEMS);

		// solve boundary value problem of the sensitivity equations for switching time j
		solveSensitivityBVP(j, slqp_->SsTimeTrajectoryStock_, MmTrajectoriesStock, SvTrajectoriesStock);

		// calculate sensitivity controller feedforward part
		calculateBVPSensitivityControllerForward(j, SvTrajectoriesStock, nominalSensitivityControllersStock_);

	} // end of j loop

	// calculate nominal rollout sensitivity to switching times
	rolloutSensitivity2SwitchingTime(nominalSensitivityControllersStock_,
			sensitivityTimeTrajectoryStock_, nablaOutputTrajectoryStock_, nablaInputTrajectoryStock_);

	// calculate the cost function derivatives w.r.t. switchingTimes
	calculateBVPCostFunctionDerivative(nominalCostFuntionDerivative_);

}


} // namespace ocs2

