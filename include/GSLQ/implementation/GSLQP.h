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

	slqp_.rollout(initState,
			controllersStock, timeTrajectoriesStock, stateTrajectoriesStock, inputTrajectoriesStock, outputTrajectoriesStock);

}


/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::rollout(const state_vector_t& initState,
		const std::vector<controller_t>& controllersStock,
		std::vector<scalar_array_t>& timeTrajectoriesStock,
		std::vector<state_vector_array_t>& stateTrajectoriesStock,
		std::vector<control_vector_array_t>& inputTrajectoriesStock)  {

	slqp_.rollout(initState,
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

	slqp_.rollout(initState,
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

	slqp_.calculateCostFunction(timeTrajectoriesStock, outputTrajectoriesStock, inputTrajectoriesStock,
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

	slqp_.calculateMeritFunction(timeTrajectoriesStock, nc1TrajectoriesStock, EvTrajectoryStock, lagrangeTrajectoriesStock, totalCost,
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

	return slqp_.calculateConstraintISE(timeTrajectoriesStock, nc1TrajectoriesStock, EvTrajectoriesStock,
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

		nominalOutputFunc.setTimeStamp( &(slqp_.nominalTimeTrajectoriesStock_[i]) );
		nominalOutputFunc.setData( &(slqp_.nominalOutputTrajectoriesStock_[i]) );

		for (int k=0; k<slqp_.SsTimeTrajectoryStock_[i].size(); k++) {

			output_vector_t nominalOutput;
			nominalOutputFunc.interpolate(slqp_.SsTimeTrajectoryStock_[i][k], nominalOutput);

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

	slqp_.getController(controllersStock);
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

	slqp_.getValueFuntion(time, output, valueFuntion);
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * calculate the cost function's derivatives at the initial time
 * 		inputs
 * 			+ initOutput: initial output
 *
 * 		output:
 * 			+ costFuntionDerivative: cost function' derivatives w.r.t. switchingTimes for given initial output vector
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getCostFuntionDerivative(const output_vector_t& initOutput,
		Eigen::Matrix<double,NUM_SUBSYSTEMS-1,1>& costFuntionDerivative)  {


	for (int j=0; j<NUM_SUBSYSTEMS-1; j++)  {

		state_matrix_t dSm  = nablaSmTrajectoryStock_[0][0][j];
		output_vector_t dSv = nablaSvTrajectoryStock_[0][0][j];
		eigen_scalar_t ds   = nablasTrajectoryStock_[0][0][j];

		costFuntionDerivative(j) = (ds + initOutput.transpose()*dSv + 0.5*initOutput.transpose()*dSm*initOutput).eval()(0);
	}
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

	slqp_.getNominalTrajectories(nominalTimeTrajectoriesStock, nominalStateTrajectoriesStock,
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
 *
 * 		modifies:
 * 			V(t,y) = y^T*Sm*y + y^T*(Sv) + s
 * 			dV(t,y) = y^T*dSm*y + y^T*(dSv) + ds
 * 			+ SsTimeTrajectoryStock_: time stamp
 * 			+ SmTrajectoryStock_: Sm matrix
 * 			+ SvTrajectoryStock_: Sv vector
 * 			+ sTrajectoryStock_: s scalar
 * 			+ nablaSmTrajectoryStock_: dSm
 * 			+ nablaSvTrajectoryStock_: dSv
 * 			+ nablasTrajectoryStock_: ds
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::solveFullSequentialRiccatiEquations(const scalar_t& learningRate)  {

	// final value for the last Riccati equations
	typename FullRiccatiEquations_t::all_s_vector_t allSsFinal;
	nabla_Sm_t nablaQmFinal;  nablaQmFinal.fill(state_matrix_t::Zero());
	nabla_Sv_t nablaQvFinal;  nablaQvFinal.fill(output_vector_t::Zero());
	nabla_s_t  nablaqFinal;   nablaqFinal.fill(eigen_scalar_t::Zero());
	for (int i=0; i<NUM_SUBSYSTEMS-1; i++) {
		nablaQvFinal[i] = nablaQvFinal_.template block<OUTPUT_DIM,1>(0,i);
		nablaqFinal[i](0) = nablaqFinal_(i);
	}
	FullRiccatiEquations_t::convert2Vector(slqp_.QmFinal_, slqp_.QvFinal_, slqp_.qFinal_, nablaQmFinal, nablaQvFinal, nablaqFinal, allSsFinal);

	for (int i=NUM_SUBSYSTEMS-1; i>=0; i--) {

		// set data for Riccati equations
		auto riccatiEquationsPtr = std::make_shared<FullRiccatiEquations_t>();
		riccatiEquationsPtr->setData(learningRate,
				i, switchingTimes_[i], switchingTimes_[i+1],
				&slqp_.nominalTimeTrajectoriesStock_[i],
				&slqp_.AmConstrainedTrajectoryStock_[i], &slqp_.BmTrajectoryStock_[i],
				&slqp_.qTrajectoryStock_[i], &slqp_.QvConstrainedTrajectoryStock_[i], &slqp_.QmConstrainedTrajectoryStock_[i],
				&slqp_.RvTrajectoryStock_[i], &slqp_.RmInverseTrajectoryStock_[i], &slqp_.RmConstrainedTrajectoryStock_[i], &slqp_.PmTrajectoryStock_[i],
				&sensitivityTimeTrajectoryStock_[i], &nablaqTrajectoryStock_[i],
				&nablaQvTrajectoryStock_[i], &nablaRvTrajectoryStock_[i]);

		// integrating the Riccati equations
		ODE45<FullRiccatiEquations_t::S_DIM_*NUM_SUBSYSTEMS> ode45(riccatiEquationsPtr);
		std::vector<double> normalizedTimeTrajectory;
		std::vector<typename FullRiccatiEquations_t::all_s_vector_t, Eigen::aligned_allocator<typename FullRiccatiEquations_t::all_s_vector_t> > allSsTrajectory;
		ode45.integrate(allSsFinal, i, i+1, allSsTrajectory, normalizedTimeTrajectory,
				1e-3, options_.AbsTolODE_, options_.RelTolODE_);

		// denormalizing time and constructing 'Sm', 'Sv', and 's'
		int N = normalizedTimeTrajectory.size();
		slqp_.SsTimeTrajectoryStock_[i].resize(N);
		slqp_.SmTrajectoryStock_[i].resize(N);
		slqp_.SvTrajectoryStock_[i].resize(N);
		slqp_.sTrajectoryStock_[i].resize(N);
		nablaSmTrajectoryStock_[i].resize(N);
		nablaSvTrajectoryStock_[i].resize(N);
		nablasTrajectoryStock_[i].resize(N);
		for (int k=0; k<normalizedTimeTrajectory.size(); k++) {

			FullRiccatiEquations_t::convert2Matrix(allSsTrajectory[N-1-k],
					slqp_.SmTrajectoryStock_[i][k], slqp_.SvTrajectoryStock_[i][k], slqp_.sTrajectoryStock_[i][k],
					nablaSmTrajectoryStock_[i][k], nablaSvTrajectoryStock_[i][k], nablasTrajectoryStock_[i][k]);
			slqp_.SsTimeTrajectoryStock_[i][k] = (switchingTimes_[i]-switchingTimes_[i+1])*(normalizedTimeTrajectory[N-1-k]-i) + switchingTimes_[i+1];
		}  // end of k loop

		// reset the final value for next Riccati equation
		allSsFinal = allSsTrajectory.back();

	}  // end of i loop

}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * calculate the sensitivity of the control input increment to switchingTimes
 * 		inputs:
 * 			+ nablaUffTimeTrajectoryStock: time stamp
 * 			+ nablaUffTrajectoryStock: sensitivity of the control input increment
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::inputIncrementSensitivity2SwitchingTime(
		std::vector<scalar_array_t>& nablaUffTimeTrajectoryStock,
		std::vector<nabla_input_matrix_array_t>& nablaUffTrajectoryStock)  {

	LinearInterpolation<control_gain_matrix_t,Eigen::aligned_allocator<control_gain_matrix_t> > BmFunc;
	LinearInterpolation<control_matrix_t,Eigen::aligned_allocator<control_matrix_t> > RmInverseFunc;

	LinearInterpolation<nabla_input_matrix_t,Eigen::aligned_allocator<nabla_input_matrix_t> > nabla_RvFunc;

	nablaUffTimeTrajectoryStock = slqp_.SsTimeTrajectoryStock_;

	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		// set data
		BmFunc.setTimeStamp(&slqp_.nominalTimeTrajectoriesStock_[i]);
		BmFunc.setData(&slqp_.BmTrajectoryStock_[i]);
		RmInverseFunc.setTimeStamp(&slqp_.nominalTimeTrajectoriesStock_[i]);
		RmInverseFunc.setData(&slqp_.RmInverseTrajectoryStock_[i]);
		nabla_RvFunc.setTimeStamp(&sensitivityTimeTrajectoryStock_[i]);
		nabla_RvFunc.setData(&nablaRvTrajectoryStock_[i]);

		// resizing the
		size_t N = nablaUffTimeTrajectoryStock[i].size();
		nablaUffTrajectoryStock[i].resize(N);

		for (size_t k=0; k<N; k++) {

			// time
			const double& t = slqp_.SsTimeTrajectoryStock_[i][k];

			// nabla_sv
			nabla_output_matrix_t nabla_Sv;
			for (size_t j=0; j<NUM_SUBSYSTEMS-1; j++)
				nabla_Sv.col(j) = nablaSvTrajectoryStock_[i][k][j];

			// Bm
			control_gain_matrix_t Bm;
			BmFunc.interpolate(t, Bm);
			int greatestLessTimeStampIndex = BmFunc.getGreatestLessTimeStampIndex();
			// RmInverse
			control_matrix_t RmInverse;
			RmInverseFunc.interpolate(t, RmInverse, greatestLessTimeStampIndex);

			// nabla_Rv
			nabla_input_matrix_t nabla_Rv;
			nabla_RvFunc.interpolate(t, nabla_Rv);

			nablaUffTrajectoryStock[i][k] = -RmInverse*(nabla_Rv+Bm.transpose()*nabla_Sv);
		}
	}
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * calculates the sensitivity of the rollout and LQ model to the switchingTimes
 * 		input:
 * 			+ nablaSvUpdated: true if dSv is already updated
 *
 * 		modifies:
 * 			+ sensitivityTimeTrajectoryStock_: time stamp
 * 			+ nablaOutputTrajectoryStock_: dy
 * 			+ nablaInputTrajectoryStock_: du
 * 			+ nablaqTrajectoryStock_: dq
 * 			+ nablaQvTrajectoryStock_: dQv
 * 			+ nablaRvTrajectoryStock_: dRv
 * 			+ nablaqFinal_: dq_f
 * 			+ nablaQvFinal_: dQv_f
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::rolloutSensitivity2SwitchingTime(bool nablaSvUpdated)  {


	std::vector<scalar_array_t> nablaUffTimeTrajectoryStock(NUM_SUBSYSTEMS);
	std::vector<nabla_input_matrix_array_t> nablaUffTrajectoryStock(NUM_SUBSYSTEMS);
	if (nablaSvUpdated==true)
		inputIncrementSensitivity2SwitchingTime(nablaUffTimeTrajectoryStock, nablaUffTrajectoryStock);

	auto rolloutSensitivityEquationsPtr = std::make_shared<RolloutSensitivityEquations_t>();

	nabla_output_vector_t nabla_YmInit;
	RolloutSensitivityEquations_t::convert2Vector(nabla_output_matrix_t::Zero(), nabla_YmInit);

	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		// initialize subsystem i
		slqp_.subsystemDynamicsPtrStock_[i]->initializeModel(slqp_.nominalTimeTrajectoriesStock_[i].front(),
				slqp_.nominalStateTrajectoriesStock_[i].front(), slqp_.nominalTimeTrajectoriesStock_[i].back(), "GSLQP");

		if (nablaSvUpdated==true)
			rolloutSensitivityEquationsPtr->setData(i, switchingTimes_, slqp_.subsystemDynamicsPtrStock_[i], &slqp_.nominalControllersStock_[i],
					&slqp_.nominalTimeTrajectoriesStock_[i], &slqp_.nominalStateTrajectoriesStock_[i], &slqp_.nominalInputTrajectoriesStock_[i],
					&slqp_.AmConstrainedTrajectoryStock_[i], &slqp_.BmTrajectoryStock_[i],
					&nablaUffTimeTrajectoryStock[i], &nablaUffTrajectoryStock[i]);
		else
			rolloutSensitivityEquationsPtr->setData(i, switchingTimes_, slqp_.subsystemDynamicsPtrStock_[i], &slqp_.nominalControllersStock_[i],
								&slqp_.nominalTimeTrajectoriesStock_[i], &slqp_.nominalStateTrajectoriesStock_[i], &slqp_.nominalInputTrajectoriesStock_[i],
								&slqp_.AmConstrainedTrajectoryStock_[i], &slqp_.BmTrajectoryStock_[i]);

		// integrating
		scalar_array_t normalizedSensitivityTimeTrajectory;
		nabla_output_vector_array_t sensitivityOutputTrajectory;
		ODE45<(NUM_SUBSYSTEMS-1)*OUTPUT_DIM> ode45(rolloutSensitivityEquationsPtr);
		ode45.integrate(nabla_YmInit, i, i+1, sensitivityOutputTrajectory, normalizedSensitivityTimeTrajectory,
				1e-3, options_.AbsTolODE_, options_.RelTolODE_);

		// denormalizing time and constructing SensitivityStateTrajectory and computing control trajectory sensitivity for subsystem i
		int N = normalizedSensitivityTimeTrajectory.size();
		sensitivityTimeTrajectoryStock_[i].resize(N);
		nablaOutputTrajectoryStock_[i].resize(N);
		nablaInputTrajectoryStock_[i].resize(N);
		for (int k=0; k<N; k++) {

			sensitivityTimeTrajectoryStock_[i][k] = switchingTimes_[i] + (switchingTimes_[i+1]-switchingTimes_[i])*(normalizedSensitivityTimeTrajectory[k]-i);
			RolloutSensitivityEquations_t::convert2Matrix(sensitivityOutputTrajectory[k], nablaOutputTrajectoryStock_[i][k]);
			rolloutSensitivityEquationsPtr->computeInputSensitivity(sensitivityTimeTrajectoryStock_[i][k], nablaOutputTrajectoryStock_[i][k],
					nablaInputTrajectoryStock_[i][k]);
		}

		// reset the initial state
		nabla_YmInit = sensitivityOutputTrajectory.back();
	}

	// calculate nabla_q, nabla_Qv, nabla_Rv
	LinearInterpolation<constraint1_state_matrix_t,Eigen::aligned_allocator<constraint1_state_matrix_t> > CmFunc;
	LinearInterpolation<constraint1_control_matrix_t,Eigen::aligned_allocator<constraint1_control_matrix_t> > DmFunc;
	LinearInterpolation<output_vector_t,Eigen::aligned_allocator<output_vector_t> > QvFunc;
	LinearInterpolation<state_matrix_t,Eigen::aligned_allocator<state_matrix_t> > QmFunc;
	LinearInterpolation<control_vector_t,Eigen::aligned_allocator<control_vector_t> > RvFunc;
	LinearInterpolation<control_matrix_t,Eigen::aligned_allocator<control_matrix_t> > RmFunc;
	LinearInterpolation<control_feedback_t,Eigen::aligned_allocator<control_feedback_t> > PmFunc;

	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		CmFunc.setTimeStamp( &(slqp_.nominalTimeTrajectoriesStock_[i]) );
		CmFunc.setData( &(slqp_.CmTrajectoryStock_[i]) );
		DmFunc.setTimeStamp( &(slqp_.nominalTimeTrajectoriesStock_[i]) );
		DmFunc.setData( &(slqp_.DmTrajectoryStock_[i]) );
		QvFunc.setTimeStamp( &(slqp_.nominalTimeTrajectoriesStock_[i]) );
		QvFunc.setData( &(slqp_.QvTrajectoryStock_[i]) );
		QmFunc.setTimeStamp( &(slqp_.nominalTimeTrajectoriesStock_[i]) );
		QmFunc.setData( &(slqp_.QmTrajectoryStock_[i]) );
		RvFunc.setTimeStamp( &(slqp_.nominalTimeTrajectoriesStock_[i]) );
		RvFunc.setData( &(slqp_.RvTrajectoryStock_[i]) );
		RmFunc.setTimeStamp( &(slqp_.nominalTimeTrajectoriesStock_[i]) );
		RmFunc.setData( &(slqp_.RmTrajectoryStock_[i]) );
		PmFunc.setTimeStamp( &(slqp_.nominalTimeTrajectoriesStock_[i]) );
		PmFunc.setData( &(slqp_.PmTrajectoryStock_[i]) );

		int N = sensitivityTimeTrajectoryStock_[i].size();
		nablaqTrajectoryStock_[i].resize(N);
		nablaQvTrajectoryStock_[i].resize(N);
		nablaRvTrajectoryStock_[i].resize(N);

		for (int k=0; k<N; k++) {

			constraint1_state_matrix_t Cm;
			CmFunc.interpolate(sensitivityTimeTrajectoryStock_[i][k], Cm);
			constraint1_control_matrix_t Dm;
			DmFunc.interpolate(sensitivityTimeTrajectoryStock_[i][k], Dm);
			output_vector_t Qv;
			QvFunc.interpolate(sensitivityTimeTrajectoryStock_[i][k], Qv);
			state_matrix_t Qm;
			QmFunc.interpolate(sensitivityTimeTrajectoryStock_[i][k], Qm);
			control_vector_t Rv;
			RvFunc.interpolate(sensitivityTimeTrajectoryStock_[i][k], Rv);
			control_matrix_t Rm;
			RmFunc.interpolate(sensitivityTimeTrajectoryStock_[i][k], Rm);
			control_feedback_t Pm;
			PmFunc.interpolate(sensitivityTimeTrajectoryStock_[i][k], Pm);

			nablaqTrajectoryStock_[i][k]  = Qv.transpose()*nablaOutputTrajectoryStock_[i][k] + Rv.transpose()*nablaInputTrajectoryStock_[i][k];
			nablaQvTrajectoryStock_[i][k] = Qm*nablaOutputTrajectoryStock_[i][k] + Pm.transpose()*nablaInputTrajectoryStock_[i][k];
			nablaRvTrajectoryStock_[i][k] = Pm*nablaOutputTrajectoryStock_[i][k] + Rm*nablaInputTrajectoryStock_[i][k];
		}

		if (i==NUM_SUBSYSTEMS-1)  {
			output_vector_t Qv = slqp_.QvTrajectoryStock_[NUM_SUBSYSTEMS-1].back();
			state_matrix_t Qm  = slqp_.QmTrajectoryStock_[NUM_SUBSYSTEMS-1].back();

			nablaqFinal_  = Qv.transpose()*nablaOutputTrajectoryStock_[NUM_SUBSYSTEMS-1].back();
			nablaQvFinal_ = Qm*nablaOutputTrajectoryStock_[NUM_SUBSYSTEMS-1].back();
		}
	}  // end of i loop

}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * run the SLQ algorithm for a given state and switching times
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::run(const state_vector_t& initState, const std::vector<scalar_t>& switchingTimes)  {

	switchingTimes_ = switchingTimes;
	initState_ = initState;

	// run the SLQ algorithm
	slqp_.run(initState, switchingTimes);

}

} // namespace ocs2

