/*
 * Implementation of GSLQP.h
 *
 *  Created on: Jan 5, 2016
 *      Author: farbod
 */

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
void GSLQP<STATE_DIM, INPUT_DIM, NUM_Subsystems>::rollout(const state_vector_t& initState,
		const std::vector<controller_t>& controllersStock,
		std::vector<scalar_array_t>& timeTrajectoriesStock,
		std::vector<state_vector_array_t>& stateTrajectoriesStock,
		std::vector<control_vector_array_t>& inputTrajectoriesStock)  {

	if (controllersStock.size() != NUM_Subsystems)
		throw std::runtime_error("controllersStock has less controllers then the number of subsystems");

	timeTrajectoriesStock.resize(NUM_Subsystems);
	stateTrajectoriesStock.resize(NUM_Subsystems);
	inputTrajectoriesStock.resize(NUM_Subsystems);

	state_vector_t x0 = initState;
	for (int i=0; i<NUM_Subsystems; i++) {

		timeTrajectoriesStock[i].clear();
		stateTrajectoriesStock[i].clear();
		inputTrajectoriesStock[i].clear();

		// set controller for subsystem i
		subsystemDynamicsPtrStock[i]->setController(controllersStock[i]);
		// simulate subsystem i
		subsystemSimulatorsStockPtr_[i].integrate(x0, switchingTimes_[i], switchingTimes_[i+1], stateTrajectoriesStock[i], timeTrajectoriesStock[i], 1e-3);

		// compute control trajectory for subsystem i
		inputTrajectoriesStock[i].resize(timeTrajectoriesStock[i].size());
		for (int k=0; k<timeTrajectoriesStock[i].size(); k++)
			subsystemDynamicsPtrStock[i]->computeInput(timeTrajectoriesStock[i][k], stateTrajectoriesStock[i][k], inputTrajectoriesStock[i][k]);

		// reset the initial state
		x0 = stateTrajectoriesStock[i].back();
	}
}



/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
void GSLQP<STATE_DIM, INPUT_DIM, NUM_Subsystems>::rolloutCost(const std::vector<scalar_array_t>& timeTrajectoriesStock,
		const std::vector<state_vector_array_t>& stateTrajectoriesStock,
		const std::vector<control_vector_array_t>& inputTrajectoriesStock,
		scalar_t& totalCost)  {

	totalCost = 0.0;
	for (int i=0; i<NUM_Subsystems; i++) {

		scalar_t currentIntermediateCost;
		scalar_t nextIntermediateCost;
		for (int k=0; k<timeTrajectoriesStock[i].size()-1; k++) {

			if (k==0) {
				subsystemCostFunctionsPtrStock_[i]->setCurrentStateAndControl(timeTrajectoriesStock[i][k], stateTrajectoriesStock[i][k], inputTrajectoriesStock[i][k]);
				subsystemCostFunctionsPtrStock_[i]->evaluate(currentIntermediateCost);
			} else {
				currentIntermediateCost = nextIntermediateCost;
			}

			// feed next state and control to cost function
			subsystemCostFunctionsPtrStock_[i]->setCurrentStateAndControl(timeTrajectoriesStock[i][k+1], stateTrajectoriesStock[i][k+1], inputTrajectoriesStock[i][k+1]);
			// evaluate intermediate cost for next time step
			subsystemCostFunctionsPtrStock_[i]->evaluate(nextIntermediateCost);

			totalCost += 0.5*(currentIntermediateCost+nextIntermediateCost)*(timeTrajectoriesStock[i][k+1]-timeTrajectoriesStock[i][k]);
		}

		// terminal cost
		if (i==NUM_Subsystems-1)  {
			scalar_t finalCost;
			subsystemCostFunctionsPtrStock_[i]->setCurrentStateAndControl(timeTrajectoriesStock[i].back(), stateTrajectoriesStock[i].back(), inputTrajectoriesStock[i].back());
			subsystemCostFunctionsPtrStock_[i]->terminalCost(finalCost);
			totalCost += finalCost;
		}
	}

}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
void GSLQP<STATE_DIM, INPUT_DIM, NUM_Subsystems>::approximateOptimalControlProblem()  {

	for (int i=0; i<NUM_Subsystems; i++) {

		for (int k=0; k<nominalTimeTrajectoriesStock_[i].size(); k++) {

			subsystemDerivativesPtrStock_[i]->setCurrentStateAndControl(nominalTimeTrajectoriesStock_[i][k],
					nominalStateTrajectoriesStock_[i][k], nominalInputTrajectoriesStock_[i][k]);
			subsystemDerivativesPtrStock_[i]->getDerivativeState(AmTrajectoryStock_[i][k]);
			subsystemDerivativesPtrStock_[i]->getDerivativesControl(BmTrajectoryStock_[i][k]);

			subsystemCostFunctionsPtrStock_[i]->setCurrentStateAndControl(nominalTimeTrajectoriesStock_[i][k],
					nominalStateTrajectoriesStock_[i][k], nominalInputTrajectoriesStock_[i][k]);
			subsystemCostFunctionsPtrStock_[i]->evaluate(qTrajectoryStock_[i][k](0));
			subsystemCostFunctionsPtrStock_[i]->stateDerivative(QvTrajectoryStock_[i][k]);
			subsystemCostFunctionsPtrStock_[i]->stateSecondDerivative(QmTrajectoryStock_[i][k]);
			subsystemCostFunctionsPtrStock_[i]->controlDerivative(RvTrajectoryStock_[i][k]);
			subsystemCostFunctionsPtrStock_[i]->controlSecondDerivative(RmTrajectoryStock_[i][k]);
			subsystemCostFunctionsPtrStock_[i]->stateControlDerivative(PmTrajectoryStock_[i][k]);
		}

		if (i==NUM_Subsystems-1)  {
			subsystemCostFunctionsPtrStock_[i]->terminalCost(qFinal_(0));
			subsystemCostFunctionsPtrStock_[i]->terminalCostStateDerivative(QvFinal_);
			subsystemCostFunctionsPtrStock_[i]->terminalCostStateSecondDerivative(QmFinal_);
		}
	}
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
void GSLQP<STATE_DIM, INPUT_DIM, NUM_Subsystems>::calculatecontroller(scalar_t& learningRate, std::vector<controller_t>& controllersStock) {

	LinearInterpolation<state_matrix_t,Eigen::aligned_allocator<state_matrix_t> > AmFunc;
	LinearInterpolation<control_gain_matrix_t,Eigen::aligned_allocator<control_gain_matrix_t> > BmFunc;

	LinearInterpolation<eigen_scalar_t,Eigen::aligned_allocator<eigen_scalar_t> > qFunc;
	LinearInterpolation<state_vector_t,Eigen::aligned_allocator<state_vector_t> > QvFunc;
	LinearInterpolation<state_matrix_t,Eigen::aligned_allocator<state_matrix_t> > QmFunc;
	LinearInterpolation<control_vector_t,Eigen::aligned_allocator<control_vector_t> > RvFunc;
	LinearInterpolation<control_matrix_t,Eigen::aligned_allocator<control_matrix_t> > RmFunc;
	LinearInterpolation<control_feedback_t,Eigen::aligned_allocator<control_feedback_t> > PmFunc;

	LinearInterpolation<state_vector_t,Eigen::aligned_allocator<state_vector_t> > nominalStateFunc;
	LinearInterpolation<control_matrix_t,Eigen::aligned_allocator<control_matrix_t> > nominalInputFunc;

	std::vector<control_vector_array_t> deltaUffStock(NUM_Subsystems);
	for (int i=0; i<NUM_Subsystems; i++) {

		AmFunc.reset();
		AmFunc.setTimeStamp(&nominalTimeTrajectoriesStock_[i]);
		AmFunc.setData(&AmTrajectoryStock_[i]);
		BmFunc.reset();
		BmFunc.setTimeStamp(&nominalTimeTrajectoriesStock_[i]);
		BmFunc.setData(&BmTrajectoryStock_[i]);

		qFunc.reset();
		qFunc.setTimeStamp(&nominalTimeTrajectoriesStock_[i]);
		qFunc.setData(&qTrajectoryStock_[i]);
		QvFunc.reset();
		QvFunc.setTimeStamp(&nominalTimeTrajectoriesStock_[i]);
		QvFunc.setData(&QvTrajectoryStock_[i]);
		QmFunc.reset();
		QmFunc.setTimeStamp(&nominalTimeTrajectoriesStock_[i]);
		QmFunc.setData(&QmTrajectoryStock_[i]);
		RvFunc.reset();
		RvFunc.setTimeStamp(&nominalTimeTrajectoriesStock_[i]);
		RvFunc.setData(&RvTrajectoryStock_[i]);
		RmFunc.reset();
		RmFunc.setTimeStamp(&nominalTimeTrajectoriesStock_[i]);
		RmFunc.setData(&RmTrajectoryStock_[i]);
		PmFunc.reset();
		PmFunc.setTimeStamp(&nominalTimeTrajectoriesStock_[i]);
		PmFunc.setData(&PmTrajectoryStock_[i]);

		nominalStateFunc.reset();
		nominalStateFunc.setTimeStamp(&nominalTimeTrajectoriesStock_[i]);
		nominalStateFunc.setData(&nominalStateTrajectoriesStock_[i]);
		nominalInputFunc.reset();
		nominalInputFunc.setTimeStamp(&nominalTimeTrajectoriesStock_[i]);
		nominalInputFunc.setData(&nominalInputTrajectoriesStock_[i]);


		int N = SsTimeTrajectoryStock_[i].size();
		controllersStock[i].time_ = SsTimeTrajectoryStock_[i];
		controllersStock[i].k_.resize(N);
		controllersStock[i].uff_.resize(N);
		deltaUffStock[i].resize(N);
		for (int k=0; k<N; k++) {

			state_matrix_t Am;
			AmFunc.interpolate(SsTimeTrajectoryStock_[i][k], Am);
			control_gain_matrix_t Bm;
			BmFunc.interpolate(SsTimeTrajectoryStock_[i][k], Bm);

			eigen_scalar_t q;
			qFunc.interpolate(SsTimeTrajectoryStock_[i][k], q);
			state_vector_t Qv;
			QvFunc.interpolate(SsTimeTrajectoryStock_[i][k], Qv);
			state_matrix_t Qm;
			QmFunc.interpolate(SsTimeTrajectoryStock_[i][k], Qm);
			control_vector_t Rv;
			RvFunc.interpolate(SsTimeTrajectoryStock_[i][k], Rv);
			control_matrix_t Rm;
			RmFunc.interpolate(SsTimeTrajectoryStock_[i][k], Rm);
			control_feedback_t Pm;
			PmFunc.interpolate(SsTimeTrajectoryStock_[i][k], Pm);

			state_vector_t nominalState;
			nominalStateFunc.interpolate(SsTimeTrajectoryStock_[i][k], nominalState);
			control_matrix_t nominalInput;
			nominalInputFunc.interpolate(SsTimeTrajectoryStock_[i][k], nominalInput);

			controllersStock[i].k_[k]   = -Rm.inverse() * (Pm + Bm.transpose()*SmTrajectoryStock_[i][k]);
			controllersStock[i].uff_[k] = nominalInput - controllersStock[i].k_[k]*nominalState;
			deltaUffStock[i][k] = -Rm.inverse() * (Rv + Bm.transpose()*SvTrajectoryStock_[i][k]);
		}  // end of k loop
	}  // end of i loop

	// finding the optimal learningRate
	lineSearch(controllersStock, deltaUffStock, learningRate);

	// updating uff by the local increment
	for (int i=0; i<NUM_Subsystems; i++)
		for (int k=0; k<SsTimeTrajectoryStock_[i].size(); k++)
			controllersStock[i].uff_[k] += learningRate*deltaUffStock[i][k];

}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
void GSLQP<STATE_DIM, INPUT_DIM, NUM_Subsystems>::lineSearch(const std::vector<controller_t>& controllersStock,
		const std::vector<control_vector_array_t>& deltaUffStock,
		scalar_t& learningRate, scalar_t& totalCost)  {

	scalar_t previousTotalCost;
	std::vector<controller_t> previousControllersStock(NUM_Subsystems);
	std::vector<scalar_array_t> previousTimeTrajectoriesStock(NUM_Subsystems);
	std::vector<state_vector_array_t> previousStateTrajectoriesStock(NUM_Subsystems);
	std::vector<control_vector_array_t> previousInputTrajectoriesStock(NUM_Subsystems);

	scalar_t currentTotalCost;
	std::vector<controller_t> currentControllersStock(controllersStock);
	std::vector<scalar_array_t> currentTimeTrajectoriesStock(NUM_Subsystems);
	std::vector<state_vector_array_t> currentStateTrajectoriesStock(NUM_Subsystems);
	std::vector<control_vector_array_t> currentInputTrajectoriesStock(NUM_Subsystems);

	// modifying uff by the local increamant
	for (int i=0; i<NUM_Subsystems; i++)
		for (int k=0; k<SsTimeTrajectoryStock_[i].size(); k++)
			currentControllersStock[i].uff_[k] += learningRate*deltaUffStock[i][k];

	// rollout
	rollout(initState_, currentControllersStock,
			currentTimeTrajectoriesStock, currentStateTrajectoriesStock, currentInputTrajectoriesStock);

	// calculate rollout cost
	rolloutCost(currentTimeTrajectoriesStock, currentStateTrajectoriesStock, currentInputTrajectoriesStock,
			currentTotalCost);

	// cashing data
	previousTotalCost = currentTotalCost;
	previousControllersStock = currentControllersStock;
	previousTimeTrajectoriesStock  = currentTimeTrajectoriesStock;
	previousStateTrajectoriesStock = currentStateTrajectoriesStock;
	previousInputTrajectoriesStock = currentInputTrajectoriesStock;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
void GSLQP<STATE_DIM, INPUT_DIM, NUM_Subsystems>::SolveSequentialRiccatiEquations()  {

	// final value for the last Riccati equations
	Eigen::Matrix<double,RiccatiEquations::S_DIM_,1> allSsFinal;
	RiccatiEquations::convert2Vector(QmFinal_, QvFinal_, qFinal_, allSsFinal);

	for (int i=NUM_Subsystems-1; i>=0; i--) {

		// set data for Riccati equations
		auto riccatiEquationsPtr = std::make_shared<RiccatiEquations>();
		riccatiEquationsPtr->setData(switchingTimes_[i], switchingTimes_[i+1],
				&AmTrajectoryStock_[i], &BmTrajectoryStock_[i],
				&qTrajectoryStock_[i], &QvTrajectoryStock_[i], &QmTrajectoryStock_[i],
				&RvTrajectoryStock_[i][i], &RmTrajectoryStock_[i][i], &PmTrajectoryStock_[i][i]);

		// integrating the Riccati equations
		ODE45<RiccatiEquations::S_DIM_> ode45(riccatiEquationsPtr);
		std::vector<double> normalizedTimeTrajectory;
		std::vector<Eigen::Matrix<double,RiccatiEquations::S_DIM_,1>, Eigen::aligned_allocator<Eigen::Matrix<double,RiccatiEquations::S_DIM_,1>> > allSsTrajectory;
		ode45.integrate(allSsFinal, i, i+1, allSsTrajectory, normalizedTimeTrajectory);

		// denormalizing time and constructing 'Sm', 'Sv', and 's'
		int N = normalizedTimeTrajectory.size();
		SsTimeTrajectoryStock_[i].resize(N);
		SmTrajectoryStock_[i].resize(N);
		SvTrajectoryStock_[i].resize(N);
		sTrajectoryStock_[i].resize(N);
		for (int k=0; k<normalizedTimeTrajectory.size(); k++) {

			RiccatiEquations::convert2Matrix(allSsTrajectory[N-1-k], SmTrajectoryStock_[i][k], SvTrajectoryStock_[i][k], sTrajectoryStock_[i][k]);
			SsTimeTrajectoryStock_[i][k] = (switchingTimes_[i]-switchingTimes_[i+1])*(normalizedTimeTrajectory[N-1-k]-i) + switchingTimes_[i+1];
		}

		// reset the final value for next Riccati equation
		allSsFinal = allSsTrajectory.back();
	}

}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
void GSLQP<STATE_DIM, INPUT_DIM, NUM_Subsystems>::run(const state_vector_t& initState, const std::vector<scalar_t>& switchingTimes)  {

	if (switchingTimes.size() != NUM_Subsystems+1)
		throw std::runtime_error("Number of switching times should be one plus the number of subsystems.");
	switchingTimes_ = switchingTimes;
	initState_ = initState;

	if (nominalRolloutIsUpdated_==false)  {
		rollout(initState, nominalControllersStock_,
				nominalTimeTrajectoriesStock_, nominalStateTrajectoriesStock_, nominalInputTrajectoriesStock_);
	}

	// linearizing the dynamics and quadratizing the cost funtion along nominal trajectories
	approximateOptimalControlProblem();

	// solve Riccati equations
	SolveSequentialRiccatiEquations();

	// calculate controller
	calculatecontroller(1.0, nominalControllersStock_);

//	// transforme the local value funtion to the global representation
//	transformeLocalValueFuntion2Global();
}

