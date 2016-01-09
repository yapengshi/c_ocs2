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
		subsystemSimulatorsStockPtr_[i]->integrate(x0, switchingTimes_[i], switchingTimes_[i+1], stateTrajectoriesStock[i], timeTrajectoriesStock[i], 1e-3);

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

		int N = nominalTimeTrajectoriesStock_[i].size();
		AmTrajectoryStock_[i].resize(N);
		BmTrajectoryStock_[i].resize(N);

		qTrajectoryStock_[i].resize(N);
		QvTrajectoryStock_[i].resize(N);
		QmTrajectoryStock_[i].resize(N);
		RvTrajectoryStock_[i].resize(N);
		RmTrajectoryStock_[i].resize(N);
		PmTrajectoryStock_[i].resize(N);

		for (int k=0; k<N; k++) {

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
void GSLQP<STATE_DIM, INPUT_DIM, NUM_Subsystems>::calculatecontroller(scalar_t& learningRateStar) {

	std::vector<controller_t> controllersStock(NUM_Subsystems);

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
	std::vector<control_vector_t> maxDeltaUffStock(NUM_Subsystems);

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

		// display
		if (options_.dispay_)  maxDeltaUffStock[i] = *std::max_element(deltaUffStock[i].begin(), deltaUffStock[i].end(),
				[] (const control_vector_t& u1, const control_vector_t& u2){ return u1.norm() < u2.norm(); });

	}  // end of i loop

	// display
	if (options_.dispay_)  {
		control_vector_t maxDeltaUff = *std::max_element(maxDeltaUffStock.begin(), maxDeltaUffStock.end(),
				[] (const control_vector_t& u1, const control_vector_t& u2){ return u1.norm() < u2.norm(); });
		std::cout << "max delta_uff norm: " << maxDeltaUff.norm() << std::endl;
	}

	// finding the optimal learningRate
	lineSearch(controllersStock, deltaUffStock, learningRateStar);

	// calculating the nominal controller
	nominalControllersStock_ = controllersStock;
	if (learningRateStar>0) {
		for (int i=0; i<NUM_Subsystems; i++)
			for (int k=0; k<SsTimeTrajectoryStock_[i].size(); k++)
				nominalControllersStock_[i].uff_[k] += learningRateStar*deltaUffStock[i][k];
	}

}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
void GSLQP<STATE_DIM, INPUT_DIM, NUM_Subsystems>::lineSearch(const std::vector<controller_t>& controllersStock,
		const std::vector<control_vector_array_t>& deltaUffStock,
		scalar_t& learningRateStar)  {

	scalar_t learningRate = learningRateStar;

	scalar_t lsTotalCost;
	std::vector<controller_t> lsControllersStock(NUM_Subsystems);
	std::vector<scalar_array_t> lsTimeTrajectoriesStock(NUM_Subsystems);
	std::vector<state_vector_array_t> lsStateTrajectoriesStock(NUM_Subsystems);
	std::vector<control_vector_array_t> lsInputTrajectoriesStock(NUM_Subsystems);

	while (learningRate > options_.minLearningRate_)  {
		// modifying uff by the local increamant
		lsControllersStock = controllersStock;
		for (int i=0; i<NUM_Subsystems; i++)
			for (int k=0; k<SsTimeTrajectoryStock_[i].size(); k++)
				lsControllersStock[i].uff_[k] += learningRate*deltaUffStock[i][k];

		// rollout
		rollout(initState_, lsControllersStock, lsTimeTrajectoriesStock, lsStateTrajectoriesStock, lsInputTrajectoriesStock);
		// calculate rollout cost
		rolloutCost(lsTimeTrajectoriesStock, lsStateTrajectoriesStock, lsInputTrajectoriesStock, lsTotalCost);

		if (options_.dispay_)  std::cout << "\t learningRate " << learningRate << " cost: " << nominalTotalCost_ << std::endl;

		// break condition 1: it exits with largest learningRate that its cost is smaller than nominal cost.
		if (lsTotalCost < nominalTotalCost_*(1-1e-3*learningRate))  {
			nominalRolloutIsUpdated_ = true;
			nominalTotalCost_ = lsTotalCost;
			nominalControllersStock_ = lsControllersStock;
			nominalTimeTrajectoriesStock_  = lsTimeTrajectoriesStock;
			nominalStateTrajectoriesStock_ = lsStateTrajectoriesStock;
			nominalInputTrajectoriesStock_ = lsInputTrajectoriesStock;
			learningRateStar = learningRate;
			break;
		} else {
			learningRate = 0.5*learningRate;
		}

	}  // end of while

	if (learningRate <= options_.minLearningRate_) {
		nominalRolloutIsUpdated_ = true;  // since the open loop input will not change, the nominal trajectories will be constatnt (no disturbance effect has been assumed)
		learningRateStar = 0.0;
	}

	// display
	if (options_.dispay_)  std::cout << "The chosen learningRate is: " << learningRateStar << std::endl;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
void GSLQP<STATE_DIM, INPUT_DIM, NUM_Subsystems>::transformeLocalValueFuntion2Global() {

	LinearInterpolation<state_vector_t,Eigen::aligned_allocator<state_vector_t> > nominalStateFunc;

	for (int i=0; i<NUM_Subsystems; i++) {

		nominalStateFunc.reset();
		nominalStateFunc.setTimeStamp(&nominalTimeTrajectoriesStock_[i]);
		nominalStateFunc.setData(&nominalStateTrajectoriesStock_[i]);

		for (int k=0; k<SsTimeTrajectoryStock_[i].size(); k++) {

			state_vector_t nominalState;
			nominalStateFunc.interpolate(SsTimeTrajectoryStock_[i][k], nominalState);

			sTrajectoryStock_[i][k] = sTrajectoryStock_[i][k] - nominalState.transpose()*SvTrajectoryStock_[i][k] +
					0.5*nominalState.transpose()*SmTrajectoryStock_[i][k]*nominalState;
			SvTrajectoryStock_[i][k] = SvTrajectoryStock_[i][k] - SmTrajectoryStock_[i][k]*nominalState;
		}  // end of k loop
	}  // enf of i loop
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
void GSLQP<STATE_DIM, INPUT_DIM, NUM_Subsystems>::getValueFuntion(const scalar_t& time, const state_vector_t& state, scalar_t& valueFuntion)  {

	int activeSubsystem = -1;
	for (int i=0; i<NUM_Subsystems; i++)  {
		activeSubsystem = i;
		if (switchingTimes_[i]<=time && time<switchingTimes_[i+1])
			break;
	}

	state_matrix_t Sm;
	LinearInterpolation<state_matrix_t,Eigen::aligned_allocator<state_matrix_t> > SmFunc(&SsTimeTrajectoryStock_[activeSubsystem], &SmTrajectoryStock_[activeSubsystem]);
	SmFunc.interpolate(time, Sm);
	state_vector_t Sv;
	LinearInterpolation<state_vector_t,Eigen::aligned_allocator<state_vector_t> > SvFunc(&SsTimeTrajectoryStock_[activeSubsystem], &SvTrajectoryStock_[activeSubsystem]);
	SvFunc.interpolate(time, Sv);
	eigen_scalar_t s;
	LinearInterpolation<eigen_scalar_t,Eigen::aligned_allocator<eigen_scalar_t> > sFunc(&SsTimeTrajectoryStock_[activeSubsystem], &sTrajectoryStock_[activeSubsystem]);
	sFunc.interpolate(time, s);

	valueFuntion = (s + state.transpose()*Sv + 0.5*state.transpose()*Sm*state).eval()(0);
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
void GSLQP<STATE_DIM, INPUT_DIM, NUM_Subsystems>::SolveSequentialRiccatiEquations(const scalar_t& learningRate)  {

	// final value for the last Riccati equations
	Eigen::Matrix<double,RiccatiEquations::S_DIM_,1> allSsFinal;
	RiccatiEquations::convert2Vector(QmFinal_, QvFinal_, qFinal_, allSsFinal);

	for (int i=NUM_Subsystems-1; i>=0; i--) {

		// set data for Riccati equations
		auto riccatiEquationsPtr = std::make_shared<RiccatiEquations>();
		riccatiEquationsPtr->setData(learningRate,
				i, switchingTimes_[i], switchingTimes_[i+1],
				&nominalTimeTrajectoriesStock_[i],
				&AmTrajectoryStock_[i], &BmTrajectoryStock_[i],
				&qTrajectoryStock_[i], &QvTrajectoryStock_[i], &QmTrajectoryStock_[i],
				&RvTrajectoryStock_[i], &RmTrajectoryStock_[i], &PmTrajectoryStock_[i]);

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

	scalar_t learningRateStar = 1.0;  // resetting learningRateStar
	size_t iteration = 0;
	while (iteration<options_.maxIteration_ && learningRateStar>0)  {

		// do a rollout if nominalRolloutIsUpdated_ is fale.
		if (nominalRolloutIsUpdated_==false)  {
			rollout(initState_, nominalControllersStock_,
					nominalTimeTrajectoriesStock_, nominalStateTrajectoriesStock_, nominalInputTrajectoriesStock_);
			rolloutCost(nominalTimeTrajectoriesStock_, nominalStateTrajectoriesStock_, nominalInputTrajectoriesStock_,
					nominalTotalCost_);
			// display
			if (options_.dispay_ && iteration==0)  std::cout << "\n#### Initial controller: \ncost: " << nominalTotalCost_ << std::endl;
		}

		// display
		if (options_.dispay_)  std::cout << "\n#### Iteration " <<  iteration << std::endl;

		// linearizing the dynamics and quadratizing the cost funtion along nominal trajectories
		approximateOptimalControlProblem();

		// solve Riccati equations
		SolveSequentialRiccatiEquations(1.0 /*nominal learningRate*/);

		// calculate controller
		nominalRolloutIsUpdated_ = false;
		learningRateStar = 1.0;  // resetting learningRateStar
		calculatecontroller(learningRateStar);

		// display
		if (options_.dispay_)  std::cout << "cost: " << nominalTotalCost_ << std::endl;

		iteration++;

	}  // end of j loop

	// display
	if (options_.dispay_)  std::cout << "\n#### Final iteration" << std::endl;

	// linearizing the dynamics and quadratizing the cost funtion along nominal trajectories
	approximateOptimalControlProblem();
	// prevents the changes in the nominal trajectories and just update the gains
	learningRateStar = 0.0;
	// solve Riccati equations
	SolveSequentialRiccatiEquations(learningRateStar);
	// calculate controller
	calculatecontroller(learningRateStar);

	transformeLocalValueFuntion2Global();

}

