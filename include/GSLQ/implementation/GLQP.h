/*
 * GLQP.h
 *
 *  Created on: Jan 5, 2016
 *      Author: farbod
 */

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM>
void GLQP<STATE_DIM, INPUT_DIM>::rollout(const state_vector_t& initState,
		const std::vector<controller_t>& controllersStock,
		std::vector<scalar_array_t>& timeTrajectoriesStock,
		std::vector<state_vector_array_t>& stateTrajectoriesStock,
		std::vector<control_vector_array_t>& controlTrajectoriesStock)  {

	state_vector_t x0 = initState;
	for (int i=0; i<numSubsystems_; i++) {

		timeTrajectoriesStock[i].clear();
		stateTrajectoriesStock[i].clear();
		controlTrajectoriesStock[i].clear();

		// set controller for subsystem i
		subsystemDynamicsPtrStock[i]->setController(controllersStock[i]);
		// simulate subsystem i
		subsystemSimulatorsStockPtr_[i]->integrate(x0, switchingTimes_[i], switchingTimes_[i+1], stateTrajectoriesStock[i], timeTrajectoriesStock[i]);

		// compute control trajectory for subsystem i
		controlTrajectoriesStock[i].resize(timeTrajectoriesStock[i].size());
		for (int k=0; k<timeTrajectoriesStock[i].size(); k++)
			subsystemDynamicsPtrStock[i]->computeInput(timeTrajectoriesStock[i][k], stateTrajectoriesStock[i][k], controlTrajectoriesStock[i][k]);

		// reset the initial state
		x0 = stateTrajectoriesStock[i].back();
	}
}



/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM>
void GLQP<STATE_DIM, INPUT_DIM>::rolloutCost(const std::vector<scalar_array_t>& timeTrajectoriesStock,
		const std::vector<state_vector_array_t>& stateTrajectoriesStock,
		const std::vector<control_vector_array_t>& controlTrajectoriesStock,
		scalar_t& totalCost)  {

	totalCost = 0.0;
	for (int i=0; i<numSubsystems_; i++) {

		scalar_t currentIntermediateCost;
		scalar_t nextIntermediateCost;
		for (int k=0; k<timeTrajectoriesStock[i].size()-1; k++) {

			if (k==0) {
				subsystemCostFunctionsPtrStock_[i]->setCurrentStateAndControl(timeTrajectoriesStock[i][k], stateTrajectoriesStock[i][k], controlTrajectoriesStock[i][k]);
				subsystemCostFunctionsPtrStock_[i]->evaluate(currentIntermediateCost);
			} else {
				currentIntermediateCost = nextIntermediateCost;
			}

			// feed next state and control to cost function
			subsystemCostFunctionsPtrStock_[i]->setCurrentStateAndControl(timeTrajectoriesStock[i][k+1], stateTrajectoriesStock[i][k+1], controlTrajectoriesStock[i][k+1]);
			// evaluate intermediate cost for next time step
			subsystemCostFunctionsPtrStock_[i]->evaluate(nextIntermediateCost);

			totalCost += 0.5*(currentIntermediateCost+nextIntermediateCost)*(timeTrajectoriesStock[i][k+1]-timeTrajectoriesStock[i][k]);
		}

		// terminal cost
		if (i==numSubsystems_-1)  {
			scalar_t finalCost;
			subsystemCostFunctionsPtrStock_[i]->terminalCost(finalCost);
			totalCost += finalCost;
		}
	}

}




/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM>
void GLQP<STATE_DIM, INPUT_DIM>::approximateOptimalControlProblem()  {

	for (int i=0; i<numSubsystems_; i++) {

		subsystemDerivativesPtrStock_[i]->setCurrentStateAndControl(0, stateOperatingPointsStock_[i], inputOperatingPointsStock_[i]);
		subsystemDerivativesPtrStock_[i]->getDerivativeState(AmStock_[i]);
		subsystemDerivativesPtrStock_[i]->getDerivativesControl(BmStock_[i]);

		subsystemCostFunctionsPtrStock_[i]->setCurrentStateAndControl(0, stateOperatingPointsStock_[i], inputOperatingPointsStock_[i]);
		subsystemCostFunctionsPtrStock_[i]->evaluate(qStock_[i](0));
		subsystemCostFunctionsPtrStock_[i]->stateDerivative(QvStock_[i]);
		subsystemCostFunctionsPtrStock_[i]->stateSecondDerivative(QmStock_[i]);
		subsystemCostFunctionsPtrStock_[i]->controlDerivative(RvStock_[i]);
		subsystemCostFunctionsPtrStock_[i]->controlSecondDerivative(RmStock_[i]);
		subsystemCostFunctionsPtrStock_[i]->stateControlDerivative(PmStock_[i]);


		if (i==numSubsystems_-1)  {
			subsystemCostFunctionsPtrStock_[i]->terminalCost(qFinal_(0));
			subsystemCostFunctionsPtrStock_[i]->terminalCostStateDerivative(QvFinal_);
			subsystemCostFunctionsPtrStock_[i]->terminalCostStateSecondDerivative(QmFinal_);
		}
	}
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM>
void GLQP<STATE_DIM, INPUT_DIM>::calculatecontroller(const scalar_t& learningRate, std::vector<controller_t>& controllersStock) {

	for (int i=0; i<numSubsystems_; i++) {

		controllersStock[i].time_ = timeTrajectoryStock_[i];

		controllersStock[i].k_.resize(timeTrajectoryStock_[i].size());
		controllersStock[i].uff_.resize(timeTrajectoryStock_[i].size());
		for (int k=0; k<timeTrajectoryStock_[i].size(); k++) {

			controllersStock[i].k_[k]    = -RmStock_[i].inverse() * (PmStock_[i] + BmStock_[i].transpose()*SmTrajectoryStock_[i][k]);
			controllersStock[i].uff_[k]  = -learningRate * RmStock_[i].inverse() * (RvStock_[i]  + BmStock_[i].transpose()*SvTrajectoryStock_[i][k])
								+ inputOperatingPointsStock_[i] - controllersStock[i].k_[k]*stateOperatingPointsStock_[i];
		}
	}
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM>
void GLQP<STATE_DIM, INPUT_DIM>::transformeLocalValueFuntion2Global() {

	for (int i=0; i<numSubsystems_; i++)
		for (int k=0; k<timeTrajectoryStock_[i].size(); k++) {

			sTrajectoryStock_[i][k] = sTrajectoryStock_[i][k] - stateOperatingPointsStock_[i].transpose()*SvTrajectoryStock_[i][k] +
					0.5*stateOperatingPointsStock_[i].transpose()*SmTrajectoryStock_[i][k]*stateOperatingPointsStock_[i];
			SvTrajectoryStock_[i][k] = SvTrajectoryStock_[i][k] - SmTrajectoryStock_[i][k]*stateOperatingPointsStock_[i];
		}
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM>
void GLQP<STATE_DIM, INPUT_DIM>::SolveRiccatiEquation(const std::vector<scalar_t>& switchingTimes)  {

	if (switchingTimes.size() != numSubsystems_+1)
		throw std::runtime_error("Number of switching times should be one plus the number of subsystems.");

	switchingTimes_ = switchingTimes;

	approximateOptimalControlProblem();

	Eigen::Matrix<double,STATE_DIM*STATE_DIM+STATE_DIM+1,1> allSsFinal;
	RiccatiEquations<STATE_DIM, INPUT_DIM>::convert2Vector(QmFinal_, QvFinal_, qFinal_, allSsFinal);
	for (int i=0; i<numSubsystems_; i++) {

		auto riccatiEquationsPtr = std::make_shared<RiccatiEquations<STATE_DIM, INPUT_DIM> >();
		riccatiEquationsPtr->setData(switchingTimes[i], switchingTimes[i+1],
				AmStock_[i], BmStock_[i],
				qStock_[i], QvStock_[i], QmStock_[i], RvStock_[i], RmStock_[i], PmStock_[i]);

		ODE45<STATE_DIM*STATE_DIM+STATE_DIM+1> ode45(riccatiEquationsPtr);

		std::vector<double> normalizedTimeTrajectory;
		std::vector<Eigen::Matrix<double,STATE_DIM*STATE_DIM+STATE_DIM+1,1>, Eigen::aligned_allocator<Eigen::Matrix<double,STATE_DIM*STATE_DIM+STATE_DIM+1,1>> > allSsTrajectory;

		bool flag = ode45.integrate(allSsFinal, i, i+1, allSsTrajectory, normalizedTimeTrajectory);

		int N = normalizedTimeTrajectory.size();
		timeTrajectoryStock_[i].resize(N);
		SmTrajectoryStock_[i].resize(N);
		SvTrajectoryStock_[i].resize(N);
		sTrajectoryStock_[i].resize(N);
		for (int k=0; k<normalizedTimeTrajectory.size(); k++) {

			RiccatiEquations<STATE_DIM, INPUT_DIM>::convert2Matrix(allSsTrajectory[N-1-k], SmTrajectoryStock_[i][k], SvTrajectoryStock_[i][k], sTrajectoryStock_[i][k]);
			timeTrajectoryStock_[i][k] = (switchingTimes[i]-switchingTimes[i+1])* + switchingTimes[i+1];
		}

		// reset the final value for next Riccati equation
		allSsFinal = allSsTrajectory.back();

	}

	// calculate controller
	calculatecontroller(1.0, controllersStock_);

	// transforme the local value funtion to the global representation
	transformeLocalValueFuntion2Global();

}

