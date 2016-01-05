/*
 * Implementation of GSLQP.h
 *
 *  Created on: Jan 5, 2016
 *      Author: farbod
 */

template <size_t STATE_DIM, size_t INPUT_DIM>
void GSLQP<STATE_DIM, INPUT_DIM>::rollout(const state_vector_t& initState,
		const std::vector<scalar_t>& switchingTimes,
		const std::vector<controller_t>& controllersStock,
		std::vector<scalar_array_t>& timeTrajectoriesStock,
		std::vector<state_vector_array_t>& stateTrajectoriesStock,
		std::vector<control_vector_array_t>& controlTrajectoriesStock)  {

	if (switchingTimes.size() != numSubsystems_+1)
		throw std::runtime_error("Number of switching times should be one plus the number of subsystems.");

	state_vector_t x0 = initState;
	for (int i=0; i<numSubsystems_; i++) {

		timeTrajectoriesStock[i].clear();
		stateTrajectoriesStock[i].clear();
		controlTrajectoriesStock[i].clear();

		// set controller for subsystem i
		subsystemDynamicsPtrStock[i]->setController(controllersStock[i]);
		// simulate subsystem i
		subsystemSimulatorsStock_[i].integrate(x0, switchingTimes[i], switchingTimes[i+1], stateTrajectoriesStock[i], timeTrajectoriesStock[i]);

		// compute control trajectory for subsystem i
		controlTrajectoriesStock[i].resize(timeTrajectoriesStock[i].size());
		for (int k=0; k<timeTrajectoriesStock[i].size(); k++)
			subsystemDynamicsPtrStock[i]->computeInput(timeTrajectoriesStock[i][k], stateTrajectoriesStock[i][k], controlTrajectoriesStock[i][k]);

		// reset the initial state
		x0 = stateTrajectoriesStock[i].back();
	}
}



template <size_t STATE_DIM, size_t INPUT_DIM>
void GSLQP<STATE_DIM, INPUT_DIM>::rolloutCost(const std::vector<scalar_array_t>& timeTrajectoriesStock,
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




template <size_t STATE_DIM, size_t INPUT_DIM>
void GSLQP<STATE_DIM, INPUT_DIM>::approximateOptimalControlProblem()  {

	for (int i=0; i<numSubsystems_; i++) {

		for (int k=0; k<nominalTimeTrajectoriesStock_[i].size(); k++) {

			subsystemDerivativesPtrStock_[i]->setCurrentStateAndControl(nominalTimeTrajectoriesStock_[i][k],
					nominalStateTrajectoriesStock_[i][k],
					nominalControlTrajectoriesStock_[i][k]);
			subsystemDerivativesPtrStock_[i]->getDerivativeState(AmTrajectoryStock_[i][k]);
			subsystemDerivativesPtrStock_[i]->getDerivativesControl(BmTrajectoryStock_[i][k]);

			subsystemCostFunctionsPtrStock_[i]->setCurrentStateAndControl(nominalTimeTrajectoriesStock_[i][k],
					nominalStateTrajectoriesStock_[i][k],
					nominalControlTrajectoriesStock_[i][k]);
			subsystemCostFunctionsPtrStock_[i]->evaluate(qTrajectoryStock_[i][k]);
			subsystemCostFunctionsPtrStock_[i]->stateDerivative(QvTrajectoryStock_[i][k]);
			subsystemCostFunctionsPtrStock_[i]->stateSecondDerivative(QmTrajectoryStock_[i][k]);
			subsystemCostFunctionsPtrStock_[i]->controlDerivative(RvTrajectoryStock_[i][k]);
			subsystemCostFunctionsPtrStock_[i]->controlSecondDerivative(RmTrajectoryStock_[i][k]);
			subsystemCostFunctionsPtrStock_[i]->stateControlDerivative(PmTrajectoryStock_[i][k]);

		}

		if (i==numSubsystems_-1)  {
			subsystemCostFunctionsPtrStock_[i]->terminalCost(qFinal_);
			subsystemCostFunctionsPtrStock_[i]->terminalCostStateDerivative(QvFinal_);
			subsystemCostFunctionsPtrStock_[i]->terminalCostStateSecondDerivative(QmFinal_);
		}
	}
}

