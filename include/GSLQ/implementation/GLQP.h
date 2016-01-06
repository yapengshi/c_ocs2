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

	timeTrajectoriesStock.resize(numSubsystems_);
	stateTrajectoriesStock.resize(numSubsystems_);
	controlTrajectoriesStock.resize(numSubsystems_);

	state_vector_t x0 = initState;
	for (int i=0; i<numSubsystems_; i++) {

		timeTrajectoriesStock[i].clear();
		stateTrajectoriesStock[i].clear();
		controlTrajectoriesStock[i].clear();

		// set controller for subsystem i
		subsystemDynamicsPtrStock[i]->setController(controllersStock[i]);
		// simulate subsystem i
		subsystemSimulatorsStockPtr_[i]->integrate(x0, switchingTimes_[i], switchingTimes_[i+1], stateTrajectoriesStock[i], timeTrajectoriesStock[i], 1e-3);

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
			subsystemCostFunctionsPtrStock_[i]->setCurrentStateAndControl(timeTrajectoriesStock[i].back(), stateTrajectoriesStock[i].back(), controlTrajectoriesStock[i].back());
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

//		std::cout << "subsystem " << i << " A and B:\n";
//		std::cout << AmStock_[i] << std::endl << BmStock_[i] << std::endl;
//
//		std::cout << "subsystem " << i << " Qm, Qv, q, Rm, Rv, Pm:\n";
//		std::cout << QmStock_[i] << std::endl << QvStock_[i] << std::endl << qStock_[i] <<
//				std::endl << RmStock_[i] << std::endl << RvStock_[i] << std::endl << PmStock_[i] << std::endl;

		if (i==numSubsystems_-1)  {
			subsystemCostFunctionsPtrStock_[i]->terminalCost(qFinal_(0));
			subsystemCostFunctionsPtrStock_[i]->terminalCostStateDerivative(QvFinal_);
			subsystemCostFunctionsPtrStock_[i]->terminalCostStateSecondDerivative(QmFinal_);

//			std::cout << "subsystem " << i << " Qm_f, Qv_f, q_f:\n";
//			std::cout << QmFinal_ << std::endl << QvFinal_ << std::endl << qFinal_ << std::endl;
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
void GLQP<STATE_DIM, INPUT_DIM>::getValueFuntion(const scalar_t& time, const state_vector_t& state, scalar_t& valueFuntion)  {

	int activeSubsystem = -1;
	for (int i=0; i<numSubsystems_; i++)  {
		activeSubsystem = i;
		if (switchingTimes_[i]<=time && time<switchingTimes_[i+1])
			break;
	}

	state_matrix_t Sm;
	LinearInterpolation<state_matrix_t,Eigen::aligned_allocator<state_matrix_t> > SmFunc(&timeTrajectoryStock_[activeSubsystem], &SmTrajectoryStock_[activeSubsystem]);
	SmFunc.interpolate(time, Sm);
	state_vector_t Sv;
	LinearInterpolation<state_vector_t,Eigen::aligned_allocator<state_vector_t> > SvFunc(&timeTrajectoryStock_[activeSubsystem], &SvTrajectoryStock_[activeSubsystem]);
	SvFunc.interpolate(time, Sv);
	eigen_scalar_t s;
	LinearInterpolation<eigen_scalar_t,Eigen::aligned_allocator<eigen_scalar_t> > sFunc(&timeTrajectoryStock_[activeSubsystem], &sTrajectoryStock_[activeSubsystem]);
	sFunc.interpolate(time, s);

	valueFuntion = (s + state.transpose()*Sv + 0.5*state.transpose()*Sm*state).eval()(0);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM>
void GLQP<STATE_DIM, INPUT_DIM>::SolveRiccatiEquation(const std::vector<scalar_t>& switchingTimes)  {

	if (switchingTimes.size() != numSubsystems_+1)
		throw std::runtime_error("Number of switching times should be one plus the number of subsystems.");

	switchingTimes_ = switchingTimes;
//	std::cout << "switchingTimes: " << switchingTimes_[0] << ", " << switchingTimes_[1] << ", " << switchingTimes_[2] << std::endl;

	approximateOptimalControlProblem();

	Eigen::Matrix<double,STATE_DIM*STATE_DIM+STATE_DIM+1,1> allSsFinal;
	RiccatiEquations<STATE_DIM, INPUT_DIM>::convert2Vector(QmFinal_, QvFinal_, qFinal_, allSsFinal);
//	std::cout << "allSsFinal: " << allSsFinal.transpose() << std::endl;

	for (int i=numSubsystems_-1; i>=0; i--) {

		auto riccatiEquationsPtr = std::make_shared<RiccatiEquations<STATE_DIM, INPUT_DIM> >();
		riccatiEquationsPtr->setData(switchingTimes[i], switchingTimes[i+1],
				AmStock_[i], BmStock_[i],
				qStock_[i], QvStock_[i], QmStock_[i], RvStock_[i], RmStock_[i], PmStock_[i]);

		ODE45<STATE_DIM*STATE_DIM+STATE_DIM+1> ode45(riccatiEquationsPtr);

		std::vector<double> normalizedTimeTrajectory;
		std::vector<Eigen::Matrix<double,STATE_DIM*STATE_DIM+STATE_DIM+1,1>, Eigen::aligned_allocator<Eigen::Matrix<double,STATE_DIM*STATE_DIM+STATE_DIM+1,1>> > allSsTrajectory;

		ode45.integrate(allSsFinal, i, i+1, allSsTrajectory, normalizedTimeTrajectory);

		int N = normalizedTimeTrajectory.size();
		timeTrajectoryStock_[i].resize(N);
		SmTrajectoryStock_[i].resize(N);
		SvTrajectoryStock_[i].resize(N);
		sTrajectoryStock_[i].resize(N);
		for (int k=0; k<normalizedTimeTrajectory.size(); k++) {

			RiccatiEquations<STATE_DIM, INPUT_DIM>::convert2Matrix(allSsTrajectory[N-1-k], SmTrajectoryStock_[i][k], SvTrajectoryStock_[i][k], sTrajectoryStock_[i][k]);
			timeTrajectoryStock_[i][k] = (switchingTimes[i]-switchingTimes[i+1])*(normalizedTimeTrajectory[N-1-k]-i) + switchingTimes[i+1];
//			std::cout << "Time " << timeTrajectoryStock_[i][k] << ",\t AllSs[" << i << "][" << k<< "]: " << allSsTrajectory[N-1-k].transpose() << std::endl;
		}
//		std::cout << "Sm[" << i << "][0]: \n" << SmTrajectoryStock_[i].front() << std::endl;
//		std::cout << "SV[" << i << "][0]: " << SvTrajectoryStock_[i].front().transpose() << std::endl;
//		std::cout << "s[" << i << "][0]: " << sTrajectoryStock_[i].front() << std::endl;

		// reset the final value for next Riccati equation
		allSsFinal = allSsTrajectory.back();

	}

	// calculate controller
	calculatecontroller(1.0, controllersStock_);

	// transforme the local value funtion to the global representation
	transformeLocalValueFuntion2Global();

}

