/*
 * GLQP.h
 *
 *  Created on: Jan 5, 2016
 *      Author: farbod
 */


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
void GLQP<STATE_DIM, INPUT_DIM, NUM_Subsystems>::rollout(const state_vector_t& initState,
		const std::vector<controller_t>& controllersStock,
		std::vector<scalar_array_t>& timeTrajectoriesStock,
		std::vector<state_vector_array_t>& stateTrajectoriesStock,
		std::vector<control_vector_array_t>& controlTrajectoriesStock)  {

	if (controllersStock.size() != NUM_Subsystems)
		throw std::runtime_error("controllersStock has less controllers then the number of subsystems");

	timeTrajectoriesStock.resize(NUM_Subsystems);
	stateTrajectoriesStock.resize(NUM_Subsystems);
	controlTrajectoriesStock.resize(NUM_Subsystems);

	state_vector_t x0 = initState;
	for (int i=0; i<NUM_Subsystems; i++) {

		timeTrajectoriesStock[i].clear();
		stateTrajectoriesStock[i].clear();
		controlTrajectoriesStock[i].clear();

		// initialize subsystem i
		subsystemDynamicsPtrStock[i]->initializeModel(switchingTimes_[i], x0, switchingTimes_[i+1], "GLQP");
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
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
void GLQP<STATE_DIM, INPUT_DIM, NUM_Subsystems>::rolloutCost(const std::vector<scalar_array_t>& timeTrajectoriesStock,
		const std::vector<state_vector_array_t>& stateTrajectoriesStock,
		const std::vector<control_vector_array_t>& controlTrajectoriesStock,
		scalar_t& totalCost)  {

	totalCost = 0.0;
	for (int i=0; i<NUM_Subsystems; i++) {

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
		if (i==NUM_Subsystems-1)  {
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
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
void GLQP<STATE_DIM, INPUT_DIM, NUM_Subsystems>::approximateOptimalControlProblem()  {

	for (int i=0; i<NUM_Subsystems; i++) {

		subsystemDerivativesPtrStock_[i]->initializeModel(switchingTimes_[i], stateOperatingPointsStock_.at(i), switchingTimes_[i+1], "GLQP");
		subsystemDerivativesPtrStock_[i]->setCurrentStateAndControl(0, stateOperatingPointsStock_.at(i), inputOperatingPointsStock_.at(i));
		subsystemDerivativesPtrStock_[i]->getDerivativeState(AmStock_.at(i));
		subsystemDerivativesPtrStock_[i]->getDerivativesControl(BmStock_.at(i));

		subsystemCostFunctionsPtrStock_[i]->setCurrentStateAndControl(0, stateOperatingPointsStock_.at(i), inputOperatingPointsStock_.at(i));
		subsystemCostFunctionsPtrStock_[i]->evaluate(qStock_.at(i)(0));
		subsystemCostFunctionsPtrStock_[i]->stateDerivative(QvStock_.at(i));
		subsystemCostFunctionsPtrStock_[i]->stateSecondDerivative(QmStock_.at(i));
		subsystemCostFunctionsPtrStock_[i]->controlDerivative(RvStock_.at(i));
		subsystemCostFunctionsPtrStock_[i]->controlSecondDerivative(RmStock_.at(i));
		subsystemCostFunctionsPtrStock_[i]->stateControlDerivative(PmStock_.at(i));

		if (INFO_ON_) {
			std::cout<< "stateOperatingPoint[" << i << "]: \n" << stateOperatingPointsStock_[i].transpose() << std::endl;
			std::cout<< "inputOperatingPoint[" << i << "]: \n" << inputOperatingPointsStock_[i].transpose() << std::endl;
			std::cout<< "A[" << i << "]: \n" << AmStock_[i] << std::endl;
			std::cout<< "B[" << i << "]: \n" << BmStock_[i] << std::endl;
			std::cout<< "q[" << i << "]: \t" << qStock_[i] << std::endl;
			std::cout<< "Qv[" << i << "]: \n" << QvStock_[i].transpose() << std::endl;
			std::cout<< "Qm[" << i << "]: \n" << QmStock_[i] << std::endl;
			std::cout<< "Rv[" << i << "]: \n" << RvStock_[i].transpose() << std::endl;
			std::cout<< "Rm[" << i << "]: \n" << RmStock_[i] << std::endl;
			std::cout<< "Pm[" << i << "]: \n" << PmStock_[i] << std::endl;
		}

		// making sure that Qm is PSD
		makePSD(QmStock_[i]);

		if (i==NUM_Subsystems-1)  {
			subsystemCostFunctionsPtrStock_[i]->terminalCost(qFinal_(0));
			subsystemCostFunctionsPtrStock_[i]->terminalCostStateDerivative(QvFinal_);
			subsystemCostFunctionsPtrStock_[i]->terminalCostStateSecondDerivative(QmFinal_);
			// making sure that Qm is PSD
			makePSD(QmFinal_);

			if (INFO_ON_) {
				std::cout<< "qFinal[" << i << "]: \t" << qFinal_ << std::endl;
				std::cout<< "QvFinal[" << i << "]: \n" << QvFinal_.transpose() << std::endl;
				std::cout<< "QmFinal[" << i << "]: \n" << QmFinal_ << std::endl;
			}
		}
	}
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
void GLQP<STATE_DIM, INPUT_DIM, NUM_Subsystems>::calculatecontroller(const scalar_t& learningRate, std::vector<controller_t>& controllersStock) {

	for (int i=0; i<NUM_Subsystems; i++) {

		controllersStock[i].time_ = timeTrajectoryStock_[i];

		controllersStock[i].k_.resize(timeTrajectoryStock_[i].size());
		controllersStock[i].uff_.resize(timeTrajectoryStock_[i].size());
		for (int k=0; k<timeTrajectoryStock_[i].size(); k++) {

			control_matrix_t RmInverse = RmStock_[i].inverse();
			controllersStock[i].k_[k]    = -RmInverse * (PmStock_[i] + BmStock_[i].transpose()*SmTrajectoryStock_[i][k]);
			controllersStock[i].uff_[k]  = -learningRate * RmInverse * (RvStock_[i]  + BmStock_[i].transpose()*SvTrajectoryStock_[i][k])
								+ inputOperatingPointsStock_[i] - controllersStock[i].k_[k]*stateOperatingPointsStock_[i];
		}
	}
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
void GLQP<STATE_DIM, INPUT_DIM, NUM_Subsystems>::transformeLocalValueFuntion2Global() {

	for (int i=0; i<NUM_Subsystems; i++)
		for (int k=0; k<timeTrajectoryStock_[i].size(); k++) {

			sTrajectoryStock_[i][k] = sTrajectoryStock_[i][k] - stateOperatingPointsStock_[i].transpose()*SvTrajectoryStock_[i][k] +
					0.5*stateOperatingPointsStock_[i].transpose()*SmTrajectoryStock_[i][k]*stateOperatingPointsStock_[i];
			SvTrajectoryStock_[i][k] = SvTrajectoryStock_[i][k] - SmTrajectoryStock_[i][k]*stateOperatingPointsStock_[i];
		}
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
template <typename Derived>
bool GLQP<STATE_DIM, INPUT_DIM, NUM_Subsystems>::makePSD(Eigen::MatrixBase<Derived>& squareMatrix) {

	if (squareMatrix.rows() != squareMatrix.cols())  throw std::runtime_error("Not a square matrix: makePSD() method is for square matrix.");

	Eigen::SelfAdjointEigenSolver<Derived> eig(squareMatrix);
	Eigen::VectorXd lambda = eig.eigenvalues();

	bool hasNegativeEigenValue = false;
	for (size_t j=0; j<lambda.size() ; j++)
		if (lambda(j) < 0.0) {
			hasNegativeEigenValue = true;
			lambda(j) = 0.0;
		}

	if (hasNegativeEigenValue)
		squareMatrix = eig.eigenvectors() * lambda.asDiagonal() * eig.eigenvectors().inverse();
//	else
//		squareMatrix = 0.5*(squareMatrix+squareMatrix.transpose()).eval();

	return hasNegativeEigenValue;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
void GLQP<STATE_DIM, INPUT_DIM, NUM_Subsystems>::getValueFuntion(const scalar_t& time, const state_vector_t& state, scalar_t& valueFuntion)  {

	int activeSubsystem = -1;
	for (int i=0; i<NUM_Subsystems; i++)  {
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
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
void GLQP<STATE_DIM, INPUT_DIM, NUM_Subsystems>::SolveRiccatiEquations()  {

	// final value for the last Riccati equations
	Eigen::Matrix<double,RiccatiEquations::S_DIM_,1> allSsFinal;
	RiccatiEquations::convert2Vector(QmFinal_, QvFinal_, qFinal_, allSsFinal);

	for (int i=NUM_Subsystems-1; i>=0; i--) {

		// set data for Riccati equations
		auto riccatiEquationsPtr = std::make_shared<RiccatiEquations>();
		riccatiEquationsPtr->setData(switchingTimes_[i], switchingTimes_[i+1],
				AmStock_[i], BmStock_[i],
				qStock_[i], QvStock_[i], QmStock_[i], RvStock_[i], RmStock_[i], PmStock_[i]);

		// integrating the Riccati equations
		ODE45<RiccatiEquations::S_DIM_> ode45(riccatiEquationsPtr);
		std::vector<double> normalizedTimeTrajectory;
		std::vector<Eigen::Matrix<double,RiccatiEquations::S_DIM_,1>, Eigen::aligned_allocator<Eigen::Matrix<double,RiccatiEquations::S_DIM_,1>> > allSsTrajectory;
		ode45.integrate(allSsFinal, i, i+1, allSsTrajectory, normalizedTimeTrajectory);

		// denormalizing time and constructing 'Sm', 'Sv', and 's'
		int N = normalizedTimeTrajectory.size();
		timeTrajectoryStock_[i].resize(N);
		SmTrajectoryStock_[i].resize(N);
		SvTrajectoryStock_[i].resize(N);
		sTrajectoryStock_[i].resize(N);
		for (int k=0; k<normalizedTimeTrajectory.size(); k++) {

			RiccatiEquations::convert2Matrix(allSsTrajectory[N-1-k], SmTrajectoryStock_[i][k], SvTrajectoryStock_[i][k], sTrajectoryStock_[i][k]);
			timeTrajectoryStock_[i][k] = (switchingTimes_[i]-switchingTimes_[i+1])*(normalizedTimeTrajectory[N-1-k]-i) + switchingTimes_[i+1];
		}

		// reset the final value for next Riccati equation
		allSsFinal = allSsTrajectory.back();
	}

}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
void GLQP<STATE_DIM, INPUT_DIM, NUM_Subsystems>::run(const std::vector<scalar_t>& switchingTimes)  {

//	std::cout << "\n#### GLQP solver starts ..." << std::endl << std::endl;

	if (switchingTimes.size() != NUM_Subsystems+1)
		throw std::runtime_error("Number of switching times should be one plus the number of subsystems.");
	switchingTimes_ = switchingTimes;

	// linearizing the dynamics and quadratizing the cost funtion along nominal trajectories
	approximateOptimalControlProblem();

	// solve Riccati equations
	SolveRiccatiEquations();

	// calculate controller
	calculatecontroller(1.0, controllersStock_);

	// transforme the local value funtion to the global representation
	transformeLocalValueFuntion2Global();

//	std::cout << "\n#### GLQP solver ends." << std::endl;
}


