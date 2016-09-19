/*
 * LQP.h
 *
 *  Created on: Aug 3, 2016
 *      Author: farbod
 */

namespace ocs2{


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
void LQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::rollout(const state_vector_t& initState,
		const std::vector<controller_t>& controllersStock,
		std::vector<scalar_array_t>& timeTrajectoriesStock,
		std::vector<state_vector_array_t>& stateTrajectoriesStock,
		std::vector<control_vector_array_t>& inputTrajectoriesStock,
		std::vector<output_vector_array_t>& outputTrajectoriesStock)  {

	if (controllersStock.size() != NUM_Subsystems)
		throw std::runtime_error("controllersStock has less controllers then the number of subsystems");

	timeTrajectoriesStock.resize(NUM_Subsystems);
	stateTrajectoriesStock.resize(NUM_Subsystems);
	inputTrajectoriesStock.resize(NUM_Subsystems);
	outputTrajectoriesStock.resize(NUM_Subsystems);

	state_vector_t x0 = initState;
	for (int i=0; i<NUM_Subsystems; i++) {

		timeTrajectoriesStock[i].clear();
		stateTrajectoriesStock[i].clear();

		// initialize subsystem i
		subsystemDynamicsPtrStock_[i]->initializeModel(switchingTimes_, x0, i, "LQP");
		// set controller for subsystem i
		subsystemDynamicsPtrStock_[i]->setController(controllersStock[i]);
		// simulate subsystem i
		subsystemSimulatorsStockPtr_[i]->integrate(x0, switchingTimes_[i], switchingTimes_[i+1], stateTrajectoriesStock[i], timeTrajectoriesStock[i], 1e-3);

		// compute control trajectory for subsystem i
		inputTrajectoriesStock[i].resize(timeTrajectoriesStock[i].size());
		outputTrajectoriesStock[i].resize(timeTrajectoriesStock[i].size());
		for (int k=0; k<timeTrajectoriesStock[i].size(); k++)   {
			subsystemDynamicsPtrStock_[i]->computeOutput(timeTrajectoriesStock[i][k], stateTrajectoriesStock[i][k], outputTrajectoriesStock[i][k]);
			subsystemDynamicsPtrStock_[i]->computeInput(timeTrajectoriesStock[i][k], outputTrajectoriesStock[i][k], inputTrajectoriesStock[i][k]);
		}

		// reset the initial state
		x0 = stateTrajectoriesStock[i].back();

		if (x0 != x0)  throw std::runtime_error("The rollout in LQP is unstable.");
	}  // end of i loop
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
void LQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::rollout(const state_vector_t& initState,
		const std::vector<controller_t>& controllersStock,
		std::vector<scalar_array_t>& timeTrajectoriesStock,
		std::vector<state_vector_array_t>& stateTrajectoriesStock,
		std::vector<control_vector_array_t>& inputTrajectoriesStock)  {

	std::vector<output_vector_array_t> outputTrajectoriesStock;

	rollout(initState, controllersStock,
			timeTrajectoriesStock, stateTrajectoriesStock, inputTrajectoriesStock,
			outputTrajectoriesStock);
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
void LQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::rolloutCost(const std::vector<scalar_array_t>& timeTrajectoriesStock,
		const std::vector<output_vector_array_t>& outputTrajectoriesStock,
		const std::vector<control_vector_array_t>& inputTrajectoriesStock,
		scalar_t& totalCost)  {

	totalCost = 0.0;
	for (int i=0; i<NUM_Subsystems; i++) {

		scalar_t currentIntermediateCost;
		scalar_t nextIntermediateCost;
		for (int k=0; k<timeTrajectoriesStock[i].size()-1; k++) {

			if (k==0) {
				subsystemCostFunctionsPtrStock_[i]->setCurrentStateAndControl(timeTrajectoriesStock[i][k], outputTrajectoriesStock[i][k], inputTrajectoriesStock[i][k]);
				subsystemCostFunctionsPtrStock_[i]->evaluate(currentIntermediateCost);
			} else {
				currentIntermediateCost = nextIntermediateCost;
			}

			// feed next state and control to cost function
			subsystemCostFunctionsPtrStock_[i]->setCurrentStateAndControl(timeTrajectoriesStock[i][k+1], outputTrajectoriesStock[i][k+1], inputTrajectoriesStock[i][k+1]);
			// evaluate intermediate cost for next time step
			subsystemCostFunctionsPtrStock_[i]->evaluate(nextIntermediateCost);

			totalCost += 0.5*(currentIntermediateCost+nextIntermediateCost)*(timeTrajectoriesStock[i][k+1]-timeTrajectoriesStock[i][k]);
		}

		// terminal cost
		if (i==NUM_Subsystems-1)  {
			scalar_t finalCost;
			subsystemCostFunctionsPtrStock_[i]->setCurrentStateAndControl(timeTrajectoriesStock[i].back(), outputTrajectoriesStock[i].back(), inputTrajectoriesStock[i].back());
			subsystemCostFunctionsPtrStock_[i]->terminalCost(finalCost);
			totalCost += finalCost;
		}
	}

}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
void LQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::approximateOptimalControlProblem()  {

	for (int i=0; i<NUM_Subsystems; i++) {

		subsystemDerivativesPtrStock_[i]->initializeModel(switchingTimes_, stateOperatingPointsStock_.at(i), i, "LQP");
		subsystemDerivativesPtrStock_[i]->setCurrentStateAndControl(0.0,
				stateOperatingPointsStock_.at(i), inputOperatingPointsStock_.at(i), outputOperatingPointsStock_.at(i));
		subsystemDerivativesPtrStock_[i]->getDerivativeState(AmStock_.at(i));
		subsystemDerivativesPtrStock_[i]->getDerivativesControl(BmStock_.at(i));

		if (runAsInitializer_==true) {
			subsystemDynamicsPtrStock_[i]->initializeModel(switchingTimes_, stateOperatingPointsStock_.at(i), i, "LQP");
			subsystemDynamicsPtrStock_[i]->computeDerivative(switchingTimes_.at(i), stateOperatingPointsStock_.at(i),
					inputOperatingPointsStock_.at(i), GvStock_.at(i));
		} else {
			GvStock_.at(i).setZero();
		}

		subsystemCostFunctionsPtrStock_[i]->setCurrentStateAndControl(0, outputOperatingPointsStock_.at(i), inputOperatingPointsStock_.at(i));
		subsystemCostFunctionsPtrStock_[i]->stateSecondDerivative(QmStock_.at(i));
		subsystemCostFunctionsPtrStock_[i]->controlSecondDerivative(RmStock_.at(i));
		subsystemCostFunctionsPtrStock_[i]->stateControlDerivative(PmStock_.at(i));
		RmInverseStock_.at(i) = RmStock_.at(i).inverse();

		if (runAsInitializer_==true) {
			QvStock_.at(i).setZero();
			RvStock_.at(i).setZero();
		} else {
			subsystemCostFunctionsPtrStock_[i]->stateDerivative(QvStock_.at(i));
			subsystemCostFunctionsPtrStock_[i]->controlDerivative(RvStock_.at(i));
		}

		if (INFO_ON_) {
			std::cout<< "stateOperatingPoint[" << i << "]: \n" << stateOperatingPointsStock_[i].transpose() << std::endl;
			std::cout<< "inputOperatingPoint[" << i << "]: \n" << inputOperatingPointsStock_[i].transpose() << std::endl;
			std::cout<< "Am[" << i << "]: \n" << AmStock_[i] << std::endl;
			std::cout<< "Bv[" << i << "]: \n" << BmStock_[i] << std::endl;
			std::cout<< "Gv[" << i << "]: \t" << GvStock_[i].transpose() << std::endl;
			std::cout<< "Qv[" << i << "]: \n" << QvStock_[i].transpose() << std::endl;
			std::cout<< "Qm[" << i << "]: \n" << QmStock_[i] << std::endl;
			std::cout<< "Rv[" << i << "]: \t" << RvStock_[i].transpose() << std::endl;
			std::cout<< "Rm[" << i << "]: \n" << RmStock_[i] << std::endl;
			std::cout<< "Pm[" << i << "]: \n" << PmStock_[i] << std::endl;
		}

		// making sure that Qm is PSD
		makePSD(QmStock_[i]);

		if (i==NUM_Subsystems-1)  {
			subsystemCostFunctionsPtrStock_[i]->terminalCostStateSecondDerivative(QmFinal_);
			// making sure that Qm is PSD
			makePSD(QmFinal_);

			if (runAsInitializer_==true)
				QvFinal_.setZero();
			else
				subsystemCostFunctionsPtrStock_[i]->terminalCostStateDerivative(QvFinal_);

			if (INFO_ON_) {
				std::cout<< "QvFinal[" << i << "]: \n" << QvFinal_.transpose() << std::endl;
				std::cout<< "QmFinal[" << i << "]: \n" << QmFinal_ << std::endl;
			}
		}
	}
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
void LQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::calculatecontroller(const scalar_t& learningRate, std::vector<controller_t>& controllersStock) {

	for (int i=0; i<NUM_Subsystems; i++) {

		controllersStock[i].time_ = timeTrajectoryStock_[i];

		size_t N = timeTrajectoryStock_[i].size();
		controllersStock[i].k_.resize(N);
		controllersStock[i].uff_.resize(N);
		for (int k=0; k<N; k++) {
			controllersStock[i].k_[k]    = -RmInverseStock_[i] * (PmStock_[i] + BmStock_[i].transpose()*SmTrajectoryStock_[i][k]);
			controllersStock[i].uff_[k]  = -learningRate * RmInverseStock_[i] * (RvStock_[i]  + BmStock_[i].transpose()*SvTrajectoryStock_[i][k])
								+ inputOperatingPointsStock_[i] - controllersStock[i].k_[k]*outputOperatingPointsStock_[i];
		}  // end of k loop

		if (INFO_ON_ ) {
			std::cout << "Controller of subsystem" << i << ":" << std::endl;
			std::cout << "learningRate " << learningRate << std::endl;
			std::cout << "time: " << controllersStock[i].time_.front() << std::endl;
			std::cout << "delta_uff: " <<  (controllersStock[i].uff_[0] + controllersStock[i].k_[0]*outputOperatingPointsStock_[i]).transpose() << std::endl;
			std::cout << "u0: " <<  inputOperatingPointsStock_[i].transpose() << std::endl;
			std::cout << "k: \n" <<  controllersStock[i].k_.front() << std::endl << std::endl;
		}
	}  // end of i loop
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
template <typename Derived>
bool LQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::makePSD(Eigen::MatrixBase<Derived>& squareMatrix) {

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
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
void LQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::getController(std::vector<controller_t>& controllersStock) {
	controllersStock = controllersStock_;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
void LQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::SolveRiccatiEquations()  {

	// general bvp solver
	SolveBVP<OUTPUT_DIM, INPUT_DIM> bvpSolver;

	// final value for the last Riccati equations
	output_vector_t SvFinal = output_vector_t::Zero();
	state_matrix_t  SmFinal = QmFinal_;

	for (int i=NUM_Subsystems-1; i>=0; i--) {

		scalar_array_t timeTrajectory{switchingTimes_[i], switchingTimes_[i+1]};

		state_matrix_array_t        AmTrajectory(2, AmStock_[i]);
		control_gain_matrix_array_t BmTrajectory(2, BmStock_[i]);
		output_vector_array_t       GvTrajectory(2, GvStock_[i]);

		output_vector_array_t  QvTrajectory(2, QvStock_[i]);
		state_matrix_array_t   QmTrajectory(2, QmStock_[i]);
		control_vector_array_t RvTrajectory(2, RvStock_[i]);
		control_matrix_array_t RmTrajectory(2, RmStock_[i]);
		control_matrix_array_t   RmInverseTrajectory(2, RmInverseStock_[i]);
		control_feedback_array_t PmTrajectory(2, PmStock_[i]);

		// set data for Riccati equations
		bvpSolver.setData(&timeTrajectory,
				&AmTrajectory, NULL,  &BmTrajectory, &GvTrajectory,
				&QvTrajectory, &QmTrajectory, &PmTrajectory,
				&RvTrajectory, &RmTrajectory, &RmInverseTrajectory);

		// solve BVP for the given time trajectory
		bvpSolver.solve(SvFinal, SmFinal,
				timeTrajectoryStock_[i], SmTrajectoryStock_[i], SvTrajectoryStock_[i],
				options_.AbsTolODE_, options_.RelTolODE_);

		//set the final value of the previous subsystem solver to the starting time value of the current subsystem's solution
		SvFinal = SvTrajectoryStock_[i].front();
		SmFinal = SmTrajectoryStock_[i].front();
	}  // end of i loop

}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
void LQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::run(const std::vector<scalar_t>& switchingTimes, const scalar_t& learningRate)  {

	if (switchingTimes.size() != NUM_Subsystems+1)
		throw std::runtime_error("Number of switching times should be one plus the number of subsystems.");
	switchingTimes_ = switchingTimes;

	// linearizing the dynamics and quadratizing the cost funtion along nominal trajectories
	approximateOptimalControlProblem();

	// solve Riccati equations
	SolveRiccatiEquations();

	// calculate controller
	calculatecontroller(learningRate, controllersStock_);
}

} // namespace ocs2
