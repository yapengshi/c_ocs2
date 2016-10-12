/*
 * OCS2Projected.h
 *
 *  Created on: Jul 21, 2016
 *      Author: farbod
 */

namespace ocs2{

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void OCS2Projected<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getCostFunction(scalar_t& costFunction) const  {
	costFunction = optimizedTotalCost_;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void OCS2Projected<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getCostFunctionDerivative(
		Eigen::Matrix<double,NUM_SUBSYSTEMS-1,1>& costFuntionDerivative) const {
	costFuntionDerivative = costFuntionDerivative_;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void OCS2Projected<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getSwitchingTimes(
		scalar_array_t& switchingTimes) const {
	switchingTimes = optimizedSwitchingTimes_;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void OCS2Projected<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getController(
		std::vector<controller_t>& optimizedControllersStock)  const  {
	optimizedControllersStock = optimizedControllersStock_;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void OCS2Projected<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getTrajectories(
		std::vector<scalar_array_t>& optimizedTimeTrajectoriesStock,
		std::vector<state_vector_array_t>& optimizedStateTrajectoriesStock,
		std::vector<control_vector_array_t>& optimizedInputTrajectoriesStock,
		std::vector<output_vector_array_t>& optimizedOutputTrajectoriesStock) const {

	optimizedTimeTrajectoriesStock   = optimizedTimeTrajectoriesStock_;
	optimizedStateTrajectoriesStock  = optimizedStateTrajectoriesStock_;
	optimizedInputTrajectoriesStock  = optimizedInputTrajectoriesStock_;
	optimizedOutputTrajectoriesStock = optimizedOutputTrajectoriesStock_;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void OCS2Projected<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getTrajectories(
		std::vector<scalar_array_t>& optimizedTimeTrajectoriesStock,
		std::vector<state_vector_array_t>& optimizedStateTrajectoriesStock,
		std::vector<control_vector_array_t>& optimizedInputTrajectoriesStock) const {

	optimizedTimeTrajectoriesStock   = optimizedTimeTrajectoriesStock_;
	optimizedStateTrajectoriesStock  = optimizedStateTrajectoriesStock_;
	optimizedInputTrajectoriesStock  = optimizedInputTrajectoriesStock_;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
size_t OCS2Projected<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::findNearestController(const Eigen::VectorXd& enquiryParameter) const  {

	if (parameterBag_.size()==0)  throw  std::runtime_error("controllerStock bag is empty.");

	// evaluating distance
	std::vector<double> distance(parameterBag_.size());
	for (size_t i=0; i<parameterBag_.size(); i++)
		distance[i] = (parameterBag_[i]-enquiryParameter).squaredNorm();

	// min index
	auto it = std::min_element(distance.begin(), distance.end());
	size_t index = std::distance(distance.begin(), it);

	return index;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void OCS2Projected<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateLinearEqualityConstraint(
		Eigen::MatrixXd& Am, Eigen::VectorXd& Bv) {

	Am = Eigen::MatrixXd::Zero(numParameters_+1, numParameters_);
	for (size_t i=0; i<numParameters_+1; i++) {
		if (i<numParameters_) 	Am(i,i)  = -1.0;
		if (i>0)				Am(i,i-1) = 1.0;
	}

	Bv = Eigen::VectorXd::Zero(numParameters_+1);
	Bv(0) = initSwitchingTimes_.front();
	Bv(numParameters_) = -initSwitchingTimes_.back();
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
bool OCS2Projected<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateGradient(const size_t& id,
		const Eigen::VectorXd& parameters, Eigen::VectorXd& gradient) {

	// switching time vector
	scalar_array_t switchingTimes(initSwitchingTimes_);
	for (size_t j=0; j<numParameters_; j++)
		switchingTimes[j+1] = parameters(j);

	gslqpSolver_->run(initState_, switchingTimes, slqpSolverPtrs_[id]);
	gslqpSolver_->getCostFuntionDerivative(costFuntionDerivative_);
	gradient = costFuntionDerivative_;

	return true;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
bool OCS2Projected<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateCost(const size_t& id,
		const Eigen::VectorXd& parameters, double& cost) {

	// switching time vector
	scalar_array_t switchingTimes(initSwitchingTimes_);
	for (size_t j=0; j<numParameters_; j++)
		switchingTimes[j+1] = parameters(j);

	// initial controller
	std::vector<controller_t> controllersStock(NUM_SUBSYSTEMS);
	calculateInitialController(initState_, switchingTimes, controllersStock);

	// run SLQP
	try {
		slqpSolverPtrs_[id]->run(initState_, switchingTimes, controllersStock);
	}
	catch (const std::exception& e){
		std::cerr << "\t     exception: " << e.what();
		return false;
	}

	double unconstraintCost, constraintISE;
	slqpSolverPtrs_[id]->getCostFuntion(unconstraintCost, constraintISE);
	cost = unconstraintCost;

	// display
	if (nlpOptions_.displayGradientDescent_)  std::cerr << "\t     constraintISE: " << constraintISE
			<< "\t#Iterations: " << slqpSolverPtrs_[id]->getNumIterations() << std::endl;

	// saving solution in the bag
	saveToBag(id, parameters);

	// status is false if the constraintISE is higher than minAbsConstraint1RMSE_
	bool status = (constraintISE <= options_.minAbsConstraint1ISE_) ? true : false;

	return status;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void OCS2Projected<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::saveToBag(size_t id, const Eigen::VectorXd& parameters)  {

	// get the nominal trajectories
	std::vector<scalar_array_t> timeTrajectoriesStock(NUM_SUBSYSTEMS);
	std::vector<state_vector_array_t>   stateTrajectoriesStock(NUM_SUBSYSTEMS);
	std::vector<control_vector_array_t> inputTrajectoriesStock(NUM_SUBSYSTEMS);
	std::vector<output_vector_array_t>  outputTrajectoriesStock(NUM_SUBSYSTEMS);
	slqpSolverPtrs_[id]->getNominalTrajectories(timeTrajectoriesStock, stateTrajectoriesStock, inputTrajectoriesStock, outputTrajectoriesStock);

	// get the optimized controller
	std::vector<controller_t> controllersStock(NUM_SUBSYSTEMS);
	slqpSolverPtrs_[id]->getController(controllersStock);

	// changing the controller structure to tracking controller
	LinearInterpolation<output_vector_t,Eigen::aligned_allocator<output_vector_t> >   nominalOutputFunc;
	LinearInterpolation<control_vector_t,Eigen::aligned_allocator<control_vector_t> > nominalInputFunc;
	for (size_t i=0; i<NUM_SUBSYSTEMS; i++) {

		nominalOutputFunc.setTimeStamp( &(timeTrajectoriesStock[i]) );
		nominalOutputFunc.setData( &(outputTrajectoriesStock[i]) );

		nominalInputFunc.setTimeStamp( &(timeTrajectoriesStock[i]) );
		nominalInputFunc.setData( &(inputTrajectoriesStock[i]) );

		controllersStock[i].deltaUff_.resize(controllersStock[i].time_.size());
		for (size_t k=0; k<controllersStock[i].time_.size(); k++) {

			const double& time = controllersStock[i].time_[k];
			output_vector_t nominalOutput;
			nominalOutputFunc.interpolate(time, nominalOutput);
			size_t greatestLessTimeStampIndex = nominalOutputFunc.getGreatestLessTimeStampIndex();
			control_vector_t nominalInput;
			nominalInputFunc.interpolate(time, nominalInput, greatestLessTimeStampIndex);

			controllersStock[i].uff_[k] = -controllersStock[i].k_[k]*nominalOutput;
			controllersStock[i].deltaUff_[k] = nominalInput;

		} // end of k loop
	}  // end of i loop

	// save the parameter and controller in the Bag
	parameterBag_.push_back(parameters);
	controllersStockBag_.push_back(controllersStock);
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void OCS2Projected<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateInitialController(const state_vector_t& initState,
		const scalar_array_t& switchingTimes,
		std::vector<controller_t>& controllersStock) {

	std::vector<controller_t> coldStartControllersStock(NUM_SUBSYSTEMS);
	std::vector<controller_t> warmStartControllersStock(NUM_SUBSYSTEMS);
	scalar_t coldStartTotalCost = std::numeric_limits<double>::max();
	scalar_t warmStartTotalCost = std::numeric_limits<double>::max();

	// using GLQP for coldStart controller
	GLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS> glqp(subsystemDynamicsPtr_, subsystemDerivativesPtr_, subsystemCostFunctionsPtr_,
			stateOperatingPoints_, inputOperatingPoints_, systemStockIndex_);
	glqp.run(switchingTimes, 0.0 /*learning rate*/);  // since learning rate is zero, the feedforward input is zero
	glqp.getController(coldStartControllersStock);

	// calculate the coldStart controllers' cost
	std::vector<scalar_array_t> timeTrajectoriesStock(NUM_SUBSYSTEMS);
	std::vector<state_vector_array_t> stateTrajectoriesStock(NUM_SUBSYSTEMS);
	std::vector<control_vector_array_t> inputTrajectoriesStock(NUM_SUBSYSTEMS);
	std::vector<output_vector_array_t> outputTrajectoriesStock(NUM_SUBSYSTEMS);
	rollout(initState_, switchingTimes, coldStartControllersStock,
			timeTrajectoriesStock, stateTrajectoriesStock, inputTrajectoriesStock, outputTrajectoriesStock);
	calculateCostFunction(timeTrajectoriesStock, outputTrajectoriesStock, inputTrajectoriesStock, coldStartTotalCost);


	size_t warmStartIndex;
	if (options_.warmStartGSLQP_==true && parameterBag_.size()>0) {

		// find most similar controller based on the parameter
		const Eigen::VectorXd parameters = Eigen::Map<const Eigen::VectorXd>(&switchingTimes[1], numParameters_);
		warmStartIndex = findNearestController(parameters);

		warmStartControllersStock = controllersStockBag_[warmStartIndex];
		// scaling the controller time to the current switchingTimes
		for (size_t i=0; i<NUM_SUBSYSTEMS; i++) {
			double scale = (switchingTimes[i+1]-switchingTimes[i]) / (warmStartControllersStock[i].time_.back()-warmStartControllersStock[i].time_.front());
			for (size_t k=0; k<warmStartControllersStock[i].time_.size(); k++) {
				warmStartControllersStock[i].time_[k] = switchingTimes[i] + scale*(warmStartControllersStock[i].time_[k]-warmStartControllersStock[i].time_.front());
				/*previously used by farbod to scale velocities with switching times: */
				//				warmStartControllersStock[i].uff_[k].head(12) += warmStartControllersStock[i].deltaUff_[k].head(12);
				//				warmStartControllersStock[i].uff_[k].tail(12) += warmStartControllersStock[i].deltaUff_[k].tail(12) / scale;
			} // end of k loop
		}  // end of i loop

		// calculate the warmStart controllers' cost
		std::vector<scalar_array_t> timeTrajectoriesStock(NUM_SUBSYSTEMS);
		std::vector<state_vector_array_t> stateTrajectoriesStock(NUM_SUBSYSTEMS);
		std::vector<control_vector_array_t> inputTrajectoriesStock(NUM_SUBSYSTEMS);
		std::vector<output_vector_array_t> outputTrajectoriesStock(NUM_SUBSYSTEMS);
		try {
			rollout(initState_, switchingTimes, warmStartControllersStock,
					timeTrajectoriesStock, stateTrajectoriesStock, inputTrajectoriesStock, outputTrajectoriesStock);
			calculateCostFunction(timeTrajectoriesStock, outputTrajectoriesStock, inputTrajectoriesStock,
					warmStartTotalCost);
		} catch (const std::exception& e) {}

	}

	// choose which controller to use based on the cost function
	if (warmStartTotalCost<=coldStartTotalCost) {
		controllersStock.swap(warmStartControllersStock);
		std::cerr << "\t     Warm start!" << std::endl;
		if (nlpOptions_.displayGradientDescent_)  std::cerr << "\t     IndexFromBack: " << static_cast<int>(warmStartIndex-(parameterBag_.size()-1))
						<< "\t#parameters: " << parameterBag_[warmStartIndex].transpose().format(CleanFmtDisplay_) << std::endl;
	} else {
		controllersStock.swap(coldStartControllersStock);
		std::cerr << "\t     Cold start!" << std::endl;
	}

}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void OCS2Projected<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::rollout(const state_vector_t& initState,
		const scalar_array_t& switchingTimes,
		const std::vector<controller_t>& controllersStock,
		std::vector<scalar_array_t>& timeTrajectoriesStock,
		std::vector<state_vector_array_t>& stateTrajectoriesStock,
		std::vector<control_vector_array_t>& inputTrajectoriesStock,
		std::vector<output_vector_array_t>& outputTrajectoriesStock)  {

	if (controllersStock.size() != NUM_SUBSYSTEMS)
		throw std::runtime_error("controllersStock has less controllers then the number of subsystems");

	timeTrajectoriesStock.resize(NUM_SUBSYSTEMS);
	stateTrajectoriesStock.resize(NUM_SUBSYSTEMS);
	inputTrajectoriesStock.resize(NUM_SUBSYSTEMS);
	outputTrajectoriesStock.resize(NUM_SUBSYSTEMS);

	state_vector_t x0 = initState;
	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		timeTrajectoriesStock[i].clear();
		stateTrajectoriesStock[i].clear();

		//		size_t maxNumSteps = options_.maxNumStepsPerSecond_*(switchingTimes_[i+1]-switchingTimes_[i]);
		//		maxNumSteps = ((1000>maxNumSteps) ? 1000 : maxNumSteps);
		size_t maxNumSteps = options_.maxNumStepsPerSecond_ * std::max( 1.0, switchingTimes[i+1]-switchingTimes[i] );

		// initialize subsystem i
		subsystemDynamicsPtrStock_[i]->initializeModel(switchingTimes, x0, i, "GSLPQ");
		// set controller for subsystem i
		subsystemDynamicsPtrStock_[i]->setController(controllersStock[i]);
		// simulate subsystem i
		subsystemSimulatorsStockPtr_[i]->integrate(x0, switchingTimes[i], switchingTimes[i+1],
				stateTrajectoriesStock[i], timeTrajectoriesStock[i],
				1e-3, options_.AbsTolODE_, options_.RelTolODE_, maxNumSteps);

		if (stateTrajectoriesStock[i].back() != stateTrajectoriesStock[i].back())
			throw std::runtime_error("System became unstable during the SLQP rollout.");

		// compute control trajectory for subsystem i
		inputTrajectoriesStock[i].resize(timeTrajectoriesStock[i].size());
		outputTrajectoriesStock[i].resize(timeTrajectoriesStock[i].size());
		for (int k=0; k<timeTrajectoriesStock[i].size(); k++) {
			subsystemDynamicsPtrStock_[i]->computeOutput(timeTrajectoriesStock[i][k], stateTrajectoriesStock[i][k], outputTrajectoriesStock[i][k]);
			subsystemDynamicsPtrStock_[i]->computeInput(timeTrajectoriesStock[i][k], outputTrajectoriesStock[i][k], inputTrajectoriesStock[i][k]);
		}

		// reset the initial state
		x0 = stateTrajectoriesStock[i].back();
	}
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void OCS2Projected<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateCostFunction(
		const std::vector<scalar_array_t>& timeTrajectoriesStock,
		const std::vector<output_vector_array_t>& outputTrajectoriesStock,
		const std::vector<control_vector_array_t>& inputTrajectoriesStock,
		scalar_t& totalCost)  {

	totalCost = 0.0;
	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		// integrates the intermediate cost using the trapezoidal approximation method
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
		}  // end of k loop

		// terminal cost
		if (i==NUM_SUBSYSTEMS-1)  {
			scalar_t finalCost;
			subsystemCostFunctionsPtrStock_[i]->setCurrentStateAndControl(timeTrajectoriesStock[i].back(), outputTrajectoriesStock[i].back(), inputTrajectoriesStock[i].back());
			subsystemCostFunctionsPtrStock_[i]->terminalCost(finalCost);
			totalCost += finalCost;
		}

	}  // end of i loop
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void OCS2Projected<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getSolution(size_t idStar)  {

	slqpSolverPtrs_[idStar]->getCostFuntion(optimizedTotalCost_, optimizedConstraintISE_);
	slqpSolverPtrs_[idStar]->getSwitchingTimes(optimizedSwitchingTimes_);
	slqpSolverPtrs_[idStar]->getController(optimizedControllersStock_);
	slqpSolverPtrs_[idStar]->getNominalTrajectories(optimizedTimeTrajectoriesStock_, optimizedStateTrajectoriesStock_,
			optimizedInputTrajectoriesStock_, optimizedOutputTrajectoriesStock_);
	slqpSolverPtrs_[idStar]->getIterationsLog(slqIterationCost_, slqIterationISE1_);
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void OCS2Projected<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::run(const state_vector_t& initState,
		const scalar_array_t& switchingTimes) {

	initState_ = initState;
	initSwitchingTimes_ = switchingTimes;

	// GLQP controller
	GLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS> glqp(subsystemDynamicsPtr_, subsystemDerivativesPtr_, subsystemCostFunctionsPtr_,
			stateOperatingPoints_, inputOperatingPoints_, systemStockIndex_);
	glqp.run(initSwitchingTimes_, 0.0 /*learning rate*/);  // since learning rate is zero, the feedforward input is zero
	glqp.getController(initControllersStock_);

	// SLQP solvers
	slqpSolverPtrs_.resize(numLineSearch_+1);
	for (size_t i=0; i<slqpSolverPtrs_.size(); i++)
		if (options_.useMultiThreading_==true)
			slqpSolverPtrs_[i] = slqp_base_ptr_t( new slqp_mp_t(subsystemDynamicsPtr_, subsystemDerivativesPtr_, subsystemCostFunctionsPtr_,
					initControllersStock_, systemStockIndex_, options_) );
		else
			slqpSolverPtrs_[i] = slqp_base_ptr_t( new slqp_t(subsystemDynamicsPtr_, subsystemDerivativesPtr_, subsystemCostFunctionsPtr_,
					initControllersStock_, systemStockIndex_, options_) );

	// GSLQP solvers
	gslqpSolver_ = gslqp_ptr_t( new gslqp_t(subsystemDynamicsPtr_, subsystemDerivativesPtr_, subsystemCostFunctionsPtr_,
			initControllersStock_, systemStockIndex_, options_) );

	// run
	Eigen::VectorXd initParameters = Eigen::Map<Eigen::VectorXd>(initSwitchingTimes_.data()+1, NUM_SUBSYSTEMS-1);
	GradientDescent::run(initParameters);

}







}  // end of ocs2 namespace
