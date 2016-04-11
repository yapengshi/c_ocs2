/*
 * Implementation of OCS2Ipopt.h
 *
 *  Created on: Jan 18, 2016
 *      Author: farbod
 */


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
bool OCS2Ipopt<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::run ()  {

	// Create a new instance of IpoptApplication
	//  (use a SmartPtr, not raw)
	// We are using the factory, since this allows us to compile this
	// example with an Ipopt Windows DLL
	SmartPtr<IpoptApplication> ipoptApplication = IpoptApplicationFactory();

	// setting some o the IPOPT options
	ipoptApplication->Options()->SetStringValue("hessian_approximation", "limited-memory");  // BFGS method
	ipoptApplication->Options()->SetNumericValue("tol", options_.tolIPOPT_);
	ipoptApplication->Options()->SetNumericValue("acceptable_tol", options_.acceptableTolIPOPT_);
	ipoptApplication->Options()->SetNumericValue("acceptable_obj_change_tol", options_.acceptableTolIPOPT_); //  This is useful for the quasi-Newton option, which has trouble to bring down the dual infeasibility.
	ipoptApplication->Options()->SetIntegerValue("max_iter", options_.maxIterationIPOPT_);
	ipoptApplication->Options()->SetIntegerValue("acceptable_iter", 2);
	ipoptApplication->Options()->SetStringValue("mu_strategy", "adaptive");
	ipoptApplication->Options()->SetStringValue("derivative_test", "none");
	ipoptApplication->Options()->SetNumericValue("derivative_test_tol", 1e-6);
	ipoptApplication->Options()->SetStringValue("jac_d_constant", "yes");  // Indicates that all inequality constraints are linear
	ipoptApplication->Options()->SetStringValue("output_file", "/home/farbod/Programs/ct_ws/src/fbcc_hyq_controller/config/CoMOptimization/result/ipopt_info.txt");
	if (options_.displayIPOPT_)
		ipoptApplication->Options()->SetStringValue("print_user_options", "yes");
	else
		ipoptApplication->Options()->SetIntegerValue("print_level", 0);

	// Intialize the IpoptApplication and process the options
	ApplicationReturnStatus status;
	status = ipoptApplication->Initialize();
	if (status != Solve_Succeeded)  throw std::runtime_error("ipopt applicatiom initialization was not successful.");

	// Ask Ipopt to solve the problem
	status = ipoptApplication->OptimizeTNLP(ocs2Nlp_);

	// extract solutions
	rolloutSolution();

	return (int)status;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
void OCS2Ipopt<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::getCost(scalar_t& optimizedTotalCost,
		Eigen::Matrix<double,NUM_Subsystems-1,1>& optimizedTotalCostDerivative,
		scalar_t& ipoptOptimizedTotalCost)  const  {

	optimizedTotalCost = optimizedTotalCost_;
	optimizedTotalCostDerivative = optimizedTotalCostDerivative_;
	ipoptOptimizedTotalCost = ipoptOptimizedTotalCost_;
}



/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
void OCS2Ipopt<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::getController(scalar_array_t& optimizedSwitchingTimes,
		std::vector<controller_t>& optimizedControllersStock)  const  {

	optimizedSwitchingTimes = ipoptOptimizedSwitchingTimes_;
	optimizedControllersStock = ipoptOptimizedControllersStock_;
}



/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
void OCS2Ipopt<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::getTrajectories(
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
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
void OCS2Ipopt<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::getTrajectories(
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
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
void OCS2Ipopt<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::rolloutSolution()  {

	dynamic_cast<IpoptCostFunntion<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>*>(GetRawPtr(ocs2Nlp_))->getSolution(
			ipoptOptimizedTotalCost_, ipoptOptimizedSwitchingTimes_, ipoptOptimizedControllersStock_);

	for (size_t i=1; i<NUM_Subsystems; i++)
		if ((ipoptOptimizedSwitchingTimes_[i+1]-ipoptOptimizedSwitchingTimes_[i]) <= options_.minAcceptedSwitchingTimeDifference_)
			ipoptOptimizedSwitchingTimes_[i+1] = ipoptOptimizedSwitchingTimes_[i]+0.001;

	// GSLQP with one iteration
	Options_t nonIteratingOptions = options_;
//		nonIteratingOptions.maxIterationGSLQP_ = 1;
		nonIteratingOptions.dispayGSLQP_ = false;
	//
	GSLQP_t gslqp(subsystemDynamicsPtr_, subsystemDerivativesPtr_, subsystemCostFunctionsPtr_,
			ipoptOptimizedControllersStock_, systemStockIndex_, nonIteratingOptions);

	gslqp.run(initState_, ipoptOptimizedSwitchingTimes_);

	// trajectories
	gslqp.getNominalTrajectories(optimizedTimeTrajectoriesStock_, optimizedStateTrajectoriesStock_,
			optimizedInputTrajectoriesStock_, optimizedOutputTrajectoriesStock_);

	// cost funtion
	gslqp.rolloutCost(optimizedTimeTrajectoriesStock_, optimizedStateTrajectoriesStock_,
			optimizedInputTrajectoriesStock_, optimizedTotalCost_);

	// cost funtion jacobian
	gslqp.getCostFuntionDerivative(initState_, optimizedTotalCostDerivative_);
}




