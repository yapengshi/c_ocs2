/*
 * Implementation of OCS2Ipopt.h
 *
 *  Created on: Jan 18, 2016
 *      Author: farbod
 */


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
bool OCS2Ipopt<STATE_DIM, INPUT_DIM, NUM_Subsystems>::run ()  {

	// Create a new instance of IpoptApplication
	//  (use a SmartPtr, not raw)
	// We are using the factory, since this allows us to compile this
	// example with an Ipopt Windows DLL
	SmartPtr<IpoptApplication> ipoptApplication = IpoptApplicationFactory();

	// Change some options
	ipoptApplication->Options()->SetStringValue("hessian_approximation", "limited-memory");  // BFGS method
	ipoptApplication->Options()->SetNumericValue("tol", options_.tolIPOPT_);
	ipoptApplication->Options()->SetNumericValue("acceptable_tol", options_.acceptableTolIPOPT_);
	ipoptApplication->Options()->SetNumericValue("acceptable_obj_change_tol", options_.acceptableTolIPOPT_); //  This is useful for the quasi-Newton option, which has trouble to bring down the dual infeasibility.
	ipoptApplication->Options()->SetIntegerValue("max_iter", options_.maxIterationIPOPT_);
	ipoptApplication->Options()->SetIntegerValue("acceptable_iter", 2);
	ipoptApplication->Options()->SetStringValue("mu_strategy", "adaptive");
	ipoptApplication->Options()->SetStringValue("derivative_test", "first-order");
	ipoptApplication->Options()->SetNumericValue("derivative_test_tol", 1e-6);
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
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
void OCS2Ipopt<STATE_DIM, INPUT_DIM, NUM_Subsystems>::getCost(scalar_t& optimizedTotalCost,
		Eigen::Matrix<double,NUM_Subsystems-1,1>& optimizedTotalCostDerivative,
		scalar_t& ipoptOptimizedTotalCost)  const  {

	optimizedTotalCost = optimizedTotalCost_;
	optimizedTotalCostDerivative = optimizedTotalCostDerivative_;
	ipoptOptimizedTotalCost = ipoptOptimizedTotalCost_;
}



/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
void OCS2Ipopt<STATE_DIM, INPUT_DIM, NUM_Subsystems>::getController(scalar_array_t optimizedSwitchingTimes,
		std::vector<controller_t>& optimizedControllersStock)  const  {

	optimizedSwitchingTimes = ipoptOptimizedSwitchingTimes_;
	optimizedControllersStock = ipoptOptimizedControllersStock_;
}



/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
void OCS2Ipopt<STATE_DIM, INPUT_DIM, NUM_Subsystems>::getTrajectories(
		std::vector<scalar_array_t>& optimizedTimeTrajectoriesStock,
		std::vector<state_vector_array_t>& optimizedStateTrajectoriesStock,
		std::vector<control_vector_array_t>& optimizedInputTrajectoriesStock) const {

	optimizedTimeTrajectoriesStock  = optimizedTimeTrajectoriesStock_;
	optimizedStateTrajectoriesStock = optimizedStateTrajectoriesStock_;
	optimizedInputTrajectoriesStock = optimizedInputTrajectoriesStock_;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
void OCS2Ipopt<STATE_DIM, INPUT_DIM, NUM_Subsystems>::rolloutSolution()  {

	static_cast<IpoptCostFunntion<STATE_DIM, INPUT_DIM, NUM_Subsystems>*>(GetRawPtr(ocs2Nlp_))->getSolution(
			ipoptOptimizedTotalCost_, ipoptOptimizedSwitchingTimes_, ipoptOptimizedControllersStock_);

	// GSLQP
	Options_t nonIteratingOptions = options_;
		nonIteratingOptions.maxIterationGSLQP_ = 1;
		nonIteratingOptions.dispayGSLQP_ = false;
	//
	GSLQP<STATE_DIM, INPUT_DIM, NUM_Subsystems> gslqp(subsystemDynamicsPtr_, subsystemDerivativesPtr_, subsystemCostFunctionsPtr_,
			ipoptOptimizedControllersStock_, systemStockIndex_, nonIteratingOptions);

	gslqp.run(initState_, ipoptOptimizedSwitchingTimes_);

	// trajectories
	gslqp.getNominalTrajectories(optimizedTimeTrajectoriesStock_, optimizedStateTrajectoriesStock_,
			optimizedInputTrajectoriesStock_);

	// cost funtion
	gslqp.rolloutCost(optimizedTimeTrajectoriesStock_, optimizedStateTrajectoriesStock_,
			optimizedInputTrajectoriesStock_, optimizedTotalCost_);

	// cost funtion jacobian
	gslqp.getCostFuntionDerivative(initState_, optimizedTotalCostDerivative_);
}




