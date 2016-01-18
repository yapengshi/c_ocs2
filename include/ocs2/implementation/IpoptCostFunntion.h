/*
 * Implementation of IpoptCostFunntion.h
 *
 *  Created on: Jan 12, 2016
 *      Author: farbod
 */

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
// returns the size of the problem
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
bool IpoptCostFunntion<STATE_DIM, INPUT_DIM, NUM_Subsystems>::get_nlp_info(Index& numParameters, Index& numConstraints,
		Index& nnz_jac_g, Index& nnz_h_lag,
		IndexStyleEnum& index_style)
{
	// number of optimization parameters
	numParameters = NumParameters_;
	// number of total equality ans inequality constraints
	numConstraints = NumConstraints_;
	// number of elements in constraint jacobian
	nnz_jac_g = 2*NumConstraints_;
	// number of elements in constraint hessian, but we
	// only need the lower left corner (since it is symmetric)
	nnz_h_lag = 0.5*(NumParameters_*(NumParameters_+1));
	// use the C style indexing (0-based)
	index_style = TNLP::C_STYLE;

	return true;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
// returns the variable bounds
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
bool IpoptCostFunntion<STATE_DIM, INPUT_DIM, NUM_Subsystems>::get_bounds_info(Index numParameters, Number* x_l, Number* x_u,
		Index numConstraints, Number* g_l, Number* g_u)
{
	// here, the numParameters and numConstraints we gave IPOPT in get_nlp_info are passed back to us.
	if (numParameters!=NumParameters_)  throw  std::runtime_error("numParameters is not correct.");
	if (numConstraints!=NumConstraints_) throw  std::runtime_error("numConstraints is not correct.");

	// parameters lower and upper bounds
	for (Index j=0; j<NumParameters_; j++)  {
		x_l[j] = -1e20;
		x_u[j] = +1e20;
	}
	x_l[0] = initSwitchingTimes_.front();
	x_u[NumParameters_-1] = initSwitchingTimes_.back();

	// constraints lower and upper bounds
	for (Index j=0; j<NumConstraints_; j++)  {
		g_l[j] = 0.0;
		g_u[j] = +1e20;
	}

	return true;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
// returns the initial point for the problem
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
bool IpoptCostFunntion<STATE_DIM, INPUT_DIM, NUM_Subsystems>::get_starting_point(Index numParameters, bool init_x, Number* x,
		bool init_z, Number* z_L, Number* z_U,
		Index m, bool init_lambda,
		Number* lambda)
{
	if (numParameters!=NumParameters_)  throw  std::runtime_error("numParameters is not correct.");
	if (init_x!=true)        throw  std::runtime_error("no request for initial parameters.");
	if (init_z!=false)       throw  std::runtime_error("no initial value for the bound multipliers is determined.");
	if (init_lambda!=false)  throw  std::runtime_error("no initial value for the constraint multipliers is determined.");

	// initialize to the given starting point
	for (Index j=0; j<NumParameters_; j++)
		x[j] = initSwitchingTimes_[j+1];

	return true;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
// returns the value of the objective function
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
bool IpoptCostFunntion<STATE_DIM, INPUT_DIM, NUM_Subsystems>::eval_f(Index numParameters, const Number* x, bool new_x, Number& obj_value)
{
	if (numParameters!=NumParameters_)  throw  std::runtime_error("numParameters is not correct.");

	if (new_x)  solveGSLQP(x);

	obj_value = currentTotalCost_;

	return true;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
// return the gradient of the objective function grad_{x} f(x)
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
bool IpoptCostFunntion<STATE_DIM, INPUT_DIM, NUM_Subsystems>::eval_grad_f(Index numParameters, const Number* x, bool new_x, Number* grad_f)
{
	if (numParameters!=NumParameters_)  throw  std::runtime_error("numParameters is not correct.");

	if (new_x)  solveGSLQP(x);

	for (size_t j=0; j<NumParameters_; j++)
		grad_f[j] = currentCostFuntionDerivative_[j];

	return true;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
// return the value of the constraints: g(x)
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
bool IpoptCostFunntion<STATE_DIM, INPUT_DIM, NUM_Subsystems>::eval_g(Index numParameters, const Number* x, bool new_x, Index numConstraints, Number* g)
{
	if (numParameters!=NumParameters_)  throw  std::runtime_error("numParameters is not correct.");
	if (numConstraints!=NumConstraints_) throw  std::runtime_error("numConstraints is not correct.");

	if (new_x)  solveGSLQP(x);

	// switching time vector
	scalar_array_t switchingTimes(initSwitchingTimes_);
	for (Index j=0; j<NumParameters_; j++)
		switchingTimes[j+1] = x[j];

	for (Index j=0; j<NumConstraints_; j++)
		g[j] = switchingTimes[j+2] - switchingTimes[j+1];

	return true;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
// return the structure or values of the jacobian
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
bool IpoptCostFunntion<STATE_DIM, INPUT_DIM, NUM_Subsystems>::eval_jac_g(Index numParameters, const Number* x, bool new_x,
		Index numConstraints, Index nele_jac, Index* iRow, Index *jCol,
		Number* values)
{
	if (numParameters!=NumParameters_)  throw  std::runtime_error("numParameters is not correct.");
	if (numConstraints!=NumConstraints_) throw  std::runtime_error("numConstraints is not correct.");

	if (new_x)  solveGSLQP(x);

	if (values == NULL) {
		// return the structure of the jacobian: (i,j) element is dg_i/dx_j
		size_t iterator = 0;
		for (size_t i=0; i<NumConstraints_; i++) {
				iRow[iterator] = i;  jCol[iterator] = i;
				iterator++;
				iRow[iterator] = i;  jCol[iterator] = i+1;
				iterator++;
			}
	}
	else {
		// return the values of the jacobian of the constraints
		size_t iterator = 0;
		for (size_t i=0; i<NumConstraints_; i++) {
			values[iterator] = -1.0;
			iterator++;
			values[iterator] = 1.0;
			iterator++;
		}
	}

	return true;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
// return the structure or hessian of the lagrangian
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
bool IpoptCostFunntion<STATE_DIM, INPUT_DIM, NUM_Subsystems>::eval_h(Index n, const Number* x, bool new_x,
			Number obj_factor, Index m, const Number* lambda,
			bool new_lambda, Index nele_hess, Index* iRow,
			Index* jCol, Number* values)  {

	if (new_x)  solveGSLQP(x);

	return false;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
void IpoptCostFunntion<STATE_DIM, INPUT_DIM, NUM_Subsystems>::finalize_solution(SolverReturn status,
		Index numParameters, const Number* x, const Number* z_L, const Number* z_U,
		Index numConstraints, const Number* g, const Number* lambda,
		Number optimizedTotalCost,
		const IpoptData* ip_data,
		IpoptCalculatedQuantities* ip_cq)
{

	// cost funtion
	optimizedTotalCost_ = optimizedTotalCost;

	// switching times
	optimizedSwitchingTimes_.front() = initSwitchingTimes_.front();
	optimizedSwitchingTimes_.back()  = initSwitchingTimes_.back();
	for (Index j=0; j<NumParameters_; j++)
		optimizedSwitchingTimes_[j+1] = x[j];

	// controller
	size_t index = findNearestController(x);
	optimizedControllersStock_ = controllersStockBag_[index];


	if (options_.displayIPOPT_) {
		std::cout << "\n## Optimal cost: " << optimizedTotalCost_ << std::endl;
		std::cout << "Number of funtion call: " << numFuntionCall_ << std::endl;

		std::cout << "Switching times are: [" << optimizedSwitchingTimes_[0] << ", ";
		for (size_t i=0; i<NumParameters_; i++)
			std::cout << optimizedSwitchingTimes_[i+1] << ", ";
		std::cout << optimizedSwitchingTimes_.back() << "]" << std::endl;
	}

}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
void IpoptCostFunntion<STATE_DIM, INPUT_DIM, NUM_Subsystems>::getSolution(scalar_t& optimizedTotalCost,
		scalar_array_t& optimizedSwitchingTimes, std::vector<controller_t>& optimizedControllersStock)  const  {

	optimizedTotalCost = optimizedTotalCost_;
	optimizedSwitchingTimes = optimizedSwitchingTimes_;
	optimizedControllersStock = optimizedControllersStock_;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
size_t IpoptCostFunntion<STATE_DIM, INPUT_DIM, NUM_Subsystems>::findNearestController(const Number* x) const  {

	if (parameterBag_.size()==0)  throw  std::runtime_error("controllerStock bag is empty.");

	Eigen::Matrix<double,NumParameters_,1> enquiry = Eigen::Map<const Eigen::Matrix<double,NumParameters_,1> >(x);

	// evaluating distance
	std::vector<double> distance(parameterBag_.size());
	for (size_t i=0; i<parameterBag_.size(); i++)
		distance[i] = (parameterBag_[i]-enquiry).squaredNorm();

	// min index
	auto it = std::min_element(distance.begin(), distance.end());
	size_t index = std::distance(distance.begin(), it);

	return index;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
void IpoptCostFunntion<STATE_DIM, INPUT_DIM, NUM_Subsystems>::solveGSLQP(const Number* x)  {

	// switching time vector
	scalar_array_t switchingTimes(initSwitchingTimes_);
	for (Index j=0; j<NumParameters_; j++)
		switchingTimes[j+1] = x[j];

//	std::cout << "switching time: ";
//	for (Index j=0; j<switchingTimes.size(); j++)
//		std::cout << switchingTimes[j] << ", ";
//	std::cout << std::endl;

	std::vector<controller_t> controllersStock(NUM_Subsystems);
	if (parameterBag_.size()==0 || options_.warmStartGSLQP_==false) {

		std::cout << "==> IPOPT cost funtion call WITHOUT Memorization" << std::endl;

		// GLQP initialization
		GLQP_t glqp(subsystemDynamicsPtr_, subsystemDerivativesPtr_, subsystemCostFunctionsPtr_,
				stateOperatingPoints_, inputOperatingPoints_, systemStockIndex_);
		glqp.run(switchingTimes);

		// GLQP controller
		glqp.getController(controllersStock);

//		////
//		std::vector<scalar_array_t> timeTrajectoriesStock;
//		std::vector<state_vector_array_t> stateTrajectoriesStock;
//		std::vector<control_vector_array_t> inputTrajectoriesStock;
//
//		// rollout
//		glqp.rollout(initState_, controllersStock, timeTrajectoriesStock, stateTrajectoriesStock, inputTrajectoriesStock);
//
//		// compute test rollout cost
//		double rolloutCost;
//		glqp.rolloutCost(timeTrajectoriesStock, stateTrajectoriesStock, inputTrajectoriesStock, rolloutCost);
//		std::cout << "GLQP rollout cost: " << rolloutCost << std::endl;
	}
	else {
		std::cout << "==> IPOPT cost funtion call WITH memorization" << std::endl;

		// find nearest controller
		size_t index = findNearestController(x);
		controllersStock = controllersStockBag_[index];
	}

	// GSLQP
	GSLQP_t gslqp(subsystemDynamicsPtr_, subsystemDerivativesPtr_, subsystemCostFunctionsPtr_,
			controllersStock, systemStockIndex_, options_);
	gslqp.run(initState_, switchingTimes);

	// cost funtion
	gslqp.getValueFuntion(0.0, initState_, currentTotalCost_);

	// cost funtion jacobian
	gslqp.getCostFuntionDerivative(initState_, currentCostFuntionDerivative_);

	numFuntionCall_++;

	// GSLQP controller
	gslqp.getController(controllersStock);

	// saving to bag
	parameterBag_.push_back(Eigen::Map<const Eigen::Matrix<double,NumParameters_,1> >(x));
	controllersStockBag_.push_back(controllersStock);

}
