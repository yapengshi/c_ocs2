/*
 * OCS2Ipopt.h
 *
 *  Created on: Jan 12, 2016
 *      Author: farbod
 */

#ifndef OCS2IPOPT_H_
#define OCS2IPOPT_H_

#include "IpIpoptApplication.hpp"

#include "ocs2/IpopotCostFunntion.h"

using namespace Ipopt;

template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
class OCS2Ipopt
{
public:
	typedef GSLQP<STATE_DIM, INPUT_DIM, NUM_Subsystems> GSLQP_t;

	typedef Dimensions<STATE_DIM, INPUT_DIM> DIMENSIONS;
	typedef typename DIMENSIONS::controller_t controller_t;
	typedef typename DIMENSIONS::scalar_t 		scalar_t;
	typedef typename DIMENSIONS::scalar_array_t scalar_array_t;
	typedef typename DIMENSIONS::eigen_scalar_t       eigen_scalar_t;
	typedef typename DIMENSIONS::eigen_scalar_array_t eigen_scalar_array_t;
	typedef typename DIMENSIONS::state_vector_t 	  state_vector_t;
	typedef typename DIMENSIONS::state_vector_array_t state_vector_array_t;
	typedef typename DIMENSIONS::control_vector_t 		control_vector_t;
	typedef typename DIMENSIONS::control_vector_array_t control_vector_array_t;
	typedef typename DIMENSIONS::control_feedback_t 	  control_feedback_t;
	typedef typename DIMENSIONS::control_feedback_array_t control_feedback_array_t;
	typedef typename DIMENSIONS::state_matrix_t 	  state_matrix_t;
	typedef typename DIMENSIONS::state_matrix_array_t state_matrix_array_t;
	typedef typename DIMENSIONS::control_matrix_t 		control_matrix_t;
	typedef typename DIMENSIONS::control_matrix_array_t control_matrix_array_t;
	typedef typename DIMENSIONS::control_gain_matrix_t 		 control_gain_matrix_t;
	typedef typename DIMENSIONS::control_gain_matrix_array_t control_gain_matrix_array_t;

	OCS2Ipopt (const std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM> > >& subsystemDynamicsPtr,
			const std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM> > >& subsystemDerivativesPtr,
			const std::vector<std::shared_ptr<CostFunctionBase<STATE_DIM, INPUT_DIM> > >& subsystemCostFunctionsPtr,
			const state_vector_array_t&   stateOperatingPoints,
			const control_vector_array_t& inputOperatingPoints,
			const std::vector<size_t>& systemStockIndex,
			const scalar_array_t& initSwitchingTimes,
			const state_vector_t& initState,
			const typename GSLQP_t::Options& options = GSLQP_t::Options() )
		: subsystemDynamicsPtr_(subsystemDynamicsPtr),
		  subsystemDerivativesPtr_(subsystemDerivativesPtr),
		  subsystemCostFunctionsPtr_(subsystemCostFunctionsPtr),
		  stateOperatingPoints_(stateOperatingPoints),
		  inputOperatingPoints_(inputOperatingPoints),
		  systemStockIndex_(systemStockIndex),
		  initSwitchingTimes_(initSwitchingTimes),
		  initState_(initState),
		  options_(options)
	{}
	virtual ~OCS2Ipopt() {}


	bool run ()  {

		// Create a new instance of your nlp
		//  (use a SmartPtr, not raw)
		SmartPtr<TNLP> ocs2Nlp = new IpopotCostFunntion<STATE_DIM, INPUT_DIM, NUM_Subsystems>(subsystemDynamicsPtr_, subsystemDerivativesPtr_, subsystemCostFunctionsPtr_,
				stateOperatingPoints_, inputOperatingPoints_, systemStockIndex_, initSwitchingTimes_, initState_, options_);

		// Create a new instance of IpoptApplication
		//  (use a SmartPtr, not raw)
		// We are using the factory, since this allows us to compile this
		// example with an Ipopt Windows DLL
		SmartPtr<IpoptApplication> ipoptApplication = IpoptApplicationFactory();

		// Change some options
		ipoptApplication->Options()->SetStringValue("hessian_approximation", "limited-memory");  // BFGS method
		ipoptApplication->Options()->SetNumericValue("tol", 1e-2);
		ipoptApplication->Options()->SetNumericValue("acceptable_tol", 1e-01);
		ipoptApplication->Options()->SetNumericValue("acceptable_obj_change_tol", 1e-01); //  This is useful for the quasi-Newton option, which has trouble to bring down the dual infeasibility.
		ipoptApplication->Options()->SetIntegerValue("max_iter", 20);
		ipoptApplication->Options()->SetIntegerValue("acceptable_iter", 2);
		ipoptApplication->Options()->SetStringValue("mu_strategy", "adaptive");
//		ipoptApplication->Options()->SetStringValue("output_file", "ipopt.out");
		ipoptApplication->Options()->SetStringValue("print_user_options", "yes");

		// Intialize the IpoptApplication and process the options
		ApplicationReturnStatus status;
		status = ipoptApplication->Initialize();
		if (status != Solve_Succeeded)  throw std::runtime_error("ipopt applicatiom initialization was not successful.");

		// Ask Ipopt to solve the problem
		status = ipoptApplication->OptimizeTNLP(ocs2Nlp);

		return (int)status;

	}


private:
	std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM> > > subsystemDynamicsPtr_;
	std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM> > > subsystemDerivativesPtr_;
	std::vector<std::shared_ptr<CostFunctionBase<STATE_DIM, INPUT_DIM> > > subsystemCostFunctionsPtr_;

	state_vector_array_t   stateOperatingPoints_;
	control_vector_array_t inputOperatingPoints_;

	std::vector<size_t> systemStockIndex_;

	scalar_array_t initSwitchingTimes_;
	state_vector_t initState_;

	typename GSLQP_t::Options options_;

};



#endif /* OCS2IPOPT_H_ */
