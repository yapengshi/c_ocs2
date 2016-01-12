/*
 * OCS2.h
 *
 *  Created on: Jan 11, 2016
 *      Author: farbod
 */

#ifndef OCS2_H_
#define OCS2_H_

#include <ceres/ceres.h>
#include <glog/logging.h>

#include "ocs2/LeastSquareCost.h"

template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
class OCS2
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

	OCS2 (const std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM> > >& subsystemDynamicsPtr,
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
	virtual ~OCS2() {}


	void run ()  {

		google::InitGoogleLogging("OCS2");

		// The variable to solve for with its initial value. It will be
		// mutated in place by the solver.
		scalar_array_t switchingTimes(initSwitchingTimes_);
		double parameters[NUM_Subsystems-1];
		for (size_t j=0; j<NUM_Subsystems-1; j++)
			parameters[j] = switchingTimes[j+1];

		// Build the problem.
		ceres::Problem problem;
		// Set up the only cost function (also known as residual).
		ceres::CostFunction* cost_function = new LeastSquareCost<STATE_DIM, INPUT_DIM, NUM_Subsystems>(subsystemDynamicsPtr_, subsystemDerivativesPtr_, subsystemCostFunctionsPtr_,
			stateOperatingPoints_,inputOperatingPoints_, systemStockIndex_, switchingTimes.front(), switchingTimes.back(), initState_, options_);
		problem.AddResidualBlock(cost_function, NULL, parameters);

		// Run the solver!
		ceres::Solver::Options ceresOptions;
		ceresOptions.minimizer_progress_to_stdout = true;
		ceres::Solver::Summary summary;
		ceres::Solve(ceresOptions, &problem, &summary);

		for (size_t j=1; j<NUM_Subsystems; j++)
			switchingTimes[j] = parameters[j-1];

		std::cout << summary.BriefReport() << std::endl;
		std::cout << "Switching times are: [" << switchingTimes[0] << ", ";
		for (size_t j=1; j<NUM_Subsystems; j++)
			std::cout << switchingTimes[j] << ", ";
		std::cout << switchingTimes.back() << "]\n";

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

#endif /* OCS2_H_ */
