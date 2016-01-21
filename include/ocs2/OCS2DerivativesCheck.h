/*
 * OCS2DerivativesCheck.h
 *
 *  Created on: Jan 19, 2016
 *      Author: farbod
 */

#ifndef OCS2DERIVATIVESCHECK_H_
#define OCS2DERIVATIVESCHECK_H_


#include "GSLQ/GSLQPSolver.h"

template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
class OCS2DerivativesCheck
{
public:
	typedef GLQP<STATE_DIM, INPUT_DIM, NUM_Subsystems> GLQP_t;
	typedef GSLQP<STATE_DIM, INPUT_DIM, NUM_Subsystems> GSLQP_t;

	typedef Dimensions<STATE_DIM, INPUT_DIM> DIMENSIONS;
	typedef typename DIMENSIONS::controller_t controller_t;
	typedef typename DIMENSIONS::Options Options_t;
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

	typedef Eigen::Matrix<double, NUM_Subsystems-1, 1> parameters_t;


	OCS2DerivativesCheck(const std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM> > >& subsystemDynamicsPtr,
			const std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM> > >& subsystemDerivativesPtr,
			const std::vector<std::shared_ptr<CostFunctionBase<STATE_DIM, INPUT_DIM> > >& subsystemCostFunctionsPtr,
			const state_vector_array_t&   stateOperatingPoints,
			const control_vector_array_t& inputOperatingPoints,
			const std::vector<size_t>& systemStockIndex)

	: gslqpSolver_(subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr,
	  				stateOperatingPoints, inputOperatingPoints, systemStockIndex, Options_t())
	{}

	virtual ~OCS2DerivativesCheck() {}


	void check(const state_vector_t& initState, const std::vector<scalar_t>& switchingTimes)  {
		Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");

		gslqpSolver_.reset();
		gslqpSolver_.options().maxIterationGSLQP_ = 30;
		gslqpSolver_.options().warmStartGSLQP_ = true;

		gslqpSolver_.run(initState, switchingTimes);
		gslqpCostDerivative_ = gslqpSolver_.costDerivative();
		std::cout << "cost0: " << gslqpSolver_.cost() << std::endl;

		// finite difference derivative
		gslqpSolver_.options().maxIterationGSLQP_ = 10;

		for (size_t i=0; i<NUM_Subsystems-1; i++) {

			double h = eps_ * std::max(fabs(switchingTimes[i+1]), 1.0);

			std::vector<scalar_t> switchingTimesPlusPerturbed(switchingTimes);
			switchingTimesPlusPerturbed[i+1] += h;
			std::cout << "x+: " << (Eigen::VectorXd::Map(switchingTimesPlusPerturbed.data(),NUM_Subsystems+1)).transpose().format(CleanFmt) << std::endl;
			//
			gslqpSolver_.run(initState, switchingTimesPlusPerturbed);
			double costPlusPerturbed = gslqpSolver_.cost();
			std::cout << "cost+: " << costPlusPerturbed << std::endl;

			std::vector<scalar_t> switchingTimesMinusPerturbed(switchingTimes);
			switchingTimesMinusPerturbed[i+1] -= h;
			std::cout << "x-: " << (Eigen::VectorXd::Map(switchingTimesMinusPerturbed.data(),NUM_Subsystems+1)).transpose().format(CleanFmt) << std::endl;
			//
			gslqpSolver_.run(initState, switchingTimesMinusPerturbed);
			double costMinusPerturbed = gslqpSolver_.cost();
			std::cout << "cost+: " << costMinusPerturbed << std::endl;

			fdCostDerivative_(i) = (costPlusPerturbed-costMinusPerturbed) / (2.0*h);
		}

		parameters_t delta = gslqpCostDerivative_-fdCostDerivative_;
		Eigen::VectorXd eigenSwitchingTimes = Eigen::VectorXd::Map(switchingTimes.data(),NUM_Subsystems+1);

		std::cout << "Switching time: \t" << eigenSwitchingTimes.transpose().format(CleanFmt) << std::endl;
		std::cout << "GSLQP derivatives: \t" << gslqpCostDerivative_.transpose().format(CleanFmt) << std::endl;
		std::cout << "FD derivatives: \t" << fdCostDerivative_.transpose().format(CleanFmt) << std::endl;
		std::cout << "Difference: \t\t" << delta.transpose().format(CleanFmt) << std::endl << std::endl;

	}




private:
	const double eps_= sqrt(1e-2);

	parameters_t gslqpCostDerivative_;
	parameters_t fdCostDerivative_;

	GSLQPSolver<STATE_DIM, INPUT_DIM, NUM_Subsystems> gslqpSolver_;

};


#endif /* OCS2DERIVATIVESCHECK_H_ */
