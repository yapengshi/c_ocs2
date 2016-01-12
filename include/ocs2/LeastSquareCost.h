/*
 * LeastSquareCost.h
 *
 *  Created on: Jan 11, 2016
 *      Author: farbod
 */

#ifndef LEASTSQUARECOST_H_
#define LEASTSQUARECOST_H_

#include <cmath>

#include <ceres/ceres.h>

#include "GSLQ/GLQP.h"
#include "GSLQ/GSLQP.h"

template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
class LeastSquareCost : public ceres::SizedCostFunction<1, NUM_Subsystems-1>
{
public:
	typedef GLQP<STATE_DIM, INPUT_DIM, NUM_Subsystems> GLQP_t;
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

	LeastSquareCost(const std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM> > >& subsystemDynamicsPtr,
			const std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM> > >& subsystemDerivativesPtr,
			const std::vector<std::shared_ptr<CostFunctionBase<STATE_DIM, INPUT_DIM> > >& subsystemCostFunctionsPtr,
			const state_vector_array_t&   stateOperatingPoints,
			const control_vector_array_t& inputOperatingPoints,
			const std::vector<size_t>& systemStockIndex,
			const double& startTime,
			const double& finalTime,
			const state_vector_t& initState,
			const typename GSLQP_t::Options& options = GSLQP_t::Options() )
		: subsystemDynamicsPtr_(subsystemDynamicsPtr),
		  subsystemDerivativesPtr_(subsystemDerivativesPtr),
		  subsystemCostFunctionsPtr_(subsystemCostFunctionsPtr),
		  stateOperatingPoints_(stateOperatingPoints),
		  inputOperatingPoints_(inputOperatingPoints),
		  systemStockIndex_(systemStockIndex),
		  startTime_(startTime),
		  finalTime_(finalTime),
		  initState_(initState),
		  options_(options)
	{}
	virtual ~LeastSquareCost() {}


	virtual bool Evaluate(double const* const* parameters, double* totalCost, double** totalCostJacobian) const  {

		scalar_array_t switchingTimes(NUM_Subsystems+1);
		switchingTimes.front() = startTime_;
		switchingTimes.back()  = finalTime_;
		for (size_t j=1; j<NUM_Subsystems; j++)
			switchingTimes[j] = parameters[0][j-1];

		// GLQP initialization
		GLQP_t glqp(subsystemDynamicsPtr_, subsystemDerivativesPtr_, subsystemCostFunctionsPtr_,
				stateOperatingPoints_, inputOperatingPoints_, systemStockIndex_);
		glqp.run(switchingTimes);

		// GLQP controller
		std::vector<controller_t> controllersStock(NUM_Subsystems);
		glqp.getController(controllersStock);

		// GSLQP
		GSLQP_t gslqp(subsystemDynamicsPtr_, subsystemDerivativesPtr_, subsystemCostFunctionsPtr_,
				controllersStock, systemStockIndex_, options_);
		gslqp.run(initState_, switchingTimes);

		// cost funtion
		gslqp.getValueFuntion(0.0, initState_, totalCost[0]);
		totalCost[0] = sqrt(totalCost[0]);

		// cost funtion jacobian
		Eigen::Matrix<double,NUM_Subsystems-1,1> costFuntionDerivative;
		gslqp.getCostFuntionDerivative(initState_, costFuntionDerivative);
		for (size_t j=0; j<NUM_Subsystems-1; j++)
			if (totalCostJacobian != NULL && totalCostJacobian[0] != NULL)
				totalCostJacobian[0][j] = costFuntionDerivative[j]/totalCost[0];

		return true;
	}


private:
	std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM> > > subsystemDynamicsPtr_;
	std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM> > > subsystemDerivativesPtr_;
	std::vector<std::shared_ptr<CostFunctionBase<STATE_DIM, INPUT_DIM> > > subsystemCostFunctionsPtr_;

	state_vector_array_t   stateOperatingPoints_;
	control_vector_array_t inputOperatingPoints_;

	std::vector<size_t> systemStockIndex_;

	double startTime_;
	double finalTime_;

	state_vector_t initState_;

	typename GSLQP_t::Options options_;

};


#endif /* LEASTSQUARECOST_H_ */
