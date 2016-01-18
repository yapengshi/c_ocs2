/*
 * GSLQPSolver.h
 *
 *  Created on: Jan 18, 2016
 *      Author: farbod
 */

#ifndef GSLQPSOLVER_H_
#define GSLQPSOLVER_H_

#include "GSLQ/GLQP.h"
#include "GSLQ/GSLQP.h"

template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
class GSLQPSolver
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


	GSLQPSolver(const std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM> > >& subsystemDynamicsPtr,
			const std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM> > >& subsystemDerivativesPtr,
			const std::vector<std::shared_ptr<CostFunctionBase<STATE_DIM, INPUT_DIM> > >& subsystemCostFunctionsPtr,
			const state_vector_array_t&   stateOperatingPoints,
			const control_vector_array_t& inputOperatingPoints,
			const std::vector<size_t>& systemStockIndex,
			const Options_t& options = Options_t::Options() )

	: subsystemDynamicsPtr_(subsystemDynamicsPtr),
	  subsystemDerivativesPtr_(subsystemDerivativesPtr),
	  subsystemCostFunctionsPtr_(subsystemCostFunctionsPtr),
	  stateOperatingPoints_(stateOperatingPoints),
	  inputOperatingPoints_(inputOperatingPoints),
	  systemStockIndex_(systemStockIndex),
	  options_(options)
	{}

	virtual ~GSLQPSolver() {}

	void run(const state_vector_t& initState, const std::vector<scalar_t>& switchingTimes);

	scalar_t cost()  { return cost_; }

	parameters_t costDerivative()  {return costDerivative_;}


protected:
	size_t findNearestController(const Eigen::Matrix<double, NUM_Subsystems-1, 1>& parameters) const;


private:
	std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM> > > subsystemDynamicsPtr_;
	std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM> > > subsystemDerivativesPtr_;
	std::vector<std::shared_ptr<CostFunctionBase<STATE_DIM, INPUT_DIM> > > subsystemCostFunctionsPtr_;

	state_vector_array_t   stateOperatingPoints_;
	control_vector_array_t inputOperatingPoints_;

	std::vector<size_t> systemStockIndex_;

	Options_t options_;

	std::vector<parameters_t, Eigen::aligned_allocator<parameters_t> > parametersBag_;
	std::vector<std::vector<controller_t> > controllersStockBag_;

	scalar_t cost_;
	parameters_t costDerivative_;

};

#include "implementation/GSLQPSolver.h"

#endif /* GSLQPSOLVER_H_ */
