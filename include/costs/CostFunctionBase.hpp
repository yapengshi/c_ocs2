/*
 * CostFunctionBase.hpp
 *
 *  Created on: 26.03.2014
 *      Author: neunertm
 */

#ifndef COSTFUNCTIONBASE_HPP_
#define COSTFUNCTIONBASE_HPP_

#include <ilqg/Dimensions.hpp>

template <size_t STATE_DIM, size_t CONTROL_DIM>
class CostFunctionBase
{
public:
	typedef Dimensions<STATE_DIM, CONTROL_DIM> DIMENSIONS;

	typedef typename DIMENSIONS::state_vector_t state_vector_t;
	typedef typename DIMENSIONS::state_matrix_t state_matrix_t;
	typedef typename DIMENSIONS::control_vector_t control_vector_t;
	typedef typename DIMENSIONS::control_matrix_t control_matrix_t;
	typedef typename DIMENSIONS::control_feedback_t control_feedback_t;
	typedef typename DIMENSIONS::scalar_t scalar_t;

	CostFunctionBase() {};
	virtual ~CostFunctionBase() {};

	virtual void setCurrentStateAndControl(const state_vector_t& x, const control_vector_t& u, const double& t) {
		x_ = x;
		u_ = u;
	}

	virtual void evaluate(scalar_t& cost) = 0;

	virtual void stateDerivative(state_vector_t& cost) = 0;
	virtual void stateSecondDerivative(state_matrix_t& cost) = 0;
	virtual void controlDerivative(control_vector_t& cost) = 0;
	virtual void controlSecondDerivative(control_matrix_t& cost) = 0;

	virtual void stateControlDerivative(control_feedback_t& cost) = 0;

	virtual void terminalCost(scalar_t& cost) = 0;
	virtual void terminalCostStateDerivative(state_vector_t& cost) = 0;
	virtual void terminalCostStateSecondDerivative(state_matrix_t& cost) = 0;

	virtual std::shared_ptr<CostFunctionBase<STATE_DIM, CONTROL_DIM> > clone() const {
		return new CostFunctionBase<STATE_DIM, CONTROL_DIM>(*this);
	}

protected:
	state_vector_t x_;
	control_vector_t u_;
};

#endif /* COSTFUNCTIONBASE_HPP_ */
