/*
 * CostFunctionBase.h
 *
 *  Created on: Jan 3, 2016
 *      Author: farbod
 */

#ifndef COSTFUNCTIONBASE_H_
#define COSTFUNCTIONBASE_H_

#include "Dimensions.h"

template <size_t STATE_DIM, size_t CONTROL_DIM>
class CostFunctionBase
{
public:
	typedef Dimensions<STATE_DIM, CONTROL_DIM> DIMENSIONS;
	typedef typename DIMENSIONS::scalar_t scalar_t;
	typedef typename DIMENSIONS::state_vector_t state_vector_t;
	typedef typename DIMENSIONS::state_matrix_t state_matrix_t;
	typedef typename DIMENSIONS::control_vector_t control_vector_t;
	typedef typename DIMENSIONS::control_matrix_t control_matrix_t;
	typedef typename DIMENSIONS::control_feedback_t control_feedback_t;

	CostFunctionBase() {};
	virtual ~CostFunctionBase() {};

	virtual void setCurrentStateAndControl(const scalar_t& t, const state_vector_t& x, const control_vector_t& u) {
		t_ = t;
		x_ = x;
		u_ = u;
	}

	virtual void evaluate(scalar_t& L) = 0;

	virtual void stateDerivative(state_vector_t& dLdx) = 0;
	virtual void stateSecondDerivative(state_matrix_t& dLdxx) = 0;
	virtual void controlDerivative(control_vector_t& dLdu) = 0;
	virtual void controlSecondDerivative(control_matrix_t& dLduu) = 0;

	virtual void stateControlDerivative(control_feedback_t& dLdxu) = 0;

	virtual void terminalCost(scalar_t& Phi) = 0;
	virtual void terminalCostStateDerivative(state_vector_t& dPhidx) = 0;
	virtual void terminalCostStateSecondDerivative(state_matrix_t& dPhidxx) = 0;

	virtual std::shared_ptr<CostFunctionBase<STATE_DIM, CONTROL_DIM> > clone() const = 0;

protected:
	scalar_t t_;
	state_vector_t x_;
	control_vector_t u_;
};



#endif /* COSTFUNCTIONBASE_H_ */
