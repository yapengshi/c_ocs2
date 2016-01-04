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
	typedef typename DIMENSIONS::scaler_t scaler_t;
	typedef typename DIMENSIONS::state_vector_t state_vector_t;
	typedef typename DIMENSIONS::state_matrix_t state_matrix_t;
	typedef typename DIMENSIONS::control_vector_t control_vector_t;
	typedef typename DIMENSIONS::control_matrix_t control_matrix_t;
	typedef typename DIMENSIONS::control_feedback_t control_feedback_t;
	typedef typename DIMENSIONS::scalar_t scalar_t;

	CostFunctionBase() {};
	virtual ~CostFunctionBase() {};

	virtual void setCurrentStateAndControl(const scaler_t& t, const state_vector_t& x, const control_vector_t& u) {
		t_ = t;
		x_ = x;
		u_ = u;
	}

	virtual void evaluate(scalar_t& J) = 0;

	virtual void stateDerivative(state_vector_t& dJdx) = 0;
	virtual void stateSecondDerivative(state_matrix_t& dJdxx) = 0;
	virtual void controlDerivative(control_vector_t& dJdu) = 0;
	virtual void controlSecondDerivative(control_matrix_t& dJduu) = 0;

	virtual void stateControlDerivative(control_feedback_t& dJdxu) = 0;

	virtual void terminalCost(scalar_t& Phi) = 0;
	virtual void terminalCostStateDerivative(state_vector_t& dPhidx) = 0;
	virtual void terminalCostStateSecondDerivative(state_matrix_t& dPhidxx) = 0;

	virtual std::shared_ptr<CostFunctionBase<STATE_DIM, CONTROL_DIM> > clone() const = 0;

protected:
	scaler_t t_;
	state_vector_t x_;
	control_vector_t u_;
};



#endif /* COSTFUNCTIONBASE_H_ */
