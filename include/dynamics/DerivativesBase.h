/*
 * DerivativesBase.h
 *
 *  Created on: Jan 3, 2016
 *      Author: farbod
 */

#ifndef DERIVATIVESBASE_H_
#define DERIVATIVESBASE_H_

#include <memory>
#include <cstring>

#include "Dimensions.h"

template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM=STATE_DIM>
class DerivativesBase
{
public:

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	typedef Dimensions<STATE_DIM, INPUT_DIM, OUTPUT_DIM> DIMENSIONS;
	typedef typename DIMENSIONS::scalar_t scalar_t;
	typedef typename DIMENSIONS::state_vector_t   state_vector_t;
	typedef typename DIMENSIONS::control_vector_t control_vector_t;
	typedef typename DIMENSIONS::output_vector_t  output_vector_t;
	typedef typename DIMENSIONS::state_matrix_t state_matrix_t;
	typedef typename DIMENSIONS::control_gain_matrix_t control_gain_matrix_t;
	typedef typename DIMENSIONS::constraint1_state_matrix_t   constraint1_state_matrix_t;
	typedef typename DIMENSIONS::constraint1_control_matrix_t constraint1_control_matrix_t;

	DerivativesBase() {}
	virtual ~DerivativesBase() {}

	virtual void setCurrentStateAndControl(const scalar_t& t, const state_vector_t& x,
			const control_vector_t& u, const output_vector_t& y) {

		setCurrentStateAndControl(t, x, u);
		y_ = y;
	}

	virtual void setCurrentStateAndControl(const scalar_t& t, const state_vector_t& x,
			const control_vector_t& u) {
		t_ = t;
		x_ = x;
		u_ = u;
		y_ = x_.template head<OUTPUT_DIM>();
	}

	virtual void initializeModel(const scalar_t& initTime, const state_vector_t& initState,
			const scalar_t& finalTime=0, const char* algorithmName=NULL)
	{}

	virtual void getDerivativeState(state_matrix_t& A) = 0;
	virtual void getDerivativesControl(control_gain_matrix_t& B) = 0;
	virtual void getConstraint1DerivativesState(size_t& numConstraint1, constraint1_state_matrix_t& C) {	numConstraint1 = 0; }
	virtual void getConstraint1DerivativesControl(size_t& numConstraint1, constraint1_control_matrix_t& D) { numConstraint1 = 0; }

	virtual std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > clone() const = 0;

protected:
	scalar_t t_;
	state_vector_t x_;
	control_vector_t u_;
	output_vector_t y_;
};

#endif /* DERIVATIVESBASE_H_ */
