/*
 * DerivativesBase.hpp
 *
 *  Created on: 26.03.2014
 *      Author: neunertm
 */

#ifndef DERIVATIVESBASE_HPP_
#define DERIVATIVESBASE_HPP_

#include <ilqg/Dimensions.hpp>

template <size_t STATE_DIM, size_t CONTROL_DIM>
class DerivativesBase
{
public:
	typedef Dimensions<STATE_DIM, CONTROL_DIM> DIMENSIONS;

	typedef typename DIMENSIONS::state_vector_t state_vector_t;
	typedef typename DIMENSIONS::state_matrix_t state_matrix_t;
	typedef typename DIMENSIONS::control_vector_t control_vector_t;
	typedef typename DIMENSIONS::control_gain_matrix_t control_gain_matrix_t;

	DerivativesBase() {};
	virtual ~DerivativesBase() {};

	virtual void setCurrentStateAndControl(const state_vector_t& x, const control_vector_t& u) {
		x_ = x;
		u_ = u;
	}

	// TODO: Return by reference?
	virtual state_matrix_t getDerivativeState() = 0;
	virtual control_gain_matrix_t getDerivativesControl() = 0;

	virtual std::shared_ptr<DerivativesBase<STATE_DIM, CONTROL_DIM> > clone() const {
		return new DerivativesBase<STATE_DIM, CONTROL_DIM>(*this);
	}

protected:
	state_vector_t x_;
	control_vector_t u_;
};


#endif /* DERIVATIVESBASE_HPP_ */
