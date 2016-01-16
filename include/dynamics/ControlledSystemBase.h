/*
 * ControlledSystemBase.h
 *
 *  Created on: Dec 29, 2015
 *      Author: farbod
 */

#ifndef CONTROLLEDSYSTEMBASE_H_
#define CONTROLLEDSYSTEMBASE_H_

#include <vector>
#include <Eigen/Dense>
#include <Eigen/StdVector>

#include "Dimensions.h"
#include "SystemBase.h"
#include "misc/LinearInterpolation.h"

template <size_t STATE_DIM, size_t INPUT_DIM>
class ControlledSystemBase : public SystemBase<STATE_DIM>
{
public:
	typedef Dimensions<STATE_DIM, INPUT_DIM> DIMENSIONS;
	typedef typename DIMENSIONS::scalar_t scalar_t;
	typedef typename DIMENSIONS::scalar_array_t scalar_array_t;
	typedef typename DIMENSIONS::state_vector_t state_vector_t;
	typedef typename DIMENSIONS::control_vector_t control_vector_t;
	typedef typename DIMENSIONS::control_vector_array_t control_vector_array_t;
	typedef typename DIMENSIONS::control_feedback_t control_feedback_t;
	typedef typename DIMENSIONS::control_feedback_array_t control_feedback_array_t;
	typedef typename DIMENSIONS::controller_t controller_t;

	ControlledSystemBase()
		: modelUpdated_(true)
	{}
	virtual ~ControlledSystemBase() {}

	void setController(const controller_t& controller) {

		controller_ = controller;

		linInterpolateUff_.setTimeStamp(&controller_.time_);
		linInterpolateUff_.setData(&controller_.uff_);

		linInterpolateK_.setTimeStamp(&controller_.time_);
		linInterpolateK_.setData(&controller_.k_);

		modelUpdated_ = true;
	}

	void setController(const scalar_array_t& controllerTime,
			const control_vector_array_t& uff,
			const control_feedback_array_t& k) {

		controller_.time_ = controllerTime;
		controller_.uff_ = uff;
		controller_.k_ = k;

		setController(controller_);
	}


	void computeInput(const scalar_t& t, const state_vector_t& x, control_vector_t& u)
	{
		control_vector_t uff;
		control_feedback_t k;

		linInterpolateUff_.interpolate(t, uff);
		linInterpolateK_.interpolate(t, k);

		u = uff + k*x;
	}


	void computeDerivative(const scalar_t& t, const state_vector_t& x, state_vector_t& dxdt)  {

		control_vector_t u;
		computeInput(t, x, u);
		computeDerivative(t, x, u, dxdt);
	}

	virtual void initializeModel(const scalar_t& initTime, const state_vector_t& initState, const scalar_t& finalTime=0)
	{}

	virtual std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM> > clone() const = 0;

	virtual void computeDerivative(
			const scalar_t& t,
			const state_vector_t& x,
			const control_vector_t& u,
			state_vector_t& dxdt) = 0;

protected:
	controller_t controller_;

	bool modelUpdated_;

	LinearInterpolation<control_vector_t, Eigen::aligned_allocator<control_vector_t> > linInterpolateUff_;
	LinearInterpolation<control_feedback_t, Eigen::aligned_allocator<control_feedback_t> > linInterpolateK_;

};


#endif /* CONTROLLEDSYSTEMBASE_H_ */
