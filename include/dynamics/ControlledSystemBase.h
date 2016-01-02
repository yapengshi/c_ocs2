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

#include "SystemBase.h"
#include "misc/LinearInterpolation.h"

template <size_t STATE_DIM, size_t INPUT_DIM>
class ControlledSystemBase : public SystemBase<STATE_DIM>
{
private:
	typedef Eigen::Matrix<double,STATE_DIM,1> STATE_T;
	typedef Eigen::Matrix<double,INPUT_DIM,1> INPUT_T;
	typedef Eigen::Matrix<double,INPUT_DIM,STATE_DIM> GAIN_T;

public:
	ControlledSystemBase() {}
	virtual ~ControlledSystemBase() {}

	void setController(const std::vector<double>& controllerTime,
			const std::vector<INPUT_T, Eigen::aligned_allocator<INPUT_T> >& uff,
			const std::vector<GAIN_T, Eigen::aligned_allocator<GAIN_T> >& k) {

		controllerTime_ = controllerTime;
		uff_ = uff;
		k_ = k;

		linInterpolateUff_.setTimeStamp(&controllerTime_);
		linInterpolateUff_.setData(&uff_);

		linInterpolateK_.setTimeStamp(&controllerTime_);
		linInterpolateK_.setData(&k_);
	}


	void computeInput(const double& t, const STATE_T& x, INPUT_T& u)
	{
		Eigen::Matrix<double,INPUT_DIM,1> uff;
		Eigen::Matrix<double,INPUT_DIM,STATE_DIM> k;

		linInterpolateUff_.interpolate(t, uff);
		linInterpolateK_.interpolate(t, k);

		u = uff + k*x;
	}


	void computeDerivative(const double& t, const STATE_T& x, STATE_T& dxdt)  {

		INPUT_T u;
		computeInput(t, x, u);
		computeDerivative(t, x, u, dxdt);
	}

	virtual void computeDerivative(
			const double& t,
			const STATE_T& x,
			const INPUT_T& u,
			STATE_T& dxdt) = 0;

private:
	std::vector<double> controllerTime_;
	std::vector<INPUT_T, Eigen::aligned_allocator<INPUT_T> > uff_;
	std::vector<GAIN_T, Eigen::aligned_allocator<GAIN_T> > k_;

	LinearInterpolation<INPUT_T, Eigen::aligned_allocator<INPUT_T> > linInterpolateUff_;
	LinearInterpolation<GAIN_T, Eigen::aligned_allocator<GAIN_T> > linInterpolateK_;

};


#endif /* CONTROLLEDSYSTEMBASE_H_ */
