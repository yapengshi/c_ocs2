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

#include "SystemBase.h"
#include "misc/LinearInterpolation.h"

template <size_t STATE_DIM, size_t INPUT_DIM>
class ControlledSystemBase : public SystemBase
{
public:
	ControlledSystemBase() {}
	virtual ~ControlledSystemBase() {}

	void setController(const std::vector<double>& controllerTime,
			const std::vector<Eigen::Matrix<double,INPUT_DIM,1> >& uff,
			const std::vector<Eigen::Matrix<double,INPUT_DIM,STATE_DIM> >& k) {

		controllerTime_ = controllerTime;
		uff_ = uff;
		k_ = k;

		linInterpolateUff_.setTimeStamp(&controllerTime_);
		linInterpolateUff_.setData(&uff_);

		linInterpolateK_.setTimeStamp(&controllerTime_);
		linInterpolateK_.setData(&k_);
	}


	void computeInput(const double& t,
			const Eigen::Matrix<double,STATE_DIM,1>& x,
			Eigen::Matrix<double,INPUT_DIM,1>& u)
	{
		Eigen::Matrix<double,INPUT_DIM,1> uff;
		Eigen::Matrix<double,INPUT_DIM,STATE_DIM> k;

		linInterpolateUff_.interpolate(t, uff);
		linInterpolateK_.interpolate(t, k);

		u = uff + k*x;
	}


	void computeDerivative(
			const Eigen::Matrix<double,STATE_DIM,1>& x,
			const double& t,
			Eigen::Matrix<double,INPUT_DIM,1>& dxdt)  {

		Eigen::Matrix<double,INPUT_DIM,1> u;
		computeInput(t, x, u);
		computeDerivative(t, x, u, dxdt);
	}

	virtual void computeDerivative(
			const double& t,
			const Eigen::Matrix<double,STATE_DIM,1>& x,
			const Eigen::Matrix<double,INPUT_DIM,1>& u,
			Eigen::Matrix<double,INPUT_DIM,1>& dxdt) = 0;

private:
	std::vector<double> controllerTime_;
	std::vector<Eigen::Matrix<double,INPUT_DIM,1> > uff_;
	std::vector<Eigen::Matrix<double,INPUT_DIM,STATE_DIM> > k_;

	LinearInterpolation<Eigen::Matrix<double,INPUT_DIM,1> > linInterpolateUff_;
	LinearInterpolation<Eigen::Matrix<double,INPUT_DIM,STATE_DIM> > linInterpolateK_;

};


#endif /* CONTROLLEDSYSTEMBASE_H_ */
