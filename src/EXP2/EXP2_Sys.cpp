/*
 * EXP2_Sys.cpp
 *
 *  Created on: Dec 27, 2015
 *      Author: farbod
 */

#include "dynamics/SystemBase.h"


class EXP2_Sys1 : public SystemBase<2>
{
public:
	EXP2_Sys1() {}
	~EXP2_Sys1() {}

	void computeDerivative(
			const Eigen::Matrix<double,2,1>& x,
			const double& t,
			Eigen::Matrix<double,2,1>& dxdt)
	{
		Eigen::Matrix<double,1,1> u;
		computeInput(x, t, u);

		Eigen::Matrix2d A;
		A << 0.6, 1.2, -0.8, 3.4;

		Eigen::Vector2d B;
		B << 1, 1;

		dxdt = A*x + B*u;
	}

	void computeInput(
			const Eigen::Matrix<double,2,1>& x,
			const double& t,
			Eigen::Matrix<double,1,1>& u)
	{

	}

	void setUff()

	void setGain()

private:
	std::vector<double> controllerTime_;

};
