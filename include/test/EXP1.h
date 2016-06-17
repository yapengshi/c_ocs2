/*
 * EXP1.h
 *
 *  Created on: Jan 11, 2016
 *      Author: farbod
 */

#ifndef EXP1_H_
#define EXP1_H_

#include <cmath>

#include "GSLQ/GSLQP.h"

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
class EXP1_Sys1 : public ControlledSystemBase<2,1>
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	EXP1_Sys1() {}
	~EXP1_Sys1() {}

	void computeDerivative( const double& t, const Eigen::Vector2d& x, const Eigen::Matrix<double,1,1>& u, Eigen::Vector2d& dxdt)  {
		dxdt(0) = x(0) + u(0)*sin(x(0));
		dxdt(1) = -x(1) - u(0)*cos(x(1));
	}

	std::shared_ptr<ControlledSystemBase<2, 1> > clone() const { return std::make_shared<EXP1_Sys1>(*this); }
};

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
class EXP1_Sys2 : public ControlledSystemBase<2,1>
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	EXP1_Sys2() {}
	~EXP1_Sys2() {}

	void computeDerivative( const double& t, const Eigen::Vector2d& x, const Eigen::Matrix<double,1,1>& u, Eigen::Vector2d& dxdt)  {
		dxdt(0) = x(1) + u(0)*sin(x(1));
		dxdt(1) = -x(0) - u(0)*cos(x(0));
	}

	std::shared_ptr<ControlledSystemBase<2, 1> > clone() const { return std::make_shared<EXP1_Sys2>(*this); }
};

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
class EXP1_Sys3 : public ControlledSystemBase<2,1>
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	EXP1_Sys3() {}
	~EXP1_Sys3() {}

	void computeDerivative( const double& t, const Eigen::Vector2d& x, const Eigen::Matrix<double,1,1>& u, Eigen::Vector2d& dxdt)  {
		dxdt(0) = -x(0) - u(0)*sin(x(0));
		dxdt(1) = x(1) + u(0)*cos(x(1));
	}

	std::shared_ptr<ControlledSystemBase<2, 1> > clone() const { return std::make_shared<EXP1_Sys3>(*this); }
};

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
class EXP1_SysDerivative1 : public DerivativesBase<2,1>
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	EXP1_SysDerivative1() {};
	~EXP1_SysDerivative1() {};

	void getDerivativeState(state_matrix_t& A)  { A << u_(0)*cos(x_(0))+1, 0, 0, u_(0)*sin(x_(1))-1; }
	void getDerivativesControl(control_gain_matrix_t& B) { B << sin(x_(0)), -cos(x_(1)); }

	std::shared_ptr<DerivativesBase<2,1> > clone() const { return std::make_shared<EXP1_SysDerivative1>(*this); }
};


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
class EXP1_SysDerivative2 : public DerivativesBase<2,1>
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	EXP1_SysDerivative2() {};
	~EXP1_SysDerivative2() {};

	void getDerivativeState(state_matrix_t& A)  { A << 0, u_(0)*cos(x_(1))+1, u_(0)*sin(x_(0))-1, 0; }
	void getDerivativesControl(control_gain_matrix_t& B) { B << sin(x_(1)), -cos(x_(0)); }

	std::shared_ptr<DerivativesBase<2,1> > clone() const { return std::make_shared<EXP1_SysDerivative2>(*this); }
};


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
class EXP1_SysDerivative3 : public DerivativesBase<2,1>
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	EXP1_SysDerivative3() {};
	~EXP1_SysDerivative3() {};

	void getDerivativeState(state_matrix_t& A)  { A << -u_(0)*cos(x_(0))-1, 0, 0, 1-u_(0)*sin(x_(1)); }
	void getDerivativesControl(control_gain_matrix_t& B) { B << -sin(x_(0)), cos(x_(1)); }

	std::shared_ptr<DerivativesBase<2,1> > clone() const { return std::make_shared<EXP1_SysDerivative3>(*this); }
};


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
class EXP1_CostFunction1 : public CostFunctionBase<2,1>
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	EXP1_CostFunction1() {};
	~EXP1_CostFunction1() {};

	void evaluate(scalar_t& L) { L = 0.5*pow(x_(0)-1.0, 2) + 0.5*pow(x_(1)+1.0, 2) + 0.5*pow(u_(0), 2); }

	void stateDerivative(state_vector_t& dLdx) { dLdx << (x_(0)-1.0), (x_(1)+1.0); }
	void stateSecondDerivative(state_matrix_t& dLdxx)  { dLdxx << 1.0, 0.0, 0.0, 1.0; }
	void controlDerivative(control_vector_t& dLdu)  { dLdu << u_; }
	void controlSecondDerivative(control_matrix_t& dLduu)  { dLduu << 1.0; }

	void stateControlDerivative(control_feedback_t& dLdxu) { dLdxu.setZero(); }

	void terminalCost(scalar_t& Phi) { Phi = 0; }
	void terminalCostStateDerivative(state_vector_t& dPhidx)  { dPhidx.setZero(); }
	void terminalCostStateSecondDerivative(state_matrix_t& dPhidxx)  { dPhidxx.setZero(); }

	std::shared_ptr<CostFunctionBase<2,1> > clone() const { return std::make_shared<EXP1_CostFunction1>(*this); };
};


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
class EXP1_CostFunction2 : public CostFunctionBase<2,1>
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	EXP1_CostFunction2() {};
	~EXP1_CostFunction2() {};

	void evaluate(scalar_t& L) { L = 0.5*pow(x_(0)-1.0, 2) + 0.5*pow(x_(1)+1.0, 2) + 0.5*pow(u_(0), 2); }

	void stateDerivative(state_vector_t& dLdx) { dLdx << (x_(0)-1.0), (x_(1)+1.0); }
	void stateSecondDerivative(state_matrix_t& dLdxx)  { dLdxx << 1.0, 0.0, 0.0, 1.0; }
	void controlDerivative(control_vector_t& dLdu)  { dLdu << u_; }
	void controlSecondDerivative(control_matrix_t& dLduu)  { dLduu << 1.0; }

	void stateControlDerivative(control_feedback_t& dLdxu) { dLdxu.setZero(); }

	void terminalCost(scalar_t& Phi) { Phi = 0; }
	void terminalCostStateDerivative(state_vector_t& dPhidx)  { dPhidx.setZero(); }
	void terminalCostStateSecondDerivative(state_matrix_t& dPhidxx)  { dPhidxx.setZero(); }

	std::shared_ptr<CostFunctionBase<2,1> > clone() const { return std::make_shared<EXP1_CostFunction2>(*this); };

};

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
class EXP1_CostFunction3 : public CostFunctionBase<2,1>
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	EXP1_CostFunction3() {};
	~EXP1_CostFunction3() {};

	void evaluate(scalar_t& L) { L = 0.5*pow(x_(0)-1.0, 2) + 0.5*pow(x_(1)+1.0, 2) + 0.5*pow(u_(0), 2); }

	void stateDerivative(state_vector_t& dLdx) { dLdx << (x_(0)-1.0), (x_(1)+1.0); }
	void stateSecondDerivative(state_matrix_t& dLdxx)  { dLdxx << 1.0, 0.0, 0.0, 1.0; }
	void controlDerivative(control_vector_t& dLdu)  { dLdu << u_; }
	void controlSecondDerivative(control_matrix_t& dLduu)  { dLduu << 1.0; }

	void stateControlDerivative(control_feedback_t& dLdxu) { dLdxu.setZero(); }

	void terminalCost(scalar_t& Phi) { Phi = 0.5*pow(x_(0)-1.0, 2) + 0.5*pow(x_(1)+1.0, 2); }
	void terminalCostStateDerivative(state_vector_t& dPhidx)  { dPhidx << (x_(0)-1.0), (x_(1)+1.0); }
	void terminalCostStateSecondDerivative(state_matrix_t& dPhidxx)  { dPhidxx << 1.0, 0.0, 0.0, 1.0; }

	std::shared_ptr<CostFunctionBase<2,1> > clone() const { return std::make_shared<EXP1_CostFunction3>(*this); };

};


#endif /* EXP1_H_ */
