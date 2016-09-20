/*
 * EXP5.h
 *
 *  Created on: Apr 15, 2016
 *      Author: markus
 */

#ifndef EXP5_OCS2_H_
#define EXP5_OCS2_H_


#include <cmath>

#include "GSLQ/GSLQP.h"


namespace ocs2{

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
class EXP5_Sys1 : public ControlledSystemBase<2,2>
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	EXP5_Sys1() {}
	~EXP5_Sys1() {}

	void computeDerivative(const double& t, const Eigen::Vector2d& x, const Eigen::Vector2d& u, Eigen::Vector2d& dxdt)  {
		dxdt(0) =u(0);
		dxdt(1) = u(1);
	}

	void computeConstriant2(const double& t, const Eigen::Vector2d& x, size_t& numConstraint1, Eigen::Vector2d& g1)  override {
		numConstraint1 = 1;
		g1(0) = x(0)*x(0) + x(1) - 1;
	}

	std::shared_ptr<ControlledSystemBase<2, 2> > clone() const { return std::make_shared<EXP5_Sys1>(*this); }
};


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
class EXP5_SysDerivative1 : public DerivativesBase<2,2>
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	EXP5_SysDerivative1() {};
	~EXP5_SysDerivative1() {};

	void getDerivativeState(state_matrix_t& A)  {
		A << 0.0, 0.0, 0.0, 0.0;
	}
	void getDerivativesControl(control_gain_matrix_t& B) {
		B.setZero(); B(0,0) = 1.0; B(1,1)=1.0;
	}
	void getConstraint2DerivativesState(constraint2_state_matrix_t& C) {
		C.topRows<1>() << 2*x_(0), 1.0;
	}

	void getFinalConstraint2DerivativesState(constraint2_state_matrix_t& F) {
		F.topRows<1>() << 2*x_(0), 1.0;
	}

	std::shared_ptr<DerivativesBase<2,2> > clone() const { return std::make_shared<EXP5_SysDerivative1>(*this); }
};



/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
class EXP5_CostFunction1 : public CostFunctionBaseOCS2<2,2>
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	EXP5_CostFunction1() {};
	~EXP5_CostFunction1() {};

	void evaluate(scalar_t& L) {
		L = 0.5 * 0.01*pow(x_(0), 2) + 0.5 * 0.01*pow(x_(1), 2) +0.5 * 0.01*u_(0)*u_(0) +0.5 * 0.01*u_(1)*u_(1);}

	void stateDerivative(state_vector_t& dLdx) {
		dLdx << 0.01*x_(0), 0.01*x_(1);
	}

	void stateSecondDerivative(state_matrix_t& dLdxx)  { dLdxx << 0.01, 0.0, 0.0, 0.01; }
	void controlDerivative(control_vector_t& dLdu)  { dLdu << 0.01* u_(0), 0.01*u_(1); }
	void controlSecondDerivative(control_matrix_t& dLduu)  { dLduu << 0.01, 0.0, 0.0, 0.01; }

	void stateControlDerivative(control_feedback_t& dLdxu) { dLdxu.setZero(); }

	void terminalCost(scalar_t& Phi) {
		Phi = 0.5* (x_(1)-1.0) * (x_(1)-1.0); }

	void terminalCostStateDerivative(state_vector_t& dPhidx)  {
		dPhidx << 0.0, (x_(1)-1.0);
	}
	void terminalCostStateSecondDerivative(state_matrix_t& dPhidxx)  {
		dPhidxx.setZero(); dPhidxx(1,1) = 1.0;}

	std::shared_ptr<CostFunctionBaseOCS2<2,2> > clone() const { return std::make_shared<EXP5_CostFunction1>(*this); };

private:
	double alpha_ = 0.01;
};


} // namespace ocs2


#endif /* EXP5_H_ */
