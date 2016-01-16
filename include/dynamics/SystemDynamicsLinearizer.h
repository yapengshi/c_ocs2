/*
 * SystemDynamicsLinearizer.h
 *
 *  Created on: 09.11.2015
 *      Author: mgiftthaler
 *
 *  Description:
 *  Within this class, we can distinguish between Linearizing first order systems or second order systems.
 *  For second order systems, we set the upper part of the derivative matrices analytically to the correct values.
 */

#ifndef SYSTEMDYNAMICSLINEARIZER_H_
#define SYSTEMDYNAMICSLINEARIZER_H_

#include <cmath>
#include <algorithm>
#include <Eigen/Dense>

#include "dynamics/ControlledSystemBase.h"
#include "dynamics/DerivativesBase.h"


template <size_t STATE_DIM, size_t INPUT_DIM>
class SystemDynamicsLinearizer : public DerivativesBase<STATE_DIM, INPUT_DIM>
{
public:
	typedef DerivativesBase<STATE_DIM, INPUT_DIM> Base;
	typedef typename Base::scalar_t scalar_t;
	typedef typename Base::state_vector_t state_vector_t;
	typedef typename Base::state_matrix_t state_matrix_t;
	typedef typename Base::control_vector_t      control_vector_t;
	typedef typename Base::control_gain_matrix_t control_gain_matrix_t;


	SystemDynamicsLinearizer(const std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM> >& nonlinearSystem,
			bool doubleSidedDerivative=true, bool isSecondOrderSystem=false)

		: nonlinearSystem_(nonlinearSystem),
		  doubleSidedDerivative_(doubleSidedDerivative),
		  isSecondOrderSystem_(isSecondOrderSystem)
	{}

	SystemDynamicsLinearizer(const SystemDynamicsLinearizer& other)

		: nonlinearSystem_(other.nonlinearSystem_->clone()),
		  doubleSidedDerivative_(other.doubleSidedDerivative_),
		  isSecondOrderSystem_(other.isSecondOrderSystem_)
	{}

	SystemDynamicsLinearizer& operator=(const SystemDynamicsLinearizer&other)  {
		nonlinearSystem_ = other.nonlinearSystem_->clone();
		doubleSidedDerivative_ = other.doubleSidedDerivative_;
		isSecondOrderSystem_ = other.isSecondOrderSystem_;
	}

	virtual ~SystemDynamicsLinearizer(){}


	virtual void setCurrentStateAndControl(const scalar_t& t, const state_vector_t& x, const control_vector_t& u) override  {

		Base::setCurrentStateAndControl(t, x, u);

		if (doubleSidedDerivative_==false)
			nonlinearSystem_->computeDerivative(Base::t_, Base::x_, Base::u_, f_);
	}

	virtual void getDerivativeState(state_matrix_t& A) override {

		for (size_t i=0; i<STATE_DIM; i++)  {

			// inspired from http://en.wikipedia.org/wiki/Numerical_differentiation#Practical_considerations_using_floating_point_arithmetic
			double h = eps_ * std::max(fabs(Base::x_(i)), 1.0);

			state_vector_t xPlusPerturbed = Base::x_;
			xPlusPerturbed(i) += h;

			// get evaluation of f(x,u)
			state_vector_t fPlusPerturbed;
			nonlinearSystem_->computeDerivative(Base::t_, xPlusPerturbed, Base::u_, fPlusPerturbed);

			if (doubleSidedDerivative_)  {

				state_vector_t xMinusPerturbed = Base::x_;
				xMinusPerturbed(i) -= h;

				state_vector_t fMinusPerturbed;
				nonlinearSystem_->computeDerivative(Base::t_, xMinusPerturbed, Base::u_, fMinusPerturbed);

				if(isSecondOrderSystem_)  {
					A.template topLeftCorner<STATE_DIM/2, STATE_DIM/2>().setZero();
					A.template topRightCorner<STATE_DIM/2, STATE_DIM/2>().setIdentity();
					A.template block<STATE_DIM/2,1>(STATE_DIM/2,i) = (fPlusPerturbed.template tail<STATE_DIM/2>() - fMinusPerturbed.template tail<STATE_DIM/2>()) / (2.0*h);
				}
				else
					A.col(i) = (fPlusPerturbed - fMinusPerturbed) / (2.0*h);
			}
			else  {
				if(isSecondOrderSystem_)  {
					A.template topLeftCorner<STATE_DIM/2, STATE_DIM/2>().setZero();
					A.template topRightCorner<STATE_DIM/2, STATE_DIM/2>().setIdentity();
					A.template block<STATE_DIM/2,1>(STATE_DIM/2,i) = (fPlusPerturbed.template tail<STATE_DIM/2>() - f_.template tail<STATE_DIM/2>()) / h;
				}
				else
					A.col(i) = (fPlusPerturbed - f_) / h;
			}
		}  // end of i loop

	}


	virtual void getDerivativesControl(control_gain_matrix_t& B) override  {

		for (size_t i=0; i<INPUT_DIM; i++) {

//			std::cin.get();

			// inspired from http://en.wikipedia.org/wiki/Numerical_differentiation#Practical_considerations_using_floating_point_arithmetic
			double h = eps_ * std::max(fabs(Base::u_(i)), 1.0);

			control_vector_t uPlusPerturbed = Base::u_;
			uPlusPerturbed(i) += h;

			// get evaluation of f(x,u)
			state_vector_t fPlusPerturbed;
			nonlinearSystem_->computeDerivative(Base::t_, Base::x_, uPlusPerturbed, fPlusPerturbed);

			if (doubleSidedDerivative_)  {

				control_vector_t uMinusPerturbed = Base::u_;
				uMinusPerturbed(i) -= h;

				state_vector_t fMinusPerturbed;
				nonlinearSystem_->computeDerivative(Base::t_, Base::x_, uMinusPerturbed, fMinusPerturbed);

				if(isSecondOrderSystem_)  {
					B.template topRows<STATE_DIM/2>().setZero();
					B.template block<STATE_DIM/2,1>(STATE_DIM/2,i) = (fPlusPerturbed.template tail<STATE_DIM/2>() - fMinusPerturbed.template tail<STATE_DIM/2>()) / (2.0*h);
				}
				else {
					B.col(i) = (fPlusPerturbed - fMinusPerturbed) / (2.0*h);

//					std::cout << ">>>> The " << i << " element out of " << INPUT_DIM << std::endl;
//					std::cout << "u+ :" << uPlusPerturbed.transpose() << std::endl;
//					std::cout << "f+ :" << fPlusPerturbed.transpose() << std::endl;
//					std::cout << "u- :" << uMinusPerturbed.transpose() << std::endl;
//					std::cout << "f- :" << fMinusPerturbed.transpose() << std::endl << std::endl;
				}
			}
			else {
				if(isSecondOrderSystem_)  {
					B.template topRows<STATE_DIM/2>().setZero();
					B.template block<STATE_DIM/2,1>(STATE_DIM/2,i) = (fPlusPerturbed.template tail<STATE_DIM/2>() - f_.template tail<STATE_DIM/2>()) / h;
				}
				else
					B.col(i) = (fPlusPerturbed - f_) / h;
			}




		}  // end of i loop
	}


	std::shared_ptr<Base> clone() const  { return std::make_shared<SystemDynamicsLinearizer<STATE_DIM, INPUT_DIM> >(*this); }


private:
	const double eps_= sqrt(Eigen::NumTraits<double>::epsilon());

	std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM> > nonlinearSystem_;
	bool doubleSidedDerivative_;
	bool isSecondOrderSystem_;

	state_vector_t f_;

};



#endif /* SYSTEMDYNAMICSLINEARIZER_H_ */
