/*
 * SolveBVP.h
 *
 *  Created on: Jun 18, 2016
 *      Author: farbod
 */

#ifndef SOLVEBVP_OCS2_H_
#define SOLVEBVP_OCS2_H_

#include <vector>
#include <memory>
#include <Eigen/Dense>
#include <Eigen/StdVector>

#include "GSLQ/BVPEquations.h"
#include "integration/Integrator.h"
#include "misc/LinearInterpolation.h"


namespace ocs2{

template <size_t STATE_DIM, size_t INPUT_DIM>
class SolveBVP
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	typedef BVPEquations<STATE_DIM, INPUT_DIM> BVP;

	typedef typename BVP::full_ode_vector_t full_ode_vector_t;
	typedef typename BVP::full_ode_vector_array_t full_ode_vector_array_t;

	typedef typename BVP::eigen_scalar_t eigen_scalar_t;
	typedef typename BVP::state_vector_t state_vector_t;
	typedef typename BVP::input_vector_t input_vector_t;
	typedef typename BVP::state_state_matrix_t state_state_matrix_t;
	typedef typename BVP::input_input_matrix_t input_input_matrix_t;
	typedef typename BVP::input_state_matrix_t input_state_matrix_t;
	typedef typename BVP::state_input_matrix_t state_input_matrix_t;

	typedef typename BVP::double_array_t double_array_t;
	typedef typename BVP::eigen_scalar_array_t eigen_scalar_array_t;
	typedef typename BVP::state_vector_array_t state_vector_array_t;
	typedef typename BVP::input_vector_array_t input_vector_array_t;
	typedef typename BVP::state_state_matrix_array_t state_state_matrix_array_t;
	typedef typename BVP::input_input_matrix_array_t input_input_matrix_array_t;
	typedef typename BVP::input_state_matrix_array_t input_state_matrix_array_t;
	typedef typename BVP::state_input_matrix_array_t state_input_matrix_array_t;


	SolveBVP()
	: bvpEquationsPtr_(new BVP()),
	  bvpOdeSolver_(bvpEquationsPtr_)
	{}

	~SolveBVP()  {}


	void setData(const double_array_t* timeStampPtr,
			const state_state_matrix_array_t* AmPtr, const state_state_matrix_array_t* OmPtr, const state_input_matrix_array_t* BmPtr, const state_vector_array_t* GvPtr,
			const state_vector_array_t* QvPtr, const state_state_matrix_array_t* QmPtr, const input_state_matrix_array_t* PmPtr,
			const input_vector_array_t* RvPtr, const input_input_matrix_array_t* RmPtr, const input_input_matrix_array_t* RmInversePtr)  {

		bvpEquationsPtr_->setData(timeStampPtr,
				AmPtr, OmPtr, BmPtr, GvPtr,
				QvPtr, QmPtr, PmPtr,
				RvPtr, RmPtr, RmInversePtr);

		startTime_ = timeStampPtr->front();
		finalTime_ = timeStampPtr->back();
	}


	/*
	 *	adaptive time step
	 */
	void solve(const state_vector_t& QvFinal, const state_state_matrix_t& QmFinal,
			double_array_t& timeTrajectory, state_state_matrix_array_t& MmTrajectory, state_vector_array_t& SvTrajectory,
			const double& absTolODE=1e-9, const double& relTolODE=1e-6)  {

		// final value matrices to vector form
		full_ode_vector_t MSvFinal;
		BVP::convert2Vector(QmFinal, QvFinal, MSvFinal);

		// integrating the Riccati equations
		double_array_t normalizedTimeTrajectory;
		full_ode_vector_array_t MSvTrajectory;
		bvpOdeSolver_.integrate(MSvFinal, 0, 1, MSvTrajectory, normalizedTimeTrajectory,
				1e-3, absTolODE, relTolODE);

		// denormalizing time and constructing 'Mm' and 'Sv'
		int N = normalizedTimeTrajectory.size();
		timeTrajectory.resize(N);
		MmTrajectory.resize(N);
		SvTrajectory.resize(N);
		for (int k=0; k<N; k++) {
			BVP::convert2Matrix(MSvTrajectory[N-1-k], MmTrajectory[k], SvTrajectory[k]);
			timeTrajectory[k] = (startTime_-finalTime_)*(timeTrajectory[N-1-k]) + finalTime_;
		}  // end of k loop

		// testing the numerical stability of the Riccati equations
		for (int k=N-1; k>=0; k--) {
			try {
				if (MmTrajectory[k] != MmTrajectory[k])  throw std::runtime_error("Mm is unstable.");
				if (SvTrajectory[k] != SvTrajectory[k])  throw std::runtime_error("Sv is unstable.");
			}
			catch(const std::exception& error)
			{
				std::cerr << "what(): " << error.what() << " at time " << timeTrajectory[k] << " [sec]." << std::endl;
				for (int kp=k; kp<k+10; kp++)  {
					if (kp >= N) continue;
					std::cerr << "Mm[" << timeTrajectory[kp] << "]:\n"<< MmTrajectory[kp].transpose() << std::endl;
					std::cerr << "Sv[" << timeTrajectory[kp] << "]:\t"<< SvTrajectory[kp].transpose() << std::endl;
				}
				exit(1);
			}

		}  // end of k loop

	}


	/*
	 * given time step
	 */
	void solve(const double_array_t& timeTrajectory, const state_vector_t& QvFinal, const state_state_matrix_t& QmFinal,
			state_state_matrix_array_t& MmTrajectory, state_vector_array_t& SvTrajectory,
			const double& absTolODE=1e-9, const double& relTolODE=1e-6)  {

		if (fabs(timeTrajectory.front()-startTime_) > 1e-6)  throw std::runtime_error("The input time vector has different startTime than the setData time.");
		if (fabs(timeTrajectory.back()-finalTime_) > 1e-6)   throw std::runtime_error("The input time vector has different finalTime than the setData time.");

		// normalizing input time trajectory
		size_t N = timeTrajectory.size();
		double_array_t normalizedTimeTrajectory(N);
		for (size_t k=0; k<N; k++)
			normalizedTimeTrajectory[k] = (finalTime_-timeTrajectory[N-1-k])/(finalTime_-startTime_);

		// final value matrices to vector form
		full_ode_vector_t MSvFinal;
		BVP::convert2Vector(QmFinal, QvFinal, MSvFinal);

		// integrating the Riccati equations
		full_ode_vector_array_t MSvTrajectory;
		bvpOdeSolver_.integrate(MSvFinal, normalizedTimeTrajectory, MSvTrajectory,
				1e-3, absTolODE, relTolODE);

		// denormalizing time and constructing 'Mm' and 'Sv'
		MmTrajectory.resize(N);
		SvTrajectory.resize(N);
		for (int k=0; k<N; k++)
			BVP::convert2Matrix(MSvTrajectory[N-1-k], MmTrajectory[k], SvTrajectory[k]);

		// testing the numerical stability of the Riccati equations
		for (int k=N-1; k>=0; k--) {
			try {
				if (MmTrajectory[k] != MmTrajectory[k])  throw std::runtime_error("Mm is unstable.");
				if (SvTrajectory[k] != SvTrajectory[k])  throw std::runtime_error("Sv is unstable.");
			}
			catch(const std::exception& error)
			{
				std::cerr << "what(): " << error.what() << " at time " << timeTrajectory[k] << " [sec]." << std::endl;
				for (int kp=k; kp<k+10; kp++)  {
					if (kp >= N) continue;
					std::cerr << "Mm[" << timeTrajectory[kp] << "]:\n"<< MmTrajectory[kp].transpose() << std::endl;
					std::cerr << "Sv[" << timeTrajectory[kp] << "]:\t"<< SvTrajectory[kp].transpose() << std::endl;
				}
				exit(1);
			}

		}  // end of k loop

	}

private:
	std::shared_ptr<BVP> bvpEquationsPtr_;
	ODE45<BVP::Full_ODE_VECTOR_DIM> bvpOdeSolver_;

	double startTime_;
	double finalTime_;

};

}  // end of OCS2 name space



#endif /* SOLVEBVP_H_ */
