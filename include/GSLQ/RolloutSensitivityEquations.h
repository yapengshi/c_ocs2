/*
 * RolloutSensitivityEquations.h
 *
 *  Created on: Jan 9, 2016
 *      Author: farbod
 */

#ifndef ROLLOUTSENSITIVITYEQUATIONS_H_
#define ROLLOUTSENSITIVITYEQUATIONS_H_

#include <functional>

#include "Dimensions.h"

#include "dynamics/ControlledSystemBase.h"

#include "misc/LinearInterpolation.h"

template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
class RolloutSensitivityEquations : public SystemBase<(NUM_Subsystems-1)*STATE_DIM>
{
public:
	typedef Eigen::Matrix<double,STATE_DIM,NUM_Subsystems-1>     nabla_state_matrix_t;
	typedef Eigen::Matrix<double,(NUM_Subsystems-1)*STATE_DIM,1> nabla_state_vector_t;
	typedef Eigen::Matrix<double,INPUT_DIM,NUM_Subsystems-1>     nabla_input_matrix_t;
	typedef Dimensions<STATE_DIM, INPUT_DIM> DIMENSIONS;
	typedef typename DIMENSIONS::controller_t controller_t;
	typedef typename DIMENSIONS::scalar_t 		scalar_t;
	typedef typename DIMENSIONS::scalar_array_t scalar_array_t;
	typedef typename DIMENSIONS::eigen_scalar_t       eigen_scalar_t;
	typedef typename DIMENSIONS::eigen_scalar_array_t eigen_scalar_array_t;
	typedef typename DIMENSIONS::state_vector_t 	  state_vector_t;
	typedef typename DIMENSIONS::state_vector_array_t state_vector_array_t;
	typedef typename DIMENSIONS::control_vector_t 		control_vector_t;
	typedef typename DIMENSIONS::control_vector_array_t control_vector_array_t;
	typedef typename DIMENSIONS::control_feedback_t 	  control_feedback_t;
	typedef typename DIMENSIONS::control_feedback_array_t control_feedback_array_t;
	typedef typename DIMENSIONS::state_matrix_t 	  state_matrix_t;
	typedef typename DIMENSIONS::state_matrix_array_t state_matrix_array_t;
	typedef typename DIMENSIONS::control_matrix_t 		control_matrix_t;
	typedef typename DIMENSIONS::control_matrix_array_t control_matrix_array_t;
	typedef typename DIMENSIONS::control_gain_matrix_t 		 control_gain_matrix_t;
	typedef typename DIMENSIONS::control_gain_matrix_array_t control_gain_matrix_array_t;

	RolloutSensitivityEquations()  {}
	~RolloutSensitivityEquations() {}

	static void convert2Vector(const nabla_state_matrix_t& nabla_Xm, nabla_state_vector_t& nabla_Xv)  {

		nabla_Xv = Eigen::Map<const nabla_state_vector_t>(nabla_Xm.data());
	}

	static void convert2Matrix(const nabla_state_vector_t& nabla_Xv, nabla_state_matrix_t& nabla_Xm)  {

		nabla_Xm = Eigen::Map<const nabla_state_matrix_t>(nabla_Xv.data());
	}

	void setData(const size_t& activeSubsystem, const scalar_array_t& switchingTimes,
			const std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM> >& subsystemDynamicsPtr, controller_t* const controllerPtr,
			scalar_array_t* const timeTrajectoryPtr, state_vector_array_t* const stateTrajectoryPtr, control_vector_array_t* const inputTrajectoryPtr,
			state_matrix_array_t* const AmTrajectoryPtr, control_gain_matrix_array_t* const BmTrajectoryPtr)  {

		activeSubsystem_ = activeSubsystem;
		switchingTimes_ = switchingTimes;

		systemFunction_ = [subsystemDynamicsPtr](const scalar_t& t, const state_vector_t& x, const control_vector_t& u, state_vector_t& dxdt) {
			subsystemDynamicsPtr->computeDerivative(t, x, u, dxdt); };

		KmFunc_.setTimeStamp(&(controllerPtr->time_));
		KmFunc_.setData(&(controllerPtr->k_));

		stateFunc_.setTimeStamp(timeTrajectoryPtr);
		stateFunc_.setData(stateTrajectoryPtr);
		inputFunc_.setTimeStamp(timeTrajectoryPtr);
		inputFunc_.setData(inputTrajectoryPtr);

		AmFunc_.setTimeStamp(timeTrajectoryPtr);
		AmFunc_.setData(AmTrajectoryPtr);
		BmFunc_.setTimeStamp(timeTrajectoryPtr);
		BmFunc_.setData(BmTrajectoryPtr);
	}

	void computeDerivative(const scalar_t& z, const nabla_state_vector_t& nabla_Xv, nabla_state_vector_t& derivatives) {

		// denormalized time
		scalar_t t = switchingTimes_[activeSubsystem_] + (switchingTimes_[activeSubsystem_+1]-switchingTimes_[activeSubsystem_])*(z-activeSubsystem_);

		nabla_state_matrix_t nabla_Xm;
		convert2Matrix(nabla_Xv, nabla_Xm);

		state_vector_t x;
		stateFunc_.interpolate(t, x);
		control_vector_t u;
		inputFunc_.interpolate(t, u);

		state_matrix_t Am;
		AmFunc_.interpolate(t, Am);
		control_gain_matrix_t Bm;
		BmFunc_.interpolate(t, Bm);

		// compute input sensitivity
		nabla_input_matrix_t nabla_Um;
		computeInputSensitivity(t, nabla_Xm, nabla_Um);

		state_vector_t dxdt;
		systemFunction_(t, x, u, dxdt);

		nabla_state_matrix_t nabla_dXmdz;
		nabla_dXmdz = (switchingTimes_[activeSubsystem_+1]-switchingTimes_[activeSubsystem_])*(Am*nabla_Xm+Bm*nabla_Um);
		for (size_t j=0; j<NUM_Subsystems-1; j++)  {
			if (j==activeSubsystem_)
				nabla_dXmdz.col(j) += dxdt;
			if (j==activeSubsystem_-1)
				nabla_dXmdz.col(j) -= dxdt;
		}

		convert2Vector(nabla_dXmdz, derivatives);
	}


	void computeInputSensitivity(const scalar_t& t, const nabla_state_matrix_t& nabla_Xm, nabla_input_matrix_t& nabla_Um) {

		control_feedback_t Km;
		KmFunc_.interpolate(t, Km);

		nabla_Um = Km*nabla_Xm;
	}


private:
	size_t activeSubsystem_;
	scalar_array_t switchingTimes_;

	std::function<void (const scalar_t& /*t*/, const state_vector_t& /*x*/, const control_vector_t& /*u*/, state_vector_t& /*dx*/)> systemFunction_;

	LinearInterpolation<control_feedback_t,Eigen::aligned_allocator<control_feedback_t> > KmFunc_;

	LinearInterpolation<state_vector_t,Eigen::aligned_allocator<state_vector_t> > stateFunc_;
	LinearInterpolation<control_vector_t,Eigen::aligned_allocator<control_vector_t> > inputFunc_;

	LinearInterpolation<state_matrix_t,Eigen::aligned_allocator<state_matrix_t> > AmFunc_;
	LinearInterpolation<control_gain_matrix_t,Eigen::aligned_allocator<control_gain_matrix_t> > BmFunc_;

};


#endif /* ROLLOUTSENSITIVITYEQUATIONS_H_ */
