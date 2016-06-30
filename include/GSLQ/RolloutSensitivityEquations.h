/*
 * RolloutSensitivityEquations.h
 *
 *  Created on: Jan 9, 2016
 *      Author: farbod
 */

#ifndef ROLLOUTSENSITIVITYEQUATIONS_OCS2_H_
#define ROLLOUTSENSITIVITYEQUATIONS_OCS2_H_

#include "Dimensions.h"

#include "dynamics/ControlledSystemBase.h"

#include "misc/LinearInterpolation.h"


namespace ocs2{

template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
class RolloutSensitivityEquations : public SystemBase<(NUM_SUBSYSTEMS-1)*OUTPUT_DIM>
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	typedef Dimensions<STATE_DIM, INPUT_DIM, OUTPUT_DIM> DIMENSIONS;
	typedef typename DIMENSIONS::template LinearFunction_t<INPUT_DIM, NUM_SUBSYSTEMS-1> sensitivity_controller_t;
	typedef typename DIMENSIONS::scalar_t 		scalar_t;
	typedef typename DIMENSIONS::scalar_array_t scalar_array_t;
	typedef typename DIMENSIONS::state_vector_t 	  state_vector_t;
	typedef typename DIMENSIONS::state_vector_array_t state_vector_array_t;
	typedef typename DIMENSIONS::control_vector_t 		control_vector_t;
	typedef typename DIMENSIONS::control_vector_array_t control_vector_array_t;
	typedef typename DIMENSIONS::output_vector_t 	   output_vector_t;
	typedef typename DIMENSIONS::output_vector_array_t output_vector_array_t;
	typedef typename DIMENSIONS::control_feedback_t 	  control_feedback_t;
	typedef typename DIMENSIONS::state_matrix_t 	  state_matrix_t;
	typedef typename DIMENSIONS::state_matrix_array_t state_matrix_array_t;
	typedef typename DIMENSIONS::control_gain_matrix_t 		 control_gain_matrix_t;
	typedef typename DIMENSIONS::control_gain_matrix_array_t control_gain_matrix_array_t;

	typedef Eigen::Matrix<double,(NUM_SUBSYSTEMS-1)*OUTPUT_DIM,1> nabla_output_vector_t;
	typedef std::vector<nabla_output_vector_t, Eigen::aligned_allocator<nabla_output_vector_t> > nabla_output_vector_array_t;
	typedef Eigen::Matrix<double,OUTPUT_DIM,NUM_SUBSYSTEMS-1>     nabla_output_matrix_t;
	typedef std::vector<nabla_output_matrix_t, Eigen::aligned_allocator<nabla_output_matrix_t> > nabla_output_matrix_array_t;
	//
	typedef Eigen::Matrix<double,INPUT_DIM,NUM_SUBSYSTEMS-1> nabla_input_matrix_t;
	typedef std::vector<nabla_input_matrix_t, Eigen::aligned_allocator<nabla_input_matrix_t> > nabla_input_matrix_array_t;
	//
	typedef Eigen::Matrix<double,1,NUM_SUBSYSTEMS-1> nabla_scalar_rowvector_t;
	typedef std::vector<nabla_scalar_rowvector_t, Eigen::aligned_allocator<nabla_scalar_rowvector_t> > nabla_scalar_rowvector_array_t;

	RolloutSensitivityEquations()  {}
	~RolloutSensitivityEquations() {}

	static void convert2Vector(const nabla_output_matrix_t& nabla_Ym, nabla_output_vector_t& nabla_Yv)  {

		nabla_Yv = Eigen::Map<const nabla_output_vector_t>(nabla_Ym.data());
	}

	static void convert2Matrix(const nabla_output_vector_t& nabla_Yv, nabla_output_matrix_t& nabla_Ym)  {

		nabla_Ym = Eigen::Map<const nabla_output_matrix_t>(nabla_Yv.data());
	}

	void setData(const size_t& activeSubsystem, const scalar_array_t& switchingTimes, const sensitivity_controller_t* sensitivityControllerPtr,
			const scalar_array_t* timeTrajectoryPtr, const output_vector_array_t*  outputTimeDerivativeTrajectoryPtr,
			const state_matrix_array_t* AmTrajectoryPtr, const control_gain_matrix_array_t* BmTrajectoryPtr)  {

		activeSubsystem_ = activeSubsystem;
		switchingTimes_ = switchingTimes;

		KmFunc_.setTimeStamp(&(sensitivityControllerPtr->time_));
		KmFunc_.setData(&(sensitivityControllerPtr->k_));

		LvFunc_.setTimeStamp(&(sensitivityControllerPtr->time_));
		LvFunc_.setData(&(sensitivityControllerPtr->uff_));

		outputTimeDevFunc_.setTimeStamp(timeTrajectoryPtr);
		outputTimeDevFunc_.setData(outputTimeDerivativeTrajectoryPtr);

		AmFunc_.setTimeStamp(timeTrajectoryPtr);
		AmFunc_.setData(AmTrajectoryPtr);
		BmFunc_.setTimeStamp(timeTrajectoryPtr);
		BmFunc_.setData(BmTrajectoryPtr);
	}

	void computeDerivative(const scalar_t& z, const nabla_output_vector_t& nabla_Yv, nabla_output_vector_t& derivatives) {

		// denormalized time
		scalar_t t = switchingTimes_[activeSubsystem_] + z*(switchingTimes_[activeSubsystem_+1]-switchingTimes_[activeSubsystem_]);

		nabla_output_matrix_t nabla_Ym;
		convert2Matrix(nabla_Yv, nabla_Ym);

		state_matrix_t Am;
		AmFunc_.interpolate(t, Am);
		size_t greatestLessTimeStampIndex = AmFunc_.getGreatestLessTimeStampIndex();
		control_gain_matrix_t Bm;
		BmFunc_.interpolate(t, Bm, greatestLessTimeStampIndex);
		output_vector_t dydt;
		outputTimeDevFunc_.interpolate(t, dydt, greatestLessTimeStampIndex);

		// compute input sensitivity
		nabla_input_matrix_t nabla_Um;
		computeInputSensitivity(t, nabla_Ym, nabla_Um);

		nabla_output_matrix_t nabla_dXmdz;
		nabla_dXmdz = (switchingTimes_[activeSubsystem_+1]-switchingTimes_[activeSubsystem_])*(Am*nabla_Ym+Bm*nabla_Um);

		for (size_t j=0; j<NUM_SUBSYSTEMS-1; j++)  {
			if (j==activeSubsystem_)
				nabla_dXmdz.col(j) += dydt;
			if (j==activeSubsystem_-1)
				nabla_dXmdz.col(j) -= dydt;
		}

		convert2Vector(nabla_dXmdz, derivatives);
	}


	void computeInputSensitivity(const scalar_t& t, const nabla_output_matrix_t& nabla_Ym, nabla_input_matrix_t& nabla_Um) {

		control_feedback_t Km;
		KmFunc_.interpolate(t, Km);
		size_t greatestLessTimeStampIndex = KmFunc_.getGreatestLessTimeStampIndex();

		nabla_input_matrix_t Lv;
		LvFunc_.interpolate(t, Lv, greatestLessTimeStampIndex);

		nabla_Um = Km*nabla_Ym + Lv;
	}


private:
	size_t activeSubsystem_;
	scalar_array_t switchingTimes_;

	LinearInterpolation<control_feedback_t,Eigen::aligned_allocator<control_feedback_t> > KmFunc_;
	LinearInterpolation<nabla_input_matrix_t,Eigen::aligned_allocator<nabla_input_matrix_t> > LvFunc_;

	LinearInterpolation<output_vector_t,Eigen::aligned_allocator<output_vector_t> > outputTimeDevFunc_;

	LinearInterpolation<state_matrix_t,Eigen::aligned_allocator<state_matrix_t> > AmFunc_;
	LinearInterpolation<control_gain_matrix_t,Eigen::aligned_allocator<control_gain_matrix_t> > BmFunc_;

};

} // namespace ocs2

#endif /* ROLLOUTSENSITIVITYEQUATIONS_OCS2_H_ */

