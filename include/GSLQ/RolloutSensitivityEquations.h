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

template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
class RolloutSensitivityEquations : public SystemBase<(NUM_SUBSYSTEMS-1)*OUTPUT_DIM>
{
public:
	typedef Dimensions<STATE_DIM, INPUT_DIM, OUTPUT_DIM> DIMENSIONS;
	typedef typename DIMENSIONS::controller_t controller_t;
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

	void setData(const size_t& activeSubsystem, const scalar_array_t& switchingTimes,
			const std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> >& subsystemDynamicsPtr, const controller_t* controllerPtr,
			const scalar_array_t* timeTrajectoryPtr, const state_vector_array_t* stateTrajectoryPtr, const control_vector_array_t* inputTrajectoryPtr,
			const state_matrix_array_t* AmTrajectoryPtr, const control_gain_matrix_array_t* BmTrajectoryPtr,
			const scalar_array_t* nablaLvTimeTrajectoryPtr=NULL, const nabla_input_matrix_array_t* nablaLvTrajectoryPtr=NULL)  {

		if (nablaLvTimeTrajectoryPtr==NULL && nablaLvTrajectoryPtr==NULL)
			nablaLvIsSet = false;
		else if (nablaLvTimeTrajectoryPtr!=NULL && nablaLvTrajectoryPtr!=NULL)
			nablaLvIsSet = true;
		else
			throw std::runtime_error("The pointers to nablaLv and its time stamp should be either set or ignored.");

		activeSubsystem_ = activeSubsystem;
		switchingTimes_ = switchingTimes;

		systemFunction_ = [&subsystemDynamicsPtr](const scalar_t& t, const state_vector_t& x, const control_vector_t& u, output_vector_t& dydt) {
			state_vector_t dxdt;
			subsystemDynamicsPtr->computeDerivative(t, x, u, dxdt);
			dydt = subsystemDynamicsPtr->computeOutputStateDerivative(t, x, u) * dxdt;
		};

		KmFunc_.setTimeStamp(&(controllerPtr->time_));
		KmFunc_.setData(&(controllerPtr->k_));

		LvFunc_.setTimeStamp(nablaLvTimeTrajectoryPtr);
		LvFunc_.setData(nablaLvTrajectoryPtr);

		stateFunc_.setTimeStamp(timeTrajectoryPtr);
		stateFunc_.setData(stateTrajectoryPtr);
		inputFunc_.setTimeStamp(timeTrajectoryPtr);
		inputFunc_.setData(inputTrajectoryPtr);

		AmFunc_.setTimeStamp(timeTrajectoryPtr);
		AmFunc_.setData(AmTrajectoryPtr);
		BmFunc_.setTimeStamp(timeTrajectoryPtr);
		BmFunc_.setData(BmTrajectoryPtr);
	}

	void computeDerivative(const scalar_t& z, const nabla_output_vector_t& nabla_Yv, nabla_output_vector_t& derivatives) {

		// denormalized time
		scalar_t t = switchingTimes_[activeSubsystem_] + (switchingTimes_[activeSubsystem_+1]-switchingTimes_[activeSubsystem_])*(z-activeSubsystem_);

		nabla_output_matrix_t nabla_Ym;
		convert2Matrix(nabla_Yv, nabla_Ym);

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
		computeInputSensitivity(t, nabla_Ym, nabla_Um);

		output_vector_t dydt;
		systemFunction_(t, x, u, dydt);

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

		nabla_input_matrix_t nabla_Lv;
		if (nablaLvIsSet==true)
			LvFunc_.interpolate(t, nabla_Lv);
		else
			nabla_Lv.setZero();

		nabla_Um = Km*nabla_Ym + nabla_Lv;
	}


private:
	size_t activeSubsystem_;
	scalar_array_t switchingTimes_;
	bool nablaLvIsSet;

	std::function<void (const scalar_t& /*t*/, const state_vector_t& /*x*/, const control_vector_t& /*u*/, output_vector_t& /*dy*/)> systemFunction_;

	LinearInterpolation<control_feedback_t,Eigen::aligned_allocator<control_feedback_t> > KmFunc_;
	LinearInterpolation<nabla_input_matrix_t,Eigen::aligned_allocator<nabla_input_matrix_t> > LvFunc_;

	LinearInterpolation<state_vector_t,Eigen::aligned_allocator<state_vector_t> > stateFunc_;
	LinearInterpolation<control_vector_t,Eigen::aligned_allocator<control_vector_t> > inputFunc_;

	LinearInterpolation<state_matrix_t,Eigen::aligned_allocator<state_matrix_t> > AmFunc_;
	LinearInterpolation<control_gain_matrix_t,Eigen::aligned_allocator<control_gain_matrix_t> > BmFunc_;

};


#endif /* ROLLOUTSENSITIVITYEQUATIONS_H_ */
