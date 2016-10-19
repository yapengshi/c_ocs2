/*
 * SensitivitySequentialRiccatiEquations.h
 *
 *  Created on: Jun 21, 2016
 *      Author: farbod
 */

#ifndef SENSITIVITYSEQUENTIALRICCATIEQUATIONS_OCS2_H_
#define SENSITIVITYSEQUENTIALRICCATIEQUATIONS_OCS2_H_

#include <array>

#include "Dimensions.h"

#include "dynamics/SystemBase.h"

#include "misc/LinearInterpolation.h"


namespace ocs2{

template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_SUBSYSTEMS>
class SensitivitySequentialRiccatiEquations : public SystemBase<(NUM_SUBSYSTEMS-1)*(STATE_DIM*STATE_DIM+STATE_DIM+1)>
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	enum { S_DIM_ = STATE_DIM*STATE_DIM+STATE_DIM+1 };
	typedef Eigen::Matrix<double,S_DIM_,1> s_vector_t;
	typedef Eigen::Matrix<double,(NUM_SUBSYSTEMS-1)*S_DIM_,1> all_s_vector_t;
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

	typedef Eigen::Matrix<double,STATE_DIM,NUM_SUBSYSTEMS-1> nabla_state_matrix_t;
	typedef Eigen::Matrix<double,INPUT_DIM,NUM_SUBSYSTEMS-1> nabla_input_matrix_t;
	typedef Eigen::Matrix<double,1,NUM_SUBSYSTEMS-1> 		 nabla_scalar_rowvector_t;
	typedef std::vector<nabla_state_matrix_t, Eigen::aligned_allocator<nabla_state_matrix_t> > nabla_state_matrix_array_t;
	typedef std::vector<nabla_input_matrix_t, Eigen::aligned_allocator<nabla_input_matrix_t> > nabla_input_matrix_array_t;
	typedef std::vector<nabla_scalar_rowvector_t, Eigen::aligned_allocator<nabla_scalar_rowvector_t> > nabla_scalar_rowvector_array_t;
	typedef std::array<state_matrix_t, NUM_SUBSYSTEMS-1> nabla_Sm_t;
	typedef std::array<state_vector_t, NUM_SUBSYSTEMS-1> nabla_Sv_t;
	typedef std::array<eigen_scalar_t, NUM_SUBSYSTEMS-1> nabla_s_t;

	SensitivitySequentialRiccatiEquations() {}
	~SensitivitySequentialRiccatiEquations() {}

	static void convert2Vector(const nabla_Sm_t& nabla_Sm, const nabla_Sv_t& nabla_Sv, const nabla_s_t& nabla_s,
			all_s_vector_t& allSs)  {

		for (size_t j=0; j<NUM_SUBSYSTEMS-1; j++)
			allSs.template segment<S_DIM_>(j*S_DIM_) << Eigen::Map<const Eigen::VectorXd>(nabla_Sm[j].data(),STATE_DIM*STATE_DIM),
					Eigen::Map<const Eigen::VectorXd>(nabla_Sv[j].data(),STATE_DIM),
					nabla_s[j];
	}

	static void convert2Matrix(const all_s_vector_t& allSs,
			nabla_Sm_t& nabla_Sm, nabla_Sv_t& nabla_Sv, nabla_s_t& nabla_s)  {

		for (size_t j=0; j<NUM_SUBSYSTEMS-1; j++) {
			nabla_Sm.at(j) = Eigen::Map<const Eigen::MatrixXd>(allSs.data()+j*S_DIM_,STATE_DIM, STATE_DIM);
			nabla_Sv.at(j) = Eigen::Map<const Eigen::VectorXd>(allSs.data()+j*S_DIM_+STATE_DIM*STATE_DIM, STATE_DIM);
			nabla_s.at(j)  = Eigen::Map<const Eigen::VectorXd>(allSs.data()+j*S_DIM_+STATE_DIM*STATE_DIM+STATE_DIM, 1);
		}
	}

	void setData(const scalar_t& learningRate,
			const size_t& activeSubsystem, const scalar_t& switchingTimeStart, const scalar_t& switchingTimeFinal,
			const scalar_array_t* SsTimePtr, const state_matrix_array_t* SmPtr, const state_vector_array_t* SvPtr,
			const scalar_array_t* timeStampPtr,
			const state_matrix_array_t* AmPtr, const control_gain_matrix_array_t* BmPtr,
			const eigen_scalar_array_t* qPtr, const state_vector_array_t* QvPtr, const state_matrix_array_t* QmPtr,
			const control_vector_array_t* RvPtr, const control_matrix_array_t* RmInversePtr, const control_matrix_array_t* RmPtr,
			const control_feedback_array_t* PmPtr,
			const scalar_array_t* sensitivityTimeStampPtr, const nabla_scalar_rowvector_array_t* nablaqPtr,
			const nabla_state_matrix_array_t* nablaQvPtr, const nabla_input_matrix_array_t* nablaRvPtr)  {

		alpha_ = learningRate;

		activeSubsystem_ = activeSubsystem;
		switchingTimeStart_ = switchingTimeStart;
		switchingTimeFinal_ = switchingTimeFinal;

		SvFunc_.setTimeStamp(SsTimePtr);
		SvFunc_.setData(SvPtr);
		SmFunc_.setTimeStamp(SsTimePtr);
		SmFunc_.setData(SmPtr);

		AmFunc_.setTimeStamp(timeStampPtr);
		AmFunc_.setData(AmPtr);
		BmFunc_.setTimeStamp(timeStampPtr);
		BmFunc_.setData(BmPtr);

		qFunc_.setTimeStamp(timeStampPtr);
		qFunc_.setData(qPtr);
		QvFunc_.setTimeStamp(timeStampPtr);
		QvFunc_.setData(QvPtr);
		QmFunc_.setTimeStamp(timeStampPtr);
		QmFunc_.setData(QmPtr);
		RvFunc_.setTimeStamp(timeStampPtr);
		RvFunc_.setData(RvPtr);
		RmInverseFunc_.setTimeStamp(timeStampPtr);
		RmInverseFunc_.setData(RmInversePtr);
		RmFunc_.setTimeStamp(timeStampPtr);
		RmFunc_.setData(RmPtr);
		PmFunc_.setTimeStamp(timeStampPtr);
		PmFunc_.setData(PmPtr);

		nablaqFunc_.setTimeStamp(sensitivityTimeStampPtr);
		nablaqFunc_.setData(nablaqPtr);
		nablaQvFunc_.setTimeStamp(sensitivityTimeStampPtr);
		nablaQvFunc_.setData(nablaQvPtr);
		nablaRvFunc_.setTimeStamp(sensitivityTimeStampPtr);
		nablaRvFunc_.setData(nablaRvPtr);
	}

	void computeDerivative(const scalar_t& z, const all_s_vector_t& allSs, all_s_vector_t& derivatives)  {

		// denormalized time
		scalar_t t = switchingTimeFinal_ - (switchingTimeFinal_-switchingTimeStart_)*z;

		convert2Matrix(allSs, __nabla_Sm, __nabla_Sv, __nabla_s);

		size_t greatestLessTimeStampIndex;

		SvFunc_.interpolate(t, __Sv);
		greatestLessTimeStampIndex = SvFunc_.getGreatestLessTimeStampIndex();
		SmFunc_.interpolate(t, __Sm, greatestLessTimeStampIndex);
		AmFunc_.interpolate(t, __Am);
		greatestLessTimeStampIndex = AmFunc_.getGreatestLessTimeStampIndex();
		BmFunc_.interpolate(t, __Bm, greatestLessTimeStampIndex);
		qFunc_.interpolate(t, __q);
		greatestLessTimeStampIndex = qFunc_.getGreatestLessTimeStampIndex();
		QvFunc_.interpolate(t, __Qv, greatestLessTimeStampIndex);
		QmFunc_.interpolate(t, __Qm, greatestLessTimeStampIndex);
		RvFunc_.interpolate(t, __Rv, greatestLessTimeStampIndex);
		RmInverseFunc_.interpolate(t, __invRm, greatestLessTimeStampIndex);
		RmFunc_.interpolate(t, __Rm, greatestLessTimeStampIndex);
		PmFunc_.interpolate(t, __Pm, greatestLessTimeStampIndex);

		nablaqFunc_.interpolate(t, __nablaq);
		greatestLessTimeStampIndex = nablaqFunc_.getGreatestLessTimeStampIndex();
		nablaQvFunc_.interpolate(t, __nablaQv, greatestLessTimeStampIndex);
		nablaRvFunc_.interpolate(t, __nablaRv, greatestLessTimeStampIndex);

		// Riccati equations for the original system
		__Lm = __invRm*(__Pm+__Bm.transpose()*__Sm);
		__Lv = __invRm*(__Rv+__Bm.transpose()*__Sv);
		__dSmdt = __Qm + __Am.transpose()*__Sm + __Sm.transpose()*__Am - __Lm.transpose()*__Rm*__Lm;
		__dSmdt = 0.5*(__dSmdt+__dSmdt.transpose()).eval();
		__dSvdt = __Qv + __Am.transpose()*__Sv - __Lm.transpose()*__Rm*__Lv;
		__dsdt  = __q - 0.5*alpha_*(2.0-alpha_)*__Lv.transpose()*__Rm*__Lv;

		// derivatives of Riccati equations
		for (size_t j=0; j<NUM_SUBSYSTEMS-1; j++) {

			// switching time gradient for the original system
			__nabla_dSmdt = __Am.transpose()*__nabla_Sm[j] + __nabla_Sm[j].transpose()*__Am - __nabla_Sm[j].transpose()*__Bm*__invRm*__Rm*__Lm
					- __Lm.transpose()*__Rm*__invRm*__Bm.transpose()*__nabla_Sm[j];
			__nabla_dSmdt = 0.5*(__nabla_dSmdt+__nabla_dSmdt.transpose()).eval();
			__nabla_dSvdt = __nablaQv.col(j) + __Am.transpose()*__nabla_Sv[j] - __nabla_Sm[j].transpose()*__Bm*__invRm*__Rm*__Lv
					- __Lm.transpose()*__Rm*__invRm*(__nablaRv.col(j) + __Bm.transpose()*__nabla_Sv[j]);
			__nabla_dsdt  = __nablaq.col(j) - 0.5*alpha_*(2-alpha_)*(__nablaRv.col(j) + __Bm.transpose()*__nabla_Sv[j]).transpose()*__invRm*__Rm*__Lv
					- 0.5*alpha_*(2-alpha_)*__Lv.transpose()*__Rm*__invRm*(__nablaRv.col(j) + __Bm.transpose()*__nabla_Sv[j]);

			// switching time gradient for the equivalent system
			__nabla_dSmdz[j] = (switchingTimeFinal_-switchingTimeStart_)*__nabla_dSmdt;
			__nabla_dSvdz[j] = (switchingTimeFinal_-switchingTimeStart_)*__nabla_dSvdt;
			__nabla_dsdz[j]  = (switchingTimeFinal_-switchingTimeStart_)*__nabla_dsdt;

			if (j==activeSubsystem_)  {
				__nabla_dSmdz[j] += __dSmdt;
				__nabla_dSvdz[j] += __dSvdt;
				__nabla_dsdz[j]  += __dsdt;
			}
			if (j==activeSubsystem_-1) {
				__nabla_dSmdz[j] -= __dSmdt;
				__nabla_dSvdz[j] -= __dSvdt;
				__nabla_dsdz[j]  -= __dsdt;
			}

		}  // end of j loop

		convert2Vector(__nabla_dSmdz, __nabla_dSvdz, __nabla_dsdz, derivatives);
	}


private:
	scalar_t alpha_;

	size_t activeSubsystem_;
	scalar_t switchingTimeStart_;
	scalar_t switchingTimeFinal_;

	LinearInterpolation<state_vector_t,Eigen::aligned_allocator<state_vector_t> > SvFunc_;
	LinearInterpolation<state_matrix_t,Eigen::aligned_allocator<state_matrix_t> > SmFunc_;

	LinearInterpolation<state_matrix_t,Eigen::aligned_allocator<state_matrix_t> > AmFunc_;
	LinearInterpolation<control_gain_matrix_t,Eigen::aligned_allocator<control_gain_matrix_t> > BmFunc_;

	LinearInterpolation<eigen_scalar_t,Eigen::aligned_allocator<eigen_scalar_t> > qFunc_;
	LinearInterpolation<state_vector_t,Eigen::aligned_allocator<state_vector_t> > QvFunc_;
	LinearInterpolation<state_matrix_t,Eigen::aligned_allocator<state_matrix_t> > QmFunc_;
	LinearInterpolation<control_vector_t,Eigen::aligned_allocator<control_vector_t> > RvFunc_;
	LinearInterpolation<control_matrix_t,Eigen::aligned_allocator<control_matrix_t> > RmInverseFunc_;
	LinearInterpolation<control_matrix_t,Eigen::aligned_allocator<control_matrix_t> > RmFunc_;
	LinearInterpolation<control_feedback_t,Eigen::aligned_allocator<control_feedback_t> > PmFunc_;

	LinearInterpolation<nabla_scalar_rowvector_t,Eigen::aligned_allocator<nabla_scalar_rowvector_t> > nablaqFunc_;
	LinearInterpolation<nabla_state_matrix_t,Eigen::aligned_allocator<nabla_state_matrix_t> > nablaQvFunc_;
	LinearInterpolation<nabla_input_matrix_t,Eigen::aligned_allocator<nabla_input_matrix_t> > nablaRvFunc_;


	// members only used in computeDerivative()
	nabla_Sm_t __nabla_Sm;
	nabla_Sv_t __nabla_Sv;
	nabla_s_t  __nabla_s;
	state_vector_t __Sv;
	state_matrix_t __Sm;
	state_matrix_t __Am;
	control_gain_matrix_t __Bm;
	eigen_scalar_t __q;
	state_vector_t __Qv;
	state_matrix_t __Qm;
	control_vector_t __Rv;
	control_matrix_t __invRm;
	control_matrix_t __Rm;
	control_feedback_t __Pm;
	nabla_scalar_rowvector_t __nablaq;
	nabla_state_matrix_t __nablaQv;
	nabla_input_matrix_t __nablaRv;
	control_feedback_t __Lm;
	control_vector_t __Lv;
	state_matrix_t __dSmdt;
	state_vector_t __dSvdt;
	eigen_scalar_t __dsdt;

	// derivatives of Riccati equations
	nabla_Sm_t 	__nabla_dSmdz;
	nabla_Sv_t 	__nabla_dSvdz;
	nabla_s_t 	__nabla_dsdz;

	// switching time gradient for the original system
	state_matrix_t __nabla_dSmdt;
	state_vector_t __nabla_dSvdt;
	eigen_scalar_t __nabla_dsdt;

};

} // namespace ocs2



#endif /* SENSITIVITYSEQUENTIALRICCATIEQUATIONS_OCS2_H_ */
