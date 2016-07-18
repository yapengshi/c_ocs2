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

		nabla_Sm_t nabla_Sm;
		nabla_Sv_t nabla_Sv;
		nabla_s_t  nabla_s;
		convert2Matrix(allSs, nabla_Sm, nabla_Sv, nabla_s);

		size_t greatestLessTimeStampIndex;

		state_vector_t Sv;
		SvFunc_.interpolate(t, Sv);
		greatestLessTimeStampIndex = SvFunc_.getGreatestLessTimeStampIndex();
		state_matrix_t Sm;
		SmFunc_.interpolate(t, Sm, greatestLessTimeStampIndex);

		state_matrix_t Am;
		AmFunc_.interpolate(t, Am);
		greatestLessTimeStampIndex = AmFunc_.getGreatestLessTimeStampIndex();
		control_gain_matrix_t Bm;
		BmFunc_.interpolate(t, Bm, greatestLessTimeStampIndex);

		eigen_scalar_t q;
		qFunc_.interpolate(t, q);
		greatestLessTimeStampIndex = qFunc_.getGreatestLessTimeStampIndex();
		state_vector_t Qv;
		QvFunc_.interpolate(t, Qv, greatestLessTimeStampIndex);
		state_matrix_t Qm;
		QmFunc_.interpolate(t, Qm, greatestLessTimeStampIndex);
		control_vector_t Rv;
		RvFunc_.interpolate(t, Rv, greatestLessTimeStampIndex);
		control_matrix_t inverseRm;
		RmInverseFunc_.interpolate(t, inverseRm, greatestLessTimeStampIndex);
		control_matrix_t Rm;
		RmFunc_.interpolate(t, Rm, greatestLessTimeStampIndex);
		control_feedback_t Pm;
		PmFunc_.interpolate(t, Pm, greatestLessTimeStampIndex);

		nabla_scalar_rowvector_t nablaq;
		nablaqFunc_.interpolate(t, nablaq);
		greatestLessTimeStampIndex = nablaqFunc_.getGreatestLessTimeStampIndex();
		nabla_state_matrix_t nablaQv;
		nablaQvFunc_.interpolate(t, nablaQv, greatestLessTimeStampIndex);
		nabla_input_matrix_t nablaRv;
		nablaRvFunc_.interpolate(t, nablaRv, greatestLessTimeStampIndex);

		// Riccati equations for the original system
		control_feedback_t Lm = inverseRm*(Pm+Bm.transpose()*Sm);
		control_vector_t   Lv = inverseRm*(Rv+Bm.transpose()*Sv);
		state_matrix_t dSmdt = Qm + Am.transpose()*Sm + Sm.transpose()*Am - Lm.transpose()*Rm*Lm;
		dSmdt = 0.5*(dSmdt+dSmdt.transpose()).eval();
		state_vector_t dSvdt = Qv + Am.transpose()*Sv - Lm.transpose()*Rm*Lv;
		eigen_scalar_t dsdt  = q - 0.5*alpha_*(2.0-alpha_)*Lv.transpose()*Rm*Lv;

		// derivatives of Riccati equations
		nabla_Sm_t nabla_dSmdz;
		nabla_Sv_t nabla_dSvdz;
		nabla_s_t nabla_dsdz;

		for (size_t j=0; j<NUM_SUBSYSTEMS-1; j++) {

			// switching time gradient for the original system
			state_matrix_t nabla_dSmdt;
			state_vector_t nabla_dSvdt;
			eigen_scalar_t nabla_dsdt;

			nabla_dSmdt = Am.transpose()*nabla_Sm[j] + nabla_Sm[j].transpose()*Am - nabla_Sm[j].transpose()*Bm*inverseRm*Rm*Lm
					- Lm.transpose()*Rm*inverseRm*Bm.transpose()*nabla_Sm[j];
			nabla_dSmdt = 0.5*(nabla_dSmdt+nabla_dSmdt.transpose()).eval();
			nabla_dSvdt = nablaQv.col(j) + Am.transpose()*nabla_Sv[j] - nabla_Sm[j].transpose()*Bm*inverseRm*Rm*Lv
					- Lm.transpose()*Rm*inverseRm*(nablaRv.col(j) + Bm.transpose()*nabla_Sv[j]);
			nabla_dsdt  = nablaq.col(j) - 0.5*alpha_*(2-alpha_)*(nablaRv.col(j) + Bm.transpose()*nabla_Sv[j]).transpose()*inverseRm*Rm*Lv
					- 0.5*alpha_*(2-alpha_)*Lv.transpose()*Rm*inverseRm*(nablaRv.col(j) + Bm.transpose()*nabla_Sv[j]);

			// switching time gradient for the equivalent system
			nabla_dSmdz[j] = (switchingTimeFinal_-switchingTimeStart_)*nabla_dSmdt;
			nabla_dSvdz[j] = (switchingTimeFinal_-switchingTimeStart_)*nabla_dSvdt;
			nabla_dsdz[j]  = (switchingTimeFinal_-switchingTimeStart_)*nabla_dsdt;

			if (j==activeSubsystem_)  {
				nabla_dSmdz[j] += dSmdt;
				nabla_dSvdz[j] += dSvdt;
				nabla_dsdz[j]  += dsdt;
			}
			if (j==activeSubsystem_-1) {
				nabla_dSmdz[j] -= dSmdt;
				nabla_dSvdz[j] -= dSvdt;
				nabla_dsdz[j]  -= dsdt;
			}

		}  // end of j loop

		convert2Vector(nabla_dSmdz, nabla_dSvdz, nabla_dsdz, derivatives);
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

};

} // namespace ocs2



#endif /* SENSITIVITYSEQUENTIALRICCATIEQUATIONS_OCS2_H_ */