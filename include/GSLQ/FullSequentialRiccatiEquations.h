/*
 * FullSequentialRiccatiEquations.h
 *
 *  Created on: Jan 9, 2016
 *      Author: farbod
 */

#ifndef FULLSEQUENTIALRICCATIEQUATIONS_H_
#define FULLSEQUENTIALRICCATIEQUATIONS_H_

#include <array>

#include "Dimensions.h"

#include "dynamics/SystemBase.h"

#include "misc/LinearInterpolation.h"

template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_SUBSYSTEMS>
class FullSequentialRiccatiEquations : public SystemBase<NUM_SUBSYSTEMS*(STATE_DIM*STATE_DIM+STATE_DIM+1)>
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	enum { S_DIM_ = STATE_DIM*STATE_DIM+STATE_DIM+1 };
	typedef Eigen::Matrix<double,S_DIM_,1> s_vector_t;
	typedef Eigen::Matrix<double,NUM_SUBSYSTEMS*S_DIM_,1> all_s_vector_t;
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

	FullSequentialRiccatiEquations() {}
	~FullSequentialRiccatiEquations() {}

	static void convert2Vector(const state_matrix_t& Sm, const state_vector_t& Sv, const eigen_scalar_t& s,
			const nabla_Sm_t& nabla_Sm, const nabla_Sv_t& nabla_Sv, const nabla_s_t& nabla_s,
			all_s_vector_t& allSs)  {

		allSs.template head<S_DIM_>() << Eigen::Map<const Eigen::VectorXd>(Sm.data(),STATE_DIM*STATE_DIM),
										 Eigen::Map<const Eigen::VectorXd>(Sv.data(),STATE_DIM),
										 s;

		for (size_t j=0; j<NUM_SUBSYSTEMS-1; j++)
			allSs.template segment<S_DIM_>(S_DIM_+j*S_DIM_) << Eigen::Map<const Eigen::VectorXd>(nabla_Sm[j].data(),STATE_DIM*STATE_DIM),
					Eigen::Map<const Eigen::VectorXd>(nabla_Sv[j].data(),STATE_DIM),
					nabla_s[j];
	}

	static void convert2Matrix(const all_s_vector_t& allSs,
			state_matrix_t& Sm, state_vector_t& Sv, eigen_scalar_t& s,
			nabla_Sm_t& nabla_Sm, nabla_Sv_t& nabla_Sv, nabla_s_t& nabla_s)  {

		Sm = Eigen::Map<const Eigen::MatrixXd>(allSs.data(),STATE_DIM, STATE_DIM);
		Sv = Eigen::Map<const Eigen::VectorXd>(allSs.data()+STATE_DIM*STATE_DIM, STATE_DIM);
		s  = Eigen::Map<const Eigen::VectorXd>(allSs.data()+STATE_DIM*STATE_DIM+STATE_DIM, 1);

		for (size_t j=0; j<NUM_SUBSYSTEMS-1; j++) {
			nabla_Sm.at(j) = Eigen::Map<const Eigen::MatrixXd>(allSs.data()+(j+1)*S_DIM_,STATE_DIM, STATE_DIM);
			nabla_Sv.at(j) = Eigen::Map<const Eigen::VectorXd>(allSs.data()+(j+1)*S_DIM_+STATE_DIM*STATE_DIM, STATE_DIM);
			nabla_s.at(j)  = Eigen::Map<const Eigen::VectorXd>(allSs.data()+(j+1)*S_DIM_+STATE_DIM*STATE_DIM+STATE_DIM, 1);
		}
	}

	void setData(const scalar_t& learningRate,
			const size_t& activeSubsystem, const scalar_t& switchingTimeStart, const scalar_t& switchingTimeFinal,
			scalar_array_t* const timeStampPtr,
			state_matrix_array_t* const AmPtr, control_gain_matrix_array_t* const BmPtr,
			eigen_scalar_array_t* const qPtr, state_vector_array_t* const QvPtr, state_matrix_array_t* const QmPtr,
			control_vector_array_t* const RvPtr, control_matrix_array_t* const RmInversePtr,
			control_feedback_array_t* const PmPtr,
			scalar_array_t* const sensitivityTimeStampPtr, nabla_scalar_rowvector_array_t* const nablaqPtr,
			nabla_state_matrix_array_t* const nablaQvPtr, nabla_input_matrix_array_t* const nablaRvPtr)  {

		alpha_ = learningRate;

		activeSubsystem_ = activeSubsystem;
		switchingTimeStart_ = switchingTimeStart;
		switchingTimeFinal_ = switchingTimeFinal;

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
		scalar_t t = switchingTimeFinal_ - (switchingTimeFinal_-switchingTimeStart_)*(z-activeSubsystem_);

		state_matrix_t Sm;
		state_vector_t Sv;
		eigen_scalar_t s;
		nabla_Sm_t nabla_Sm;
		nabla_Sv_t nabla_Sv;
		nabla_s_t  nabla_s;
		convert2Matrix(allSs, Sm, Sv, s, nabla_Sm, nabla_Sv, nabla_s);

		size_t greatestLessTimeStampIndex;

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
		state_matrix_t dSmdt = Qm + Am.transpose()*Sm + Sm.transpose()*Am - (Pm+Bm.transpose()*Sm).transpose()*inverseRm*(Pm+Bm.transpose()*Sm);
		dSmdt = 0.5*(dSmdt+dSmdt.transpose()).eval();
		state_vector_t dSvdt = Qv + Am.transpose()*Sv - (Pm+Bm.transpose()*Sm).transpose()*inverseRm*(Rv+Bm.transpose()*Sv);
		eigen_scalar_t dsdt  = q - 0.5*alpha_*(2.0-alpha_)*(Rv+Bm.transpose()*Sv).transpose()*inverseRm*(Rv+Bm.transpose()*Sv);
		// Riccati equations for the equivalent system
		state_matrix_t dSmdz = (switchingTimeFinal_-switchingTimeStart_)*dSmdt;
		state_vector_t dSvdz = (switchingTimeFinal_-switchingTimeStart_)*dSvdt;
		eigen_scalar_t dsdz  = (switchingTimeFinal_-switchingTimeStart_)*dsdt;

		// derivatives of Riccati equations
		nabla_Sm_t nabla_dSmdz;
		nabla_Sv_t nabla_dSvdz;
		nabla_s_t nabla_dsdz;

		for (size_t j=0; j<NUM_SUBSYSTEMS-1; j++) {

			// switching time gradient for the original system
			state_matrix_t nabla_dSmdt;
			state_vector_t nabla_dSvdt;
			eigen_scalar_t nabla_dsdt;

			nabla_dSmdt = Am.transpose()*nabla_Sm[j] + nabla_Sm[j].transpose()*Am - nabla_Sm[j].transpose()*Bm*inverseRm*(Pm+Bm.transpose()*Sm)
					- (Pm+Bm.transpose()*Sm).transpose()*inverseRm*Bm.transpose()*nabla_Sm[j];
			nabla_dSmdt = 0.5*(nabla_dSmdt+nabla_dSmdt.transpose()).eval();
			nabla_dSvdt = nablaQv.col(j) + Am.transpose()*nabla_Sv[j] - nabla_Sm[j].transpose()*Bm*inverseRm*(Rv + Bm.transpose()*Sv)
					- (Pm+Bm.transpose()*Sm).transpose()*inverseRm*(nablaRv.col(j) + Bm.transpose()*nabla_Sv[j]);
			nabla_dsdt  = nablaq.col(j) - 0.5*alpha_*(2-alpha_)*(nablaRv.col(j) + Bm.transpose()*nabla_Sv[j]).transpose()*inverseRm*(Rv+Bm.transpose()*Sv)
					- 0.5*alpha_*(2-alpha_)*(Rv+Bm.transpose()*Sv).transpose()*inverseRm*(nablaRv.col(j) + Bm.transpose()*nabla_Sv[j]);

			// switching time gradient for the equvalent system
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

		convert2Vector(dSmdz, dSvdz, dsdz, nabla_dSmdz, nabla_dSvdz, nabla_dsdz, derivatives);
	}


private:
	scalar_t alpha_;

	size_t activeSubsystem_;
	scalar_t switchingTimeStart_;
	scalar_t switchingTimeFinal_;

	LinearInterpolation<state_matrix_t,Eigen::aligned_allocator<state_matrix_t> > AmFunc_;
	LinearInterpolation<control_gain_matrix_t,Eigen::aligned_allocator<control_gain_matrix_t> > BmFunc_;

	LinearInterpolation<eigen_scalar_t,Eigen::aligned_allocator<eigen_scalar_t> > qFunc_;
	LinearInterpolation<state_vector_t,Eigen::aligned_allocator<state_vector_t> > QvFunc_;
	LinearInterpolation<state_matrix_t,Eigen::aligned_allocator<state_matrix_t> > QmFunc_;
	LinearInterpolation<control_vector_t,Eigen::aligned_allocator<control_vector_t> > RvFunc_;
	LinearInterpolation<control_matrix_t,Eigen::aligned_allocator<control_matrix_t> > RmInverseFunc_;
	LinearInterpolation<control_feedback_t,Eigen::aligned_allocator<control_feedback_t> > PmFunc_;

	LinearInterpolation<nabla_scalar_rowvector_t,Eigen::aligned_allocator<nabla_scalar_rowvector_t> > nablaqFunc_;
	LinearInterpolation<nabla_state_matrix_t,Eigen::aligned_allocator<nabla_state_matrix_t> > nablaQvFunc_;
	LinearInterpolation<nabla_input_matrix_t,Eigen::aligned_allocator<nabla_input_matrix_t> > nablaRvFunc_;

};


#endif /* FULLSEQUENTIALRICCATIEQUATIONS_H_ */
