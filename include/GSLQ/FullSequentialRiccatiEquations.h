/*
 * FullSequentialRiccatiEquations.h
 *
 *  Created on: Jan 9, 2016
 *      Author: farbod
 */

#ifndef FULLSEQUENTIALRICCATIEQUATIONS_H_
#define FULLSEQUENTIALRICCATIEQUATIONS_H_

#include "Dimensions.h"

#include "dynamics/SystemBase.h"

#include "misc/LinearInterpolation.h"

template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
class FullSequentialRiccatiEquations : public SystemBase<NUM_Subsystems*(STATE_DIM*STATE_DIM+STATE_DIM+1)>
{
public:
	enum { S_DIM_ = STATE_DIM*STATE_DIM+STATE_DIM+1 };
	typedef Eigen::Matrix<double,S_DIM_,1> s_vector_t;
	typedef Eigen::Matrix<double,NUM_Subsystems*S_DIM_,1> all_s_vector_t;
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

	typedef Eigen::Matrix<double,STATE_DIM,NUM_Subsystems-1> nabla_state_matrix_t;
	typedef Eigen::Matrix<double,INPUT_DIM,NUM_Subsystems-1> nabla_input_matrix_t;
	typedef Eigen::Matrix<double,1,NUM_Subsystems-1> 		 nabla_scalar_rowvector_t;
	typedef std::vector<nabla_state_matrix_t, Eigen::aligned_allocator<nabla_state_matrix_t> > nabla_state_matrix_array_t;
	typedef std::vector<nabla_input_matrix_t, Eigen::aligned_allocator<nabla_input_matrix_t> > nabla_input_matrix_array_t;
	typedef std::vector<nabla_scalar_rowvector_t, Eigen::aligned_allocator<nabla_scalar_rowvector_t> > nabla_scalar_rowvector_array_t;

	FullSequentialRiccatiEquations() {}
	~FullSequentialRiccatiEquations() {}

	static void convert2Vector(const state_matrix_t& Sm, const state_vector_t& Sv, const eigen_scalar_t& s,
			const state_matrix_array_t& nabla_Sm, const state_vector_array_t& nabla_Sv, const eigen_scalar_array_t& nabla_s,
			all_s_vector_t& allSs)  {

		if (nabla_Sm.size() != NUM_Subsystems-1)  throw std::runtime_error("nabla_Sm std::vector size should be the number of subsystems minus one.");
		if (nabla_Sv.size() != NUM_Subsystems-1)  throw std::runtime_error("nabla_Sv std::vector size should be the number of subsystems minus one.");
		if (nabla_s.size() != NUM_Subsystems-1)   throw std::runtime_error("nabla_s std::vector should be the number of subsystems minus one.");

		allSs.template head<S_DIM_>() = (s_vector_t() << Eigen::Map<const Eigen::VectorXd>(Sm.data(),STATE_DIM*STATE_DIM),
														 Eigen::Map<const Eigen::VectorXd>(Sv.data(),STATE_DIM),
														 s).finished();

		for (size_t i=0; i<NUM_Subsystems-1; i++) {

			allSs.template segment<S_DIM_>(S_DIM_+i*S_DIM_) = (s_vector_t() << Eigen::Map<const Eigen::VectorXd>(nabla_Sm[i].data(),STATE_DIM*STATE_DIM),
					Eigen::Map<const Eigen::VectorXd>(nabla_Sv[i].data(),STATE_DIM),
					nabla_s[i]).finished();
		}
	}

	static void convert2Matrix(const s_vector_t& allSs,
			state_matrix_t& Sm, state_vector_t& Sv, eigen_scalar_t& s,
			state_matrix_array_t& nabla_Sm, state_vector_array_t& nabla_Sv, eigen_scalar_array_t& nabla_s)  {

		Sm = Eigen::Map<const Eigen::MatrixXd>(allSs.data(),STATE_DIM, STATE_DIM);
		Sv = Eigen::Map<const Eigen::VectorXd>(allSs.data()+STATE_DIM*STATE_DIM, STATE_DIM);
		s  = Eigen::Map<const Eigen::VectorXd>(allSs.data()+STATE_DIM*STATE_DIM+STATE_DIM, 1);

		for (size_t i=0; i<NUM_Subsystems-1; i++) {

			nabla_Sm.at(i) = Eigen::Map<const Eigen::MatrixXd>(allSs.data()+(i+1)*S_DIM_,STATE_DIM, STATE_DIM);
			nabla_Sv.at(i) = Eigen::Map<const Eigen::VectorXd>(allSs.data()+(i+1)*S_DIM_+STATE_DIM*STATE_DIM, STATE_DIM);
			nabla_s.at(i)  = Eigen::Map<const Eigen::VectorXd>(allSs.data()+(i+1)*S_DIM_+STATE_DIM*STATE_DIM+STATE_DIM, 1);
		}
	}

	void setData(const scalar_t& learningRate,
			const size_t& activeSubsystem, const scalar_t& switchingTimeStart, const scalar_t& switchingTimeFinal,
			scalar_array_t* const timeStampPtr,
			state_matrix_array_t* const AmPtr, control_gain_matrix_array_t* const BmPtr,
			eigen_scalar_array_t* const qPtr, state_vector_array_t* const QvPtr, state_matrix_array_t* const QmPtr,
			control_vector_array_t* const RvPtr, control_matrix_array_t* const RmPtr,
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

	void computeDerivative(const scalar_t& z, const all_s_vector_t& allSs, s_vector_t& derivatives) {

		// denormalized time
		scalar_t t = switchingTimeFinal_ - (switchingTimeFinal_-switchingTimeStart_)*(z-activeSubsystem_);

		state_matrix_t Sm;
		state_vector_t Sv;
		eigen_scalar_t s;
		state_matrix_array_t nabla_Sm(NUM_Subsystems-1);
		state_vector_array_t nabla_Sv(NUM_Subsystems-1);
		eigen_scalar_array_t nabla_s(NUM_Subsystems-1);
		convert2Matrix(allSs, Sm, Sv, s, nabla_Sm, nabla_Sv, nabla_s);

		state_matrix_t Am;
		AmFunc_.interpolate(t, Am);
		control_gain_matrix_t Bm;
		BmFunc_.interpolate(t, Bm);

		eigen_scalar_t q;
		qFunc_.interpolate(t, q);
		state_vector_t Qv;
		QvFunc_.interpolate(t, Qv);
		state_matrix_t Qm;
		QmFunc_.interpolate(t, Qm);
		control_vector_t Rv;
		RvFunc_.interpolate(t, Rv);
		control_matrix_t Rm;
		RmFunc_.interpolate(t, Rm);
		control_feedback_t Pm;
		PmFunc_.interpolate(t, Pm);

		nabla_scalar_rowvector_t nablaq;
		nablaqFunc_.interpolate(t, nablaq);
		nabla_state_matrix_t nablaQv;
		nablaQvFunc_.interpolate(t, nablaQv);
		nabla_input_matrix_t nablaRv;
		nablaRvFunc_.interpolate(t, nablaRv);

		state_matrix_t dSmdt, dSmdz;
		state_vector_t dSvdt, dSvdz;
		eigen_scalar_t dsdt, dsdz;

		state_matrix_array_t dnabla_Smdt(NUM_Subsystems-1);
		state_vector_array_t dnabla_Svdt(NUM_Subsystems-1);
		eigen_scalar_array_t dnabla_sdt(NUM_Subsystems-1);

		state_matrix_array_t dnabla_Smdz(NUM_Subsystems-1);
		state_vector_array_t dnabla_Svdz(NUM_Subsystems-1);
		eigen_scalar_array_t dnabla_sdz(NUM_Subsystems-1);

		// Riccati equations for the original system
		dSmdt = Qm + Am.transpose()*Sm + Sm.transpose()*Am - (Pm+Bm.transpose()*Sm).transpose()*Rm.inverse()*(Pm+Bm.transpose()*Sm);
		dSmdt = 0.5*(dSmdt+dSmdt.transpose()).eval();
		dSvdt = Qv + Am.transpose()*Sv - (Pm+Bm.transpose()*Sm).transpose()*Rm.inverse()*(Rv+Bm.transpose()*Sv);
		dsdt  = q - 0.5*alpha_*(2.0-alpha_)*(Rv+Bm.transpose()*Sv).transpose()*Rm.inverse()*(Rv+Bm.transpose()*Sv);
		// Riccati equations for the equivalent system
		dSmdz = (switchingTimeFinal_-switchingTimeStart_)*dSmdt;
		dSvdz = (switchingTimeFinal_-switchingTimeStart_)*dSvdt;
		dsdz  = (switchingTimeFinal_-switchingTimeStart_)*dsdt;

		for (size_t i=0; i<NUM_Subsystems-1; i++) {


		}

		convert2Vector(dSmdz, dSvdz, dsdz, dnabla_Smdz, dnabla_Svdz, dnabla_sdz, derivatives);
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
	LinearInterpolation<control_matrix_t,Eigen::aligned_allocator<control_matrix_t> > RmFunc_;
	LinearInterpolation<control_feedback_t,Eigen::aligned_allocator<control_feedback_t> > PmFunc_;

	LinearInterpolation<nabla_scalar_rowvector_t,Eigen::aligned_allocator<nabla_scalar_rowvector_t> > nablaqFunc_;
	LinearInterpolation<nabla_state_matrix_t,Eigen::aligned_allocator<nabla_state_matrix_t> > nablaQvFunc_;
	LinearInterpolation<nabla_input_matrix_t,Eigen::aligned_allocator<nabla_input_matrix_t> > nablaRvFunc_;

};


#endif /* FULLSEQUENTIALRICCATIEQUATIONS_H_ */
