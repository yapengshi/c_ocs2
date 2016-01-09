/*
 * SequentialRiccatiEquations.h
 *
 *  Created on: Jan 7, 2016
 *      Author: farbod
 */

#ifndef SEQUENTIALRICCATIEQUATIONS_H_
#define SEQUENTIALRICCATIEQUATIONS_H_


#include "dynamics/SystemBase.h"

#include "misc/LinearInterpolation.h"

template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
class SequentialRiccatiEquations : public SystemBase<STATE_DIM*STATE_DIM+STATE_DIM+1>
{
public:
	enum { S_DIM_ = STATE_DIM*STATE_DIM+STATE_DIM+1 };
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

	SequentialRiccatiEquations() {}
	~SequentialRiccatiEquations() {}

	static void convert2Vector(const state_matrix_t& Sm, const state_vector_t& Sv, const eigen_scalar_t& s,
			Eigen::Matrix<double,S_DIM_,1>& allSs)  {

		allSs << Eigen::Map<const Eigen::VectorXd>(Sm.data(),STATE_DIM*STATE_DIM),
				Eigen::Map<const Eigen::VectorXd>(Sv.data(),STATE_DIM),
				s;
	}

	static void convert2Matrix(const Eigen::Matrix<double,S_DIM_,1>& allSs,
			state_matrix_t& Sm, state_vector_t& Sv, eigen_scalar_t& s)  {

		Sm = Eigen::Map<const Eigen::MatrixXd>(allSs.data(),STATE_DIM,STATE_DIM);
		Sv = Eigen::Map<const Eigen::VectorXd>(allSs.data()+STATE_DIM*STATE_DIM, STATE_DIM);
		s  = allSs.template tail<1>();
	}

	void setData(const scalar_t& learningRate,
			const size_t& activeSubsystem, const scalar_t& switchingTimeStart, const scalar_t& switchingTimeFinal,
			scalar_array_t* const timeStampPtr,
			state_matrix_array_t* const AmPtr, control_gain_matrix_array_t* const BmPtr,
			eigen_scalar_array_t* const qPtr, state_vector_array_t* const QvPtr, state_matrix_array_t* const QmPtr,
			control_vector_array_t* const RvPtr, control_matrix_array_t* const RmPtr,
			control_feedback_array_t* const PmPtr)  {

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
	}

	void computeDerivative(const scalar_t& z,
			const Eigen::Matrix<double,S_DIM_,1>& allSs,
			Eigen::Matrix<double,S_DIM_,1>& derivatives) {

		// denormalized time
		scalar_t t = switchingTimeFinal_ - (switchingTimeFinal_-switchingTimeStart_)*(z-activeSubsystem_);

		state_matrix_t Sm;
		state_vector_t Sv;
		eigen_scalar_t s;
		convert2Matrix(allSs, Sm, Sv, s);

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

		state_matrix_t dSmdt, dSmdz;
		state_vector_t dSvdt, dSvdz;
		eigen_scalar_t dsdt, dsdz;

		// Riccati equations for the original system
		dSmdt = Qm + Am.transpose()*Sm + Sm.transpose()*Am - (Pm+Bm.transpose()*Sm).transpose()*Rm.inverse()*(Pm+Bm.transpose()*Sm);
		dSmdt = 0.5*(dSmdt+dSmdt.transpose()).eval();
		dSvdt = Qv + Am.transpose()*Sv - (Pm+Bm.transpose()*Sm).transpose()*Rm.inverse()*(Rv+Bm.transpose()*Sv);
		dsdt  = q - 0.5*alpha_*(2.0-alpha_)*(Rv+Bm.transpose()*Sv).transpose()*Rm.inverse()*(Rv+Bm.transpose()*Sv);
		// Riccati equations for the equivalent system
		dSmdz = (switchingTimeFinal_-switchingTimeStart_)*dSmdt;
		dSvdz = (switchingTimeFinal_-switchingTimeStart_)*dSvdt;
		dsdz  = (switchingTimeFinal_-switchingTimeStart_)*dsdt;

		convert2Vector(dSmdz, dSvdz, dsdz, derivatives);
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

};



#endif /* SEQUENTIALRICCATIEQUATIONS_H_ */
