/*
 * SequentialRiccatiEquations.h
 *
 *  Created on: Jan 7, 2016
 *      Author: farbod
 */

#ifndef SEQUENTIALRICCATIEQUATIONS_OCS2_H_
#define SEQUENTIALRICCATIEQUATIONS_OCS2_H_

#include "Dimensions.h"

#include "dynamics/SystemBase.h"

#include "misc/LinearInterpolation.h"


namespace ocs2{

template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
class SequentialRiccatiEquations : public SystemBase<OUTPUT_DIM*OUTPUT_DIM+OUTPUT_DIM+1>
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	enum { S_DIM_ = OUTPUT_DIM*OUTPUT_DIM+OUTPUT_DIM+1 };
	typedef Eigen::Matrix<double,S_DIM_,1> s_vector_t;
	typedef Dimensions<OUTPUT_DIM, INPUT_DIM> DIMENSIONS;
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

	static void convert2Vector(const state_matrix_t& Sm, const state_vector_t& Sv, const eigen_scalar_t& s, s_vector_t& allSs)  {

		allSs << Eigen::Map<const Eigen::VectorXd>(Sm.data(),OUTPUT_DIM*OUTPUT_DIM),
				Eigen::Map<const Eigen::VectorXd>(Sv.data(),OUTPUT_DIM),
				s;
	}

	static void convert2Matrix(const s_vector_t& allSs, state_matrix_t& Sm, state_vector_t& Sv, eigen_scalar_t& s)  {

		Sm = Eigen::Map<const Eigen::MatrixXd>(allSs.data(),OUTPUT_DIM,OUTPUT_DIM);
		Sv = Eigen::Map<const Eigen::VectorXd>(allSs.data()+OUTPUT_DIM*OUTPUT_DIM, OUTPUT_DIM);
		s  = allSs.template tail<1>();
	}

	void setData(const scalar_t& learningRate,
			const size_t& activeSubsystem, const scalar_t& switchingTimeStart, const scalar_t& switchingTimeFinal,
			const scalar_array_t* timeStampPtr,
			const state_matrix_array_t* AmPtr, const control_gain_matrix_array_t* BmPtr,
			const eigen_scalar_array_t* qPtr, const state_vector_array_t* QvPtr, const state_matrix_array_t* QmPtr,
			const control_vector_array_t* RvPtr, const control_matrix_array_t* RmInversePtr, const control_matrix_array_t* RmPtr,
			const control_feedback_array_t* PmPtr)  {

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
		RmFunc_.setTimeStamp(timeStampPtr);
		RmFunc_.setData(RmPtr);
		PmFunc_.setTimeStamp(timeStampPtr);
		PmFunc_.setData(PmPtr);
	}

	void computeDerivative(const scalar_t& z, const s_vector_t& allSs, s_vector_t& derivatives) {

		// denormalized time
		scalar_t t = switchingTimeFinal_ - (switchingTimeFinal_-switchingTimeStart_)*(z-activeSubsystem_);

		state_matrix_t Sm;
		state_vector_t Sv;
		eigen_scalar_t s;
		convert2Matrix(allSs, Sm, Sv, s);

		state_matrix_t Am;
		AmFunc_.interpolate(t, Am);
		size_t greatestLessTimeStampIndex = AmFunc_.getGreatestLessTimeStampIndex();
		control_gain_matrix_t Bm;
		BmFunc_.interpolate(t, Bm, greatestLessTimeStampIndex);

		eigen_scalar_t q;
		qFunc_.interpolate(t, q, greatestLessTimeStampIndex);
		state_vector_t Qv;
		QvFunc_.interpolate(t, Qv, greatestLessTimeStampIndex);
		state_matrix_t Qm;
		QmFunc_.interpolate(t, Qm, greatestLessTimeStampIndex);
		control_vector_t Rv;
		RvFunc_.interpolate(t, Rv, greatestLessTimeStampIndex);
		control_matrix_t RmInverse;
		RmInverseFunc_.interpolate(t, RmInverse, greatestLessTimeStampIndex);
		control_matrix_t Rm;
		RmFunc_.interpolate(t, Rm, greatestLessTimeStampIndex);
		control_feedback_t Pm;
		PmFunc_.interpolate(t, Pm, greatestLessTimeStampIndex);

		state_matrix_t dSmdt, dSmdz;
		state_vector_t dSvdt, dSvdz;
		eigen_scalar_t dsdt, dsdz;

		// Riccati equations for the original system
		control_feedback_t Lm = RmInverse*(Pm+Bm.transpose()*Sm);
		control_vector_t   Lv = RmInverse*(Rv+Bm.transpose()*Sv);
		dSmdt = Qm + Am.transpose()*Sm + Sm.transpose()*Am - Lm.transpose()*Rm*Lm;
		dSmdt = 0.5*(dSmdt+dSmdt.transpose()).eval();
		dSvdt = Qv + Am.transpose()*Sv - Lm.transpose()*Rm*Lv;
		dsdt  = q - 0.5*alpha_*(2.0-alpha_)*Lv.transpose()*Rm*Lv;

		// Riccati equations for the equivalent system
		dSmdz = (switchingTimeFinal_-switchingTimeStart_)*dSmdt;
		dSvdz = (switchingTimeFinal_-switchingTimeStart_)*dSvdt;
		dsdz  = (switchingTimeFinal_-switchingTimeStart_)*dsdt;

		convert2Vector(dSmdz, dSvdz, dsdz, derivatives);
	}

protected:
	template <typename Derived>
	bool makePSD(Eigen::MatrixBase<Derived>& squareMatrix) {

		if (squareMatrix.rows() != squareMatrix.cols())  throw std::runtime_error("Not a square matrix: makePSD() method is for square matrix.");

		Eigen::SelfAdjointEigenSolver<Derived> eig(squareMatrix);
		Eigen::VectorXd lambda = eig.eigenvalues();

		bool hasNegativeEigenValue = false;
		for (size_t j=0; j<lambda.size() ; j++)
			if (lambda(j) < 0.0) {
				hasNegativeEigenValue = true;
				lambda(j) = 0.0;
			}

		if (hasNegativeEigenValue)
			squareMatrix = eig.eigenvectors() * lambda.asDiagonal() * eig.eigenvectors().inverse();
		//	else
		//		squareMatrix = 0.5*(squareMatrix+squareMatrix.transpose()).eval();

		return hasNegativeEigenValue;
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
	LinearInterpolation<control_matrix_t,Eigen::aligned_allocator<control_matrix_t> > RmFunc_;
	LinearInterpolation<control_feedback_t,Eigen::aligned_allocator<control_feedback_t> > PmFunc_;

};

}

#endif /* SEQUENTIALRICCATIEQUATIONS_H_ */
