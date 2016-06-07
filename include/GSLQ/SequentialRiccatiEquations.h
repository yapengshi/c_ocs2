/*
 * SequentialRiccatiEquations.h
 *
 *  Created on: Jan 7, 2016
 *      Author: farbod
 */

#ifndef SEQUENTIALRICCATIEQUATIONS_H_
#define SEQUENTIALRICCATIEQUATIONS_H_

#include "Dimensions.h"

#include "dynamics/SystemBase.h"

#include "misc/LinearInterpolation.h"

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
			scalar_array_t* const timeStampPtr,
			state_matrix_array_t* const AmPtr, control_gain_matrix_array_t* const BmPtr,
			eigen_scalar_array_t* const qPtr, state_vector_array_t* const QvPtr, state_matrix_array_t* const QmPtr,
			control_vector_array_t* const RvPtr, control_matrix_array_t* const RmInversePtr,
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
		RmRmInverseFunc_.setTimeStamp(timeStampPtr);
		RmRmInverseFunc_.setData(RmInversePtr);
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
		RmRmInverseFunc_.interpolate(t, RmInverse, greatestLessTimeStampIndex);
		control_feedback_t Pm;
		PmFunc_.interpolate(t, Pm, greatestLessTimeStampIndex);

		state_matrix_t dSmdt, dSmdz;
		state_vector_t dSvdt, dSvdz;
		eigen_scalar_t dsdt, dsdz;

		// Riccati equations for the original system
		dSmdt = Qm + Am.transpose()*Sm + Sm.transpose()*Am - (Pm+Bm.transpose()*Sm).transpose()*RmInverse*(Pm+Bm.transpose()*Sm);
		dSmdt = 0.5*(dSmdt+dSmdt.transpose()).eval();
		dSvdt = Qv + Am.transpose()*Sv - (Pm+Bm.transpose()*Sm).transpose()*RmInverse*(Rv+Bm.transpose()*Sv);
		dsdt  = q - 0.5*alpha_*(2.0-alpha_)*(Rv+Bm.transpose()*Sv).transpose()*RmInverse*(Rv+Bm.transpose()*Sv);

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
	LinearInterpolation<control_matrix_t,Eigen::aligned_allocator<control_matrix_t> > RmRmInverseFunc_;
	LinearInterpolation<control_feedback_t,Eigen::aligned_allocator<control_feedback_t> > PmFunc_;

};


#endif /* SEQUENTIALRICCATIEQUATIONS_H_ */
