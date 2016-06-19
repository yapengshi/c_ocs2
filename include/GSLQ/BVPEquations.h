/*
 * BVPEquations.h
 *
 *  Created on: Jun 18, 2016
 *      Author: farbod
 */

#ifndef BVPEQUATIONS_OCS2_H_
#define BVPEQUATIONS_OCS2_H_

#include <vector>
#include <Eigen/Dense>
#include <Eigen/StdVector>

#include "dynamics/SystemBase.h"
#include "misc/LinearInterpolation.h"


namespace ocs2{

template <size_t STATE_DIM, size_t INPUT_DIM>
class BVPEquations : public SystemBase<STATE_DIM*STATE_DIM+STATE_DIM>
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	enum {
		Full_ODE_VECTOR_DIM = STATE_DIM*STATE_DIM+STATE_DIM
	};

	typedef Eigen::Matrix<double, Full_ODE_VECTOR_DIM, 1> full_ode_vector_t;
	typedef std::vector<full_ode_vector_t, Eigen::aligned_allocator<full_ode_vector_t> > full_ode_vector_array_t;

	typedef Eigen::Matrix<double, 1, 1> eigen_scalar_t;
	typedef Eigen::Matrix<double, STATE_DIM, 1> state_vector_t;
	typedef Eigen::Matrix<double, INPUT_DIM, 1> input_vector_t;
	typedef Eigen::Matrix<double, STATE_DIM, STATE_DIM> state_state_matrix_t;
	typedef Eigen::Matrix<double, INPUT_DIM, INPUT_DIM> input_input_matrix_t;
	typedef Eigen::Matrix<double, INPUT_DIM, STATE_DIM> input_state_matrix_t;
	typedef Eigen::Matrix<double, STATE_DIM, INPUT_DIM> state_input_matrix_t;

	typedef std::vector<double> double_array_t;
	typedef std::vector<eigen_scalar_t, Eigen::aligned_allocator<eigen_scalar_t> > eigen_scalar_array_t;
	typedef std::vector<state_vector_t, Eigen::aligned_allocator<state_vector_t> > state_vector_array_t;
	typedef std::vector<input_vector_t, Eigen::aligned_allocator<input_vector_t> > input_vector_array_t;
	typedef std::vector<state_state_matrix_t, Eigen::aligned_allocator<state_state_matrix_t> > state_state_matrix_array_t;
	typedef std::vector<input_input_matrix_t, Eigen::aligned_allocator<input_input_matrix_t> > input_input_matrix_array_t;
	typedef std::vector<input_state_matrix_t, Eigen::aligned_allocator<input_state_matrix_t> > input_state_matrix_array_t;
	typedef std::vector<state_input_matrix_t, Eigen::aligned_allocator<state_input_matrix_t> > state_input_matrix_array_t;


	BVPEquations() {}
	~BVPEquations() {}

	static void convert2Vector(const state_state_matrix_t& Mm, const state_vector_t& Sv, full_ode_vector_t& MSv)  {

		MSv << Eigen::Map<const Eigen::VectorXd>(Mm.data(),STATE_DIM*STATE_DIM),
				Eigen::Map<const Eigen::VectorXd>(Sv.data(),STATE_DIM);
	}

	static void convert2Matrix(const full_ode_vector_t& MSv, state_state_matrix_t& Mm, state_vector_t& Sv)  {

		Mm = Eigen::Map<const Eigen::MatrixXd>(MSv.data(),STATE_DIM,STATE_DIM);
		Sv = Eigen::Map<const Eigen::VectorXd>(MSv.data()+STATE_DIM*STATE_DIM, STATE_DIM);
	}

	void setData(const double_array_t* timeStampPtr,
			const state_state_matrix_array_t* AmPtr, const state_state_matrix_array_t* OmPtr, const state_input_matrix_array_t* BmPtr, const state_vector_array_t* GvPtr,
			const state_vector_array_t* QvPtr, const state_state_matrix_array_t* QmPtr, const input_state_matrix_array_t* PmPtr,
			const input_vector_array_t* RvPtr, const input_input_matrix_array_t* RmPtr, const input_input_matrix_array_t* RmInversePtr)  {
		AmFunc_.setTimeStamp(timeStampPtr);
		OmFunc_.setTimeStamp(timeStampPtr);
		BmFunc_.setTimeStamp(timeStampPtr);
		GvFunc_.setTimeStamp(timeStampPtr);

		QvFunc_.setTimeStamp(timeStampPtr);
		QmFunc_.setTimeStamp(timeStampPtr);
		PmFunc_.setTimeStamp(timeStampPtr);

		RvFunc_.setTimeStamp(timeStampPtr);
		RmFunc_.setTimeStamp(timeStampPtr);
		RmInverseFunc_.setTimeStamp(timeStampPtr);

		if (AmPtr) AmFunc_.setData(AmPtr); else  throw std::runtime_error("Am is not set.");
		if (OmPtr) OmFunc_.setData(OmPtr); else  OmFunc_.setZero();
		if (BmPtr) BmFunc_.setData(BmPtr); else  throw std::runtime_error("Bm is not set.");
		if (GvPtr) GvFunc_.setData(GvPtr); else  GvFunc_.setZero();

		if (QvPtr) QvFunc_.setData(QvPtr); else  QvFunc_.setZero();
		if (QmPtr) QmFunc_.setData(QmPtr); else  QmFunc_.setZero();
		if (PmPtr) PmFunc_.setData(PmPtr); else  PmFunc_.setZero();

		if (RvPtr) RvFunc_.setData(RvPtr); else  RvFunc_.setZero();
		if (RmPtr) RmFunc_.setData(RmPtr); else  throw std::runtime_error("Rm is not set.");
		if (RmInversePtr) RmInverseFunc_.setData(RmInversePtr); else  throw std::runtime_error("RmInverse is not set.");

		startTime_ = timeStampPtr->front();
		finalTime_ = timeStampPtr->back();
	}

	void computeDerivative(const double& z, const full_ode_vector_t& MSv, full_ode_vector_t& derivatives) override {

		// denormalized time
		if (z>1 || z<0)  throw std::runtime_error("The normalized time should be between zero and one.");

		double t = finalTime_ - (finalTime_-startTime_)*z;

		state_state_matrix_t Mm;
		state_vector_t Sv;
		convert2Matrix(MSv, Mm, Sv);

		state_state_matrix_t Am;
		AmFunc_.interpolate(t, Am);
		size_t greatestLessTimeStampIndex = AmFunc_.getGreatestLessTimeStampIndex();
		state_state_matrix_t Om;
		OmFunc_.interpolate(t, Om, greatestLessTimeStampIndex);
		state_input_matrix_t Bm;
		BmFunc_.interpolate(t, Bm, greatestLessTimeStampIndex);
		state_vector_t Gv;
		GvFunc_.interpolate(t, Gv, greatestLessTimeStampIndex);

		state_vector_t Qv;
		QvFunc_.interpolate(t, Qv, greatestLessTimeStampIndex);
		state_state_matrix_t Qm;
		QmFunc_.interpolate(t, Qm, greatestLessTimeStampIndex);
		input_state_matrix_t Pm;
		PmFunc_.interpolate(t, Pm, greatestLessTimeStampIndex);
		input_vector_t Rv;
		RvFunc_.interpolate(t, Rv, greatestLessTimeStampIndex);
		input_input_matrix_t Rm;
		RmFunc_.interpolate(t, Rm, greatestLessTimeStampIndex);
		input_input_matrix_t RmInverse;
		RmInverseFunc_.interpolate(t, RmInverse, greatestLessTimeStampIndex);

		state_state_matrix_t dMmdt, dMmdz;
		state_vector_t dSvdt, dSvdz;

		// Uv = -Lv - Km*x
		input_vector_t       Lv = RmInverse*(Rv+Bm.transpose()*Sv);
		input_state_matrix_t Km = RmInverse*(Pm+Bm.transpose()*Mm);

		// Riccati equations for the original system
		dMmdt = Qm + Am.transpose()*Mm + Mm.transpose()*Am + Mm.transpose()*Om*Mm - Km.transpose()*Rm*Km;
		dMmdt = 0.5*(dMmdt+dMmdt.transpose()).eval();
		dSvdt = (Qv+Mm*Gv) + Am.transpose()*Sv + Mm.transpose()*Om*Sv - Km.transpose()*Rm*Lv;

		// Riccati equations for the equivalent system
		dMmdz = (finalTime_-startTime_)*dMmdt;
		dSvdz = (finalTime_-startTime_)*dSvdt;

		convert2Vector(dMmdz, dSvdz, derivatives);
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
	double startTime_;
	double finalTime_;

	LinearInterpolation<state_state_matrix_t, Eigen::aligned_allocator<state_state_matrix_t> > AmFunc_;
	LinearInterpolation<state_state_matrix_t, Eigen::aligned_allocator<state_state_matrix_t> > OmFunc_;
	LinearInterpolation<state_input_matrix_t, Eigen::aligned_allocator<state_input_matrix_t> > BmFunc_;
	LinearInterpolation<state_vector_t, Eigen::aligned_allocator<state_vector_t> >             GvFunc_;

	LinearInterpolation<state_vector_t, Eigen::aligned_allocator<state_vector_t> > QvFunc_;
	LinearInterpolation<input_vector_t, Eigen::aligned_allocator<input_vector_t> > RvFunc_;

	LinearInterpolation<state_state_matrix_t, Eigen::aligned_allocator<state_state_matrix_t> > QmFunc_;
	LinearInterpolation<input_state_matrix_t, Eigen::aligned_allocator<input_state_matrix_t> > PmFunc_;
	LinearInterpolation<input_input_matrix_t, Eigen::aligned_allocator<input_input_matrix_t> > RmFunc_;
	LinearInterpolation<input_input_matrix_t, Eigen::aligned_allocator<input_input_matrix_t> > RmInverseFunc_;

};

}  // end of OCS2 name space



#endif /* BVPEQUATIONS_H_OCS2_ */
