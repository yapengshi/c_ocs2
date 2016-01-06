/*
 * PartialRiccatiEquations.h
 *
 *  Created on: Jan 5, 2016
 *      Author: farbod
 */

#ifndef RICCATIEQUATIONS_H_
#define RICCATIEQUATIONS_H_

#include "dynamics/SystemBase.h"

template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
class PartialRiccatiEquations : public SystemBase<STATE_DIM*STATE_DIM+STATE_DIM+1>
{
public:
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

	PartialRiccatiEquations() {}
	~PartialRiccatiEquations() {}

	static void convert2Vector(const state_matrix_t& Sm, const state_vector_t& Sv, const eigen_scalar_t& s,
			Eigen::Matrix<double,STATE_DIM*STATE_DIM+STATE_DIM+1,1>& allSs)  {

		allSs << Eigen::Map<const Eigen::VectorXd>(Sm.data(),STATE_DIM*STATE_DIM),
				Eigen::Map<const Eigen::VectorXd>(Sv.data(),STATE_DIM),
				s;
	}

	static void convert2Matrix(const Eigen::Matrix<double,STATE_DIM*STATE_DIM+STATE_DIM+1,1>& allSs,
			state_matrix_t& Sm, state_vector_t& Sv, eigen_scalar_t& s)  {

		Sm = Eigen::Map<const Eigen::MatrixXd>(allSs.data(),STATE_DIM,STATE_DIM);
		Sv = Eigen::Map<const Eigen::VectorXd>(allSs.data()+STATE_DIM*STATE_DIM, STATE_DIM);
		s  = allSs.template tail<1>();
	}

	void setData(const scalar_t& timeStart, const scalar_t& timeFinal,
			const state_matrix_t& Am, const control_gain_matrix_t& Bm,
			const eigen_scalar_t& q, const state_vector_t& Qv, const state_matrix_t& Qm,
			const control_vector_t& Rv, const control_matrix_t& Rm,
			const control_feedback_t& Pm)  {

		timeStart_ = timeStart;
		timeFinal_ = timeFinal;
		Am_ = Am;
		Bm_ = Bm;
		q_ = q;
		Qv_ = Qv;
		Qm_ = Qm;
		Rv_ = Rv;
		Rm_ = Rm;
		Pm_ = Pm;
	}

	void computeDerivative(const scalar_t& t,
			const Eigen::Matrix<double,STATE_DIM*STATE_DIM+STATE_DIM+1,1>& state,
			Eigen::Matrix<double,STATE_DIM*STATE_DIM+STATE_DIM+1,1>& derivative) {

		state_matrix_t Sm;
		state_vector_t Sv;
		eigen_scalar_t s;
		convert2Matrix(state, Sm, Sv, s);

		state_matrix_t dSmdt, dSmdz;
		state_vector_t dSvdt, dSvdz;
		eigen_scalar_t dsdt, dsdz;

		// Riccati equations for the original system
		dSmdt = Qm_ + Am_.transpose()*Sm + Sm.transpose()*Am_ - (Pm_+Bm_.transpose()*Sm).transpose()*Rm_.inverse()*(Pm_+Bm_.transpose()*Sm);
		dSmdt = (dSmdt+dSmdt.transpose()).eval()*0.5;
		dSvdt = Qv_ + Am_.transpose()*Sv - (Pm_+Bm_.transpose()*Sm).transpose()*Rm_.inverse()*(Rv_+Bm_.transpose()*Sv);
		dsdt  = q_ - 0.5*(Rv_+Bm_.transpose()*Sv).transpose()*Rm_.inverse()*(Rv_+Bm_.transpose()*Sv);
		// Riccati equations for the equivalent system
		dSmdz = (timeFinal_-timeStart_)*dSmdt;
		dSvdz = (timeFinal_-timeStart_)*dSvdt;
		dsdz  = (timeFinal_-timeStart_)*dsdt;

		convert2Vector(dSmdz, dSvdz, dsdz, derivative);
	}


private:
	scalar_t timeStart_;
	scalar_t timeFinal_;

	state_matrix_t Am_;
	control_gain_matrix_t Bm_;

	eigen_scalar_t q_;
	state_vector_t Qv_;
	state_matrix_t Qm_;
	control_vector_t Rv_;
	control_matrix_t Rm_;
	control_feedback_t Pm_;

};



#endif /* RICCATIEQUATIONS_H_ */
