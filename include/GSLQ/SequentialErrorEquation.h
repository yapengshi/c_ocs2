/*
 * SequentialErrorEquation.h
 *
 *  Created on: Apr 13, 2016
 *      Author: farbod
 */

#ifndef SEQUENTIALERROREQUATION_H_
#define SEQUENTIALERROREQUATION_H_

#include "Dimensions.h"

#include "dynamics/SystemBase.h"

#include "misc/LinearInterpolation.h"

template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
class SequentialErrorEquation : public SystemBase<OUTPUT_DIM>
{
public:
	typedef Dimensions<OUTPUT_DIM, INPUT_DIM> DIMENSIONS;
	typedef typename DIMENSIONS::scalar_t 		scalar_t;
	typedef typename DIMENSIONS::scalar_array_t scalar_array_t;
	typedef typename DIMENSIONS::state_vector_t 	  state_vector_t;
	typedef typename DIMENSIONS::state_vector_array_t state_vector_array_t;
	typedef typename DIMENSIONS::state_matrix_t 	  state_matrix_t;
	typedef typename DIMENSIONS::state_matrix_array_t state_matrix_array_t;

	SequentialErrorEquation() {}
	~SequentialErrorEquation() {}

	void setData(const size_t& activeSubsystem, const scalar_t& switchingTimeStart, const scalar_t& switchingTimeFinal,
			scalar_array_t* const timeStampPtr, state_vector_array_t* const GvPtr, state_matrix_array_t* const GmPtr)  {

		activeSubsystem_ = activeSubsystem;
		switchingTimeStart_ = switchingTimeStart;
		switchingTimeFinal_ = switchingTimeFinal;

		GvFunc_.setTimeStamp(timeStampPtr);
		GvFunc_.setData(GvPtr);
		GmFunc_.setTimeStamp(timeStampPtr);
		GmFunc_.setData(GmPtr);
	}

	void computeDerivative(const scalar_t& z, const state_vector_t& Sve, state_vector_t& derivatives) {

		// denormalized time
		scalar_t t = switchingTimeFinal_ - (switchingTimeFinal_-switchingTimeStart_)*(z-activeSubsystem_);

		state_vector_t Gv;
		GvFunc_.interpolate(t, Gv);
		size_t greatestLessTimeStampIndex = GvFunc_.getGreatestLessTimeStampIndex();
		state_matrix_t Gm;
		GmFunc_.interpolate(t, Gm, greatestLessTimeStampIndex);

		// Error equation for the equivalent system
		derivatives = (switchingTimeFinal_-switchingTimeStart_)*(Gm.transpose()*Sve+Gv);
	}


private:
	size_t activeSubsystem_;
	scalar_t switchingTimeStart_;
	scalar_t switchingTimeFinal_;

	LinearInterpolation<state_vector_t,Eigen::aligned_allocator<state_vector_t> > GvFunc_;
	LinearInterpolation<state_matrix_t,Eigen::aligned_allocator<state_matrix_t> > GmFunc_;

};



#endif /* SEQUENTIALERROREQUATION_H_ */
