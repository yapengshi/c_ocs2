/*
 * Dimensions.h
 *
 *  Created on: Jan 3, 2016
 *      Author: farbod
 */

#ifndef DIMENSIONS_H_
#define DIMENSIONS_H_

#include <vector>
#include <Eigen/Dense>
#include <Eigen/StdVector>


template <size_t STATE_DIM, size_t INPUT_DIM>
class Dimensions
{

public:
	typedef Eigen::Matrix<double, STATE_DIM, 1> state_vector_t;
	typedef std::vector<state_vector_t, Eigen::aligned_allocator<state_vector_t> > state_vector_array_t;

	typedef Eigen::Matrix<double, STATE_DIM, STATE_DIM> state_matrix_t;
	typedef std::vector<state_matrix_t, Eigen::aligned_allocator<state_matrix_t> > state_matrix_array_t;

	typedef Eigen::Matrix<double, STATE_DIM, INPUT_DIM> control_gain_matrix_t;
	typedef std::vector<control_gain_matrix_t, Eigen::aligned_allocator<control_gain_matrix_t> > control_gain_matrix_array_t;

	typedef Eigen::Matrix<double, INPUT_DIM, STATE_DIM> control_feedback_t;
	typedef std::vector<control_feedback_t, Eigen::aligned_allocator<control_feedback_t> > control_feedback_array_t;

	typedef Eigen::Matrix<double, INPUT_DIM, 1> control_vector_t;
	typedef std::vector<control_vector_t, Eigen::aligned_allocator<control_vector_t> > control_vector_array_t;

	typedef Eigen::Matrix<double, INPUT_DIM, INPUT_DIM> control_matrix_t;
	typedef std::vector<control_matrix_t, Eigen::aligned_allocator<control_matrix_t> > control_matrix_array_t;

    typedef Eigen::Matrix<double, STATE_DIM*STATE_DIM , 1 > state_matrix_vectorized_t;

	typedef double scalar_t;
	typedef std::vector<scalar_t> scalar_array_t;

	struct controller_t {
		scalar_array_t time_;
		control_vector_array_t uff_;
		control_feedback_array_t k_;
	};

	enum DIMS {
		STATE_DIM_ = STATE_DIM,
		INPUT_DIM_ = INPUT_DIM
	};

private:

};



#endif /* DIMENSIONS_H_ */
