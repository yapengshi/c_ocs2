/*
 * dimensions.hpp
 *
 *  Created on: 26.03.2014
 *      Author: neunertm
 */

#ifndef DIMENSIONS_HPP_
#define DIMENSIONS_HPP_

#include <vector>
#include <Eigen/Dense>


template <size_t STATE_DIM, size_t CONTROL_DIM>
class Dimensions {

public:
	enum DIMS {
		STATE_SIZE = STATE_DIM,
		CONTROL_SIZE = CONTROL_DIM
	};

	typedef Eigen::Matrix<double, STATE_DIM, 1> state_vector_t;
	typedef std::vector<state_vector_t, Eigen::aligned_allocator<state_vector_t> > state_vector_array_t;

	typedef Eigen::Matrix<double, STATE_DIM, STATE_DIM> state_matrix_t;
	typedef std::vector<state_matrix_t, Eigen::aligned_allocator<state_matrix_t> > state_matrix_array_t;

	typedef Eigen::Matrix<double, STATE_DIM, CONTROL_DIM> control_gain_matrix_t;
	typedef std::vector<control_gain_matrix_t, Eigen::aligned_allocator<control_gain_matrix_t> > control_gain_matrix_array_t;

	typedef Eigen::Matrix<double, CONTROL_DIM, STATE_DIM> control_feedback_t;
	typedef std::vector<control_feedback_t, Eigen::aligned_allocator<control_feedback_t> > control_feedback_array_t;

	typedef Eigen::Matrix<double, CONTROL_DIM, 1> control_vector_t;
	typedef std::vector<control_vector_t, Eigen::aligned_allocator<control_vector_t> > control_vector_array_t;

	typedef Eigen::Matrix<double, CONTROL_DIM, CONTROL_DIM> control_matrix_t;
	typedef std::vector<control_matrix_t, Eigen::aligned_allocator<control_matrix_t> > control_matrix_array_t;

    typedef Eigen::Matrix<double, STATE_DIM*STATE_DIM , 1 > state_matrix_vectorized_t;

	typedef double scalar_t;
	typedef std::vector<scalar_t> scalar_array_t;

private:

};
#endif /* DIMENSIONS_HPP_ */
