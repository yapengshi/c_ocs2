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


template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM=STATE_DIM>
class Dimensions
{

public:
	enum DIMS {
		STATE_DIM_  = STATE_DIM,
		INPUT_DIM_  = INPUT_DIM,
		OUTPUT_DIM_ = OUTPUT_DIM,
		MAX_CONSTRAINT1_DIM_ = INPUT_DIM,
	};

	typedef Eigen::Matrix<double, STATE_DIM, 1> state_vector_t;
	typedef std::vector<state_vector_t, Eigen::aligned_allocator<state_vector_t> > state_vector_array_t;

	typedef Eigen::Matrix<double, OUTPUT_DIM, 1> output_vector_t;
	typedef std::vector<output_vector_t, Eigen::aligned_allocator<output_vector_t> > output_vector_array_t;

	typedef Eigen::Matrix<double, OUTPUT_DIM, OUTPUT_DIM> state_matrix_t;
	typedef std::vector<state_matrix_t, Eigen::aligned_allocator<state_matrix_t> > state_matrix_array_t;

	typedef Eigen::Matrix<double, OUTPUT_DIM, INPUT_DIM> control_gain_matrix_t;
	typedef std::vector<control_gain_matrix_t, Eigen::aligned_allocator<control_gain_matrix_t> > control_gain_matrix_array_t;

	typedef Eigen::Matrix<double, INPUT_DIM, OUTPUT_DIM> control_feedback_t;
	typedef std::vector<control_feedback_t, Eigen::aligned_allocator<control_feedback_t> > control_feedback_array_t;

	typedef Eigen::Matrix<double, INPUT_DIM, 1> control_vector_t;
	typedef std::vector<control_vector_t, Eigen::aligned_allocator<control_vector_t> > control_vector_array_t;

	typedef Eigen::Matrix<double, INPUT_DIM, INPUT_DIM> control_matrix_t;
	typedef std::vector<control_matrix_t, Eigen::aligned_allocator<control_matrix_t> > control_matrix_array_t;

    typedef Eigen::Matrix<double, OUTPUT_DIM*OUTPUT_DIM , 1 > state_matrix_vectorized_t;

    typedef Eigen::Matrix<double, MAX_CONSTRAINT1_DIM_, 1> constraint1_vector_t;
    typedef std::vector<constraint1_vector_t, Eigen::aligned_allocator<constraint1_vector_t> > constraint1_vector_array_t;

    typedef Eigen::Matrix<double, MAX_CONSTRAINT1_DIM_, 1> constraint1_matrix_t;
    typedef std::vector<constraint1_matrix_t, Eigen::aligned_allocator<constraint1_matrix_t> > constraint1_matrix_array_t;

    typedef Eigen::Matrix<double, MAX_CONSTRAINT1_DIM_, OUTPUT_DIM> constraint1_state_matrix_t;
    typedef std::vector<constraint1_state_matrix_t, Eigen::aligned_allocator<constraint1_state_matrix_t> > constraint1_state_matrix_array_t;

    typedef Eigen::Matrix<double, MAX_CONSTRAINT1_DIM_, INPUT_DIM> constraint1_control_matrix_t;
    typedef std::vector<constraint1_control_matrix_t, Eigen::aligned_allocator<constraint1_control_matrix_t> > constraint1_control_matrix_array_t;

    typedef Eigen::Matrix<double, MAX_CONSTRAINT1_DIM_, INPUT_DIM> control_constraint1_matrix_t;
    typedef std::vector<control_constraint1_matrix_t, Eigen::aligned_allocator<control_constraint1_matrix_t> > control_constraint1_matrix_array_t;

	typedef double scalar_t;
	typedef std::vector<scalar_t> scalar_array_t;

	typedef Eigen::Matrix<double, 1, 1> eigen_scalar_t;
	typedef std::vector<eigen_scalar_t, Eigen::aligned_allocator<eigen_scalar_t> > eigen_scalar_array_t;

	template <int DIM>
	class LinearFunction_t {
	public:
		scalar_array_t time_;
		std::vector<Eigen::Matrix<double, DIM, 1>, Eigen::aligned_allocator<Eigen::Matrix<double, DIM, 1>> > uff_;
		std::vector<Eigen::Matrix<double, DIM, 1>, Eigen::aligned_allocator<Eigen::Matrix<double, DIM, 1>> > deltaUff_;
		std::vector<Eigen::Matrix<double, DIM, OUTPUT_DIM>, Eigen::aligned_allocator<Eigen::Matrix<double, DIM, OUTPUT_DIM>> > k_;
	};
	typedef LinearFunction_t<INPUT_DIM> controller_t;

	struct Options {
		Options() :
			maxIterationGSLQP_(10),
			minLearningRateGSLQP_(0.05),
			minRelCostGSLQP_(1e-3),
			meritFunctionRho_(1.0),
			dispayGSLQP_(false),
			warmStartGSLQP_(false),

			AbsTolODE_(1e-9),
			RelTolODE_(1e-6),
			simulationIsConstrained_(false),
			minAbsConstraint1RMSE_(1e-3),
			minRelConstraint1RMSE_(1e-3),

			displayIPOPT_(true),
			tolIPOPT_(1e-2),
			acceptableTolIPOPT_(1e-1),
			maxIterationIPOPT_(20),
			minAcceptedSwitchingTimeDifference_(0.0)
		{}

		size_t maxIterationGSLQP_;
		double minLearningRateGSLQP_;
		double minRelCostGSLQP_;
		double meritFunctionRho_;
		bool dispayGSLQP_;
		bool warmStartGSLQP_;

		double AbsTolODE_;
		double RelTolODE_;
		bool simulationIsConstrained_;
		double minAbsConstraint1RMSE_;
		double minRelConstraint1RMSE_;

		bool displayIPOPT_;
		double tolIPOPT_;
		double acceptableTolIPOPT_;
		size_t maxIterationIPOPT_;
		double minAcceptedSwitchingTimeDifference_;
	};


private:

};

using HyQDimensions = Dimensions<36,12>;

#endif /* DIMENSIONS_H_ */
