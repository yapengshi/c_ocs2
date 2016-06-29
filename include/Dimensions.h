/*
 * Dimensions.h
 *
 *  Created on: Jan 3, 2016
 *      Author: farbod
 */

#ifndef DIMENSIONS_OCS2_H_
#define DIMENSIONS_OCS2_H_

#include <vector>
#include <Eigen/Dense>
#include <Eigen/StdVector>


namespace ocs2{

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
	struct LinearFunction_t {
		scalar_array_t time_;
		std::vector<Eigen::Matrix<double, DIM, 1>, Eigen::aligned_allocator<Eigen::Matrix<double, DIM, 1>> > uff_;
		std::vector<Eigen::Matrix<double, DIM, 1>, Eigen::aligned_allocator<Eigen::Matrix<double, DIM, 1>> > deltaUff_;
		std::vector<Eigen::Matrix<double, DIM, OUTPUT_DIM>, Eigen::aligned_allocator<Eigen::Matrix<double, DIM, OUTPUT_DIM>> > k_;
	};
	typedef LinearFunction_t<INPUT_DIM> controller_t;

	struct Options {
	public:
		Options() :
			maxIterationGSLQP_(10),
			minLearningRateGSLQP_(0.05),
			maxLearningRateGSLQP_(1.0),
			minRelCostGSLQP_(1e-3),
			meritFunctionRho_(1.0),
			constraintStepSize_(1.0),
			lineSearchByMeritFuntion_(false),
			dispayGSLQP_(false),
			warmStartGSLQP_(false),

			AbsTolODE_(1e-9),
			RelTolODE_(1e-6),
			simulationIsConstrained_(false),
			minAbsConstraint1ISE_(1e-3),
			minRelConstraint1ISE_(1e-3),

			displayIPOPT_(true),
			tolIPOPT_(1e-2),
			acceptableTolIPOPT_(1e-1),
			maxIterationIPOPT_(20),
			minAcceptedSwitchingTimeDifference_(0.0)
		{}

		void print()
		{
			std::cout << " #### ========================== Options ============================ ####" << std::endl;
			std::cout << "maxIterationGSLQP_                 " << maxIterationGSLQP_                  << std::endl;
			std::cout << "minLearningRateGSLQP_              " << minLearningRateGSLQP_               << std::endl;
			std::cout << "maxLearningRateGSLQP_              " << maxLearningRateGSLQP_               << std::endl;
			std::cout << "minRelCostGSLQP_                   " << minRelCostGSLQP_                    << std::endl;
			std::cout << "meritFunctionRho_                  " << meritFunctionRho_                   << std::endl;
			std::cout << "constraintStepSize_                " << constraintStepSize_                 << std::endl;
			std::cout << "lineSearchByMeritFuntion_          " << lineSearchByMeritFuntion_           << std::endl;
			std::cout << "dispayGSLQP_                       " << dispayGSLQP_                        << std::endl;
			std::cout << "warmStartGSLQP_                    " << warmStartGSLQP_                     << std::endl;
            std::cout << "                                   " << ""                                  << std::endl;
			std::cout << "AbsTolODE_                         " << AbsTolODE_                          << std::endl;
			std::cout << "RelTolODE_                         " << RelTolODE_                          << std::endl;
			std::cout << "simulationIsConstrained_           " << simulationIsConstrained_            << std::endl;
			std::cout << "minAbsConstraint1ISE_              " << minAbsConstraint1ISE_               << std::endl;
			std::cout << "minRelConstraint1ISE_              " << minRelConstraint1ISE_               << std::endl;
            std::cout << "                                   " << ""                                  << std::endl;
			std::cout << "displayIPOPT_                      " << displayIPOPT_                       << std::endl;
			std::cout << "tolIPOPT_                          " << tolIPOPT_                           << std::endl;
			std::cout << "acceptableTolIPOPT_                " << acceptableTolIPOPT_                 << std::endl;
			std::cout << "maxIterationIPOPT_                 " << maxIterationIPOPT_                  << std::endl;
			std::cout << "minAcceptedSwitchingTimeDifference_" << minAcceptedSwitchingTimeDifference_ << std::endl;
			std::cout << " #### ============================ end ============================== ####" << std::endl;
			std::cout << std::endl;
		}

		size_t maxIterationGSLQP_;
		double minLearningRateGSLQP_;
		double maxLearningRateGSLQP_;
		double minRelCostGSLQP_;
		double meritFunctionRho_;
		double constraintStepSize_;
		bool lineSearchByMeritFuntion_;
		bool dispayGSLQP_;
		bool warmStartGSLQP_;

		double AbsTolODE_;
		double RelTolODE_;
		bool simulationIsConstrained_;
		double minAbsConstraint1ISE_;
		double minRelConstraint1ISE_;

		bool displayIPOPT_;
		double tolIPOPT_;
		double acceptableTolIPOPT_;
		size_t maxIterationIPOPT_;
		double minAcceptedSwitchingTimeDifference_;
	};

	struct MP_Options{
		MP_Options()
		{
			nThreads_ = 1;
			debugPrintMP_ = false;
		}
		size_t nThreads_;
		bool debugPrintMP_;
	};

private:

};

using HyQDimensions = Dimensions<36,12>;

} // namespace ocs2

#endif /* DIMENSIONS_H_ */
