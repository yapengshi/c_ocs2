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
		MAX_CONSTRAINT2_DIM_ = INPUT_DIM
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

    typedef Eigen::Matrix<double, MAX_CONSTRAINT2_DIM_, 1> constraint2_vector_t;
    typedef std::vector<constraint2_vector_t, Eigen::aligned_allocator<constraint2_vector_t> > constraint2_vector_array_t;

    typedef Eigen::Matrix<double, MAX_CONSTRAINT2_DIM_, OUTPUT_DIM> constraint2_state_matrix_t;
    typedef std::vector<constraint2_state_matrix_t, Eigen::aligned_allocator<constraint2_state_matrix_t> > constraint2_state_matrix_array_t;

	typedef double scalar_t;
	typedef std::vector<scalar_t> scalar_array_t;

	typedef Eigen::Matrix<double, 1, 1> eigen_scalar_t;
	typedef std::vector<eigen_scalar_t, Eigen::aligned_allocator<eigen_scalar_t> > eigen_scalar_array_t;

	template <int DIM1, int DIM2=1>
	struct LinearFunction_t {
	public:
		scalar_array_t time_;
		std::vector<Eigen::Matrix<double, DIM1, DIM2>, Eigen::aligned_allocator<Eigen::Matrix<double, DIM1, DIM2>> > uff_;
		std::vector<Eigen::Matrix<double, DIM1, DIM2>, Eigen::aligned_allocator<Eigen::Matrix<double, DIM1, DIM2>> > deltaUff_;
		std::vector<Eigen::Matrix<double, DIM1, OUTPUT_DIM>, Eigen::aligned_allocator<Eigen::Matrix<double, DIM1, OUTPUT_DIM>> > k_;

		void swap(LinearFunction_t& arg){
			uff_.swap(arg.uff_);
			deltaUff_.swap(arg.deltaUff_);
			k_.swap(arg.k_);
		}

		void setZero(){
			std::fill(uff_.begin(), uff_.end(), Eigen::Matrix<double, DIM1, DIM2>::Zero());
			std::fill(deltaUff_.begin(), deltaUff_.end(), Eigen::Matrix<double, DIM1, DIM2>::Zero());
			std::fill(k_.begin(), k_.end(), Eigen::Matrix<double, DIM1, OUTPUT_DIM>::Zero());
		}
	};
	typedef LinearFunction_t<INPUT_DIM> controller_t;

	enum RICCATI_INTEGRATOR_TYPE{
		ODE45 = 1,
		ADAMS_BASHFORTH = 2,
		BULIRSCH_STOER = 3
	};

	struct Options {
	public:
		Options() :
			maxIterationGSLQP_(10),
			minLearningRateGSLQP_(0.05),
			maxLearningRateGSLQP_(1.0),
			lineSearchContractionRate_(0.5),
			minRelCostGSLQP_(1e-3),
			stateConstraintPenaltyCoeff_(0.0),
			stateConstraintPenaltyBase_(1.0),
			meritFunctionRho_(1.0),
			constraintStepSize_(1.0),
			lineSearchByMeritFuntion_(false),
			dispayGSLQP_(false),
			displayShortSummary_(false),
			warmStartGSLQP_(false),
			useLQForDerivatives_(false),

			AbsTolODE_(1e-9),
			RelTolODE_(1e-6),
			maxNumStepsPerSecond_(5000),
			simulationIsConstrained_(false),
			minSimulationTimeDuration_(1e-3),
			minAbsConstraint1ISE_(1e-3),
			minRelConstraint1ISE_(1e-3),

			displayGradientDescent_(true),
			tolGradientDescent_(1e-2),
			acceptableTolGradientDescent_(1e-1),
			maxIterationGradientDescent_(20),
			minLearningRateNLP_(0.05),
		    maxLearningRateNLP_(1.0),
		    useAscendingLineSearchNLP_(true),
			minAcceptedSwitchingTimeDifference_(0.0),

			RiccatiIntegratorType_(ODE45),
			adams_integrator_dt_(0.001),

			useMultiThreading_(false),
			nThreads_(4),
			debugPrintMP_(false),
			lsStepsizeGreedy_(true)
		{}

		void print()
		{
			std::cout << " #### ========================== Options ============================ ####" << std::endl;
			std::cout << "maxIterationGSLQP_                 " << maxIterationGSLQP_                  << std::endl;
			std::cout << "minLearningRateGSLQP_              " << minLearningRateGSLQP_               << std::endl;
			std::cout << "maxLearningRateGSLQP_              " << maxLearningRateGSLQP_               << std::endl;
			std::cout << "lineSearchContractionRate_         " << lineSearchContractionRate_          << std::endl;
			std::cout << "minRelCostGSLQP_                   " << minRelCostGSLQP_                    << std::endl;
			std::cout << "meritFunctionRho_                  " << meritFunctionRho_                   << std::endl;
			std::cout << "constraintStepSize_                " << constraintStepSize_                 << std::endl;
			std::cout << "lineSearchByMeritFuntion_          " << lineSearchByMeritFuntion_           << std::endl;
			std::cout << "dispayGSLQP_                       " << dispayGSLQP_                        << std::endl;
			std::cout << "displayShortSummary_               " << displayShortSummary_                << std::endl;
			std::cout << "warmStartGSLQP_                    " << warmStartGSLQP_                     << std::endl;
            std::cout << "                                   " << ""                                  << std::endl;
			std::cout << "AbsTolODE_                         " << AbsTolODE_                          << std::endl;
			std::cout << "RelTolODE_                         " << RelTolODE_                          << std::endl;
			std::cout << "simulationIsConstrained_           " << simulationIsConstrained_            << std::endl;
			std::cout << "minAbsConstraint1ISE_              " << minAbsConstraint1ISE_               << std::endl;
			std::cout << "minRelConstraint1ISE_              " << minRelConstraint1ISE_               << std::endl;
            std::cout << "                                   " << ""                                  << std::endl;
			std::cout << "displayGradientDescent_            " << displayGradientDescent_             << std::endl;
			std::cout << "tolGradientDescent_                " << tolGradientDescent_                 << std::endl;
			std::cout << "acceptableTolGradientDescent_      " << acceptableTolGradientDescent_       << std::endl;
			std::cout << "maxIterationGradientDescent_       " << maxIterationGradientDescent_        << std::endl;
			std::cout << "minAcceptedSwitchingTimeDifference_" << minAcceptedSwitchingTimeDifference_ << std::endl;
            std::cout << "                                   " << ""                                  << std::endl;
			std::cout << "RiccatiIntegratorType_             " << RiccatiIntegratorType_              << std::endl;
			std::cout << "adams_integrator_dt_               " << adams_integrator_dt_                << std::endl;
            std::cout << "                                   " << ""                                  << std::endl;
			std::cout << "useMultiThreading_                 " << useMultiThreading_                  << std::endl;
			std::cout << "nThreads_                 		 " << nThreads_                			  << std::endl;
			std::cout << "debugPrintMP_		                 " << debugPrintMP_		                  << std::endl;
			std::cout << "lsStepsizeGreedy_                  " << lsStepsizeGreedy_                   << std::endl;
			std::cout << " #### ============================ end ============================== ####" << std::endl;
			std::cout << std::endl;
		}

		size_t maxIterationGSLQP_;
		double minLearningRateGSLQP_;
		double maxLearningRateGSLQP_;
		double lineSearchContractionRate_;
		double minRelCostGSLQP_;
		double stateConstraintPenaltyCoeff_;
		double stateConstraintPenaltyBase_;
		double meritFunctionRho_;
		double constraintStepSize_;
		bool lineSearchByMeritFuntion_;
		bool dispayGSLQP_;
		bool displayShortSummary_;
		bool warmStartGSLQP_;
		bool useLQForDerivatives_;

		double AbsTolODE_;
		double RelTolODE_;
		size_t maxNumStepsPerSecond_;
		bool simulationIsConstrained_;
		double minSimulationTimeDuration_;
		double minAbsConstraint1ISE_;
		double minRelConstraint1ISE_;

		bool displayGradientDescent_;
		double tolGradientDescent_;
		double acceptableTolGradientDescent_;
		size_t maxIterationGradientDescent_;
		double minLearningRateNLP_;
		double maxLearningRateNLP_;
		bool useAscendingLineSearchNLP_;
		double minAcceptedSwitchingTimeDifference_;

		size_t RiccatiIntegratorType_;
		double adams_integrator_dt_;

		bool useMultiThreading_;
		size_t nThreads_;
		bool debugPrintMP_;
		//mp line search options
		bool lsStepsizeGreedy_;	// otherwise it's merit-greedy
	};
};

using HyQDimensions = Dimensions<36,12>;

} // namespace ocs2

#endif /* DIMENSIONS_H_ */
