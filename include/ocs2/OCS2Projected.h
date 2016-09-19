/*
 * OCS2Projected.h
 *
 *  Created on: Jul 21, 2016
 *      Author: farbod
 */

#ifndef OCS2_OCS2PROJECTED_H_
#define OCS2_OCS2PROJECTED_H_

#include <vector>
#include <array>
#include <memory>
#include <iterator>
#include <algorithm>

#include <c_gradient_descent/GradientDescent.h>

#include "GSLQ/GLQP.h"
#include "GSLQ/GSLQP.h"


namespace ocs2{


template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
class OCS2Projected : private nlp::GradientDescent
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	typedef std::shared_ptr<OCS2Projected<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS> > Ptr;

	enum {NumConstraints_=NUM_SUBSYSTEMS-2};

	typedef SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>  slqp_t;
	typedef GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS> gslqp_t;
	typedef typename slqp_t::Ptr  slqp_ptr_t;
	typedef typename gslqp_t::Ptr gslqp_ptr_t;

	typedef Dimensions<STATE_DIM, INPUT_DIM, OUTPUT_DIM> DIMENSIONS;
	typedef typename DIMENSIONS::controller_t controller_t;
	typedef typename DIMENSIONS::Options Options_t;
	typedef typename DIMENSIONS::scalar_t 		scalar_t;
	typedef typename DIMENSIONS::scalar_array_t scalar_array_t;
	typedef typename DIMENSIONS::eigen_scalar_t       eigen_scalar_t;
	typedef typename DIMENSIONS::eigen_scalar_array_t eigen_scalar_array_t;
	typedef typename DIMENSIONS::state_vector_t 	  state_vector_t;
	typedef typename DIMENSIONS::state_vector_array_t state_vector_array_t;
	typedef typename DIMENSIONS::control_vector_t 		control_vector_t;
	typedef typename DIMENSIONS::control_vector_array_t control_vector_array_t;
	typedef typename DIMENSIONS::output_vector_t 	   output_vector_t;
	typedef typename DIMENSIONS::output_vector_array_t output_vector_array_t;
	typedef typename DIMENSIONS::control_feedback_t 	  control_feedback_t;
	typedef typename DIMENSIONS::control_feedback_array_t control_feedback_array_t;
	typedef typename DIMENSIONS::state_matrix_t 	  state_matrix_t;
	typedef typename DIMENSIONS::state_matrix_array_t state_matrix_array_t;
	typedef typename DIMENSIONS::control_matrix_t 		control_matrix_t;
	typedef typename DIMENSIONS::control_matrix_array_t control_matrix_array_t;
	typedef typename DIMENSIONS::control_gain_matrix_t 		 control_gain_matrix_t;
	typedef typename DIMENSIONS::control_gain_matrix_array_t control_gain_matrix_array_t;

	OCS2Projected(const std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > >& subsystemDynamicsPtr,
			const std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > >& subsystemDerivativesPtr,
			const std::vector<std::shared_ptr<CostFunctionBaseOCS2<STATE_DIM, INPUT_DIM> > >& subsystemCostFunctionsPtr,
			const state_vector_array_t&   stateOperatingPoints,
			const control_vector_array_t& inputOperatingPoints,
			const std::vector<size_t>& systemStockIndex,
			const Options_t& options = Options_t::Options() )

	: subsystemDynamicsPtr_(subsystemDynamicsPtr),
	  subsystemDerivativesPtr_(subsystemDerivativesPtr),
	  subsystemCostFunctionsPtr_(subsystemCostFunctionsPtr),
	  stateOperatingPoints_(stateOperatingPoints),
	  inputOperatingPoints_(inputOperatingPoints),
	  systemStockIndex_(systemStockIndex),
	  options_(options),
	  subsystemDynamicsPtrStock_(NUM_SUBSYSTEMS),
	  subsystemCostFunctionsPtrStock_(NUM_SUBSYSTEMS),
	  subsystemSimulatorsStockPtr_(NUM_SUBSYSTEMS),
	  initSwitchingTimes_(NUM_SUBSYSTEMS+1),
	  optimizedSwitchingTimes_(NUM_SUBSYSTEMS+1),
	  optimizedControllersStock_(NUM_SUBSYSTEMS)
	{
		// NLP optimizer options
		nlpOptions_.displayGradientDescent_ = options_.displayIPOPT_;
		nlpOptions_.maxIterations_ 	  = options_.maxIterationIPOPT_;
		nlpOptions_.minRelCost_    	  = options_.acceptableTolIPOPT_;
		nlpOptions_.maxLearningRate_  = options_.maxLearningRateNLP_;
		nlpOptions_.minLearningRate_  = options_.minLearningRateNLP_;
		nlpOptions_.minDisToBoundary_ = options_.minAcceptedSwitchingTimeDifference_;
		nlpOptions_.useAscendingLineSearchNLP_ = options_.useAscendingLineSearchNLP_;
		adjustOptions();

		// setting up subsystemSimulatorsStockPtr
		if (subsystemDynamicsPtr.size() != subsystemCostFunctionsPtr.size())
			throw std::runtime_error("Number of cost functions is not equal to the number of subsystems.");
		if (subsystemDynamicsPtr.size()-1 < *std::max_element(systemStockIndex.begin(), systemStockIndex.end()))
			throw std::runtime_error("systemStockIndex points to non-existing subsystem");
		if (systemStockIndex.size() != NUM_SUBSYSTEMS)
			throw std::runtime_error("systemStockIndex has less elements than the number of subsystems");
		for (int i=0; i<NUM_SUBSYSTEMS; i++) {
			subsystemDynamicsPtrStock_[i] = subsystemDynamicsPtr[systemStockIndex[i]]->clone();
			subsystemCostFunctionsPtrStock_[i] = subsystemCostFunctionsPtr[systemStockIndex[i]]->clone();
			subsystemSimulatorsStockPtr_[i] = std::shared_ptr<ODE45<STATE_DIM>>( new ODE45<STATE_DIM>(subsystemDynamicsPtrStock_[i]) );
		}  // end of i loop

	}

	~OCS2Projected() {}

	void getCostFunction(scalar_t& costFunction) const;

	void getCostFunctionDerivative(Eigen::Matrix<double,NUM_SUBSYSTEMS-1,1>& costFuntionDerivative) const;

	void getController(std::vector<controller_t>& controllersStock) const;

	void getSwitchingTimes(scalar_array_t& switchingTimes) const;

	void getTrajectories(std::vector<scalar_array_t>& nominalTimeTrajectoriesStock,
			std::vector<state_vector_array_t>& nominalStateTrajectoriesStock,
			std::vector<control_vector_array_t>& nominalInputTrajectoriesStock,
			std::vector<output_vector_array_t>& nominalOutputTrajectoriesStock) const;

	void getTrajectories(std::vector<scalar_array_t>& nominalTimeTrajectoriesStock,
			std::vector<state_vector_array_t>& nominalStateTrajectoriesStock,
			std::vector<control_vector_array_t>& nominalInputTrajectoriesStock) const;

	void getSLQIterationsLog(eigen_scalar_array_t& slqIterationCost, eigen_scalar_array_t& slqIterationISE1) const {
		slqIterationCost = slqIterationCost_;
		slqIterationISE1 = slqIterationISE1_;
	}

	void getOCS2IterationsLog(eigen_scalar_array_t& iterationCost) const { /*iterationCost = iterationCost_; */ // ISSUE FIXME TODO
	}

	void rollout(const state_vector_t& initState,
			const scalar_array_t& switchingTimes,
			const std::vector<controller_t>& controllersStock,
			std::vector<scalar_array_t>& timeTrajectoriesStock,
			std::vector<state_vector_array_t>& stateTrajectoriesStock,
			std::vector<control_vector_array_t>& inputTrajectoriesStock,
			std::vector<output_vector_array_t>& outputTrajectoriesStock);

	void calculateCostFunction(const std::vector<scalar_array_t>& timeTrajectoriesStock,
			const std::vector<output_vector_array_t>& stateTrajectoriesStock,
			const std::vector<control_vector_array_t>& inputTrajectoriesStock,
			scalar_t& totalCost);

	void run(const state_vector_t& initState, const scalar_array_t& switchingTimes);

private:
	size_t findNearestController(const Eigen::VectorXd& enquiryParameter) const;

	void calculateLinearEqualityConstraint(Eigen::MatrixXd& Am, Eigen::VectorXd& Bv) override;

	bool calculateGradient(const size_t& id, const Eigen::VectorXd& parameters, Eigen::VectorXd& gradient) override;

	bool calculateCost(const size_t& id, const Eigen::VectorXd& parameters, double& cost) override;

	void getSolution(size_t idStar) override;

	void saveToBag(size_t id, const Eigen::VectorXd& parameters);

	void calculateInitialController(const state_vector_t& initState,
			const scalar_array_t& switchingTimes,
			std::vector<controller_t>&  controllersStock);


private:
	std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > > subsystemDynamicsPtr_;
	std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > > subsystemDerivativesPtr_;
	std::vector<std::shared_ptr<CostFunctionBaseOCS2<STATE_DIM, INPUT_DIM> > > subsystemCostFunctionsPtr_;

	state_vector_array_t   stateOperatingPoints_;
	control_vector_array_t inputOperatingPoints_;

	std::vector<size_t> systemStockIndex_;

	scalar_array_t initSwitchingTimes_;
	state_vector_t initState_;
	std::vector<controller_t> initControllersStock_;

	Options_t options_;

	// for rollout function
	std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > > subsystemDynamicsPtrStock_;
	std::vector<std::shared_ptr<CostFunctionBaseOCS2<OUTPUT_DIM, INPUT_DIM> > > subsystemCostFunctionsPtrStock_;
	std::vector<std::shared_ptr<ODE45<STATE_DIM> > > subsystemSimulatorsStockPtr_;

	// optimized solution variables
	scalar_t       optimizedTotalCost_;
	scalar_t       optimizedConstraintISE_;
	scalar_array_t optimizedSwitchingTimes_;
	std::vector<controller_t>  optimizedControllersStock_;
	std::vector<scalar_array_t> optimizedTimeTrajectoriesStock_;
	std::vector<state_vector_array_t>   optimizedStateTrajectoriesStock_;
	std::vector<control_vector_array_t> optimizedInputTrajectoriesStock_;
	std::vector<output_vector_array_t>  optimizedOutputTrajectoriesStock_;
	Eigen::Matrix<double,NUM_SUBSYSTEMS-1,1> costFuntionDerivative_;

	scalar_t currentTotalCost_;
	Eigen::VectorXd currentCostFuntionDerivative_;

	std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd> > parameterBag_;
	std::vector<std::vector<controller_t> > controllersStockBag_;

	gslqp_ptr_t gslqpSolver_;
	std::vector<slqp_ptr_t> slqpSolverPtrs_;

	eigen_scalar_array_t slqIterationCost_;
	eigen_scalar_array_t slqIterationISE1_;
};

}  // end of ocs2 namespace

#include "implementation/OCS2Projected.h"

#endif /* OCS2_OCS2PROJECTED_H_ */
