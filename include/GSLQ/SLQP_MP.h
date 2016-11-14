/*
 * SLQP_MP.h
 *
 * Multicore version of SLQP
 *
 *  Created on: Jun 20, 2016
 *      Author: mgiftthaler
 */


#ifndef SLQP_MP_OCS2_H_
#define SLQP_MP_OCS2_H_

#include <memory>
#include <thread>
#include <atomic>
#include <mutex>
#include <condition_variable>

#include "Dimensions.h"

#include "integration/KillIntegrationEventHandler.h"

#include "GSLQ/SLQP_BASE.h"

namespace ocs2{

template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
class SLQP_MP : public SLQP_BASE<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>
{
public:

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	typedef SLQP_BASE<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS> BASE;

	typedef SequentialRiccatiEquations<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS> RiccatiEquations_t;
	typedef SequentialErrorEquation<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS> ErrorEquation_t;

	typedef Dimensions<STATE_DIM, INPUT_DIM, OUTPUT_DIM> DIMENSIONS;

	typedef typename DIMENSIONS::template LinearFunction_t<Eigen::Dynamic> lagrange_t;

	typedef typename DIMENSIONS::controller_t controller_t;
	typedef typename DIMENSIONS::Options Options_t;
	typedef typename DIMENSIONS::scalar_t scalar_t;
	typedef typename DIMENSIONS::scalar_array_t scalar_array_t;
	typedef typename DIMENSIONS::eigen_scalar_t eigen_scalar_t;
	typedef typename DIMENSIONS::eigen_scalar_array_t eigen_scalar_array_t;
	typedef typename DIMENSIONS::state_vector_t state_vector_t;
	typedef typename DIMENSIONS::state_vector_array_t state_vector_array_t;
	typedef typename DIMENSIONS::output_vector_t output_vector_t;
	typedef typename DIMENSIONS::output_vector_array_t output_vector_array_t;
	typedef typename DIMENSIONS::control_vector_t control_vector_t;
	typedef typename DIMENSIONS::control_vector_array_t control_vector_array_t;
	typedef typename DIMENSIONS::control_feedback_t control_feedback_t;
	typedef typename DIMENSIONS::control_feedback_array_t control_feedback_array_t;
	typedef typename DIMENSIONS::state_matrix_t state_matrix_t;
	typedef typename DIMENSIONS::state_matrix_array_t state_matrix_array_t;
	typedef typename DIMENSIONS::control_matrix_t control_matrix_t;
	typedef typename DIMENSIONS::control_matrix_array_t control_matrix_array_t;
	typedef typename DIMENSIONS::control_gain_matrix_t control_gain_matrix_t;
	typedef typename DIMENSIONS::control_gain_matrix_array_t control_gain_matrix_array_t;
	typedef typename DIMENSIONS::constraint1_vector_t constraint1_vector_t;
	typedef typename DIMENSIONS::constraint1_vector_array_t constraint1_vector_array_t;
	typedef typename DIMENSIONS::constraint1_state_matrix_t constraint1_state_matrix_t;
	typedef typename DIMENSIONS::constraint1_state_matrix_array_t constraint1_state_matrix_array_t;
	typedef typename DIMENSIONS::constraint1_control_matrix_t constraint1_control_matrix_t;
	typedef typename DIMENSIONS::constraint1_control_matrix_array_t constraint1_control_matrix_array_t;
	typedef typename DIMENSIONS::control_constraint1_matrix_t control_constraint1_matrix_t;
	typedef typename DIMENSIONS::control_constraint1_matrix_array_t control_constraint1_matrix_array_t;
	typedef typename DIMENSIONS::constraint2_vector_t       constraint2_vector_t;
	typedef typename DIMENSIONS::constraint2_vector_array_t constraint2_vector_array_t;
	typedef typename DIMENSIONS::constraint2_state_matrix_t       constraint2_state_matrix_t;
	typedef typename DIMENSIONS::constraint2_state_matrix_array_t constraint2_state_matrix_array_t;


	enum WORKER_STATE {
		IDLE = 0,
		LINE_SEARCH,
		APPROXIMATE_LQ,
		CALCULATE_CONTROLLER_AND_LAGRANGIAN,
		SHUTDOWN
	};

	SLQP_MP(const std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > >& subsystemDynamicsPtr,
			const std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > >& subsystemDerivativesPtr,
			const std::vector<std::shared_ptr<CostFunctionBaseOCS2<OUTPUT_DIM, INPUT_DIM> > >& subsystemCostFunctionsPtr,
			const std::vector<controller_t>& initialControllersStock,
			const std::vector<size_t>& systemStockIndex,
			const Options_t& options = Options_t())

	:
		SLQP_BASE<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>(options),
		feedForwardConstraintInputStock_(NUM_SUBSYSTEMS),
		workerTask_(IDLE),
		subsystemProcessed_(0),
		learningRateStar_(1.0),
		KMax_subsystem_approx_(NUM_SUBSYSTEMS),
		KMax_subsystem_ctrl_(NUM_SUBSYSTEMS),
		kTaken_approx_(NUM_SUBSYSTEMS),
		kCompleted_approx_(NUM_SUBSYSTEMS),
		kTaken_ctrl_(NUM_SUBSYSTEMS),
		kCompleted_ctrl_(NUM_SUBSYSTEMS)
	{
		Eigen::initParallel();
//		std::cout << "initialized Eigen Parallel with " << Eigen::nbThreads( ) << " threads. " << std::endl;

		// resize instances to correct number of threads + 1
		dynamics_.resize(options.nThreads_+1);
		linearizedSystems_.resize(options.nThreads_+1);
		costFunctions_.resize(options.nThreads_+1);
		integratorsODE45_.resize(options.nThreads_+1);
		killIntegrationEventHandlers_.resize(options.nThreads_+1);
		controllers_.resize(options.nThreads_+1);

		killIntegrationEventHandler_ = std::make_shared<KillIntegrationEventHandler<STATE_DIM>> ();

		// for all threads + 1
		for (size_t i=0; i<options.nThreads_+1; i++)
		{
			// .. initialize all subsystems, etc.
			for(size_t j = 0; j<NUM_SUBSYSTEMS; j++)
			{
				// initialize dynamics
				dynamics_[i].push_back(subsystemDynamicsPtr[j]->clone());

				// initialize linearized systems
				linearizedSystems_[i].push_back(subsystemDerivativesPtr[j]->clone());

				// initialize cost functions
				costFunctions_[i].push_back(subsystemCostFunctionsPtr[j]->clone());

				// initialize controllers
				controllers_[i].push_back(initialControllersStock[j]);

				// initialize integrators
				killIntegrationEventHandlers_[i].push_back(std::make_shared<KillIntegrationEventHandler<STATE_DIM>> ());
				integratorsODE45_[i].push_back(std::shared_ptr<ODE45<STATE_DIM>>(new ODE45<STATE_DIM> ((dynamics_[i]).back(), killIntegrationEventHandler_)));
			}
		}

		// check first elements of all instances for correct number of subsystems
		assert(dynamics_[0].size() == subsystemDynamicsPtr.size() && "wrong number of subsystems in dynamics");
		assert(linearizedSystems_[0].size() == subsystemDerivativesPtr.size() && "wrong number of subsystems in derivatives");
		assert(costFunctions_[0].size() == subsystemCostFunctionsPtr.size() && "wrong number of subsystems in cost functions");
		assert(controllers_[0].size() == initialControllersStock.size() && "wrong number of subsystems in controller");

		// additional checks on stock index at runtime
		if (subsystemDynamicsPtr.size()-1 < *std::max_element(systemStockIndex.begin(), systemStockIndex.end()))
			throw std::runtime_error("systemStockIndex points to non-existing subsystem");
		if (systemStockIndex.size() != NUM_SUBSYSTEMS)
			throw std::runtime_error("systemStockIndex has less elements then the number of subsystems");

		// for controller design
		nominalOutputFunc_.resize(options.nThreads_+1);
		nominalInputFunc_.resize(options.nThreads_+1);
		BmFunc_.resize(options.nThreads_+1);
		PmFunc_.resize(options.nThreads_+1);
		RmInverseFunc_.resize(options.nThreads_+1);
		RvFunc_.resize(options.nThreads_+1);
		EvProjectedFunc_.resize(options.nThreads_+1);
		CmProjectedFunc_.resize(options.nThreads_+1);
		DmProjectedFunc_.resize(options.nThreads_+1);
		nominalLagrangeMultiplierFunc_.resize(options.nThreads_+1);
		DmDagerFunc_.resize(options.nThreads_+1);
		RmFunc_.resize(options.nThreads_+1);

		// initialize threads
		launchWorkerThreads();
	}

	~SLQP_MP();

	// only for interfacing with GSLQP
	void rollout(const state_vector_t& initState,
			const std::vector<controller_t>& controllersStock,
			std::vector<scalar_array_t>& timeTrajectoriesStock,
			std::vector<state_vector_array_t>& stateTrajectoriesStock,
			std::vector<control_vector_array_t>& inputTrajectoriesStock,
			std::vector<output_vector_array_t>& outputTrajectoriesStock,
			std::vector<std::vector<size_t> >& nc1TrajectoriesStock,
			std::vector<constraint1_vector_array_t>& EvTrajectoryStock,
			std::vector<std::vector<size_t> >& nc2TrajectoriesStock,
			std::vector<constraint2_vector_array_t>& HvTrajectoryStock,
			std::vector<size_t>& nc2FinalStock,
			std::vector<constraint2_vector_t>& HvFinalStock) override;

	// only for interfacing with GSLQP
	void rollout(const state_vector_t& initState,
			const std::vector<controller_t>& controllersStock,
			std::vector<scalar_array_t>& timeTrajectoriesStock,
			std::vector<state_vector_array_t>& stateTrajectoriesStock,
			std::vector<control_vector_array_t>& inputTrajectoriesStock,
			std::vector<output_vector_array_t>& outputTrajectoriesStock) override;

	// only for interfacing with GSLQP
	void rollout(const state_vector_t& initState,
			const std::vector<controller_t>& controllersStock,
			std::vector<scalar_array_t>& timeTrajectoriesStock,
			std::vector<state_vector_array_t>& stateTrajectoriesStock,
			std::vector<control_vector_array_t>& inputTrajectoriesStock) override;

	void rollout(
			const size_t& threadId,
			const state_vector_t& initState,
			const std::vector<controller_t>& controllersStock,
			std::vector<scalar_array_t>& timeTrajectoriesStock,
			std::vector<state_vector_array_t>& stateTrajectoriesStock,
			std::vector<control_vector_array_t>& inputTrajectoriesStock,
			std::vector<output_vector_array_t>& outputTrajectoriesStock,
			std::vector<std::vector<size_t> >& nc1TrajectoriesStock,
			std::vector<constraint1_vector_array_t>& EvTrajectoryStock,
			std::vector<std::vector<size_t> >& nc2TrajectoriesStock,
			std::vector<constraint2_vector_array_t>& HvTrajectoryStock,
			std::vector<size_t>& nc2FinalStock,
			std::vector<constraint2_vector_t>& HvFinalStock);

	void rollout(
			const size_t& threadId,
			const state_vector_t& initState,
			const std::vector<controller_t>& controllersStock,
			std::vector<scalar_array_t>& timeTrajectoriesStock,
			std::vector<state_vector_array_t>& stateTrajectoriesStock,
			std::vector<control_vector_array_t>& inputTrajectoriesStock,
			std::vector<output_vector_array_t>& outputTrajectoriesStock);

	void rollout(
			const size_t& threadId,
			const state_vector_t& initState,
			const std::vector<controller_t>& controllersStock,
			std::vector<scalar_array_t>& timeTrajectoriesStock,
			std::vector<state_vector_array_t>& stateTrajectoriesStock,
			std::vector<control_vector_array_t>& inputTrajectoriesStock);

	void rollout(const state_vector_t& initState,
			const std::vector<controller_t>& controllersStock,
			const double& stoppingTime,
			state_vector_t& stateVectorWhereStopped,
			control_vector_t& controlInputWhereStopped,
			output_vector_t& outputWhereStopped,
			size_t& numSubsystemWhereStopped) override;

	void calculateCostFunction(
			const std::vector<scalar_array_t>& timeTrajectoriesStock,
			const std::vector<output_vector_array_t>& stateTrajectoriesStock,
			const std::vector<control_vector_array_t>& inputTrajectoriesStock,
			scalar_t& totalCost,
			size_t threadId);

	// for interfacing with GSLQP
	void calculateCostFunction(
			const std::vector<scalar_array_t>& timeTrajectoriesStock,
			const std::vector<output_vector_array_t>& stateTrajectoriesStock,
			const std::vector<control_vector_array_t>& inputTrajectoriesStock,
			scalar_t& totalCost) override;

	void calculateCostFunction(const std::vector<scalar_array_t>& timeTrajectoriesStock,
			const std::vector<output_vector_array_t>& stateTrajectoriesStock,
			const std::vector<control_vector_array_t>& inputTrajectoriesStock,
			const std::vector<std::vector<size_t> >& nc2TrajectoriesStock,
			const std::vector<constraint2_vector_array_t>& HvTrajectoryStock,
			const std::vector<size_t>& nc2FinalStock,
			const std::vector<constraint2_vector_t>& HvFinalStock,
			scalar_t& totalCost) override;

	void calculateCostFunction(const std::vector<scalar_array_t>& timeTrajectoriesStock,
				const std::vector<output_vector_array_t>& stateTrajectoriesStock,
				const std::vector<control_vector_array_t>& inputTrajectoriesStock,
				const std::vector<std::vector<size_t> >& nc2TrajectoriesStock,
				const std::vector<constraint2_vector_array_t>& HvTrajectoryStock,
				const std::vector<size_t>& nc2FinalStock,
				const std::vector<constraint2_vector_t>& HvFinalStock,
				scalar_t& totalCost,
				size_t threadId);

	void calculateMeritFunction(const std::vector<scalar_array_t>& timeTrajectoriesStock,
			const std::vector<std::vector<size_t> >& nc1TrajectoriesStock,
			const std::vector<constraint1_vector_array_t>& EvTrajectoryStock,
			const std::vector<std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd> > >&  lagrangeTrajectoriesStock,
			const scalar_t& totalCost,
			scalar_t& meritFuntionValue,
			scalar_t& constraintISE) override;

	double calculateConstraintISE(const std::vector<scalar_array_t>& timeTrajectoriesStock,
			const std::vector<std::vector<size_t>>& nc1TrajectoriesStock,
			const std::vector<constraint1_vector_array_t>& EvTrajectoriesStock,
			scalar_t& constraintISE) override;

	void getController(std::vector<controller_t>& controllersStock) override;

	void setController(const std::vector<controller_t>& controllersStock) override;

	void getValueFuntion(const scalar_t& time, const output_vector_t& output, scalar_t& valueFuntion) override;

	void getCostFuntion(scalar_t& costFunction, scalar_t& constriantISE) override;

	void getNominalTrajectories(
			std::vector<scalar_array_t>& nominalTimeTrajectoriesStock,
			std::vector<state_vector_array_t>& nominalStateTrajectoriesStock,
			std::vector<control_vector_array_t>& nominalInputTrajectoriesStock,
			std::vector<output_vector_array_t>& nominalOutputTrajectoriesStock) override;

	void run(const state_vector_t& initState, const std::vector<scalar_t>& switchingTimes,
			const std::vector<controller_t>& initialControllersStock=std::vector<controller_t>()) override;


	// get subsystem dynamics on main thread
	std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM>>>& getSubsystemDynamicsPtrStock() override {
		return dynamics_[BASE::options_.nThreads_];
	}


	// allows to externally update the cost functions reference states
	void setNewCostReferenceState(const output_vector_t& newReference);


protected:
	void approximateOptimalControlProblem();

	void calculateControllerAndLagrangian(std::vector<controller_t>& controllersStock,
			std::vector<lagrange_t>& lagrangeMultiplierFunctionsStock,
			std::vector<control_vector_array_t>& feedForwardConstraintInputStock,
			bool firstCall=true);

	void calculateRolloutLagrangeMultiplier(const std::vector<scalar_array_t>& timeTrajectoriesStock,
			const std::vector<output_vector_array_t>& outputTrajectoriesStock,
			const std::vector<lagrange_t>& lagrangeMultiplierFunctionsStock,
			std::vector<std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd> > >&  lagrangeTrajectoriesStock);

	void calculateRolloutCostate(
			const std::vector<scalar_array_t>& timeTrajectoriesStock,
			const std::vector<output_vector_array_t>& outputTrajectoriesStock,
			std::vector<output_vector_array_t>& costateTrajectoriesStock);

	void lineSearch();

private:

	void launchWorkerThreads();
	void threadWork(size_t threadId);

	// main function for sub-tasks
	void approximateSubsystemLQ(const size_t sysIndex); // computes the linearized dynamics for a particular subsystem
	void calculateControllerAndLagrangian();

	// worker functions
	size_t approximateSubsystemLQWorker(size_t threadId, size_t subsystemProcessed);		//Worker functions
	size_t calculateControllerAndLagrangianWorker(size_t threadId, size_t subsystemProcessed);
	void lineSearchWorker(size_t threadId);

	// execute methods
	size_t executeApproximateSubsystemLQ(size_t threadId, size_t k, size_t subsystemProcessed);	// Computes the linearized dynamics
	size_t executeCalculateControllerAndLagrangian(size_t threadId, size_t k, size_t subsystemProcessed);
	void executeLineSearch(
			size_t threadId, double learningRate,
			scalar_t& lsTotalCost,
			scalar_t& lsTotalMerit,
			scalar_t& lsConstraint1ISE,
			std::vector<controller_t>& lsControllersStock,
			std::vector<lagrange_t>& lsLagrangeControllersStock,
			std::vector<scalar_array_t>& lsTimeTrajectoriesStock,
			std::vector<state_vector_array_t>& lsStateTrajectoriesStock,
			std::vector<control_vector_array_t>& lsInputTrajectoriesStock,
			std::vector<output_vector_array_t>& lsOutputTrajectoriesStock,
			std::vector<std::vector<size_t> >& lsNc1TrajectoriesStock,
			std::vector<constraint1_vector_array_t>& lsEvTrajectoryStock,
			std::vector<std::vector<size_t> >&   lsNc2TrajectoriesStock,
			std::vector<constraint2_vector_array_t>& lsHvTrajectoryStock,
			std::vector<size_t>&               	lsNc2FinalStock,
			std::vector<constraint2_vector_t>& 	lsHvFinalStock,
			std::vector<std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd> >>& lsLagrangeTrajectoriesStock);



	/*a heuristic that generates a unique id for a process, such that we can manage the tasks.
	 * Generates a unique identifiers for subsystem, task, iteration:
	 * */
	size_t generateUniqueProcessID (const size_t& iterateNo, const int workerState, const int subsystemId)
	{
		return (10e9*(workerState +1) + 10e6 * (subsystemId +1) + iterateNo+1);
	}

	// for nice debug printing
	inline void printString(const std::string& text){std::cout << text << std::endl;}

	std::vector<std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > > > dynamics_;
	std::vector<std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > > > linearizedSystems_;
	std::vector<std::vector<std::shared_ptr<CostFunctionBaseOCS2<OUTPUT_DIM, INPUT_DIM> > > > costFunctions_;

	std::vector<std::vector<std::shared_ptr<ODE45<STATE_DIM> > > > integratorsODE45_;
	std::vector<std::vector<std::shared_ptr<KillIntegrationEventHandler<STATE_DIM>>>> killIntegrationEventHandlers_;
	std::shared_ptr<KillIntegrationEventHandler<STATE_DIM>> killIntegrationEventHandler_;

	std::vector<std::vector<controller_t> > controllers_;

	std::vector<control_vector_array_t> feedForwardConstraintInputStock_;


	// multithreading helper variables
	std::vector<std::thread> workerThreads_;
	std::atomic_bool workersActive_;
	std::atomic_int workerTask_;
	std::atomic_int subsystemProcessed_;		// the subsystem the threads are currently working on

	std::mutex workerWakeUpMutex_;
	std::condition_variable workerWakeUpCondition_;
	std::mutex kCompletedMutex_;
	std::condition_variable kCompletedCondition_;
	std::mutex lineSearchResultMutex_;
	std::mutex alphaBestFoundMutex_;
	std::condition_variable alphaBestFoundCondition_;

	double learningRateStar_;

	std::vector<size_t> KMax_subsystem_approx_;		// denotes the number of integration steps for a particular subsystem i
	std::vector<size_t> KMax_subsystem_ctrl_;
	std::vector<std::atomic_size_t> kTaken_approx_;
	std::vector<std::atomic_size_t> kCompleted_approx_;
	std::vector<std::atomic_size_t> kTaken_ctrl_;
	std::vector<std::atomic_size_t> kCompleted_ctrl_;

	std::atomic_size_t alphaTaken_;
	size_t alphaMax_;
	size_t alphaExpBest_;
	size_t alphaExpMax_;
	std::atomic_bool alphaBestFound_;
	std::vector<size_t> alphaProcessed_;
	double lowestTotalMerit_;
	double lowestTotalCost_;
	double lowestConstraint1ISE_;
	std::atomic_size_t lsWorkerCompleted_;

	// for controller design
	// functions for controller and lagrange multiplier
	std::vector<LinearInterpolation<output_vector_t,Eigen::aligned_allocator<output_vector_t> >>   	nominalOutputFunc_;
	std::vector<LinearInterpolation<control_vector_t,Eigen::aligned_allocator<control_vector_t> >> 	nominalInputFunc_;
	std::vector<LinearInterpolation<control_gain_matrix_t,Eigen::aligned_allocator<control_gain_matrix_t> >> BmFunc_;
	std::vector<LinearInterpolation<control_feedback_t,Eigen::aligned_allocator<control_feedback_t> >> PmFunc_;
	std::vector<LinearInterpolation<control_matrix_t,Eigen::aligned_allocator<control_matrix_t> >>     RmInverseFunc_;
	std::vector<LinearInterpolation<control_vector_t,Eigen::aligned_allocator<control_vector_t> >>     RvFunc_;
	std::vector<LinearInterpolation<control_vector_t,Eigen::aligned_allocator<control_vector_t> >>     EvProjectedFunc_;
	std::vector<LinearInterpolation<control_feedback_t,Eigen::aligned_allocator<control_feedback_t> >> CmProjectedFunc_;
	std::vector<LinearInterpolation<control_matrix_t,Eigen::aligned_allocator<control_matrix_t> >>     DmProjectedFunc_;

	// functions for lagrange multiplier only
	std::vector<LinearInterpolation<Eigen::VectorXd,Eigen::aligned_allocator<Eigen::VectorXd> >> nominalLagrangeMultiplierFunc_;
	std::vector<LinearInterpolation<control_constraint1_matrix_t,Eigen::aligned_allocator<control_constraint1_matrix_t> >> DmDagerFunc_;
	std::vector<LinearInterpolation<control_matrix_t,Eigen::aligned_allocator<control_matrix_t> >> RmFunc_;

	bool nominalLagrangeMultiplierUpdated_;

	// needed for lineSearch
	std::vector<controller_t> initLScontrollersStock_;
	std::vector<lagrange_t> initLSlagrangeMultiplierFunctionsStock_;

public:
	template <size_t GSLQP_STATE_DIM, size_t GSLQP_INPUT_DIM, size_t GSLQP_OUTPUT_DIM, size_t GSLQP_NUM_SUBSYSTEMS>
	friend class GSLQP;
};

} // namespace ocs2

#include "implementation/SLQP_MP.h"


#endif /* SLQP_MP_H_ */
