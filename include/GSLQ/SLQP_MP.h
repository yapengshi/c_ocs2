/*
 * SLQP_MP.h
 *
 * Multicore version of SLQP
 *
 *  Created on: Jun 20, 2016
 *      Author: markus, farbod
 */



#ifndef SLQP_MP_OCS2_H_
#define SLQP_MP_OCS2_H_

#include <iostream>
#include <memory>
#include <thread>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <vector>
#include <array>
#include <algorithm>
#include <cstddef>
#include <Eigen/Dense>
#include <Eigen/StdVector>

#include "Dimensions.h"

#include "dynamics/ControlledSystemBase.h"
#include "dynamics/DerivativesBase.h"
#include "costs/CostFunctionBaseOCS2.h"

#include "integration/Integrator.h"
#include "misc/LinearInterpolation.h"

#include "GSLQ/SequentialRiccatiEquations.h"
#include "GSLQ/SequentialErrorEquation.h"

namespace ocs2{

template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
class SLQP_MP
{
public:

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	typedef SequentialRiccatiEquations<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS> RiccatiEquations_t;
	typedef SequentialErrorEquation<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS> ErrorEquation_t;

	typedef Dimensions<STATE_DIM, INPUT_DIM, OUTPUT_DIM> DIMENSIONS;

	typedef typename DIMENSIONS::template LinearFunction_t<Eigen::Dynamic> lagrange_t;

	typedef typename DIMENSIONS::controller_t controller_t;
	typedef typename DIMENSIONS::Options Options_t;
	typedef typename DIMENSIONS::MP_Options MP_Options_t;
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
	typedef typename DIMENSIONS::constraint1_matrix_t constraint1_matrix_t;
	typedef typename DIMENSIONS::constraint1_matrix_array_t constraint1_matrix_array_t;
	typedef typename DIMENSIONS::constraint1_state_matrix_t constraint1_state_matrix_t;
	typedef typename DIMENSIONS::constraint1_state_matrix_array_t constraint1_state_matrix_array_t;
	typedef typename DIMENSIONS::constraint1_control_matrix_t constraint1_control_matrix_t;
	typedef typename DIMENSIONS::constraint1_control_matrix_array_t constraint1_control_matrix_array_t;
	typedef typename DIMENSIONS::control_constraint1_matrix_t control_constraint1_matrix_t;
	typedef typename DIMENSIONS::control_constraint1_matrix_array_t control_constraint1_matrix_array_t;

	enum WORKER_STATE {
		IDLE,
		LINE_SEARCH,
		APPROXIMATE_LQ,
		CALCULATE_CONTROLLER_AND_LAGRANGIAN,
		SHUTDOWN
	};

	SLQP_MP(const std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > >& subsystemDynamicsPtr,
			const std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > >& subsystemDerivativesPtr,
			std::vector<std::shared_ptr<CostFunctionBaseOCS2<OUTPUT_DIM, INPUT_DIM> > >& subsystemCostFunctionsPtr,
			const std::vector<controller_t>& initialControllersStock,
			const std::vector<size_t>& systemStockIndex,
			const Options_t& options = Options_t::Options(),
			const MP_Options_t& mp_options = MP_Options_t::MP_Options())

	: nominalTimeTrajectoriesStock_(NUM_SUBSYSTEMS),
	  nominalStateTrajectoriesStock_(NUM_SUBSYSTEMS),
	  nominalInputTrajectoriesStock_(NUM_SUBSYSTEMS),
	  nominalOutputTrajectoriesStock_(NUM_SUBSYSTEMS),
	  nominalcostateTrajectoriesStock_(NUM_SUBSYSTEMS),
	  nominalLagrangeTrajectoriesStock_(NUM_SUBSYSTEMS),
	  lagrangeControllerStock_(NUM_SUBSYSTEMS),
	  AmTrajectoryStock_(NUM_SUBSYSTEMS),
	  BmTrajectoryStock_(NUM_SUBSYSTEMS),
	  nc1TrajectoriesStock_(NUM_SUBSYSTEMS),
	  EvTrajectoryStock_(NUM_SUBSYSTEMS),
	  CmTrajectoryStock_(NUM_SUBSYSTEMS),
	  DmTrajectoryStock_(NUM_SUBSYSTEMS),
	  qTrajectoryStock_(NUM_SUBSYSTEMS),
	  QvTrajectoryStock_(NUM_SUBSYSTEMS),
	  QmTrajectoryStock_(NUM_SUBSYSTEMS),
	  RvTrajectoryStock_(NUM_SUBSYSTEMS),
	  RmTrajectoryStock_(NUM_SUBSYSTEMS),
	  PmTrajectoryStock_(NUM_SUBSYSTEMS),
	  RmInverseTrajectoryStock_(NUM_SUBSYSTEMS),
	  AmConstrainedTrajectoryStock_(NUM_SUBSYSTEMS),
	  QmConstrainedTrajectoryStock_(NUM_SUBSYSTEMS),
	  QvConstrainedTrajectoryStock_(NUM_SUBSYSTEMS),
	  RmConstrainedTrajectoryStock_(NUM_SUBSYSTEMS),
	  DmDagerTrajectoryStock_(NUM_SUBSYSTEMS),
	  EvProjectedTrajectoryStock_(NUM_SUBSYSTEMS),
	  CmProjectedTrajectoryStock_(NUM_SUBSYSTEMS),
	  DmProjectedTrajectoryStock_(NUM_SUBSYSTEMS),
	  SsTimeTrajectoryStock_(NUM_SUBSYSTEMS),
	  sTrajectoryStock_(NUM_SUBSYSTEMS),
	  SvTrajectoryStock_(NUM_SUBSYSTEMS),
	  SveTrajectoryStock_(NUM_SUBSYSTEMS),
	  SmTrajectoryStock_(NUM_SUBSYSTEMS),
	  switchingTimes_(NUM_SUBSYSTEMS+1),
	  iteration_(0),
	  options_(options),
	  mp_options_(mp_options),
	  feedForwardConstraintInputStock_(NUM_SUBSYSTEMS)
	{
		// resize instances to correct number of threads + 1
		dynamics_.resize(mp_options_.nThreads_+1);
		linearizedSystems_.resize(mp_options_.nThreads_+1);
		costFunctions_.resize(mp_options_.nThreads_+1);
		integratorsODE45_.resize(mp_options_.nThreads_+1);
		controllers_.resize(mp_options_.nThreads_+1);

		// for all threads + 1
		for (size_t i=0; i<mp_options.nThreads_+1; i++)
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
				integratorsODE45_[i].push_back(std::shared_ptr<ODE45<STATE_DIM>>(new ODE45<STATE_DIM> ((dynamics_[i]).back())));
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
		nominalOutputFunc_.resize(mp_options_.nThreads_+1);
		nominalInputFunc_.resize(mp_options_.nThreads_+1);
		BmFunc_.resize(mp_options_.nThreads_+1);
		PmFunc_.resize(mp_options_.nThreads_+1);
		RmInverseFunc_.resize(mp_options_.nThreads_+1);
		RvFunc_.resize(mp_options_.nThreads_+1);
		EvProjectedFunc_.resize(mp_options_.nThreads_+1);
		CmProjectedFunc_.resize(mp_options_.nThreads_+1);
		DmProjectedFunc_.resize(mp_options_.nThreads_+1);
		nominalLagrangeMultiplierFunc_.resize(mp_options_.nThreads_+1);
		DmDagerFunc_.resize(mp_options_.nThreads_+1);
		RmFunc_.resize(mp_options_.nThreads_+1);

		// initialize threads
		launchWorkerThreads();
	}

	~SLQP_MP();

	void rollout(
			const size_t threadId,
			const state_vector_t& initState,
			const std::vector<controller_t>& controllersStock,
			std::vector<scalar_array_t>& timeTrajectoriesStock,
			std::vector<state_vector_array_t>& stateTrajectoriesStock,
			std::vector<control_vector_array_t>& inputTrajectoriesStock,
			std::vector<output_vector_array_t>& outputTrajectoriesStock,
			std::vector<std::vector<size_t> >& nc1TrajectoriesStock,
			std::vector<constraint1_vector_array_t>& EvTrajectoryStock);

	void rollout(
			const size_t threadId,
			const state_vector_t& initState,
			const std::vector<controller_t>& controllersStock,
			std::vector<scalar_array_t>& timeTrajectoriesStock,
			std::vector<state_vector_array_t>& stateTrajectoriesStock,
			std::vector<control_vector_array_t>& inputTrajectoriesStock,
			std::vector<output_vector_array_t>& outputTrajectoriesStock);

	void rollout(
			const size_t threadId,
			const state_vector_t& initState,
			const std::vector<controller_t>& controllersStock,
			std::vector<scalar_array_t>& timeTrajectoriesStock,
			std::vector<state_vector_array_t>& stateTrajectoriesStock,
			std::vector<control_vector_array_t>& inputTrajectoriesStock);

	void calculateCostFunction(const std::vector<scalar_array_t>& timeTrajectoriesStock,
			const std::vector<output_vector_array_t>& stateTrajectoriesStock,
			const std::vector<control_vector_array_t>& inputTrajectoriesStock,
			scalar_t& totalCost);

	void calculateMeritFunction(const std::vector<scalar_array_t>& timeTrajectoriesStock,
			const std::vector<std::vector<size_t> >& nc1TrajectoriesStock,
			const std::vector<constraint1_vector_array_t>& EvTrajectoryStock,
			const std::vector<std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd> > >&  lagrangeTrajectoriesStock,
			const scalar_t& totalCost,
			scalar_t& meritFuntionValue,
			scalar_t& constraintISE);

	double calculateConstraintISE(const std::vector<scalar_array_t>& timeTrajectoriesStock,
			const std::vector<std::vector<size_t>>& nc1TrajectoriesStock,
			const std::vector<constraint1_vector_array_t>& EvTrajectoriesStock,
			scalar_t& constraintISE);

	void getController(std::vector<controller_t>& controllersStock);

	void getValueFuntion(const scalar_t& time, const output_vector_t& output, scalar_t& valueFuntion); // todo: getValueFunction

	void getCostFuntion(const output_vector_t& initOutput, scalar_t& costFunction, scalar_t& constriantCostFunction);

	void getNominalTrajectories(std::vector<scalar_array_t>& nominalTimeTrajectoriesStock,
			std::vector<state_vector_array_t>& nominalStateTrajectoriesStock,
			std::vector<control_vector_array_t>& nominalInputTrajectoriesStock,
			std::vector<output_vector_array_t>& nominalOutputTrajectoriesStock);

	void run(const state_vector_t& initState, const std::vector<scalar_t>& switchingTimes);


protected:
	void solveSequentialRiccatiEquations(const scalar_t& learningRate);

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

	void lineSearch(const std::vector<control_vector_array_t>& feedForwardConstraintInputStock,
			scalar_t& learningRateStar,
			scalar_t maxLearningRateStar=1.0);

	template <typename Derived>
	bool makePSD(Eigen::MatrixBase<Derived>& squareMatrix);


private:

	void launchWorkerThreads();
	void threadWork(size_t threadId);

	// main function for sub-tasks
	void approximateSubsystemLQ(); // computes the linearized dynamics for a particular subsystem
	void calculateControllerAndLagrangian();

	// worker functions
	void approximateSubsystemLQWorker(size_t threadId);		//Worker function for linearized dynamics
	void calculateControllerAndLagrangianWorker(size_t threadId);

	// execute methods
	void executeApproximateSubsystemLQ(size_t threadId, size_t k);	// Computes the linearized dynamics
	void executeCalculateControllerAndLagrangian(size_t threadId, size_t k);

	// for generating unique identifiers for subsystem, task, iteration:
	// note: arguments must not be passed by value here
	size_t generateUniqueProcessID (const size_t& iterateNo, const std::atomic_int& workerState, const std::atomic_int& subsystemId)
	{
		return (10e9*(workerState +1) + 10e6 * (subsystemId +1) + iterateNo);
	}

	std::vector<std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > > > dynamics_;
	std::vector<std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > > > linearizedSystems_;
	std::vector<std::vector<std::shared_ptr<CostFunctionBaseOCS2<OUTPUT_DIM, INPUT_DIM> > > > costFunctions_;

	std::vector<std::vector<std::shared_ptr<ODE45<STATE_DIM> > > > integratorsODE45_;

	std::vector<std::vector<controller_t> > controllers_;

	scalar_t nominalTotalCost_;
	scalar_t nominalTotalMerit_;
	scalar_t nominalConstraint1ISE_;
	std::vector<scalar_array_t> 		nominalTimeTrajectoriesStock_;
	std::vector<state_vector_array_t> 	nominalStateTrajectoriesStock_;
	std::vector<control_vector_array_t> nominalInputTrajectoriesStock_;
	std::vector<output_vector_array_t> 	nominalOutputTrajectoriesStock_;
	std::vector<output_vector_array_t> 	nominalcostateTrajectoriesStock_;
	std::vector<std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd> > >  nominalLagrangeTrajectoriesStock_;

	std::vector<lagrange_t> lagrangeControllerStock_;

	std::vector<state_matrix_array_t> 			AmTrajectoryStock_;
	std::vector<control_gain_matrix_array_t> 	BmTrajectoryStock_;

	std::vector<std::vector<size_t> > 				nc1TrajectoriesStock_;  // nc1: Number of the active constraints
	std::vector<constraint1_vector_array_t> 		EvTrajectoryStock_;
	std::vector<constraint1_state_matrix_array_t>   CmTrajectoryStock_;
	std::vector<constraint1_control_matrix_array_t> DmTrajectoryStock_;

	eigen_scalar_t  qFinal_;
	output_vector_t QvFinal_;
	state_matrix_t  QmFinal_;
	std::vector<eigen_scalar_array_t> 	qTrajectoryStock_;
	std::vector<output_vector_array_t> 	QvTrajectoryStock_;
	std::vector<state_matrix_array_t> 	QmTrajectoryStock_;
	std::vector<control_vector_array_t> RvTrajectoryStock_;
	std::vector<control_matrix_array_t> RmTrajectoryStock_;
	std::vector<control_feedback_array_t> PmTrajectoryStock_;

	std::vector<control_matrix_array_t> RmInverseTrajectoryStock_;
	std::vector<state_matrix_array_t>   AmConstrainedTrajectoryStock_;
	std::vector<state_matrix_array_t>   QmConstrainedTrajectoryStock_;
	std::vector<output_vector_array_t>  QvConstrainedTrajectoryStock_;
	std::vector<control_matrix_array_t> RmConstrainedTrajectoryStock_;
	std::vector<control_constraint1_matrix_array_t> DmDagerTrajectoryStock_;
	std::vector<control_vector_array_t>   EvProjectedTrajectoryStock_;  // DmDager * Ev
	std::vector<control_feedback_array_t> CmProjectedTrajectoryStock_;  // DmDager * Cm
	std::vector<control_matrix_array_t>   DmProjectedTrajectoryStock_;  // DmDager * Dm


	std::vector<scalar_array_t> 	  	SsTimeTrajectoryStock_;
	std::vector<eigen_scalar_array_t> 	sTrajectoryStock_;
	std::vector<output_vector_array_t> 	SvTrajectoryStock_;
	std::vector<output_vector_array_t> 	SveTrajectoryStock_;
	std::vector<state_matrix_array_t> 	SmTrajectoryStock_;

	scalar_array_t switchingTimes_;
	state_vector_t initState_;
	size_t iteration_;
	Options_t options_;

	MP_Options_t mp_options_;
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

	std::atomic_size_t KMax_subsystem_;		// denotes the number of integration steps for a particular subsystem i

	std::atomic_size_t alphaTaken_;
	size_t alphaMax_;
	size_t alphaExpBest_;
	size_t alphaExpMax_;
	std::atomic_bool alphaBestFound_;
	std::vector<size_t> alphaProcessed_;

	std::atomic_size_t kTaken_;
	std::atomic_size_t kCompleted_;

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

	std::vector<control_vector_array_t> feedForwardConstraintInputStock_;
	bool nominalLagrangeMultiplierUpdated_;

public:
	template <size_t GSLQP_STATE_DIM, size_t GSLQP_INPUT_DIM, size_t GSLQP_OUTPUT_DIM, size_t GSLQP_NUM_SUBSYSTEMS>
	friend class GSLQP;
};

} // namespace ocs2

#include "implementation/SLQP_MP.h"


#endif /* SLQP_MP_H_ */
