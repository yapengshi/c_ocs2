/*
 * GSLQ.h
 *
 *  Created on: Dec 18, 2015
 *      Author: farbod
 */

#ifndef GSLQ_H_
#define GSLQ_H_

#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/StdVector>

#include "Dimensions.h"

#include "dynamics/ControlledSystemBase.h"
#include "dynamics/DerivativesBase.hpp"
#include "costs/CostFunctionBase.hpp"

#include "integration/Integrator.h"


template <size_t STATE_DIM, size_t INPUT_DIM>
class GSLQ
{
public:
	typedef Dimensions<STATE_DIM, INPUT_DIM> DIMENSIONS;
	typedef typename DIMENSIONS::controller_t controller_t;
	typedef typename DIMENSIONS::scalar_t 		scalar_t;
	typedef typename DIMENSIONS::scalar_array_t scalar_array_t;
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


	GSLQ( std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM> > > subsystemDynamicsPtr,
			std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM> > > subsystemDerivativesPtr,
			std::vector<std::shared_ptr<CostFunctionBase<STATE_DIM, INPUT_DIM> > > subsystemCostFunctionsPtr,
			std::vector<size_t> systemStockIndex)
		: numSubsystems_(systemStockIndex.size()),
		  subsystemDynamicsPtrStock(numSubsystems_),
		  subsystemDerivativesPtrStock_(numSubsystems_),
		  subsystemCostFunctionsPtrStock_(numSubsystems_),
		  subsystemSimulatorsStock_(numSubsystems_),
		  nominalControllersStock_(numSubsystems_),
		  nominalTimeTrajectoriesStock_(numSubsystems_),
		  nominalStateTrajectoriesStock_(numSubsystems_),
		  nominalControlTrajectoriesStock_(numSubsystems_),
		  AmTrajectoryStock_(numSubsystems_),
		  BmTrajectoryStock_(numSubsystems_),
		  qTrajectoryStock_(numSubsystems_),
		  QvTrajectoryStock_(numSubsystems_),
		  QmTrajectoryStock_(numSubsystems_),
		  RvTrajectoryStock_(numSubsystems_),
		  RmTrajectoryStock_(numSubsystems_),
		  PmTrajectoryStock_(numSubsystems_),
		  maxIteration_(10)
	{

		if (subsystemDynamicsPtr.size() != subsystemDerivativesPtr.size())
			throw std::runtime_error("Number of subsystem derivaties is not equal to the number of subsystems.");
		if (subsystemDynamicsPtr.size() != subsystemCostFunctionsPtr.size())
			throw std::runtime_error("Number of cost functions is not equal to the number of subsystems.");
		if (subsystemDynamicsPtr.size() != *std::max_element(systemStockIndex.begin(), systemStockIndex.end()))
			throw std::runtime_error("systemStockIndex points to non-existing subsystem");

		for (int i=0; i<numSubsystems_; i++) {

			subsystemDynamicsPtrStock[i] = subsystemDynamicsPtr[systemStockIndex[i]]->clone();
			subsystemDerivativesPtrStock_[i] = subsystemDerivativesPtr[systemStockIndex[i]]->clone();
			subsystemCostFunctionsPtrStock_[i] = subsystemCostFunctionsPtr[systemStockIndex[i]]->clone();

			subsystemSimulatorsStock_[i] = ODE45<STATE_DIM>(subsystemDynamicsPtrStock[i]);
		}
	}

	~GSLQ() {}


	void rollout(const state_vector_t& initState,
			const std::vector<scalar_t>& switchingTimes,
			const std::vector<controller_t>& controllersStock,
			std::vector<scalar_array_t>& timeTrajectoriesStock,
			std::vector<state_vector_array_t>& stateTrajectoriesStock,
			std::vector<control_vector_array_t>& controlTrajectoriesStock);

	void rolloutCost(const std::vector<scalar_array_t>& timeTrajectoriesStock,
			const std::vector<state_vector_array_t>& stateTrajectoriesStock,
			const std::vector<control_vector_array_t>& controlTrajectoriesStock,
			scalar_t& totalCost);

	void approximateOptimalControlProblem();





private:
	std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM> > > subsystemDynamicsPtrStock;
	std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM> > > subsystemDerivativesPtrStock_;
	std::vector<std::shared_ptr<CostFunctionBase<STATE_DIM, INPUT_DIM> > > subsystemCostFunctionsPtrStock_;

	std::vector<ODE45<STATE_DIM> > subsystemSimulatorsStock_;

	std::vector<controller_t> nominalControllersStock_;
	std::vector<scalar_array_t> nominalTimeTrajectoriesStock_;
	std::vector<state_vector_array_t> nominalStateTrajectoriesStock_;
	std::vector<control_vector_array_t> nominalControlTrajectoriesStock_;

	size_t numSubsystems_;
	size_t maxIteration_;


	std::vector<state_matrix_array_t> AmTrajectoryStock_;
	std::vector<control_gain_matrix_array_t> BmTrajectoryStock_;

	scalar_t       qFinal_;
	state_vector_t QvFinal_;
	state_matrix_t QmFinal_;
	std::vector<scalar_array_t> 	  qTrajectoryStock_;
	std::vector<state_vector_array_t> QvTrajectoryStock_;
	std::vector<state_matrix_array_t> QmTrajectoryStock_;
	std::vector<control_vector_array_t> RvTrajectoryStock_;
	std::vector<control_matrix_array_t> RmTrajectoryStock_;
	std::vector<control_feedback_array_t> PmTrajectoryStock_;

};

#include "implementation/GSLQ.h"

#endif /* GSLQ_H_ */
