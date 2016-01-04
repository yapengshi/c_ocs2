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
	typedef typename DIMENSIONS::scalar_t scalar_t;
	typedef typename DIMENSIONS::scalar_array_t scalar_array_t;
	typedef typename DIMENSIONS::state_vector_t state_vector_t;
	typedef typename DIMENSIONS::state_vector_array_t state_vector_array_t;
	typedef typename DIMENSIONS::control_vector_t control_vector_t;
	typedef typename DIMENSIONS::control_vector_array_t control_vector_array_t;
	typedef typename DIMENSIONS::control_feedback_t control_feedback_t;
	typedef typename DIMENSIONS::control_feedback_array_t control_feedback_array_t;


	GSLQ( std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM> > > subsystemsDynamicsPtr,
			std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM> > > subsystemsDerivativesPtr,
			std::vector<std::shared_ptr<CostFunctionBase<STATE_DIM, INPUT_DIM> > > subsystemsCostFunctionsPtr,
			std::vector<size_t> systemStockIndex)
		: numSubsystems_(systemStockIndex.size()),
		  subsystemDynamicsPtrStock(numSubsystems_),
		  subsystemDerivativesPtrStock_(numSubsystems_),
		  subsystemCostFunctionsPtrStock_(numSubsystems_),
		  subsystemSimulatorsStock_(numSubsystems_),
		  controllersStock_(numSubsystems_),
		  timeTrajectoriesStock_(numSubsystems_),
		  stateTrajectoriesStock_(numSubsystems_),
		  controlTrajectoriesStock_(numSubsystems_),
		  maxIteration_(10)
	{

		if (subsystemsDynamicsPtr.size() != subsystemsDerivativesPtr.size())
			throw std::runtime_error("Number of subsystem derivaties is not equal to the number of subsystems.");
		if (subsystemsDynamicsPtr.size() != subsystemsCostFunctionsPtr.size())
			throw std::runtime_error("Number of cost functions is not equal to the number of subsystems.");
		if (subsystemsDynamicsPtr.size() != *std::max_element(systemStockIndex.begin(), systemStockIndex.end()))
			throw std::runtime_error("systemStockIndex points to non-existing subsystem");

		for (int i=0; i<numSubsystems_; i++) {

			subsystemDynamicsPtrStock[i] = subsystemsDynamicsPtr[systemStockIndex[i]]->clone();
			subsystemDerivativesPtrStock_[i] = subsystemsDerivativesPtr[systemStockIndex[i]]->clone();
			subsystemCostFunctionsPtrStock_[i] = subsystemsCostFunctionsPtr[systemStockIndex[i]]->clone();

			subsystemSimulatorsStock_[i] = ODE45<STATE_DIM>(subsystemDynamicsPtrStock[i]);
		}

	}




	void rollout(const state_vector_t& initState, const std::vector<scalar_t>& switchingTimes)  {

		if (switchingTimes.size() != numSubsystems_+1)
			throw std::runtime_error("Number of switching times should be one plus the number of subsystems.");

		state_vector_t x0 = initState;
		for (int i=0; i<numSubsystems_; i++) {

			timeTrajectoriesStock_[i].clear();
			stateTrajectoriesStock_[i].clear();
			controlTrajectoriesStock_[i].clear();

			// set controller for subsystem i
			subsystemDynamicsPtrStock[i]->setController(controllersStock_[i]);
			// simulate subsystem i
			subsystemSimulatorsStock_[i].integrate(x0, switchingTimes[i], switchingTimes[i+1], stateTrajectoriesStock_[i], timeTrajectoriesStock_[i]);

			// compute control trajectory for subsystem i
			controlTrajectoriesStock_[i].resize(timeTrajectoriesStock_[i].size());
			for (int j=0; j<timeTrajectoriesStock_[i].size(); j++)
				subsystemDynamicsPtrStock[i]->computeInput(timeTrajectoriesStock_[i](j), stateTrajectoriesStock_[i](j), controlTrajectoriesStock_[i](j));

			// reset the initial state
			x0 = stateTrajectoriesStock_[i].back();
		}
	}

	void rolloutCost(scalar_t& J) {


	}



private:
	std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM> > > subsystemDynamicsPtrStock;
	std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM> > > subsystemDerivativesPtrStock_;
	std::vector<std::shared_ptr<CostFunctionBase<STATE_DIM, INPUT_DIM> > > subsystemCostFunctionsPtrStock_;

	std::vector<ODE45<STATE_DIM> > subsystemSimulatorsStock_;

	std::vector<controller_t> controllersStock_;
	std::vector<scalar_array_t> timeTrajectoriesStock_;
	std::vector<state_vector_array_t> stateTrajectoriesStock_;
	std::vector<control_vector_array_t> controlTrajectoriesStock_;

	size_t numSubsystems_;
	size_t maxIteration_;

};


#endif /* GSLQ_H_ */
