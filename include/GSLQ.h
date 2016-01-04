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

#include "dynamics/ControlledSystemBase.h"
#include "dynamics/DerivativesBase.hpp"
#include "costs/CostFunctionBase.hpp"

template <size_t STATE_DIM, size_t CONTROL_DIM>
class GSLQ
{
public:
	GSLQ( std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, CONTROL_DIM> > > subsystemsDynamics,
			std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, CONTROL_DIM> > > subsystemsDerivatives,
			std::vector<std::shared_ptr<CostFunctionBase<STATE_DIM, CONTROL_DIM> > > subsystemsCostFunctions,
			std::vector<size_t> systemStockIndex)
		: numSubsystems_(systemStockIndex.size()),
		  subsystemsDynamicsStock_(numSubsystems_),
		  subsystemsDerivativesStock_(numSubsystems_),
		  subsystemsCostFunctionsStock_(numSubsystems_),
		  maxIteration_(10)
	{

		if (subsystemsDynamics.size() != subsystemsDerivatives.size())
			throw std::runtime_error("Number of subsystem derivaties is not equal to the number of subsystems.");
		if (subsystemsDynamics.size() != subsystemsCostFunctions.size())
			throw std::runtime_error("Number of cost functions is not equal to the number of subsystems.");
		if (subsystemsDynamics.size() != *std::max_element(systemStockIndex.begin(), systemStockIndex.end()))
			throw std::runtime_error("systemStockIndex points to non-existing subsystem");

		for (size_t i; i<numSubsystems_; i++) {

			subsystemsDynamicsStock_[i] = subsystemsDynamics[systemStockIndex[i]]->clone();
			subsystemsDerivativesStock_[i] = subsystemsDerivatives[systemStockIndex[i]]->clone();
			subsystemsCostFunctionsStock_[i] = subsystemsCostFunctions[systemStockIndex[i]]->clone();
		}

	}




	bool rollout(std::vector<int> switchingTimes)
	{
		if (switchingTimes.size() != numSubsystems_)
			throw std::runtime_error("Number of switching times should be one plus the number of subsystems.");


	}



private:
	std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, CONTROL_DIM> > > subsystemsDynamicsStock_;
	std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, CONTROL_DIM> > > subsystemsDerivativesStock_;
	std::vector<std::shared_ptr<CostFunctionBase<STATE_DIM, CONTROL_DIM> > > subsystemsCostFunctionsStock_;

	size_t numSubsystems_;
	size_t maxIteration_;

};


#endif /* GSLQ_H_ */
