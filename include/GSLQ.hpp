/*
 * GSLQ.hpp
 *
 *  Created on: Dec 18, 2015
 *      Author: farbod
 */

#ifndef GSLQ_HPP_
#define GSLQ_HPP_

#include <vector>
#include <algorithm>
#include <Eigen/Dense>

#include <ilqg/iLQGMP.hpp>

#include <dynamics/DynamicsBase.hpp>
#include <dynamics/DerivativesBase.hpp>

#include <costs/CostFunctionBase.hpp>

template <size_t STATE_DIM, size_t CONTROL_DIM>
class GSLQ
{

public:
	GSLQ( std::vector<std::shared_ptr<DynamicsBase<STATE_DIM, CONTROL_DIM> > > subsystemsDynamics,
			std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, CONTROL_DIM> > > subsystemsDerivatives,
			std::vector<std::shared_ptr<CostFunctionBase<STATE_DIM, CONTROL_DIM> > > subsystemsCostFunctions,
			size_t maxIteration,
			std::vector<size_t> systemStock)
	: numSubsystems_(systemStock.size()),
	  maxIteration_(maxIteration)

	{
		if (subsystemsDynamics.size()==subsystemsDerivatives.size()) throw std::runtime_error("Number of subsystem derivaties is not equal to the number of subsystems.");
		if (subsystemsDynamics.size()==subsystemsCostFunctions.size()) throw std::runtime_error("Number of cost functions is not equal to the number of subsystems.");
		if (subsystemsDynamics.size()==*std::max_element(systemStock.begin(), systemStock.end())) throw std::runtime_error("systemStock points to non-existing subsystem");



		for (size_t i; i<numSubsystems_; i++) {
			iLQGMP<STATE_DIM, CONTROL_DIM>::multi_processing_settings_t mpSettings;
			mpSettings.nThreads = nThreads_;
			std::vector<std::shared_ptr<DynamicsBase<STATE_DIM, CONTROL_DIM> > > dynamics(nThreads_);
			std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, CONTROL_DIM> > > derivatives(nThreads_);
			std::vector<std::shared_ptr<CostFunctionBase<STATE_DIM, CONTROL_DIM> > > costFunction(nThreads_);

			for (size_t j; j<nThreads_; j++) {
				dynamics[j]     = subsystemsDynamics[systemStock[i]]->clone();
				derivatives[j]  = subsystemsDerivatives[systemStock[i]]->clone();
				costFunction[j] = subsystemsCostFunctions[systemStock[i]]->clone();
			}


			ilqc.push_back(iLQGMP<STATE_DIM, CONTROL_DIM>(dynamics, derivatives, costFunction, maxIteration_, mpSettings));
		}
	}




	bool run(double dt, double dt_sim, std::vector<int> switchingTimes)
	{

		if (switchingTimes.size()==numSubsystems_) throw std::runtime_error("Number of switching times should be one plus the number of subsystems.");

		bool foundBetter;

		for (size_t j=0; j<maxIteration_; j++)  {
			foundBetter = ilqg.runIteration(K, dt, dt_sim);


			if (foundBetter)  {
				ilqg.retrieveLastRollout(x_rollout);
				ilqg.retrieveController(u_rollout);
				ilqg.retrieveLastLinearizedModel(A, B);
			}
			else  {
				break;
			}
		}


	}








private:
	size_t numSubsystems_;
	size_t maxIteration_;

	const size_t nThreads_ = 4;

	std::vector<iLQGMP<STATE_DIM, CONTROL_DIM> > ilqc;



};





#endif /* GSLQ_HPP_ */
