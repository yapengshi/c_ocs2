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
		  controllersStock_(numSubsystems_),
		  timeTrajectoriesStock_(numSubsystems_),
		  stateTrajectoriesStock_(numSubsystems_),
		  controlTrajectoriesStock_(numSubsystems_),
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




	void rollout(const state_vector_t& initState,
			const std::vector<scalar_t>& switchingTimes,
			const std::vector<controller_t>& controllersStock,
			std::vector<scalar_array_t>& timeTrajectoriesStock,
			std::vector<state_vector_array_t>& stateTrajectoriesStock,
			std::vector<control_vector_array_t>& controlTrajectoriesStock)  {

		if (switchingTimes.size() != numSubsystems_+1)
			throw std::runtime_error("Number of switching times should be one plus the number of subsystems.");

		state_vector_t x0 = initState;
		for (int i=0; i<numSubsystems_; i++) {

			timeTrajectoriesStock[i].clear();
			stateTrajectoriesStock[i].clear();
			controlTrajectoriesStock[i].clear();

			// set controller for subsystem i
			subsystemDynamicsPtrStock[i]->setController(controllersStock[i]);
			// simulate subsystem i
			subsystemSimulatorsStock_[i].integrate(x0, switchingTimes[i], switchingTimes[i+1], stateTrajectoriesStock[i], timeTrajectoriesStock[i]);

			// compute control trajectory for subsystem i
			controlTrajectoriesStock[i].resize(timeTrajectoriesStock[i].size());
			for (int k=0; k<timeTrajectoriesStock[i].size(); k++)
				subsystemDynamicsPtrStock[i]->computeInput(timeTrajectoriesStock[i][k], stateTrajectoriesStock[i][k], controlTrajectoriesStock[i][k]);

			// reset the initial state
			x0 = stateTrajectoriesStock[i].back();
		}
	}

	void rolloutCost(const std::vector<scalar_array_t>& timeTrajectoriesStock,
			const std::vector<state_vector_array_t>& stateTrajectoriesStock,
			const std::vector<control_vector_array_t>& controlTrajectoriesStock,
			scalar_t& totalCost)  {

		totalCost = 0.0;
		for (int i=0; i<numSubsystems_; i++) {

			scalar_t currentIntermediateCost;
			scalar_t nextIntermediateCost;
			for (int k=0; k<timeTrajectoriesStock[i].size()-1; k++) {

				if (k==0) {
					subsystemCostFunctionsPtrStock_[i]->setCurrentStateAndControl(timeTrajectoriesStock[i][k], stateTrajectoriesStock[i][k], controlTrajectoriesStock[i][k]);
					subsystemCostFunctionsPtrStock_[i]->evaluate(currentIntermediateCost);
				} else {
					currentIntermediateCost = nextIntermediateCost;
				}

				// feed next state and control to cost function
				subsystemCostFunctionsPtrStock_[i]->setCurrentStateAndControl(timeTrajectoriesStock[i][k+1], stateTrajectoriesStock[i][k+1], controlTrajectoriesStock[i][k+1]);
				// evaluate intermediate cost for next time step
				subsystemCostFunctionsPtrStock_[i]->evaluate(nextIntermediateCost);

				totalCost += 0.5*(currentIntermediateCost+nextIntermediateCost)*(timeTrajectoriesStock[i][k+1]-timeTrajectoriesStock[i][k]);
			}

			// terminal cost
			if (i==numSubsystems_-1)  {
				scalar_t finalCost;
				subsystemCostFunctionsPtrStock_[i]->terminalCost(finalCost);
				totalCost += finalCost;
			}
		}

	}

	void approximateOptimalControlProblem()  {

		for (int i=0; i<numSubsystems_; i++) {

			for (int k=0; k<nominalTimeTrajectoriesStock_[i].size(); k++) {

				subsystemDerivativesPtrStock_[i]->setCurrentStateAndControl(nominalTimeTrajectoriesStock_[i][k],
						nominalStateTrajectoriesStock_[i][k],
						nominalControlTrajectoriesStock_[i][k]);
				subsystemDerivativesPtrStock_[i]->getDerivativeState(AmTrajectoryStock_[i][k]);
				subsystemDerivativesPtrStock_[i]->getDerivativesControl(BmTrajectoryStock_[i][k]);

				subsystemCostFunctionsPtrStock_[i]->setCurrentStateAndControl(nominalTimeTrajectoriesStock_[i][k],
						nominalStateTrajectoriesStock_[i][k],
						nominalControlTrajectoriesStock_[i][k]);
				subsystemCostFunctionsPtrStock_[i]->evaluate(qTrajectoryStock_[i][k]);
				subsystemCostFunctionsPtrStock_[i]->stateDerivative(QvTrajectoryStock_[i][k]);
				subsystemCostFunctionsPtrStock_[i]->stateSecondDerivative(QmTrajectoryStock_[i][k]);
				subsystemCostFunctionsPtrStock_[i]->controlDerivative(RvTrajectoryStock_[i][k]);
				subsystemCostFunctionsPtrStock_[i]->controlSecondDerivative(RmTrajectoryStock_[i][k]);
				subsystemCostFunctionsPtrStock_[i]->stateControlDerivative(PmTrajectoryStock_[i][k]);

			}

			if (i==numSubsystems_-1)  {
				subsystemCostFunctionsPtrStock_[i]->terminalCost(qFinal_);
				subsystemCostFunctionsPtrStock_[i]->terminalCostStateDerivative(QvFinal_);
				subsystemCostFunctionsPtrStock_[i]->terminalCostStateSecondDerivative(QmFinal_);
			}
		}
	}





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


#endif /* GSLQ_H_ */
