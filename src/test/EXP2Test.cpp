/*
 * EXP2TEST.cpp
 *
 *  Created on: Dec 27, 2015
 *      Author: farbod
 */

#include "test/EXP2.h"


int main (int argc, char* argv[])
{

	// subsystem dynamics
	std::vector<std::shared_ptr<ControlledSystemBase<2,1> > > subsystemDynamicsPtr {std::make_shared<EXP2_Sys1>(), std::make_shared<EXP2_Sys2>()};

	// subsystem derivatives
	std::vector<std::shared_ptr<DerivativesBase<2,1> > > subsystemDerivativesPtr {std::make_shared<EXP2_SysDerivative1>(), std::make_shared<EXP2_SysDerivative2>()};

	// subsystem cost functions
	std::vector<std::shared_ptr<CostFunctionBase<2,1> > > subsystemCostFunctionsPtr {std::make_shared<EXP2_CostFunction1>(), std::make_shared<EXP2_CostFunction2>()};

	GLQP<2,1>::state_vector_array_t   stateOperatingPoints(2, Eigen::Matrix<double,2,1>::Zero());;
	GLQP<2,1>::control_vector_array_t inputOperatingPoints(2, Eigen::Matrix<double,1,1>::Zero());;
	std::vector<size_t> systemStockIndex {0, 1};

	// GLQP
	GLQP<2,1> glqp(subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr, stateOperatingPoints, inputOperatingPoints, systemStockIndex);

	std::vector<double> switchingTimes {0, 0.8, 2};
	glqp.SolveRiccatiEquation(switchingTimes);

	// get controller
	std::vector<GLQP<2,1>::controller_t> controllersStock(2);
	glqp.getController(controllersStock);

	// rollout
	Eigen::Vector2d initState(0.0, 2.0);
	std::vector<GLQP<2,1>::scalar_array_t> timeTrajectoriesStock;
	std::vector<GLQP<2,1>::state_vector_array_t> stateTrajectoriesStock;
	std::vector<GLQP<2,1>::control_vector_array_t> controlTrajectoriesStock;
	glqp.rollout(initState, controllersStock,
				timeTrajectoriesStock, stateTrajectoriesStock, controlTrajectoriesStock);

	// compute cost
	double totalCost;
	glqp.rolloutCost(timeTrajectoriesStock, stateTrajectoriesStock, controlTrajectoriesStock, totalCost);

}
