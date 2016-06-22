/*
 * EXP1_SLQP_MP_Test.cpp
 *
 *  Created on: June 20, 2016
 *      Author: markus
 */


#include <iostream>
#include <cstdlib>

#include <fstream>
#include <cereal/archives/xml.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/Eigen.hpp>

#include "test/EXP1.h"
#include "GSLQ/GLQP.h"
#include "GSLQ/SLQP.h"
#include "GSLQ/SLQP_MP.h"

#include <PathTweaker.h>

using namespace ocs2;

int main (int argc, char* argv[])
{
	// subsystem dynamics
	std::vector<std::shared_ptr<ControlledSystemBase<2,1> > > subsystemDynamicsPtr {std::make_shared<EXP1_Sys1>(), std::make_shared<EXP1_Sys2>(), std::make_shared<EXP1_Sys3>()};

	// subsystem derivatives
	std::vector<std::shared_ptr<DerivativesBase<2,1> > > subsystemDerivativesPtr {std::make_shared<EXP1_SysDerivative1>(), std::make_shared<EXP1_SysDerivative2>(), std::make_shared<EXP1_SysDerivative3>()};

	// subsystem cost functions
	std::vector<std::shared_ptr<CostFunctionBaseOCS2<2,1> > > subsystemCostFunctionsPtr {std::make_shared<EXP1_CostFunction1>(), std::make_shared<EXP1_CostFunction2>(), std::make_shared<EXP1_CostFunction3>()};


	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/
	GSLQP<2,1,2,3>::state_vector_array_t   stateOperatingPoints(3, GSLQP<2,1,2,3>::state_vector_t::Zero());
	GSLQP<2,1,2,3>::control_vector_array_t inputOperatingPoints(3, GSLQP<2,1,2,3>::control_vector_t::Zero());
	std::vector<size_t> systemStockIndex {0, 1, 2};

	Eigen::Vector2d initState(2.0, 3.0);

	std::vector<double> switchingTimes {0, 0.2262, 1.0176, 3};
	if (argc>1)  switchingTimes[1] = std::atof(argv[1]);
	if (argc>2)  switchingTimes[2] = std::atof(argv[2]);

	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/
	// GLQP
	GLQP<2,1,2,3> glqp(subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr, stateOperatingPoints, inputOperatingPoints, systemStockIndex);

	glqp.run(switchingTimes);

	// get controller
	std::vector<GLQP<2,1,2,3>::controller_t> controllersStock(3);
	std::vector<GLQP<2,1,2,3>::controller_t> controllersStock_mp(3);

	glqp.getController(controllersStock);
	glqp.getController(controllersStock_mp);


	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/
	SLQP<2,1,2,3>::Options_t slqpOptions;
	slqpOptions.dispayGSLQP_ = 1;
	slqpOptions.lineSearchByMeritFuntion_ = false;

	SLQP_MP<2,1,2,3>::MP_Options_t mpOptions;
	mpOptions.nThreads_ = 4;
	mpOptions.debugPrintMP_ = true;

	// SLQP
	SLQP<2,1,2,3> slqp(subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr, controllersStock, systemStockIndex, slqpOptions);
	SLQP_MP<2,1,2,3> slqp_mp(subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr,
			controllersStock_mp, systemStockIndex, slqpOptions, mpOptions);

	std::cout << "Starting SLQP_MP" << std::endl;
	slqp_mp.run(initState, switchingTimes);

	/*
	std::cout << "Starting SLQP" << std::endl;
	slqp.run(initState, switchingTimes);

	// get controller
	slqp_mp.getController(controllersStock_mp);
	slqp.getController(controllersStock);

	// rollout
	std::vector<SLQP<2,1,2,3>::scalar_array_t> timeTrajectoriesStock;
	std::vector<SLQP<2,1,2,3>::state_vector_array_t> stateTrajectoriesStock;
	std::vector<SLQP<2,1,2,3>::control_vector_array_t> controlTrajectoriesStock;
	std::vector<SLQP<2,1,2,3>::scalar_array_t> timeTrajectoriesStock_mp;
	std::vector<SLQP<2,1,2,3>::state_vector_array_t> stateTrajectoriesStock_mp;
	std::vector<SLQP<2,1,2,3>::control_vector_array_t> controlTrajectoriesStock_mp;
	slqp_mp.rollout(0, initState, controllersStock_mp, timeTrajectoriesStock_mp, stateTrajectoriesStock_mp, controlTrajectoriesStock_mp);
	slqp.rollout(initState, controllersStock, timeTrajectoriesStock, stateTrajectoriesStock, controlTrajectoriesStock);

	// compute cost
	double rolloutCost_mp;
	double rolloutCost;
	slqp_mp.calculateCostFunction(timeTrajectoriesStock_mp, stateTrajectoriesStock_mp, controlTrajectoriesStock_mp, rolloutCost_mp);
	slqp.calculateCostFunction(timeTrajectoriesStock, stateTrajectoriesStock, controlTrajectoriesStock, rolloutCost);

	// value funtion
	double totalCost;
	slqp.getValueFuntion(0.0, initState, totalCost);

	std::cout << "Total cost comparison: " << std::endl;
	std::cout << "SLQ   	" << totalCost << std::endl;
	*/

	double totalCost_mp;
	slqp_mp.getValueFuntion(0.0, initState, totalCost_mp);

	std::cout << "SLQ_MP	" << totalCost_mp << std::endl;


}

