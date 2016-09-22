/*
 * EXP1GSLQPTest.cpp
 *
 *  Created on: Sept 20, 2016
 *      Author: markus
 */


#include <iostream>
#include <cstdlib>

#include <fstream>
#include <cereal/archives/xml.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/Eigen.hpp>

#include "test/EXP5.h"
#include "GSLQ/GLQP.h"
#include "GSLQ/SLQP.h"
#include "GSLQ/SLQP_MP.h"

#include <PathTweaker.h>

using namespace ocs2;

int main (int argc, char* argv[])
{
	// subsystem dynamics
	std::vector<std::shared_ptr<ControlledSystemBase<4,2> > > subsystemDynamicsPtr {std::make_shared<EXP5_Sys1>()};

	// subsystem derivatives
	std::vector<std::shared_ptr<DerivativesBase<4,2> > > subsystemDerivativesPtr {std::make_shared<EXP5_SysDerivative1>()};

	// subsystem cost functions
	std::vector<std::shared_ptr<CostFunctionBaseOCS2<4,2> > > subsystemCostFunctionsPtr {std::make_shared<EXP5_CostFunction1>()};


	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/
	GSLQP<4,2,4,1>::state_vector_array_t   stateOperatingPoints(1, GSLQP<4,2,4,1>::state_vector_t::Zero());
	GSLQP<4,2,4,1>::control_vector_array_t inputOperatingPoints(1, GSLQP<4,2,4,1>::control_vector_t::Zero());
	std::vector<size_t> systemStockIndex {0};

	Eigen::Matrix<double, 4 ,1 > initState;
	initState << 0.0, 0.0, 0.1, 0.0;

	std::vector<double> switchingTimes {0, 3};
	std::vector<double> switchingTimes_mp {0, 3};


	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/
	// GLQP
	GLQP<4,2,4,1> glqp(subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr,
			stateOperatingPoints, inputOperatingPoints, systemStockIndex);

	glqp.run(switchingTimes);

	// get controller
	std::vector<GLQP<4,2,4,1>::controller_t> controllersStock(1);
	std::vector<GLQP<4,2,4,1>::controller_t> controllersStock_mp(1);

	glqp.getController(controllersStock);
	glqp.getController(controllersStock_mp);


	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/
	GSLQP<4,2,4,1>::Options_t gslqpOptions;
	gslqpOptions.dispayGSLQP_ = 1;
	gslqpOptions.useMultiThreading_ = false;
	gslqpOptions.minLearningRateGSLQP_ = 0.01;
	gslqpOptions.minRelCostGSLQP_ = 1e-4;
	gslqpOptions.stateConstraintPenaltyCoeff_ = 1.0;
	gslqpOptions.stateConstraintPenaltyBase_ = 2.0;
	gslqpOptions.lineSearchByMeritFuntion_ = false;

	// GSLQ - single core version
	GSLQP<4,2,4,1> gslqp(subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr,
			controllersStock, systemStockIndex, gslqpOptions);

	// GSLQ MP version
	GSLQP<4,2,4,1>::Options_t gslqpOptions_mp = gslqpOptions;
	gslqpOptions_mp.useMultiThreading_ = true;
	GSLQP<4,2,4,1>::MP_Options_t mpOptions;
	mpOptions.nThreads_ = 4;
	mpOptions.debugPrintMP_ = 0;
	mpOptions.lsStepsizeGreedy_ = 1;
	GSLQP<4,2,4,1> gslqp_mp(subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr,
			controllersStock_mp, systemStockIndex, gslqpOptions_mp, mpOptions);

	SLQP<4,2,4,1> slqp(subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr, controllersStock, systemStockIndex, gslqpOptions);
	SLQP_MP<4,2,4,1> slqp_mp(subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr, controllersStock, systemStockIndex, gslqpOptions, mpOptions);


	// run both the mp and the single core versions of slqp
	slqp.run(initState, switchingTimes);
	slqp_mp.run(initState, switchingTimes);

	// try to run gslqp
	gslqp.run(initState, switchingTimes);
	gslqp_mp.run(initState, switchingTimes_mp);


	// get controller
	gslqp.getController(controllersStock);
	gslqp_mp.getController(controllersStock_mp);

	// rollout both versions
	std::vector<GSLQP<4,2,4,1>::scalar_array_t> timeTrajectoriesStock, timeTrajectoriesStock_mp;
	std::vector<GSLQP<4,2,4,1>::state_vector_array_t> stateTrajectoriesStock, stateTrajectoriesStock_mp;
	std::vector<GSLQP<4,2,4,1>::control_vector_array_t> controlTrajectoriesStock, controlTrajectoriesStock_mp;
	gslqp.rollout(initState, controllersStock, timeTrajectoriesStock, stateTrajectoriesStock, controlTrajectoriesStock);
	gslqp_mp.rollout(initState, controllersStock_mp, timeTrajectoriesStock_mp, stateTrajectoriesStock_mp, controlTrajectoriesStock_mp);

	// compute cost for both versions
	double rolloutCost, rolloutCost_mp;
	gslqp.calculateCostFunction(timeTrajectoriesStock, stateTrajectoriesStock, controlTrajectoriesStock, rolloutCost);
	gslqp_mp.calculateCostFunction(timeTrajectoriesStock_mp, stateTrajectoriesStock_mp, controlTrajectoriesStock_mp, rolloutCost_mp);

	// value function for both versions
	double totalCost;
	double totalCost_mp;
	gslqp.getValueFuntion(0.0, initState, totalCost);
	gslqp_mp.getValueFuntion(0.0, initState, totalCost_mp);


	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/
	std::cout << std::endl;

	std::cout << "Single core switching times are: [" << switchingTimes[0] << ", ";
	for (size_t i=1; i<switchingTimes.size()-1; i++)  std::cout << switchingTimes[i] << ", ";
	std::cout << switchingTimes.back() << "]\n";

	std::cout << "MP switching times are: [" << switchingTimes_mp[0] << ", ";
	for (size_t i=1; i<switchingTimes_mp.size()-1; i++)  std::cout << switchingTimes_mp[i] << ", ";
	std::cout << switchingTimes_mp.back() << "]\n";

	std::cout << "The single core total cost: " << totalCost << std::endl;
	std::cout << "The MP total cost: " << totalCost_mp << std::endl;

	std::cout << "The single core total cost in the test rollout: " << rolloutCost << std::endl;
	std::cout << "The MP total cost in the test rollout: " << rolloutCost_mp << std::endl;


}

