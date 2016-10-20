/*
 * A unit test
 *
 *  Created on: Sept 20, 2016
 *      Author: mgiftthaler
 */

#include <iostream>
#include <cstdlib>
#include <ctime>

#include <fstream>
#include <cereal/archives/xml.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/Eigen.hpp>

#include "test/EXP1.h"
#include "GSLQ/GLQP.h"
#include "GSLQ/SLQP.h"
#include "GSLQ/SLQP_MP.h"

#include <gtest/gtest.h>

using namespace ocs2;

TEST(Exp1_gslqp_test, Exp1_gslqp_test)
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
	std::vector<double> switchingTimes_mp {0, 0.2262, 1.0176, 3};

	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/
	// GLQP
	GLQP<2,1,2,3> glqp(subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr,
			stateOperatingPoints, inputOperatingPoints, systemStockIndex);

	glqp.run(switchingTimes);

	// get controller
	std::vector<GLQP<2,1,2,3>::controller_t> controllersStock(3);
	std::vector<GLQP<2,1,2,3>::controller_t> controllersStock_mp(3);

	glqp.getController(controllersStock);
	glqp.getController(controllersStock_mp);


	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/
	GSLQP<2,1,2,3>::Options_t gslqpOptions;
	gslqpOptions.dispayGSLQP_ = 1;
	gslqpOptions.lineSearchByMeritFuntion_ = false;
	gslqpOptions.useMultiThreading_ = false;

	// GSLQ - single core version
	GSLQP<2,1,2,3> gslqp(subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr,
			controllersStock, systemStockIndex, gslqpOptions);

	// GSLQ MP version
	GSLQP<2,1,2,3>::Options_t gslqpOptions_mp = gslqpOptions;
	gslqpOptions_mp.useMultiThreading_ = true;
	gslqpOptions_mp.nThreads_ = 4;
	gslqpOptions_mp.debugPrintMP_ = 0;
	gslqpOptions_mp.lsStepsizeGreedy_ = 1;
	GSLQP<2,1,2,3> gslqp_mp(subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr,
			controllersStock_mp, systemStockIndex, gslqpOptions_mp);

	// run both the mp and the single core version
	gslqp.run(initState, switchingTimes);
	gslqp_mp.run(initState, switchingTimes_mp);

	// get controller
	gslqp.getController(controllersStock);
	gslqp_mp.getController(controllersStock_mp);

	// rollout both versions
	std::vector<GSLQP<2,1,2,3>::scalar_array_t> timeTrajectoriesStock, timeTrajectoriesStock_mp;
	std::vector<GSLQP<2,1,2,3>::state_vector_array_t> stateTrajectoriesStock, stateTrajectoriesStock_mp;
	std::vector<GSLQP<2,1,2,3>::control_vector_array_t> controlTrajectoriesStock, controlTrajectoriesStock_mp;
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

	// value function derivative for both versions
	Eigen::Matrix<double,2,1> costFunctionDerivative, costFunctionDerivative_mp;
	gslqp.getCostFuntionDerivative(costFunctionDerivative);
	gslqp_mp.getCostFuntionDerivative(costFunctionDerivative_mp);


	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/
	std::cout << std::endl;

	std::cout << "Single core switching times are: [" << switchingTimes[0] << ", ";
	for (size_t i=1; i<switchingTimes.size()-1; i++)
		std::cout << switchingTimes[i] << ", ";
	std::cout << switchingTimes.back() << "]\n";

	std::cout << "MP switching times are: [" << switchingTimes_mp[0] << ", ";
	for (size_t i=1; i<switchingTimes_mp.size()-1; i++)
		std::cout << switchingTimes_mp[i] << ", ";
	std::cout << switchingTimes_mp.back() << "]\n";

	for (size_t i=0; i<switchingTimes_mp.size(); i++)
			ASSERT_LT(fabs(switchingTimes_mp[i]-switchingTimes[i]), 1e-5);

	std::cout << "The single core total cost in the test rollout: " << rolloutCost << std::endl;
	std::cout << "The MP total cost in the test rollout: " << rolloutCost_mp << std::endl;
	ASSERT_LT(fabs(rolloutCost_mp - rolloutCost), 1e-5);

	std::cout << "The single core total cost derivative: " << costFunctionDerivative.transpose() << std::endl;
	std::cout << "The MP total cost derivative: " << costFunctionDerivative_mp.transpose() << std::endl;
	for (size_t i=0; i<costFunctionDerivative_mp.size(); i++)
				ASSERT_LT(fabs(costFunctionDerivative[i]-costFunctionDerivative_mp[i]), 1e-5);

	std::cout << "SLQ_MP cost			" << totalCost_mp << std::endl;
	std::cout << "single core cost 		" << totalCost << std::endl;
	ASSERT_LT(fabs(totalCost_mp - totalCost), 1e-5);
}


int main(int argc, char** argv)
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
