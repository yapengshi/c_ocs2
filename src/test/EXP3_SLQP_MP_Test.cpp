/*
 * EXP3_SLQP_MP_Test.cpp
 *
 *  Created on: July 17, 2016
 *      Author: markus
 */

#include <iostream>
#include <cstdlib>

#include <PathTweaker.h>

#include <fstream>
#include <cereal/archives/xml.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/Eigen.hpp>

#include "test/EXP3.h"
#include "GSLQ/GLQP.h"
#include "GSLQ/SLQP.h"
#include "GSLQ/SLQP_MP.h"


using namespace ocs2;

int main (int argc, char* argv[])
{
	// subsystem dynamics
	std::vector<std::shared_ptr<ControlledSystemBase<2,2> > > subsystemDynamicsPtr {std::make_shared<EXP3_Sys1>(), std::make_shared<EXP3_Sys2>(), std::make_shared<EXP3_Sys3>()};

	// subsystem derivatives
	std::vector<std::shared_ptr<DerivativesBase<2,2> > > subsystemDerivativesPtr {std::make_shared<EXP3_SysDerivative1>(), std::make_shared<EXP3_SysDerivative2>(), std::make_shared<EXP3_SysDerivative3>()};

	// subsystem cost functions
	std::vector<std::shared_ptr<CostFunctionBaseOCS2<2,2> > > subsystemCostFunctionsPtr {std::make_shared<EXP3_CostFunction1>(), std::make_shared<EXP3_CostFunction2>(), std::make_shared<EXP3_CostFunction3>()};


	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/
	SLQP_MP<2,2,2,3>::state_vector_array_t   stateOperatingPoints(3, SLQP_MP<2,2,2,3>::state_vector_t::Zero());
	SLQP_MP<2,2,2,3>::control_vector_array_t inputOperatingPoints(3, SLQP_MP<2,2,2,3>::control_vector_t::Zero());
	std::vector<size_t> systemStockIndex {0, 1, 2};

	Eigen::Vector2d initState(2.0, 3.0);

	std::vector<double> switchingTimes {0, 0.2262, 1.0176, 3};
//	if (argc>1)  switchingTimes[1] = std::atof(argv[1]);
//	if (argc>2)  switchingTimes[2] = std::atof(argv[2]);

	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/
	// GLQP
	GLQP<2,2,2,3> glqp(subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr, stateOperatingPoints, inputOperatingPoints, systemStockIndex);

	glqp.run(switchingTimes);

	// get controller
	std::vector<GLQP<2,2,2,3>::controller_t> controllersStock(3);
	glqp.getController(controllersStock);


	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/
	SLQP_MP<2,2,2,3>::Options_t slqpOptions;
	slqpOptions.dispayGSLQP_ = 1;
	slqpOptions.maxIterationGSLQP_ = 100;
	slqpOptions.meritFunctionRho_ = 2000.0;
	slqpOptions.constraintStepSize_ = 0.2;
	slqpOptions.lineSearchByMeritFuntion_ = false;
	if (argc>1) slqpOptions.meritFunctionRho_ = std::atof(argv[1]);

	SLQP_MP<2,2,2,3>::MP_Options_t mpOptions;
	mpOptions.nThreads_ = 4;
	mpOptions.debugPrintMP_ = 0;

	// Single core SLQP
	SLQP<2,2,2,3> slqp(subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr, controllersStock, systemStockIndex, slqpOptions);

	// SLQP_MP
	SLQP_MP<2,2,2,3> slqp_mp(subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr, controllersStock, systemStockIndex, slqpOptions, mpOptions);

	// run slqp for reference
//	std::cout << " =========================== Starting single core SLQP ===============================" << std::endl;
//	slqp.run(initState, switchingTimes);
//	std::cout << " =========================== End of single core SLQP ===============================" << std::endl;

	std::cout << " =========================== Starting multi core SLQP ===============================" << std::endl;
	slqp_mp.run(initState, switchingTimes);
	std::cout << " =========================== End of multi core SLQP ===============================" << std::endl;

	exit(0);

	// get controller
	slqp.getController(controllersStock);

	// rollout
	std::vector<SLQP_MP<2,2,2,3>::scalar_array_t> timeTrajectoriesStock;
	std::vector<SLQP_MP<2,2,2,3>::state_vector_array_t> stateTrajectoriesStock;
	std::vector<SLQP_MP<2,2,2,3>::control_vector_array_t> controlTrajectoriesStock;
//	slqp.rollout(initState, controllersStock, timeTrajectoriesStock, stateTrajectoriesStock, controlTrajectoriesStock);

	// compute cost
	double rolloutCost;
//	slqp.calculateCostFunction(timeTrajectoriesStock, stateTrajectoriesStock, controlTrajectoriesStock, rolloutCost);

	// value funtion
	double totalCost;
	slqp.getValueFuntion(0.0, initState, totalCost);

//	// value funtion derivative
//	Eigen::Matrix<double,2,1> costFuntionDerivative;
//	slqp.getCostFuntionDerivative(costFuntionDerivative);

	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/
	std::cout << std::endl;

	std::cout << "Switching times are: [" << switchingTimes[0] << ", ";
	for (size_t i=1; i<switchingTimes.size()-1; i++)  std::cout << switchingTimes[i] << ", ";
	std::cout << switchingTimes.back() << "]\n";

	std::cout << "The total cost: " << totalCost << std::endl;
//	std::cout << "The total cost in the test rollout: " << rolloutCost << std::endl;
//	std::cout << "The total cost derivative: " << costFuntionDerivative.transpose() << std::endl;

}
