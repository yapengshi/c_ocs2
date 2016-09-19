/*
 * EXP1_SLQP_MP_Test.cpp
 *
 *  Created on: June 20, 2016
 *      Author: markus
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

#include <PathTweaker.h>

using namespace ocs2;

int main (int argc, char* argv[])
{
	const size_t nSys = 3; 	// number of subsystems
	const size_t stateDim = 2;
	const size_t controlDim = 1;
	const size_t outputDim = stateDim;

	// subsystem dynamics
	std::vector<std::shared_ptr<ControlledSystemBase<stateDim, controlDim, outputDim> > > subsystemDynamicsPtr {
		std::make_shared<EXP1_Sys1>(), std::make_shared<EXP1_Sys2>(), std::make_shared<EXP1_Sys3>()
	};

	// subsystem derivatives
	std::vector<std::shared_ptr<DerivativesBase<stateDim, controlDim> > > subsystemDerivativesPtr {
		std::make_shared<EXP1_SysDerivative1>(), std::make_shared<EXP1_SysDerivative2>(), std::make_shared<EXP1_SysDerivative3>()
	};

	// subsystem cost functions
	std::vector<std::shared_ptr<CostFunctionBaseOCS2<stateDim, controlDim> > > subsystemCostFunctionsPtr {
		std::make_shared<EXP1_CostFunction1>(), std::make_shared<EXP1_CostFunction2>(), std::make_shared<EXP1_CostFunction3>()
	};


	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/
	GSLQP<stateDim,controlDim, outputDim, nSys>::state_vector_array_t   stateOperatingPoints(nSys, GSLQP<stateDim,controlDim, outputDim, nSys>::state_vector_t::Zero());
	GSLQP<stateDim,controlDim, outputDim, nSys>::control_vector_array_t inputOperatingPoints(nSys, GSLQP<stateDim,controlDim, outputDim, nSys>::control_vector_t::Zero());
	std::vector<size_t> systemStockIndex {0, 1, 2};

	Eigen::Vector2d initState(2.0, 3.0);

	std::vector<double> switchingTimes {0, 0.2262, 1.0176, 3};
	if (argc>1)  switchingTimes[1] = std::atof(argv[1]);
	if (argc>2)  switchingTimes[2] = std::atof(argv[2]);

	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/
	// GLQP
	GLQP<stateDim,controlDim, outputDim, nSys> glqp(
			subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr, stateOperatingPoints, inputOperatingPoints, systemStockIndex);

	glqp.run(switchingTimes);

	// get controller
	std::vector<GLQP<stateDim,controlDim, outputDim, nSys>::controller_t> controllersStock(nSys);
	std::vector<GLQP<stateDim,controlDim, outputDim, nSys>::controller_t> controllersStock_mp(nSys);

	glqp.getController(controllersStock);
	glqp.getController(controllersStock_mp);


	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/
	SLQP<stateDim, controlDim, outputDim, nSys>::Options_t slqpOptions;
	slqpOptions.dispayGSLQP_ = 0;
	slqpOptions.lineSearchByMeritFuntion_ = false;

	SLQP_MP<stateDim, controlDim, outputDim, nSys>::Options_t slqpOptions_mp;
	slqpOptions_mp.dispayGSLQP_ = 0;
	slqpOptions.lineSearchByMeritFuntion_ = false;

	SLQP_MP<stateDim, controlDim, outputDim, nSys>::MP_Options_t mpOptions;
	mpOptions.nThreads_ = 4;
	mpOptions.debugPrintMP_ = false;

	// SLQP
	SLQP	<stateDim,controlDim, outputDim, nSys> slqp		(subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr, controllersStock, 	systemStockIndex, slqpOptions);
	SLQP_MP <stateDim,controlDim, outputDim, nSys> slqp_mp	(subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr, controllersStock_mp, systemStockIndex, slqpOptions_mp, mpOptions);

	std::cout << "Starting SLQP_MP" << std::endl;
	slqp_mp.run(initState, switchingTimes);

	std::cout << "Starting SLQP single core" << std::endl;
	slqp.run(initState, switchingTimes);


	double totalCost_mp;
	double constrCost_mp;
	double totalValue_mp;
	slqp_mp.getValueFuntion(0.0, initState, totalValue_mp);
	slqp_mp.getCostFuntion(totalCost_mp, constrCost_mp);

	double  totalCost;
	double constrCost;
	double totalValue;
	slqp.getValueFuntion(0.0, initState, totalValue);
	slqp.getCostFuntion(totalCost, constrCost);

	std::cout << "SLQ_MP cost			" << totalCost_mp << std::endl;
	std::cout << "single core cost 		" << totalCost << std::endl;

}

