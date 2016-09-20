/*
 * EXP2TEST.cpp
 *
 *  Created on: Dec 27, 2015
 *      Author: farbod
 */

#include <iostream>
#include <cstdlib>

#include <fstream>
#include <cereal/archives/xml.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/Eigen.hpp>

#include "test/EXP2.h"
#include "GSLQ/GLQP.h"

#include "ocs2/OCS2Projected.h"


using namespace ocs2;

int main (int argc, char* argv[])
{
	// subsystem dynamics
	std::vector<std::shared_ptr<ControlledSystemBase<2,1> > > subsystemDynamicsPtr {std::make_shared<EXP2_Sys1>(), std::make_shared<EXP2_Sys2>()};

	// subsystem derivatives
	std::vector<std::shared_ptr<DerivativesBase<2,1> > > subsystemDerivativesPtr {std::make_shared<EXP2_SysDerivative1>(), std::make_shared<EXP2_SysDerivative2>()};

	// subsystem cost functions
	std::vector<std::shared_ptr<CostFunctionBaseOCS2<2,1> > > subsystemCostFunctionsPtr {std::make_shared<EXP2_CostFunction1>(), std::make_shared<EXP2_CostFunction2>()};

	GSLQP<2,1,2,2>::state_vector_array_t   stateOperatingPoints(2, GSLQP<2,1,2,2>::state_vector_t::Zero());
	GSLQP<2,1,2,2>::control_vector_array_t inputOperatingPoints(2, GSLQP<2,1,2,2>::control_vector_t::Zero());
	std::vector<size_t> systemStockIndex {0, 1};

	std::vector<double> initSwitchingTimes {0, 1, 2};

	if (argc>1)
		initSwitchingTimes[1] = std::atof(argv[1]);

	Eigen::Vector2d initState(0.0, 2.0);

	OCS2Projected<2,1,2,2>::Options_t gslqpOptions;
	gslqpOptions.maxIterationGSLQP_ = 50;
	gslqpOptions.maxIterationGradientDescent_ = 5;
	gslqpOptions.warmStartGSLQP_ = false;
	gslqpOptions.dispayGSLQP_ = false;
	gslqpOptions.displayGradientDescent_ = true;
	gslqpOptions.useLQForDerivatives_ = false;
	gslqpOptions.minLearningRateNLP_ = 0.01;
	gslqpOptions.acceptableTolGradientDescent_ = 1e-3;
	gslqpOptions.useAscendingLineSearchNLP_ = true;
	gslqpOptions.minAcceptedSwitchingTimeDifference_ = 0.01;
	gslqpOptions.print();



	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/

	// setup single core version
	OCS2Projected<2,1,2,2> ocs2 (subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr,
			stateOperatingPoints, inputOperatingPoints, systemStockIndex, gslqpOptions);

	// setup multi core version
	gslqpOptions.useMultiThreading_ = true;
	OCS2Projected<2,1,2,2> ocs2_mp (subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr,
				stateOperatingPoints, inputOperatingPoints, systemStockIndex, gslqpOptions);

	ocs2.run(initState, initSwitchingTimes);
	ocs2_mp.run(initState, initSwitchingTimes);

	std::vector<double> resultingSwitchingTimes,resultingSwitchingTimes_mp;
	ocs2.getSwitchingTimes(resultingSwitchingTimes);
	ocs2_mp.getSwitchingTimes(resultingSwitchingTimes_mp);

	double cost, cost_mp;
	ocs2.getCostFunction(cost);
	ocs2_mp.getCostFunction(cost_mp);

	std::cout << "resulting costs: single core: " << cost << ", mp " << cost_mp << std::endl;
}



