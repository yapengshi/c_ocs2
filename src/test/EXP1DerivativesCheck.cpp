/*
 * EXP1DerivativesCheck.cpp
 *
 *  Created on: Jan 19, 2016
 *      Author: farbod
 */

#include <iostream>

#include "test/EXP1.h"
#include "ocs2/OCS2DerivativesCheck.h"

int main (int argc, char* argv[])
{
	// subsystem dynamics
	std::vector<std::shared_ptr<ControlledSystemBase<2,1> > > subsystemDynamicsPtr {std::make_shared<EXP1_Sys1>(), std::make_shared<EXP1_Sys2>(), std::make_shared<EXP1_Sys3>()};

	// subsystem derivatives
	std::vector<std::shared_ptr<DerivativesBase<2,1> > > subsystemDerivativesPtr {std::make_shared<EXP1_SysDerivative1>(), std::make_shared<EXP1_SysDerivative2>(), std::make_shared<EXP1_SysDerivative3>()};

	// subsystem cost functions
	std::vector<std::shared_ptr<CostFunctionBase<2,1> > > subsystemCostFunctionsPtr {std::make_shared<EXP1_CostFunction1>(), std::make_shared<EXP1_CostFunction2>(), std::make_shared<EXP1_CostFunction3>()};

	GSLQP<2,1,2,3>::state_vector_array_t   stateOperatingPoints(3, GSLQP<2,1,2,3>::state_vector_t::Zero());
	GSLQP<2,1,2,3>::control_vector_array_t inputOperatingPoints(3, GSLQP<2,1,2,3>::control_vector_t::Zero());
	std::vector<size_t> systemStockIndex {0, 1, 2};

//	std::vector<double> initSwitchingTimes {0, 1.0, 2.0, 3};
	std::vector<double> initSwitchingTimes {0, 0.2262, 1.0176, 3};
	if (argc>1)  initSwitchingTimes[1] = std::atof(argv[1]);
	if (argc>2)  initSwitchingTimes[2] = std::atof(argv[2]);

	Eigen::Vector2d initState(2.0, 3.0);


	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/
	OCS2DerivativesCheck<2,1,2,3> derivativesCheck(subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr,
				stateOperatingPoints, inputOperatingPoints, systemStockIndex);

	derivativesCheck.check(initState, initSwitchingTimes);


}
