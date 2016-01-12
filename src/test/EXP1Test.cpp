/*
 * EXP1Test.cpp
 *
 *  Created on: Jan 11, 2016
 *      Author: farbod
 */


#include <iostream>
#include <cstdlib>

#include <fstream>
#include <cereal/archives/xml.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/Eigen.hpp>

#include "test/EXP1.h"
#include "GSLQ/GLQP.h"

//#include "ocs2/OCS2.h"
#include "ocs2/OCS2Ipopt.h"

int main (int argc, char* argv[])
{
	// subsystem dynamics
	std::vector<std::shared_ptr<ControlledSystemBase<2,1> > > subsystemDynamicsPtr {std::make_shared<EXP1_Sys1>(), std::make_shared<EXP1_Sys2>(), std::make_shared<EXP1_Sys3>()};

	// subsystem derivatives
	std::vector<std::shared_ptr<DerivativesBase<2,1> > > subsystemDerivativesPtr {std::make_shared<EXP1_SysDerivative1>(), std::make_shared<EXP1_SysDerivative2>(), std::make_shared<EXP1_SysDerivative3>()};

	// subsystem cost functions
	std::vector<std::shared_ptr<CostFunctionBase<2,1> > > subsystemCostFunctionsPtr {std::make_shared<EXP1_CostFunction1>(), std::make_shared<EXP1_CostFunction2>(), std::make_shared<EXP1_CostFunction3>()};

	GSLQP<2,1,3>::state_vector_array_t   stateOperatingPoints(3, GSLQP<2,1,3>::state_vector_t::Zero());
	GSLQP<2,1,3>::control_vector_array_t inputOperatingPoints(3, GSLQP<2,1,3>::control_vector_t::Zero());
	std::vector<size_t> systemStockIndex {0, 1, 2};

	std::vector<double> initSwitchingTimes {0, 1.0, 2.0, 3};
//	std::vector<double> initSwitchingTimes {0, 0.2262, 1.0176, 3};
	if (argc>1)  initSwitchingTimes[1] = std::atof(argv[1]);
	if (argc>2)  initSwitchingTimes[2] = std::atof(argv[2]);

	Eigen::Vector2d initState(2.0, 3.0);

	GSLQP<2,1,3>::Options gslqpOptions;
	gslqpOptions.maxIteration_ = 50;
	gslqpOptions.warmStart_ = true;
	gslqpOptions.dispay_ = 1;


	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/
//	OCS2<2,1,3> ocs2(subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr,
//			stateOperatingPoints, inputOperatingPoints, systemStockIndex, initSwitchingTimes, initState, gslqpOptions);

	OCS2Ipopt<2,1,3> ocs2 (subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr,
			stateOperatingPoints, inputOperatingPoints, systemStockIndex, initSwitchingTimes, initState, gslqpOptions);

	ocs2.run();

	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/

//	std::string resultDir = "/home/farbod/Programs/ct_ws/src/c_ocs2/cereal/test/exp1_test";
//	std::string stateFile = resultDir + "/exp1State.xml";
//	std::string timeFile = resultDir + "/exp1Time.xml";
//	std::string inputFile = resultDir + "/exp1Input.xml";
//	std::string stateSensitivityFile = resultDir + "/exp1StateSensitivity.xml";
//	std::string timeSensitivityFile = resultDir + "/exp1TimeSensitivity.xml";
//
//	{ // we need these brackets to make sure the archive goes out of scope and flushes
//		std::ofstream xmlState(stateFile);
//		cereal::XMLOutputArchive archive_state_xml(xmlState);
//		archive_state_xml(CEREAL_NVP(stateTrajectory));
//
//		std::ofstream xmlTime(timeFile);
//		cereal::XMLOutputArchive archive_time_xml(xmlTime);
//		archive_time_xml(CEREAL_NVP(timeEigenTrajectory));
//
//		std::ofstream xmlInput(inputFile);
//		cereal::XMLOutputArchive archive_input_xml(xmlInput);
//		archive_input_xml(CEREAL_NVP(inputTrajectory));
//
//		std::ofstream xmlStateSensitivity(stateSensitivityFile);
//		cereal::XMLOutputArchive archive_stateSensitivity_xml(xmlStateSensitivity);
//		archive_stateSensitivity_xml(CEREAL_NVP(sensitivityStateTrajectory));
//
//		std::ofstream xmlTimeSensitivity(timeSensitivityFile);
//		cereal::XMLOutputArchive archive_timeSensitivity_xml(xmlTimeSensitivity);
//		archive_timeSensitivity_xml(CEREAL_NVP(sensitivityTimeTrajectory));
//	}

}

