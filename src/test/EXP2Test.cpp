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

#include "ocs2/OCS2Ipopt.h"


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

//	std::vector<double> initSwitchingTimes {0, 0.184, 2};
	std::vector<double> initSwitchingTimes {0, 1, 2};
	if (argc>1)  initSwitchingTimes[1] = std::atof(argv[1]);

	Eigen::Vector2d initState(0.0, 2.0);

	OCS2Ipopt<2,1,2,2>::Options_t gslqpOptions;
	gslqpOptions.maxIterationIPOPT_ = 5;
	gslqpOptions.warmStartGSLQP_ = true;
//	gslqpOptions.dispayGSLQP_ = true;
//	gslqpOptions.displayIPOPT_ = false;


	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/
	OCS2Ipopt<2,1,2,2> ocs2 (subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr,
			stateOperatingPoints, inputOperatingPoints, systemStockIndex, initSwitchingTimes, initState, gslqpOptions);

	ocs2.run();

	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/

//	std::string resultDir = "/home/farbod/Programs/ct_ws/src/c_ocs2/cereal/test/exp2_test";
//	std::string stateFile = resultDir + "/exp2State.xml";
//	std::string timeFile = resultDir + "/exp2Time.xml";
//	std::string inputFile = resultDir + "/exp2Input.xml";
//	std::string stateSensitivityFile = resultDir + "/exp2StateSensitivity.xml";
//	std::string timeSensitivityFile = resultDir + "/exp2TimeSensitivity.xml";
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



