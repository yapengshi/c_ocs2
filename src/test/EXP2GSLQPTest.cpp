/*
 * EXP2GSLQPTest.cpp
 *
 *  Created on: Jan 13, 2016
 *      Author: farbod
 */

#include <iostream>
#include <cstdlib>

#include <fstream>
#include <cereal/archives/xml.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/Eigen.hpp>

#include "test/EXP2.h"
#include "GSLQ/GSLQP.h"

#include <PathTweaker.h>

//#define LOG_DATA

using namespace ocs2;

int main (int argc, char* argv[])
{
	// subsystem dynamics
	std::vector<std::shared_ptr<ControlledSystemBase<2,1> > > subsystemDynamicsPtr {std::make_shared<EXP2_Sys1>(), std::make_shared<EXP2_Sys2>()};

	// subsystem derivatives
	std::vector<std::shared_ptr<DerivativesBase<2,1> > > subsystemDerivativesPtr {std::make_shared<EXP2_SysDerivative1>(), std::make_shared<EXP2_SysDerivative2>()};

	// subsystem cost functions
	std::vector<std::shared_ptr<CostFunctionBaseOCS2<2,1> > > subsystemCostFunctionsPtr {std::make_shared<EXP2_CostFunction1>(), std::make_shared<EXP2_CostFunction2>()};


	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/
	GLQP<2,1,2,2>::state_vector_array_t   stateOperatingPoints(2, GLQP<2,1,2,2>::state_vector_t::Zero());
	GLQP<2,1,2,2>::control_vector_array_t inputOperatingPoints(2, GLQP<2,1,2,2>::control_vector_t::Zero());
	std::vector<size_t> systemStockIndex {0, 1};

	Eigen::Vector2d initState(0.0, 2.0);

	std::vector<double> switchingTimes {0, 0.184, 2};
	if (argc>1)  switchingTimes[1] = std::atof(argv[1]);

	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/
	// GLQP
	GLQP<2,1,2,2> glqp(subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr, stateOperatingPoints, inputOperatingPoints, systemStockIndex);

	glqp.run(switchingTimes);

	// get controller
	std::vector<GLQP<2,1,2,2>::controller_t> controllersStock(2);
	glqp.getController(controllersStock);

//	// rollout
//	std::vector<GLQP<2,1,2,2>::scalar_array_t> timeTrajectoriesStock;
//	std::vector<GLQP<2,1,2,2>::state_vector_array_t> stateTrajectoriesStock;
//	std::vector<GLQP<2,1,2,2>::control_vector_array_t> controlTrajectoriesStock;
//	glqp.rollout(initState, controllersStock, timeTrajectoriesStock, stateTrajectoriesStock, controlTrajectoriesStock);
//
//	// compute cost
//	double rolloutCost;
//	glqp.rolloutCost(timeTrajectoriesStock, stateTrajectoriesStock, controlTrajectoriesStock, rolloutCost);
//
//	// value funtion
//	double totalCost;
//	glqp.getValueFuntion(0.0, initState, totalCost);


	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/
	GSLQP<2,1,2,2>::Options_t gslqpOptions;
	gslqpOptions.dispayGSLQP_ = 1;
	gslqpOptions.lineSearchByMeritFuntion_ = false;

	GSLQP<2,1,2,2>::MP_Options_t mpOptions;
	mpOptions.nThreads_ = 4;
	mpOptions.debugPrintMP_ = 0;
	mpOptions.lsStepsizeGreedy_ = 1;


	// GSLQ
	GSLQP<2,1,2,2> gslqp(subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr, controllersStock, systemStockIndex, gslqpOptions, mpOptions);

	gslqp.run(initState, switchingTimes);

	// get controller
	gslqp.getController(controllersStock);

	// rollout
	std::vector<GSLQP<2,1,2,2>::scalar_array_t> timeTrajectoriesStock;
	std::vector<GSLQP<2,1,2,2>::state_vector_array_t> stateTrajectoriesStock;
	std::vector<GSLQP<2,1,2,2>::control_vector_array_t> controlTrajectoriesStock;
	gslqp.rollout(initState, controllersStock, timeTrajectoriesStock, stateTrajectoriesStock, controlTrajectoriesStock);

	// compute cost
	double rolloutCost;
	gslqp.calculateCostFunction(timeTrajectoriesStock, stateTrajectoriesStock, controlTrajectoriesStock, rolloutCost);

	// value funtion
	double totalCost;
	gslqp.getValueFuntion(0.0, initState, totalCost);

	// value funtion derivative
	Eigen::Matrix<double,1,1> costFuntionDerivative;
	gslqp.getCostFuntionDerivative(costFuntionDerivative);

	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/
	std::cout << std::endl;

	std::cout << "Switching times are: [" << switchingTimes[0] << ", ";
	for (size_t i=1; i<switchingTimes.size()-1; i++)  std::cout << switchingTimes[i] << ", ";
	std::cout << switchingTimes[2] << "]\n";

	std::cout << "The total cost: " << totalCost << std::endl;
	std::cout << "The total cost in the test rollout: " << rolloutCost << std::endl;
	std::cout << "The total cost derivative: " << costFuntionDerivative.transpose() << std::endl;

	GSLQP<2,1,2,2>::eigen_scalar_array_t timeEigenTrajectory;
	GSLQP<2,1,2,2>::state_vector_array_t stateTrajectory;
	GSLQP<2,1,2,2>::control_vector_array_t inputTrajectory;

	for (size_t i=0; i<switchingTimes.size()-1; i++)  {

		for (size_t k=0; k<timeTrajectoriesStock[i].size(); k++)  {

			timeEigenTrajectory.push_back((Eigen::MatrixXd(1,1) << timeTrajectoriesStock[i][k]).finished());
			stateTrajectory.push_back(stateTrajectoriesStock[i][k]);
			inputTrajectory.push_back(controlTrajectoriesStock[i][k]);
		}
	}


	// Sensitivity2SwitchingTime
	std::vector<GSLQP<2,1,2,2>::scalar_array_t> sensitivityTimeTrajectoriesStock;
	std::vector<GSLQP<2,1,2,2>::nabla_output_matrix_array_t> sensitivityStateTrajectoriesStock;
	std::vector<GSLQP<2,1,2,2>::nabla_input_matrix_array_t> sensitivityInputTrajectoriesStock;
	gslqp.getRolloutSensitivity2SwitchingTime(sensitivityTimeTrajectoriesStock, sensitivityStateTrajectoriesStock, sensitivityInputTrajectoriesStock);

	GSLQP<2,1,2,2>::eigen_scalar_array_t sensitivityTimeTrajectory;
	GSLQP<2,1,2,2>::nabla_output_matrix_array_t sensitivityStateTrajectory;
	for (size_t i=0; i<switchingTimes.size()-1; i++)  {
		for (size_t k=0; k<sensitivityTimeTrajectoriesStock[i].size(); k++)  {
			sensitivityTimeTrajectory.push_back((Eigen::MatrixXd(1,1) << sensitivityTimeTrajectoriesStock[i][k]).finished());
			sensitivityStateTrajectory.push_back(sensitivityStateTrajectoriesStock[i][k]);
		}
	}


#ifdef LOG_DATA
	PathTweaker pathTweaker(argv);

	std::string resultDir = pathTweaker.getDirectory() +"/src/c_ocs2/cereal/test/exp2_test";

	std::cout << "Saving to directory " << resultDir << std::endl;

	std::string stateFile = resultDir + "/exp2State.xml";
	std::string timeFile = resultDir + "/exp2Time.xml";
	std::string inputFile = resultDir + "/exp2Input.xml";
	std::string stateSensitivityFile = resultDir + "/exp2StateSensitivity.xml";
	std::string timeSensitivityFile = resultDir + "/exp2TimeSensitivity.xml";

	{ // we need these brackets to make sure the archive goes out of scope and flushes
		std::ofstream xmlState(stateFile);
		cereal::XMLOutputArchive archive_state_xml(xmlState);
		archive_state_xml(CEREAL_NVP(stateTrajectory));

		std::ofstream xmlTime(timeFile);
		cereal::XMLOutputArchive archive_time_xml(xmlTime);
		archive_time_xml(CEREAL_NVP(timeEigenTrajectory));

		std::ofstream xmlInput(inputFile);
		cereal::XMLOutputArchive archive_input_xml(xmlInput);
		archive_input_xml(CEREAL_NVP(inputTrajectory));

		std::ofstream xmlStateSensitivity(stateSensitivityFile);
		cereal::XMLOutputArchive archive_stateSensitivity_xml(xmlStateSensitivity);
		archive_stateSensitivity_xml(CEREAL_NVP(sensitivityStateTrajectory));

		std::ofstream xmlTimeSensitivity(timeSensitivityFile);
		cereal::XMLOutputArchive archive_timeSensitivity_xml(xmlTimeSensitivity);
		archive_timeSensitivity_xml(CEREAL_NVP(sensitivityTimeTrajectory));
	}
#endif

}



