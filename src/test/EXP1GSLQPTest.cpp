/*
 * EXP1GSLQPTest.cpp
 *
 *  Created on: Jan 12, 2016
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
	std::vector<double> switchingTimes_mp {0, 0.2262, 1.0176, 3};

	if (argc>1)  switchingTimes[1] = std::atof(argv[1]);
	if (argc>2)  switchingTimes[2] = std::atof(argv[2]);
	if (argc>1)  switchingTimes_mp[1] = std::atof(argv[1]);
	if (argc>2)  switchingTimes_mp[2] = std::atof(argv[2]);

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

	//	// rollout
	//	std::vector<GLQP<2,1,2,3>::scalar_array_t> timeTrajectoriesStock;
	//	std::vector<GLQP<2,1,2,3>::state_vector_array_t> stateTrajectoriesStock;
	//	std::vector<GLQP<2,1,2,3>::control_vector_array_t> controlTrajectoriesStock;
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
	//	SLQP<2,1,2,3>::Options_t slqpOptions;
	//	slqpOptions.dispayGSLQP_ = 1;
	//	slqpOptions.lineSearchByMeritFuntion_ = false;
	//
	//	// SLQP
	//	SLQP<2,1,2,3> slqp(subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr, controllersStock, systemStockIndex, slqpOptions);
	//
	//	slqp.run(initState, switchingTimes);
	//
	//	// get controller
	//	slqp.getController(controllersStock);
	//
	//	// rollout
	//	std::vector<SLQP<2,1,2,3>::scalar_array_t> timeTrajectoriesStock;
	//	std::vector<SLQP<2,1,2,3>::state_vector_array_t> stateTrajectoriesStock;
	//	std::vector<SLQP<2,1,2,3>::control_vector_array_t> controlTrajectoriesStock;
	//	slqp.rollout(initState, controllersStock, timeTrajectoriesStock, stateTrajectoriesStock, controlTrajectoriesStock);
	//
	//	// compute cost
	//	double rolloutCost;
	//	slqp.calculateCostFunction(timeTrajectoriesStock, stateTrajectoriesStock, controlTrajectoriesStock, rolloutCost);
	//
	//	// value funtion
	//	double totalCost;
	//	slqp.getValueFuntion(0.0, initState, totalCost);


	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/
	GSLQP<2,1,2,3>::Options_t gslqpOptions;
	gslqpOptions.dispayGSLQP_ = 1;
	gslqpOptions.lineSearchByMeritFuntion_ = false;
	gslqpOptions.useMultiThreading_ = false;

	GSLQP<2,1,2,3>::MP_Options_t mpOptions;
	mpOptions.nThreads_ = 4;
	mpOptions.debugPrintMP_ = 0;
	mpOptions.lsStepsizeGreedy_ = 1;

	// GSLQ
	GSLQP<2,1,2,3> gslqp(subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr,
			controllersStock, systemStockIndex, gslqpOptions, mpOptions);

	// GSLQ MP
	gslqpOptions.useMultiThreading_ = true;
	GSLQP<2,1,2,3> gslqp_mp(subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr,
			controllersStock_mp, systemStockIndex, gslqpOptions, mpOptions);

	gslqp.run(initState, switchingTimes);
	gslqp_mp.run(initState, switchingTimes_mp);

	// get controller
	gslqp.getController(controllersStock);
	gslqp_mp.getController(controllersStock_mp);

	// rollout
	std::vector<GSLQP<2,1,2,3>::scalar_array_t> timeTrajectoriesStock, timeTrajectoriesStock_mp;
	std::vector<GSLQP<2,1,2,3>::state_vector_array_t> stateTrajectoriesStock, stateTrajectoriesStock_mp;
	std::vector<GSLQP<2,1,2,3>::control_vector_array_t> controlTrajectoriesStock, controlTrajectoriesStock_mp;
	gslqp.rollout(initState, controllersStock, timeTrajectoriesStock, stateTrajectoriesStock, controlTrajectoriesStock);
	gslqp_mp.rollout(initState, controllersStock_mp, timeTrajectoriesStock_mp, stateTrajectoriesStock_mp, controlTrajectoriesStock_mp);

	// compute cost
	double rolloutCost, rolloutCost_mp;
	gslqp.calculateCostFunction(timeTrajectoriesStock, stateTrajectoriesStock, controlTrajectoriesStock, rolloutCost);
	gslqp_mp.calculateCostFunction(timeTrajectoriesStock_mp, stateTrajectoriesStock_mp, controlTrajectoriesStock_mp, rolloutCost_mp);

	// value function
	double totalCost;
	double totalCost_mp;
	gslqp.getValueFuntion(0.0, initState, totalCost);
	gslqp_mp.getValueFuntion(0.0, initState, totalCost_mp);

	// value function derivative
	Eigen::Matrix<double,2,1> costFunctionDerivative, costFunctionDerivative_mp;
	gslqp.getCostFuntionDerivative(costFunctionDerivative);
	gslqp_mp.getCostFuntionDerivative(costFunctionDerivative_mp);


	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/
	std::cout << std::endl;

	std::cout << "Switching times are: [" << switchingTimes[0] << ", ";
	for (size_t i=1; i<switchingTimes.size()-1; i++)  std::cout << switchingTimes[i] << ", ";
	std::cout << switchingTimes.back() << "]\n";

	std::cout << "MP Switching times are: [" << switchingTimes_mp[0] << ", ";
	for (size_t i=1; i<switchingTimes_mp.size()-1; i++)  std::cout << switchingTimes_mp[i] << ", ";
	std::cout << switchingTimes_mp.back() << "]\n";

	std::cout << "The total cost: " << totalCost << std::endl;
	std::cout << "The total cost in the test rollout: " << rolloutCost << std::endl;
	std::cout << "The total cost derivative: " << costFunctionDerivative.transpose() << std::endl;
	std::cout << "The MP total cost: " << totalCost_mp << std::endl;
	std::cout << "The MP total cost in the test rollout: " << rolloutCost_mp << std::endl;
	std::cout << "The MP total cost derivative: " << costFunctionDerivative_mp.transpose() << std::endl;

	GSLQP<2,1,2,3>::eigen_scalar_array_t timeEigenTrajectory;
	GSLQP<2,1,2,3>::state_vector_array_t stateTrajectory;
	GSLQP<2,1,2,3>::control_vector_array_t inputTrajectory;

	for (size_t i=0; i<switchingTimes.size()-1; i++)  {

		for (size_t k=0; k<timeTrajectoriesStock[i].size(); k++)  {

			timeEigenTrajectory.push_back((Eigen::MatrixXd(1,1) << timeTrajectoriesStock[i][k]).finished());
			stateTrajectory.push_back(stateTrajectoriesStock[i][k]);
			inputTrajectory.push_back(controlTrajectoriesStock[i][k]);
		}
	}


	// Sensitivity2SwitchingTime
	std::vector<GSLQP<2,1,2,3>::scalar_array_t> sensitivityTimeTrajectoriesStock;
	std::vector<GSLQP<2,1,2,3>::nabla_output_matrix_array_t> sensitivityStateTrajectoriesStock;
	std::vector<GSLQP<2,1,2,3>::nabla_input_matrix_array_t> sensitivityInputTrajectoriesStock;
	gslqp.getRolloutSensitivity2SwitchingTime(sensitivityTimeTrajectoriesStock, sensitivityStateTrajectoriesStock, sensitivityInputTrajectoriesStock);

	GSLQP<2,1,2,3>::eigen_scalar_array_t sensitivityTimeTrajectory;
	GSLQP<2,1,2,3>::nabla_output_matrix_array_t sensitivityStateTrajectory;
	for (size_t i=0; i<switchingTimes.size()-1; i++)  {
		for (size_t k=0; k<sensitivityTimeTrajectoriesStock[i].size(); k++)  {
			sensitivityTimeTrajectory.push_back((Eigen::MatrixXd(1,1) << sensitivityTimeTrajectoriesStock[i][k]).finished());
			sensitivityStateTrajectory.push_back(sensitivityStateTrajectoriesStock[i][k]);
		}
	}


	PathTweaker pathTweaker(argv);

	std::string resultDir = pathTweaker.getDirectory() +"/src/c_ocs2/cereal/test/exp1_test";

	std::cout << "Saving to directory " << resultDir << std::endl;

	std::string stateFile = resultDir + "/exp1State.xml";
	std::string timeFile = resultDir + "/exp1Time.xml";
	std::string inputFile = resultDir + "/exp1Input.xml";
	std::string stateSensitivityFile = resultDir + "/exp1StateSensitivity.xml";
	std::string timeSensitivityFile = resultDir + "/exp1TimeSensitivity.xml";

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

}

