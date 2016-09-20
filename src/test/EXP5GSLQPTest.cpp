/*
 * EXP1GSLQPTest.cpp
 *
 *  Created on: Sept 20, 2016
 *      Author: markus (based on original version by farbod)
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
	gslqpOptions.stateConstraintPenaltyBase_ = 1.0;
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
//	gslqp.run(initState, switchingTimes);
//	gslqp_mp.run(initState, switchingTimes_mp);


//	// get controller
//	gslqp.getController(controllersStock);
//	gslqp_mp.getController(controllersStock_mp);
//
//	// rollout both versions
//	std::vector<GSLQP<2,2,2,1>::scalar_array_t> timeTrajectoriesStock, timeTrajectoriesStock_mp;
//	std::vector<GSLQP<2,2,2,1>::state_vector_array_t> stateTrajectoriesStock, stateTrajectoriesStock_mp;
//	std::vector<GSLQP<2,2,2,1>::control_vector_array_t> controlTrajectoriesStock, controlTrajectoriesStock_mp;
//	gslqp.rollout(initState, controllersStock, timeTrajectoriesStock, stateTrajectoriesStock, controlTrajectoriesStock);
//	gslqp_mp.rollout(initState, controllersStock_mp, timeTrajectoriesStock_mp, stateTrajectoriesStock_mp, controlTrajectoriesStock_mp);
//
//	// compute cost for both versions
//	double rolloutCost, rolloutCost_mp;
//	gslqp.calculateCostFunction(timeTrajectoriesStock, stateTrajectoriesStock, controlTrajectoriesStock, rolloutCost);
//	gslqp_mp.calculateCostFunction(timeTrajectoriesStock_mp, stateTrajectoriesStock_mp, controlTrajectoriesStock_mp, rolloutCost_mp);
//
//	// value function for both versions
//	double totalCost;
//	double totalCost_mp;
//	gslqp.getValueFuntion(0.0, initState, totalCost);
//	gslqp_mp.getValueFuntion(0.0, initState, totalCost_mp);
//
//	// value function derivative for both versions
//	Eigen::Matrix<double,2,1> costFunctionDerivative, costFunctionDerivative_mp;
////	gslqp.getCostFuntionDerivative(costFunctionDerivative);
////	gslqp_mp.getCostFuntionDerivative(costFunctionDerivative_mp);
//
//
//	/******************************************************************************************************/
//	/******************************************************************************************************/
//	/******************************************************************************************************/
//	std::cout << std::endl;
//
//	std::cout << "Single core switching times are: [" << switchingTimes[0] << ", ";
//	for (size_t i=1; i<switchingTimes.size()-1; i++)  std::cout << switchingTimes[i] << ", ";
//	std::cout << switchingTimes.back() << "]\n";
//
//	std::cout << "MP switching times are: [" << switchingTimes_mp[0] << ", ";
//	for (size_t i=1; i<switchingTimes_mp.size()-1; i++)  std::cout << switchingTimes_mp[i] << ", ";
//	std::cout << switchingTimes_mp.back() << "]\n";
//
//	std::cout << "The single core total cost: " << totalCost << std::endl;
//	std::cout << "The MP total cost: " << totalCost_mp << std::endl;
//
//	std::cout << "The single core total cost in the test rollout: " << rolloutCost << std::endl;
//	std::cout << "The MP total cost in the test rollout: " << rolloutCost_mp << std::endl;
//
//	std::cout << "The single core total cost derivative: " << costFunctionDerivative.transpose() << std::endl;
//	std::cout << "The MP total cost derivative: " << costFunctionDerivative_mp.transpose() << std::endl;
//
//	GSLQP<2,2,2,1>::eigen_scalar_array_t timeEigenTrajectory, timeEigenTrajectory_mp;
//	GSLQP<2,2,2,1>::state_vector_array_t stateTrajectory, stateTrajectory_mp;
//	GSLQP<2,2,2,1>::control_vector_array_t inputTrajectory, inputTrajectory_mp;
//
//	for (size_t i=0; i<switchingTimes.size()-1; i++)  {
//
//		for (size_t k=0; k<timeTrajectoriesStock[i].size(); k++)  {
//			timeEigenTrajectory.push_back((Eigen::MatrixXd(1,1) << timeTrajectoriesStock[i][k]).finished());
//			stateTrajectory.push_back(stateTrajectoriesStock[i][k]);
//			inputTrajectory.push_back(controlTrajectoriesStock[i][k]);
//		}
//
//		for (size_t k=0; k<timeTrajectoriesStock_mp[i].size(); k++)  {
//			timeEigenTrajectory_mp.push_back((Eigen::MatrixXd(1,1) << timeTrajectoriesStock_mp[i][k]).finished());
//			stateTrajectory_mp.push_back(stateTrajectoriesStock_mp[i][k]);
//			inputTrajectory_mp.push_back(controlTrajectoriesStock_mp[i][k]);
//		}
//	}
//
//
//	// Sensitivity2SwitchingTime
//	std::vector<GSLQP<2,2,2,1>::scalar_array_t> sensitivityTimeTrajectoriesStock, sensitivityTimeTrajectoriesStock_mp;
//	std::vector<GSLQP<2,2,2,1>::nabla_output_matrix_array_t> sensitivityStateTrajectoriesStock, sensitivityStateTrajectoriesStock_mp;
//	std::vector<GSLQP<2,2,2,1>::nabla_input_matrix_array_t> sensitivityInputTrajectoriesStock, sensitivityInputTrajectoriesStock_mp;
//	gslqp.getRolloutSensitivity2SwitchingTime(sensitivityTimeTrajectoriesStock, sensitivityStateTrajectoriesStock, sensitivityInputTrajectoriesStock);
//	gslqp_mp.getRolloutSensitivity2SwitchingTime(sensitivityTimeTrajectoriesStock_mp, sensitivityStateTrajectoriesStock_mp, sensitivityInputTrajectoriesStock_mp);
//
//	GSLQP<2,2,2,1>::eigen_scalar_array_t sensitivityTimeTrajectory, sensitivityTimeTrajectory_mp;
//	GSLQP<2,2,2,1>::nabla_output_matrix_array_t sensitivityStateTrajectory, sensitivityStateTrajectory_mp;
//	for (size_t i=0; i<switchingTimes.size()-1; i++)  {
//		for (size_t k=0; k<sensitivityTimeTrajectoriesStock[i].size(); k++)  {
//			sensitivityTimeTrajectory.push_back((Eigen::MatrixXd(1,1) << sensitivityTimeTrajectoriesStock[i][k]).finished());
//			sensitivityStateTrajectory.push_back(sensitivityStateTrajectoriesStock[i][k]);
//		}
//		for (size_t k=0; k<sensitivityTimeTrajectoriesStock_mp[i].size(); k++)  {
//			sensitivityTimeTrajectory_mp.push_back((Eigen::MatrixXd(1,1) << sensitivityTimeTrajectoriesStock_mp[i][k]).finished());
//			sensitivityStateTrajectory_mp.push_back(sensitivityStateTrajectoriesStock_mp[i][k]);
//		}
//	}
//
//
//	PathTweaker pathTweaker(argv);
//
//	std::string resultDir = pathTweaker.getDirectory() +"/src/c_ocs2/cereal/test/exp1_test";
//
//	std::cout << "Saving to directory " << resultDir << std::endl;
//
//	std::string stateFile 	= resultDir + "/exp1State.xml";
//	std::string timeFile 	= resultDir + "/exp1Time.xml";
//	std::string inputFile 	= resultDir + "/exp1Input.xml";
//	std::string stateSensitivityFile 	= resultDir + "/exp1StateSensitivity.xml";
//	std::string timeSensitivityFile 	= resultDir + "/exp1TimeSensitivity.xml";
//
//	std::string stateFile_mp 	= resultDir + "/exp1State_mp.xml";
//	std::string timeFile_mp 	= resultDir + "/exp1Time_mp.xml";
//	std::string inputFile_mp 	= resultDir + "/exp1Input_mp.xml";
//	std::string stateSensitivityFile_mp = resultDir + "/exp1StateSensitivity_mp.xml";
//	std::string timeSensitivityFile_mp 	= resultDir + "/exp1TimeSensitivity_mp.xml";
//
//	{ // we need these brackets to make sure the archive goes out of scope and flushes
//
//		// serialize single core data
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
//
//		// serialize multi core data
//		std::ofstream xmlState_mp(stateFile_mp);
//		cereal::XMLOutputArchive archive_state_xml_mp(xmlState_mp);
//		archive_state_xml_mp(CEREAL_NVP(stateTrajectory_mp));
//
//		std::ofstream xmlTime_mp(timeFile_mp);
//		cereal::XMLOutputArchive archive_time_xml_mp(xmlTime_mp);
//		archive_time_xml_mp(CEREAL_NVP(timeEigenTrajectory_mp));
//
//		std::ofstream xmlInput_mp(inputFile_mp);
//		cereal::XMLOutputArchive archive_input_xml_mp(xmlInput_mp);
//		archive_input_xml_mp(CEREAL_NVP(inputTrajectory_mp));
//
//		std::ofstream xmlStateSensitivity_mp(stateSensitivityFile_mp);
//		cereal::XMLOutputArchive archive_stateSensitivity_xml_mp(xmlStateSensitivity_mp);
//		archive_stateSensitivity_xml_mp(CEREAL_NVP(sensitivityStateTrajectory_mp));
//
//		std::ofstream xmlTimeSensitivity_mp(timeSensitivityFile_mp);
//		cereal::XMLOutputArchive archive_timeSensitivity_xml_mp(xmlTimeSensitivity_mp);
//		archive_timeSensitivity_xml_mp(CEREAL_NVP(sensitivityTimeTrajectory_mp));
//	}

}

