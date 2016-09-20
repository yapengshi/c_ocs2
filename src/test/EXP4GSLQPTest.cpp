/*
 * EXP4GSLQPTest.cpp
 *
 *  Created on: Apr 25, 2016
 *      Author: farbod
 */


#include <iostream>
#include <cstdlib>

#include <fstream>
#include <cereal/archives/xml.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/Eigen.hpp>

#include "test/EXP4.h"
#include "GSLQ/GLQP.h"

#include <PathTweaker.h>

//#define LOG_DATA


using namespace ocs2;

int main (int argc, char* argv[])
{
	// subsystem dynamics
	std::vector<std::shared_ptr<ControlledSystemBase<6,4> > > subsystemDynamicsPtr {std::make_shared<EXP4_Sys>(), std::make_shared<EXP4_Sys>()};

	// subsystem derivatives
	std::vector<std::shared_ptr<DerivativesBase<6,4> > > subsystemDerivativesPtr {std::make_shared<EXP4_SysDerivative>(), std::make_shared<EXP4_SysDerivative>()};

	// subsystem cost functions
	std::vector<std::shared_ptr<CostFunctionBaseOCS2<6,4> > > subsystemCostFunctionsPtr {std::make_shared<EXP4_CostFunction>(), std::make_shared<EXP4_CostFunction>()};


	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/
	GSLQP<6,4,6,2>::state_vector_t initState;
	initState.setZero();
	initState.tail<2>() << 2.0, sqrt(5.0);

	GSLQP<6,4,6,2>::state_vector_array_t   stateOperatingPoints(2, initState);
	GSLQP<6,4,6,2>::control_vector_array_t inputOperatingPoints(2, GSLQP<6,4,6,2>::control_vector_t::Zero());
	std::vector<size_t> systemStockIndex {0, 1};

	std::vector<double> switchingTimes {0, 2, 4};
	if (argc>1)  switchingTimes[1] = std::atof(argv[1]);

	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/
	// GLQP
	GLQP<6,4,6,2> glqp(subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr, stateOperatingPoints, inputOperatingPoints, systemStockIndex);

	glqp.run(switchingTimes);

	// get controller
	std::vector<GLQP<6,4,6,2>::controller_t> controllersStock(3), controllersStock_mp(3);
	glqp.getController(controllersStock);


	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/
	GSLQP<6,4,6,2>::Options_t gslqpOptions;
	gslqpOptions.dispayGSLQP_ = 1;
	gslqpOptions.maxIterationGSLQP_ = 31;
	gslqpOptions.minRelCostGSLQP_ = 0.1;
	gslqpOptions.lineSearchByMeritFuntion_ = true;

	// GSLQ
	GSLQP<6,4,6,2> gslqp(subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr, controllersStock, systemStockIndex, gslqpOptions);

	// GSLQ MP Version
	gslqpOptions.useMultiThreading_ = true;
	GSLQP<6,4,6,2> gslqp_mp(subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr, controllersStock, systemStockIndex, gslqpOptions);

	// run two different versions
	gslqp.run(initState, switchingTimes);
	gslqp_mp.run(initState, switchingTimes);

	// get controllers
	gslqp.getController(controllersStock);
	gslqp_mp.getController(controllersStock_mp);

	// rollout
	std::vector<GSLQP<6,4,6,2>::scalar_array_t> timeTrajectoriesStock, timeTrajectoriesStock_mp;
	std::vector<GSLQP<6,4,6,2>::state_vector_array_t> stateTrajectoriesStock, stateTrajectoriesStock_mp;
	std::vector<GSLQP<6,4,6,2>::control_vector_array_t> controlTrajectoriesStock, controlTrajectoriesStock_mp;
	gslqp.rollout(initState, controllersStock, timeTrajectoriesStock, stateTrajectoriesStock, controlTrajectoriesStock);
	gslqp_mp.rollout(initState, controllersStock_mp, timeTrajectoriesStock_mp, stateTrajectoriesStock_mp, controlTrajectoriesStock_mp);

	// compute cost
	double rolloutCost, rolloutCost_mp;
	gslqp.calculateCostFunction(timeTrajectoriesStock, stateTrajectoriesStock, controlTrajectoriesStock, rolloutCost);
	gslqp_mp.calculateCostFunction(timeTrajectoriesStock_mp, stateTrajectoriesStock_mp, controlTrajectoriesStock_mp, rolloutCost_mp);

	// value funtion
	double totalCost, totalCost_mp;
	gslqp.getValueFuntion(0.0, initState, totalCost);
	gslqp_mp.getValueFuntion(0.0, initState, totalCost_mp);


	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/
	std::cout << std::endl;

	std::cout << "Switching times are: [" << switchingTimes[0] << ", ";
	for (size_t i=1; i<switchingTimes.size()-1; i++)  std::cout << switchingTimes[i] << ", ";
	std::cout << switchingTimes.back() << "]\n";

	std::cout << "The total cost: " << totalCost << std::endl;
	std::cout << "The total cost in the test rollout: " << rolloutCost << std::endl;
	std::cout << "The total cost mp: " << totalCost_mp << std::endl;
	std::cout << "The total cost mp in the test rollout: " << rolloutCost_mp << std::endl;

	GSLQP<6,4,6,2>::eigen_scalar_array_t timeEigenTrajectory;
	GSLQP<6,4,6,2>::state_vector_array_t stateTrajectory;
	GSLQP<6,4,6,2>::control_vector_array_t inputTrajectory;

	for (size_t i=0; i<switchingTimes.size()-1; i++)  {

		for (size_t k=0; k<timeTrajectoriesStock[i].size(); k++)  {

			timeEigenTrajectory.push_back((Eigen::MatrixXd(1,1) << timeTrajectoriesStock[i][k]).finished());
			stateTrajectory.push_back(stateTrajectoriesStock[i][k]);
			inputTrajectory.push_back(controlTrajectoriesStock[i][k]);
		}
	}

#ifdef LOG_DATA
	PathTweaker pathTweaker(argv);

	std::string resultDir = pathTweaker.getDirectory() +"/src/c_ocs2/cereal/test/exp4_test";

	std::cout << "Saving to directory " << resultDir << std::endl;
	std::string stateFile = resultDir + "/exp4State.xml";
	std::string timeFile = resultDir + "/exp4Time.xml";
	std::string inputFile = resultDir + "/exp4Input.xml";

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
	}

#endif

}




