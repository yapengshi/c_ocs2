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


int main (int argc, char* argv[])
{

	// subsystem dynamics
	std::vector<std::shared_ptr<ControlledSystemBase<2,1> > > subsystemDynamicsPtr {std::make_shared<EXP2_Sys1>(), std::make_shared<EXP2_Sys2>()};

	// subsystem derivatives
	std::vector<std::shared_ptr<DerivativesBase<2,1> > > subsystemDerivativesPtr {std::make_shared<EXP2_SysDerivative1>(), std::make_shared<EXP2_SysDerivative2>()};

	// subsystem cost functions
	std::vector<std::shared_ptr<CostFunctionBase<2,1> > > subsystemCostFunctionsPtr {std::make_shared<EXP2_CostFunction1>(), std::make_shared<EXP2_CostFunction2>()};

	GLQP<2,1>::state_vector_array_t   stateOperatingPoints(2, GLQP<2,1>::state_vector_t::Zero());
	GLQP<2,1>::control_vector_array_t inputOperatingPoints(2, GLQP<2,1>::control_vector_t::Zero());
	std::vector<size_t> systemStockIndex {0, 1};

	// GLQP
	GLQP<2,1> glqp(subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr, stateOperatingPoints, inputOperatingPoints, systemStockIndex);

	std::vector<double> switchingTimes {0, 0.184, 2};
	if (argc>1)  switchingTimes[1] = std::atof(argv[1]);
	glqp.SolveRiccatiEquation(switchingTimes);

	// get controller
	std::vector<GLQP<2,1>::controller_t> controllersStock(2);
	glqp.getController(controllersStock);

	// rollout
	Eigen::Vector2d initState(0.0, 2.0);
	std::vector<GLQP<2,1>::scalar_array_t> timeTrajectoriesStock;
	std::vector<GLQP<2,1>::state_vector_array_t> stateTrajectoriesStock;
	std::vector<GLQP<2,1>::control_vector_array_t> controlTrajectoriesStock;
	glqp.rollout(initState, controllersStock,
				timeTrajectoriesStock, stateTrajectoriesStock, controlTrajectoriesStock);

	// compute cost
	double rolloutCost;
	glqp.rolloutCost(timeTrajectoriesStock, stateTrajectoriesStock, controlTrajectoriesStock, rolloutCost);

	// value funtion
	double totalCost;
	glqp.getValueFuntion(0.0, initState, totalCost);


	std::cout << "Switching times are: [" << switchingTimes[0] << ", " << switchingTimes[1] << ", " << switchingTimes[2] << "]\n";
	std::cout << "The total cost: " << totalCost << std::endl;
	std::cout << "The total cost in the test rollout: " << rolloutCost << std::endl;


	GLQP<2,1>::eigen_scalar_array_t timeEigenTrajectory;
	GLQP<2,1>::state_vector_array_t stateTrajectory;
	GLQP<2,1>::control_vector_array_t inputTrajectory;

	for (size_t i=0; i<2; i++)  {

		for (size_t k=0; k<timeTrajectoriesStock[i].size(); k++)  {

			timeEigenTrajectory.push_back((Eigen::MatrixXd(1,1) << timeTrajectoriesStock[i][k]).finished());
			stateTrajectory.push_back(stateTrajectoriesStock[i][k]);
			inputTrajectory.push_back(controlTrajectoriesStock[i][k]);
		}
	}

	std::string resultDir = "/home/farbod/Programs/ct_ws/src/c_ocs2/cereal/test/exp2_test";
	std::string stateFile = resultDir + "/exp2State.xml";
	std::string timeFile = resultDir + "/exp2Time.xml";
	std::string inputFile = resultDir + "/exp2Input.xml";


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

}


