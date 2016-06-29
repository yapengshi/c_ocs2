/*
 * loadTask.h
 *
 *  Created on: 29.06.2016
 *      Author: mgiftthaler
 */

#ifndef OCS2_LOADTASK_H_
#define OCS2_LOADTASK_H_

#include <Eigen/Dense>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>

#include "Dimensions.h"
#include "GSLQ/GSLQP.h"

#include "GSLQ/GSLQP.h"
#include "GSLQ/GLQP.h"
#include "GSLQ/SLQP.h"

namespace ocs2{

template <typename Derived>
void loadMatrix(const std::string& filename, const std::string& matrixName, Eigen::MatrixBase<Derived>& matrix)
{
	size_t rows = matrix.rows();
	size_t cols = matrix.cols();

	boost::property_tree::ptree pt;
	boost::property_tree::read_info(filename, pt);

	double scaling = pt.get<double>(matrixName + ".scaling", 1);

	for (size_t i=0; i<rows; i++)
	{
		for (size_t j=0; j<cols; j++)
		{
			matrix(i,j) = scaling*pt.get<double>(matrixName + "." + "(" +std::to_string(i) + "," + std::to_string(j) + ")" , 0.0);
		}
	}
}

template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM>
void loadOptions(const std::string& filename, typename ocs2::GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM,1>::Options_t& opt, bool verbose = true)
{
	boost::property_tree::ptree pt;
	boost::property_tree::read_info(filename, pt);

	try	{opt.maxIterationGSLQP_ = pt.get<int>("ocs2.maxIterationGSLQP");}
	catch (const std::exception& e){
		if (verbose)
			std::cout << " #### Option loader : option 'maxIterationGSLQP' not specified, using default " << opt.maxIterationGSLQP_ << std::endl;
	}

	try	{opt.minLearningRateGSLQP_ = pt.get<double>("ocs2.minLearningRateGSLQP");}
	catch (const std::exception& e){
		if (verbose)
			std::cout << " #### Option loader : option 'minLearningRateGSLQP' not specified, using default " << opt.minLearningRateGSLQP_ << std::endl;
	}

	try	{opt.maxLearningRateGSLQP_ = pt.get<double>("ocs2.maxLearningRateGSLQP");}
	catch (const std::exception& e){
		if (verbose)
			std::cout << " #### Option loader : option 'maxLearningRateGSLQP' not specified, using default " << opt.maxLearningRateGSLQP_ << std::endl;
	}

	try	{opt.minRelCostGSLQP_ = pt.get<double>("ocs2.minRelCostGSLQP");}
	catch (const std::exception& e){
		if (verbose)
			std::cout << " #### Option loader : option 'minRelCostGSLQP' not specified, using default " << opt.minRelCostGSLQP_ << std::endl;
	}

	try	{opt.meritFunctionRho_ = pt.get<double>("ocs2.meritFunctionRho");}
	catch (const std::exception& e){
		if (verbose)
			std::cout << " #### Option loader : option 'meritFunctionRho' not specified, using default " << opt.meritFunctionRho_ << std::endl;
	}

	try	{opt.constraintStepSize_ = pt.get<double>("ocs2.constraintStepSize");}
	catch (const std::exception& e){
		if (verbose)
			std::cout << " #### Option loader : option 'constraintStepSize' not specified, using default " << opt.constraintStepSize_ << std::endl;
	}

	try	{opt.lineSearchByMeritFuntion_ = (bool) pt.get<int>("ocs2.lineSearchByMeritFunction");}
	catch (const std::exception& e){
		if (verbose)
			std::cout << " #### Option loader : option 'lineSearchByMeritFunction' not specified, using default " << opt.lineSearchByMeritFuntion_ << std::endl;
	}

	try	{opt.dispayGSLQP_ = (bool) pt.get<int>("ocs2.dispayGSLQP");}
	catch (const std::exception& e){
		if (verbose)
			std::cout << " #### Option loader : option 'dispayGSLQP' not specified, using default " << opt.dispayGSLQP_ << std::endl;
	}

	try	{opt.warmStartGSLQP_ = (bool) pt.get<int>("ocs2.warmStartGSLQP");}
	catch (const std::exception& e){
		if (verbose)
			std::cout << " #### Option loader : option 'warmStartGSLQP' not specified, using default " << opt.warmStartGSLQP_ << std::endl;
	}

	try	{opt.AbsTolODE_ = pt.get<double>("ocs2.AbsTolODE");}
	catch (const std::exception& e){
		if (verbose)
			std::cout << " #### Option loader : option 'AbsTolODE' not specified, using default " << opt.AbsTolODE_ << std::endl;
	}

	try	{opt.RelTolODE_ = pt.get<double>("ocs2.RelTolODE");}
	catch (const std::exception& e){
		if (verbose)
			std::cout << " #### Option loader : option 'RelTolODE' not specified, using default " << opt.RelTolODE_ << std::endl;
	}

	try	{opt.simulationIsConstrained_ =(bool) pt.get<int>("ocs2.simulationIsConstrained");}
	catch (const std::exception& e){
		if (verbose)
			std::cout << " #### Option loader : option 'simulationIsConstrained' not specified, using default " << opt.simulationIsConstrained_ << std::endl;
	}

	try	{opt.minAbsConstraint1ISE_ = pt.get<double>("ocs2.minAbsConstraint1ISE");}
	catch (const std::exception& e){
		if (verbose)
			std::cout << " #### Option loader : option 'minAbsConstraint1ISE' not specified, using default " << opt.minAbsConstraint1ISE_ << std::endl;
	}

	try	{opt.minRelConstraint1ISE_ = pt.get<double>("ocs2.minRelConstraint1ISE");}
	catch (const std::exception& e){
		if (verbose)
			std::cout << " #### Option loader : option 'minRelConstraint1ISE' not specified, using default " << opt.minRelConstraint1ISE_ << std::endl;
	}

	try	{opt.displayIPOPT_ = (bool) pt.get<int>("ocs2.displayIPOPT");}
	catch (const std::exception& e){
		if (verbose)
			std::cout << " #### Option loader : option 'displayIPOPT' not specified, using default " << opt.displayIPOPT_ << std::endl;
	}

	try	{opt.tolIPOPT_ = pt.get<double>("ocs2.tolIPOPT");}
	catch (const std::exception& e){
		if (verbose)
			std::cout << " #### Option loader : option 'tolIPOPT' not specified, using default " << opt.tolIPOPT_ << std::endl;
	}

	try	{opt.acceptableTolIPOPT_ = pt.get<double>("ocs2.acceptableTolIPOPT");}
	catch (const std::exception& e){
		if (verbose)
			std::cout << " #### Option loader : option 'acceptableTolIPOPT' not specified, using default " << opt.acceptableTolIPOPT_ << std::endl;
	}

	try	{opt.maxIterationIPOPT_ = pt.get<int>("ocs2.maxIterationIPOPT");}
	catch (const std::exception& e){
		if (verbose)
			std::cout << " #### Option loader : option 'maxIterationIPOPT' not specified, using default " << opt.maxIterationIPOPT_ << std::endl;
	}

	try	{opt.minAcceptedSwitchingTimeDifference_ = pt.get<double>("ocs2.minAcceptedSwitchingTimeDifference");}
	catch (const std::exception& e){
		if (verbose)
			std::cout << " #### Option loader : option 'minAcceptedSwitchingTimeDifference' not specified, using default " << opt.minAcceptedSwitchingTimeDifference_ << std::endl;
	}

}


}

#endif /* OCS2_LOADTASK_H_ */
