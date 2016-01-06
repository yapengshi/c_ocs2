/*
 * GLQP.h
 *
 *  Created on: Jan 5, 2016
 *      Author: farbod
 */

#ifndef GLQP_H_
#define GLQP_H_

#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/StdVector>

#include "Dimensions.h"

#include "dynamics/ControlledSystemBase.h"
#include "dynamics/DerivativesBase.hpp"
#include "costs/CostFunctionBase.hpp"

#include "integration/Integrator.h"


template <size_t STATE_DIM, size_t INPUT_DIM>
class GLQP
{
public:
	typedef Dimensions<STATE_DIM, INPUT_DIM> DIMENSIONS;
	typedef typename DIMENSIONS::controller_t controller_t;
	typedef typename DIMENSIONS::scalar_t 		scalar_t;
	typedef typename DIMENSIONS::scalar_array_t scalar_array_t;
	typedef typename DIMENSIONS::state_vector_t 	  state_vector_t;
	typedef typename DIMENSIONS::state_vector_array_t state_vector_array_t;
	typedef typename DIMENSIONS::control_vector_t 		control_vector_t;
	typedef typename DIMENSIONS::control_vector_array_t control_vector_array_t;
	typedef typename DIMENSIONS::control_feedback_t 	  control_feedback_t;
	typedef typename DIMENSIONS::control_feedback_array_t control_feedback_array_t;
	typedef typename DIMENSIONS::state_matrix_t 	  state_matrix_t;
	typedef typename DIMENSIONS::state_matrix_array_t state_matrix_array_t;
	typedef typename DIMENSIONS::control_matrix_t 		control_matrix_t;
	typedef typename DIMENSIONS::control_matrix_array_t control_matrix_array_t;
	typedef typename DIMENSIONS::control_gain_matrix_t 		 control_gain_matrix_t;
	typedef typename DIMENSIONS::control_gain_matrix_array_t control_gain_matrix_array_t;


	GLQP( std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM> > > subsystemDynamicsPtr,
			std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM> > > subsystemDerivativesPtr,
			std::vector<std::shared_ptr<CostFunctionBase<STATE_DIM, INPUT_DIM> > > subsystemCostFunctionsPtr,
			state_vector_array_t   stateOperatingPoints,
			control_vector_array_t inputOperatingPoints,
			std::vector<size_t> systemStockIndex)
		: numSubsystems_(systemStockIndex.size()),
		  subsystemDynamicsPtrStock(numSubsystems_),
		  subsystemDerivativesPtrStock_(numSubsystems_),
		  subsystemCostFunctionsPtrStock_(numSubsystems_),
		  subsystemSimulatorsStock_(numSubsystems_),
		  controllersStock_(numSubsystems_),
		  stateOperatingPointsStock_(numSubsystems_),
		  inputOperatingPointsStock_(numSubsystems_),
		  AmStock_(numSubsystems_),
		  BmStock_(numSubsystems_),
		  qStock_(numSubsystems_),
		  QvStock_(numSubsystems_),
		  QmStock_(numSubsystems_),
		  RvStock_(numSubsystems_),
		  RmStock_(numSubsystems_),
		  PmStock_(numSubsystems_),
		  timeTrajectoryStock_(numSubsystems_),
		  sTrajectoryStock_(numSubsystems_),
		  SvTrajectoryStock_(numSubsystems_),
		  SmTrajectoryStock_(numSubsystems_),
		  maxIteration_(10)
	{

		if (subsystemDynamicsPtr.size() != subsystemDerivativesPtr.size())
			throw std::runtime_error("Number of subsystem derivaties is not equal to the number of subsystems.");
		if (subsystemDynamicsPtr.size() != subsystemCostFunctionsPtr.size())
			throw std::runtime_error("Number of cost functions is not equal to the number of subsystems.");
		if (subsystemDynamicsPtr.size() != stateOperatingPoints.size())
			throw std::runtime_error("Number of state operating points is not equal to the number of subsystems.");
		if (subsystemDynamicsPtr.size() != inputOperatingPoints.size())
			throw std::runtime_error("Number of input operating points is not equal to the number of subsystems.");
		if (subsystemDynamicsPtr.size() != *std::max_element(systemStockIndex.begin(), systemStockIndex.end()))
			throw std::runtime_error("systemStockIndex points to non-existing subsystem");


		for (int i=0; i<numSubsystems_; i++) {

			subsystemDynamicsPtrStock[i] = subsystemDynamicsPtr[systemStockIndex[i]]->clone();
			subsystemDerivativesPtrStock_[i] = subsystemDerivativesPtr[systemStockIndex[i]]->clone();
			subsystemCostFunctionsPtrStock_[i] = subsystemCostFunctionsPtr[systemStockIndex[i]]->clone();

			stateOperatingPointsStock_[i] = stateOperatingPoints[systemStockIndex[i]];
			inputOperatingPointsStock_[i] = inputOperatingPoints[systemStockIndex[i]];

			subsystemSimulatorsStock_[i] = ODE45<STATE_DIM>(subsystemDynamicsPtrStock[i]);
		}
	}

	~GLQP() {}


	void rollout(const state_vector_t& initState,
			const std::vector<scalar_t>& switchingTimes,
			const std::vector<controller_t>& controllersStock,
			std::vector<scalar_array_t>& timeTrajectoriesStock,
			std::vector<state_vector_array_t>& stateTrajectoriesStock,
			std::vector<control_vector_array_t>& controlTrajectoriesStock);

	void rolloutCost(const std::vector<scalar_array_t>& timeTrajectoriesStock,
			const std::vector<state_vector_array_t>& stateTrajectoriesStock,
			const std::vector<control_vector_array_t>& controlTrajectoriesStock,
			scalar_t& totalCost);

	void approximateOptimalControlProblem();

	void calculatecontroller(const scalar_t& learningRate, std::vector<controller_t>& controllersStock);

	void transformeLocalValueFuntion2Global();

	void getcontroller(std::vector<controller_t>& controllersStock) { controllersStock = controllersStock_;}

	void SolveRiccatiEquation(const std::vector<scalar_t>& switchingTimes);


private:
	std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM> > > subsystemDynamicsPtrStock;
	std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM> > > subsystemDerivativesPtrStock_;
	std::vector<std::shared_ptr<CostFunctionBase<STATE_DIM, INPUT_DIM> > > subsystemCostFunctionsPtrStock_;

	state_vector_array_t   stateOperatingPointsStock_;
	control_vector_array_t inputOperatingPointsStock_;

	std::vector<ODE45<STATE_DIM> > subsystemSimulatorsStock_;

	std::vector<controller_t> controllersStock_;

	size_t numSubsystems_;
	size_t maxIteration_;

	state_matrix_array_t        AmStock_;
	control_gain_matrix_array_t BmStock_;

	scalar_t       qFinal_;
	state_vector_t QvFinal_;
	state_matrix_t QmFinal_;
	scalar_array_t 	  	 qStock_;
	state_vector_array_t QvStock_;
	state_matrix_array_t QmStock_;
	control_vector_array_t RvStock_;
	control_matrix_array_t RmStock_;
	control_feedback_array_t PmStock_;

	std::vector<scalar_array_t> 	  timeTrajectoryStock_;
	std::vector<scalar_array_t> 	  sTrajectoryStock_;
	std::vector<state_vector_array_t> SvTrajectoryStock_;
	std::vector<state_matrix_array_t> SmTrajectoryStock_;


};

#include "implementation/GLQP.h"



#endif /* GLQP_H_ */