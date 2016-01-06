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
#include "dynamics/DerivativesBase.h"
#include "costs/CostFunctionBase.h"

#include "integration/Integrator.h"

#include "GSLQ/RiccatiEquations.h"


template <size_t STATE_DIM, size_t INPUT_DIM>
class GLQP
{
public:
	typedef Dimensions<STATE_DIM, INPUT_DIM> DIMENSIONS;
	typedef typename DIMENSIONS::controller_t controller_t;
	typedef typename DIMENSIONS::scalar_t 		scalar_t;
	typedef typename DIMENSIONS::eigen_scalar_t       eigen_scalar_t;
	typedef typename DIMENSIONS::eigen_scalar_array_t eigen_scalar_array_t;
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


	GLQP(const std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM> > >& subsystemDynamicsPtr,
			const std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM> > >& subsystemDerivativesPtr,
			const std::vector<std::shared_ptr<CostFunctionBase<STATE_DIM, INPUT_DIM> > >& subsystemCostFunctionsPtr,
			const state_vector_array_t&   stateOperatingPoints,
			const control_vector_array_t& inputOperatingPoints,
			const std::vector<size_t>& systemStockIndex)
		: numSubsystems_(systemStockIndex.size()),
		  subsystemDynamicsPtrStock(numSubsystems_),
		  subsystemDerivativesPtrStock_(numSubsystems_),
		  subsystemCostFunctionsPtrStock_(numSubsystems_),
		  stateOperatingPointsStock_(numSubsystems_),
		  inputOperatingPointsStock_(numSubsystems_),
		  subsystemSimulatorsStockPtr_(numSubsystems_),
		  controllersStock_(numSubsystems_),
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
		  switchingTimes_(numSubsystems_+1)
	{

		if (subsystemDynamicsPtr.size() != subsystemDerivativesPtr.size())
			throw std::runtime_error("Number of subsystem derivaties is not equal to the number of subsystems.");
		if (subsystemDynamicsPtr.size() != subsystemCostFunctionsPtr.size())
			throw std::runtime_error("Number of cost functions is not equal to the number of subsystems.");
		if (subsystemDynamicsPtr.size() != stateOperatingPoints.size())
			throw std::runtime_error("Number of state operating points is not equal to the number of subsystems.");
		if (subsystemDynamicsPtr.size() != inputOperatingPoints.size())
			throw std::runtime_error("Number of input operating points is not equal to the number of subsystems.");
		if (subsystemDynamicsPtr.size()-1 != *std::max_element(systemStockIndex.begin(), systemStockIndex.end()))
			throw std::runtime_error("systemStockIndex points to non-existing subsystem");

		for (int i=0; i<numSubsystems_; i++) {

			subsystemDynamicsPtrStock[i] = subsystemDynamicsPtr[systemStockIndex[i]]->clone();
			subsystemDerivativesPtrStock_[i] = subsystemDerivativesPtr[systemStockIndex[i]]->clone();
			subsystemCostFunctionsPtrStock_[i] = subsystemCostFunctionsPtr[systemStockIndex[i]]->clone();

			stateOperatingPointsStock_[i] = stateOperatingPoints[systemStockIndex[i]];
			inputOperatingPointsStock_[i] = inputOperatingPoints[systemStockIndex[i]];

			subsystemSimulatorsStockPtr_[i] = std::make_shared<ODE45<STATE_DIM> >(subsystemDynamicsPtrStock[i]);
		}
	}

	~GLQP() {}

	void rollout(const state_vector_t& initState,
			const std::vector<controller_t>& controllersStock,
			std::vector<scalar_array_t>& timeTrajectoriesStock,
			std::vector<state_vector_array_t>& stateTrajectoriesStock,
			std::vector<control_vector_array_t>& controlTrajectoriesStock);

	void rolloutCost(const std::vector<scalar_array_t>& timeTrajectoriesStock,
			const std::vector<state_vector_array_t>& stateTrajectoriesStock,
			const std::vector<control_vector_array_t>& controlTrajectoriesStock,
			scalar_t& totalCost);

	void getController(std::vector<controller_t>& controllersStock) { controllersStock = controllersStock_;}

	void SolveRiccatiEquation(const std::vector<scalar_t>& switchingTimes);


protected:
	void approximateOptimalControlProblem();

	void calculatecontroller(const scalar_t& learningRate, std::vector<controller_t>& controllersStock);

	void transformeLocalValueFuntion2Global();


private:
	size_t numSubsystems_;

	std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM> > > subsystemDynamicsPtrStock;
	std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM> > > subsystemDerivativesPtrStock_;
	std::vector<std::shared_ptr<CostFunctionBase<STATE_DIM, INPUT_DIM> > > subsystemCostFunctionsPtrStock_;

	state_vector_array_t   stateOperatingPointsStock_;
	control_vector_array_t inputOperatingPointsStock_;

	std::vector<std::shared_ptr<ODE45<STATE_DIM> > > subsystemSimulatorsStockPtr_;

	std::vector<controller_t> controllersStock_;

	state_matrix_array_t        AmStock_;
	control_gain_matrix_array_t BmStock_;

	eigen_scalar_t qFinal_;
	state_vector_t QvFinal_;
	state_matrix_t QmFinal_;
	eigen_scalar_array_t qStock_;
	state_vector_array_t QvStock_;
	state_matrix_array_t QmStock_;
	control_vector_array_t RvStock_;
	control_matrix_array_t RmStock_;
	control_feedback_array_t PmStock_;

	std::vector<scalar_array_t> 	  timeTrajectoryStock_;
	std::vector<eigen_scalar_array_t> sTrajectoryStock_;
	std::vector<state_vector_array_t> SvTrajectoryStock_;
	std::vector<state_matrix_array_t> SmTrajectoryStock_;

	scalar_array_t switchingTimes_;


};

#include "implementation/GLQP.h"


#endif /* GLQP_H_ */
