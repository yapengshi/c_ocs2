/*
 * GLQP.h
 *
 *  Created on: Jan 5, 2016
 *      Author: farbod
 */

#ifndef GLQP_H_
#define GLQP_H_

#include <iostream>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/StdVector>

#include "Dimensions.h"

#include "dynamics/ControlledSystemBase.h"
#include "dynamics/DerivativesBase.h"
#include "costs/CostFunctionBase.h"

#include "integration/Integrator.h"
#include "misc/LinearInterpolation.h"

#include "GSLQ/PartialRiccatiEquations.h"


template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
class GLQP
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	const bool INFO_ON_ = false;
	typedef PartialRiccatiEquations<OUTPUT_DIM, INPUT_DIM, NUM_Subsystems> RiccatiEquations_t;
	//
	typedef Dimensions<STATE_DIM, INPUT_DIM, OUTPUT_DIM> DIMENSIONS;
	typedef typename DIMENSIONS::controller_t controller_t;
	typedef typename DIMENSIONS::scalar_t 		scalar_t;
	typedef typename DIMENSIONS::eigen_scalar_t       eigen_scalar_t;
	typedef typename DIMENSIONS::eigen_scalar_array_t eigen_scalar_array_t;
	typedef typename DIMENSIONS::scalar_array_t scalar_array_t;
	typedef typename DIMENSIONS::state_vector_t 	  state_vector_t;
	typedef typename DIMENSIONS::state_vector_array_t state_vector_array_t;
	typedef typename DIMENSIONS::control_vector_t 		control_vector_t;
	typedef typename DIMENSIONS::control_vector_array_t control_vector_array_t;
	typedef typename DIMENSIONS::output_vector_t 	  output_vector_t;
	typedef typename DIMENSIONS::output_vector_array_t output_vector_array_t;
	typedef typename DIMENSIONS::control_feedback_t 	  control_feedback_t;
	typedef typename DIMENSIONS::control_feedback_array_t control_feedback_array_t;
	typedef typename DIMENSIONS::state_matrix_t 	  state_matrix_t;
	typedef typename DIMENSIONS::state_matrix_array_t state_matrix_array_t;
	typedef typename DIMENSIONS::control_matrix_t 		control_matrix_t;
	typedef typename DIMENSIONS::control_matrix_array_t control_matrix_array_t;
	typedef typename DIMENSIONS::control_gain_matrix_t 		 control_gain_matrix_t;
	typedef typename DIMENSIONS::control_gain_matrix_array_t control_gain_matrix_array_t;


	GLQP(const std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > >& subsystemDynamicsPtr,
			const std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > >& subsystemDerivativesPtr,
			const std::vector<std::shared_ptr<CostFunctionBase<OUTPUT_DIM, INPUT_DIM> > >& subsystemCostFunctionsPtr,
			const state_vector_array_t&   stateOperatingPoints,
			const control_vector_array_t& inputOperatingPoints,
			const std::vector<size_t>& systemStockIndex)

		: subsystemDynamicsPtrStock_(NUM_Subsystems),
		  subsystemDerivativesPtrStock_(NUM_Subsystems),
		  subsystemCostFunctionsPtrStock_(NUM_Subsystems),
		  stateOperatingPointsStock_(NUM_Subsystems),
		  inputOperatingPointsStock_(NUM_Subsystems),
		  outputOperatingPointsStock_(NUM_Subsystems),
		  subsystemSimulatorsStockPtr_(NUM_Subsystems),
		  controllersStock_(NUM_Subsystems),
		  AmStock_(NUM_Subsystems),
		  BmStock_(NUM_Subsystems),
		  qStock_(NUM_Subsystems),
		  QvStock_(NUM_Subsystems),
		  QmStock_(NUM_Subsystems),
		  RvStock_(NUM_Subsystems),
		  RmStock_(NUM_Subsystems),
		  PmStock_(NUM_Subsystems),
		  timeTrajectoryStock_(NUM_Subsystems),
		  sTrajectoryStock_(NUM_Subsystems),
		  SvTrajectoryStock_(NUM_Subsystems),
		  SmTrajectoryStock_(NUM_Subsystems),
		  switchingTimes_(NUM_Subsystems+1)
	{

		if (subsystemDynamicsPtr.size() != subsystemDerivativesPtr.size())
			throw std::runtime_error("Number of subsystem derivaties is not equal to the number of subsystems.");
		if (subsystemDynamicsPtr.size() != subsystemCostFunctionsPtr.size())
			throw std::runtime_error("Number of cost functions is not equal to the number of subsystems.");
		if (subsystemDynamicsPtr.size() != stateOperatingPoints.size())
			throw std::runtime_error("Number of state operating points is not equal to the number of subsystems.");
		if (subsystemDynamicsPtr.size() != inputOperatingPoints.size())
			throw std::runtime_error("Number of input operating points is not equal to the number of subsystems.");
		if (subsystemDynamicsPtr.size()-1 < *std::max_element(systemStockIndex.begin(), systemStockIndex.end()))
			throw std::runtime_error("systemStockIndex points to non-existing subsystem");
		if (systemStockIndex.size() != NUM_Subsystems)
			throw std::runtime_error("systemStockIndex has less elements then the number of subsystems");

		for (int i=0; i<NUM_Subsystems; i++) {

			subsystemDynamicsPtrStock_[i] = subsystemDynamicsPtr[systemStockIndex[i]]->clone();
			subsystemDerivativesPtrStock_[i] = subsystemDerivativesPtr[systemStockIndex[i]]->clone();
			subsystemCostFunctionsPtrStock_[i] = subsystemCostFunctionsPtr[systemStockIndex[i]]->clone();

			stateOperatingPointsStock_[i] = stateOperatingPoints[systemStockIndex[i]];
			inputOperatingPointsStock_[i] = inputOperatingPoints[systemStockIndex[i]];
			subsystemDynamicsPtrStock_[i]->computeOutput(0.0 /*time*/, stateOperatingPointsStock_[i], outputOperatingPointsStock_[i]);

			subsystemSimulatorsStockPtr_[i] = std::make_shared<ODE45<STATE_DIM> >(subsystemDynamicsPtrStock_[i]);

		}
	}

	~GLQP() {}

	void rollout(const state_vector_t& initState,
			const std::vector<controller_t>& controllersStock,
			std::vector<scalar_array_t>& timeTrajectoriesStock,
			std::vector<state_vector_array_t>& stateTrajectoriesStock,
			std::vector<control_vector_array_t>& controlTrajectoriesStock,
			std::vector<output_vector_array_t>& outputTrajectoriesStock);

	void rollout(const state_vector_t& initState,
			const std::vector<controller_t>& controllersStock,
			std::vector<scalar_array_t>& timeTrajectoriesStock,
			std::vector<state_vector_array_t>& stateTrajectoriesStock,
			std::vector<control_vector_array_t>& controlTrajectoriesStock);

	void rolloutCost(const std::vector<scalar_array_t>& timeTrajectoriesStock,
			const std::vector<output_vector_array_t>& outputTrajectoriesStock,
			const std::vector<control_vector_array_t>& controlTrajectoriesStock,
			scalar_t& totalCost);

	void getController(std::vector<controller_t>& controllersStock);

	void getValueFuntion(const scalar_t& time, const output_vector_t& output, scalar_t& valueFuntion);

	void run(const std::vector<scalar_t>& switchingTimes, const scalar_t& learningRate=1.0);


protected:
	void SolveRiccatiEquations();

	void approximateOptimalControlProblem();

	void calculatecontroller(const scalar_t& learningRate, std::vector<controller_t>& controllersStock);

	void transformeLocalValueFuntion2Global();

	template <typename Derived>
	bool makePSD(Eigen::MatrixBase<Derived>& squareMatrix);

private:
	std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > > subsystemDynamicsPtrStock_;
	std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > > subsystemDerivativesPtrStock_;
	std::vector<std::shared_ptr<CostFunctionBase<OUTPUT_DIM, INPUT_DIM> > > subsystemCostFunctionsPtrStock_;

	state_vector_array_t   stateOperatingPointsStock_;
	control_vector_array_t inputOperatingPointsStock_;
	output_vector_array_t  outputOperatingPointsStock_;

	std::vector<std::shared_ptr<ODE45<STATE_DIM> > > subsystemSimulatorsStockPtr_;

	std::vector<controller_t> controllersStock_;

	state_matrix_array_t        AmStock_;
	control_gain_matrix_array_t BmStock_;

	eigen_scalar_t qFinal_;
	output_vector_t QvFinal_;
	state_matrix_t QmFinal_;
	eigen_scalar_array_t qStock_;
	output_vector_array_t QvStock_;
	state_matrix_array_t QmStock_;
	control_vector_array_t RvStock_;
	control_matrix_array_t RmStock_;
	control_feedback_array_t PmStock_;

	std::vector<scalar_array_t> 	  timeTrajectoryStock_;
	std::vector<eigen_scalar_array_t> sTrajectoryStock_;
	std::vector<output_vector_array_t> SvTrajectoryStock_;
	std::vector<state_matrix_array_t> SmTrajectoryStock_;

	scalar_array_t switchingTimes_;


};

#include "implementation/GLQP.h"


#endif /* GLQP_H_ */
