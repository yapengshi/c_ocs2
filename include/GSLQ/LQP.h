/*
 * LQP.h
 *
 *  Created on: Aug 3, 2016
 *      Author: farbod
 */

#ifndef LQP_OCS2_H_
#define LQP_OCS2_H_

#include <vector>
#include <algorithm>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/StdVector>

#include "Dimensions.h"
#include "dynamics/ControlledSystemBase.h"
#include "dynamics/DerivativesBase.h"
#include "costs/CostFunctionBaseOCS2.h"

#include "integration/Integrator.h"
#include "GSLQ/SolveBVP.h"


namespace ocs2{

template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
class LQP
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	const bool INFO_ON_ = false;
	//
	typedef Dimensions<STATE_DIM, INPUT_DIM, OUTPUT_DIM> DIMENSIONS;
	typedef typename DIMENSIONS::controller_t controller_t;
	typedef typename DIMENSIONS::Options Options_t;
	typedef typename DIMENSIONS::scalar_t 		 scalar_t;
	typedef typename DIMENSIONS::scalar_array_t  scalar_array_t;
	typedef typename DIMENSIONS::eigen_scalar_t       eigen_scalar_t;
	typedef typename DIMENSIONS::eigen_scalar_array_t eigen_scalar_array_t;
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


	LQP(const std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > >& subsystemDynamicsPtr,
			const std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > >& subsystemDerivativesPtr,
			const std::vector<std::shared_ptr<CostFunctionBaseOCS2<OUTPUT_DIM, INPUT_DIM> > >& subsystemCostFunctionsPtr,
			const state_vector_array_t&   stateOperatingPoints,
			const control_vector_array_t& inputOperatingPoints,
			const std::vector<size_t>& systemStockIndex,
			const Options_t& options = Options_t::Options(),
			const bool& runAsInitializer = true)

		: subsystemDynamicsPtrStock_(NUM_Subsystems),
		  subsystemDerivativesPtrStock_(NUM_Subsystems),
		  subsystemCostFunctionsPtrStock_(NUM_Subsystems),
		  stateOperatingPointsStock_(NUM_Subsystems),
		  inputOperatingPointsStock_(NUM_Subsystems),
		  outputOperatingPointsStock_(NUM_Subsystems),
		  options_(options),
		  runAsInitializer_(runAsInitializer),
		  subsystemSimulatorsStockPtr_(NUM_Subsystems),
		  controllersStock_(NUM_Subsystems),
		  AmStock_(NUM_Subsystems),
		  BmStock_(NUM_Subsystems),
		  GvStock_(NUM_Subsystems),
		  QvStock_(NUM_Subsystems),
		  QmStock_(NUM_Subsystems),
		  RvStock_(NUM_Subsystems),
		  RmStock_(NUM_Subsystems),
		  RmInverseStock_(NUM_Subsystems),
		  PmStock_(NUM_Subsystems),
		  timeTrajectoryStock_(NUM_Subsystems),
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

			subsystemSimulatorsStockPtr_[i] = std::shared_ptr<ODE45<STATE_DIM>>( new ODE45<STATE_DIM>(subsystemDynamicsPtrStock_[i]) );
		}
	}

	~LQP() {}

	void rollout(const state_vector_t& initState,
			const std::vector<controller_t>& controllersStock,
			std::vector<scalar_array_t>& timeTrajectoriesStock,
			std::vector<state_vector_array_t>& stateTrajectoriesStock,
			std::vector<control_vector_array_t>& inputTrajectoriesStock,
			std::vector<output_vector_array_t>& outputTrajectoriesStock);

	void rollout(const state_vector_t& initState,
			const std::vector<controller_t>& controllersStock,
			std::vector<scalar_array_t>& timeTrajectoriesStock,
			std::vector<state_vector_array_t>& stateTrajectoriesStock,
			std::vector<control_vector_array_t>& inputTrajectoriesStock);

	void rolloutCost(const std::vector<scalar_array_t>& timeTrajectoriesStock,
			const std::vector<output_vector_array_t>& outputTrajectoriesStock,
			const std::vector<control_vector_array_t>& inputTrajectoriesStock,
			scalar_t& totalCost);

	void getController(std::vector<controller_t>& controllersStock);

	void run(const std::vector<scalar_t>& switchingTimes, const scalar_t& learningRate=1.0);

protected:
	void SolveRiccatiEquations();

	void approximateOptimalControlProblem();

	void calculatecontroller(const scalar_t& learningRate, std::vector<controller_t>& controllersStock);

	template <typename Derived>
	bool makePSD(Eigen::MatrixBase<Derived>& squareMatrix);

private:
	std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > > subsystemDynamicsPtrStock_;
	std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > > subsystemDerivativesPtrStock_;
	std::vector<std::shared_ptr<CostFunctionBaseOCS2<OUTPUT_DIM, INPUT_DIM> > > subsystemCostFunctionsPtrStock_;

	state_vector_array_t   stateOperatingPointsStock_;
	control_vector_array_t inputOperatingPointsStock_;
	output_vector_array_t  outputOperatingPointsStock_;

	Options_t options_;

	bool runAsInitializer_;

	std::vector<std::shared_ptr<ODE45<STATE_DIM> > > subsystemSimulatorsStockPtr_;

	std::vector<controller_t> controllersStock_;

	state_matrix_array_t        AmStock_;
	control_gain_matrix_array_t BmStock_;
	output_vector_array_t       GvStock_;

	output_vector_t QvFinal_;
	state_matrix_t QmFinal_;
	output_vector_array_t QvStock_;
	state_matrix_array_t QmStock_;
	control_vector_array_t RvStock_;
	control_matrix_array_t RmStock_;
	control_matrix_array_t RmInverseStock_;
	control_feedback_array_t PmStock_;

	std::vector<scalar_array_t> 	   timeTrajectoryStock_;
	std::vector<output_vector_array_t> SvTrajectoryStock_;
	std::vector<state_matrix_array_t>  SmTrajectoryStock_;

	scalar_array_t switchingTimes_;
};

} // namespace ocs2

#include "implementation/LQP.h"


#endif /* LQP_OCS2_H_ */
