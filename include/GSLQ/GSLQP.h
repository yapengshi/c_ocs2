/*
 * GSLQP.h
 *
 *  Created on: Dec 18, 2015
 *      Author: farbod
 */

#ifndef GSLQP_H_
#define GSLQP_H_

#include <vector>
#include <array>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/StdVector>

#include "Dimensions.h"

#include "dynamics/ControlledSystemBase.h"
#include "dynamics/DerivativesBase.h"
#include "costs/CostFunctionBase.h"

#include "integration/Integrator.h"
#include "misc/LinearInterpolation.h"

#include "GSLQ/SequentialRiccatiEquations.h"
#include "GSLQ/SequentialErrorEquation.h"
#include "GSLQ/FullSequentialRiccatiEquations.h"
#include "GSLQ/RolloutSensitivityEquations.h"

template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
class GSLQP
{
public:
	typedef SequentialRiccatiEquations<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS> RiccatiEquations_t;
	typedef SequentialErrorEquation<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS> ErrorEquation_t;
	typedef FullSequentialRiccatiEquations<OUTPUT_DIM, INPUT_DIM, NUM_SUBSYSTEMS> FullRiccatiEquations_t;

	typedef Dimensions<STATE_DIM, INPUT_DIM, OUTPUT_DIM> DIMENSIONS;
	typedef typename DIMENSIONS::controller_t controller_t;
	typedef typename DIMENSIONS::Options Options_t;
	typedef typename DIMENSIONS::scalar_t 		scalar_t;
	typedef typename DIMENSIONS::scalar_array_t scalar_array_t;
	typedef typename DIMENSIONS::eigen_scalar_t       eigen_scalar_t;
	typedef typename DIMENSIONS::eigen_scalar_array_t eigen_scalar_array_t;
	typedef typename DIMENSIONS::state_vector_t 	  state_vector_t;
	typedef typename DIMENSIONS::state_vector_array_t state_vector_array_t;
	typedef typename DIMENSIONS::output_vector_t 	   output_vector_t;
	typedef typename DIMENSIONS::output_vector_array_t output_vector_array_t;
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
	typedef typename DIMENSIONS::constraint1_vector_t       constraint1_vector_t;
	typedef typename DIMENSIONS::constraint1_vector_array_t constraint1_vector_array_t;
	typedef typename DIMENSIONS::constraint1_matrix_t       constraint1_matrix_t;
	typedef typename DIMENSIONS::constraint1_matrix_array_t constraint1_matrix_array_t;
	typedef typename DIMENSIONS::constraint1_state_matrix_t       constraint1_state_matrix_t;
	typedef typename DIMENSIONS::constraint1_state_matrix_array_t constraint1_state_matrix_array_t;
	typedef typename DIMENSIONS::constraint1_control_matrix_t       constraint1_control_matrix_t;
	typedef typename DIMENSIONS::constraint1_control_matrix_array_t constraint1_control_matrix_array_t;
	typedef typename DIMENSIONS::control_constraint1_matrix_t       control_constraint1_matrix_t;
	typedef typename DIMENSIONS::control_constraint1_matrix_array_t control_constraint1_matrix_array_t;

	typedef RolloutSensitivityEquations<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS> RolloutSensitivityEquations_t;
	typedef typename RolloutSensitivityEquations_t::nabla_output_vector_t       nabla_output_vector_t;
	typedef typename RolloutSensitivityEquations_t::nabla_output_vector_array_t nabla_output_vector_array_t;
	typedef typename RolloutSensitivityEquations_t::nabla_output_matrix_t       nabla_output_matrix_t;
	typedef typename RolloutSensitivityEquations_t::nabla_output_matrix_array_t nabla_output_matrix_array_t;
	typedef typename RolloutSensitivityEquations_t::nabla_input_matrix_t       nabla_input_matrix_t;
	typedef typename RolloutSensitivityEquations_t::nabla_input_matrix_array_t nabla_input_matrix_array_t;
	typedef typename RolloutSensitivityEquations_t::nabla_scalar_rowvector_t       nabla_scalar_rowvector_t;
	typedef typename RolloutSensitivityEquations_t::nabla_scalar_rowvector_array_t nabla_scalar_rowvector_array_t;
	typedef Eigen::Matrix<double,DIMENSIONS::MAX_CONSTRAINT1_DIM_,NUM_SUBSYSTEMS-1>     nabla_constraint1_matrix_t;
	typedef std::vector<nabla_constraint1_matrix_t, Eigen::aligned_allocator<nabla_constraint1_matrix_t> > nabla_constraint1_matrix_array_t;

	typedef std::array<state_matrix_t, NUM_SUBSYSTEMS-1> nabla_Sm_t;
	typedef std::array<output_vector_t, NUM_SUBSYSTEMS-1> nabla_Sv_t;
	typedef std::array<eigen_scalar_t, NUM_SUBSYSTEMS-1> nabla_s_t;
	typedef std::vector<nabla_Sm_t> nabla_Sm_array_t;
	typedef std::vector<nabla_Sv_t> nabla_Sv_array_t;
    typedef std::vector<nabla_s_t>  nabla_s_array_t;


	GSLQP(const std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > >& subsystemDynamicsPtr,
			const std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > >& subsystemDerivativesPtr,
			const std::vector<std::shared_ptr<CostFunctionBase<OUTPUT_DIM, INPUT_DIM> > >& subsystemCostFunctionsPtr,
			const std::vector<controller_t>& initialControllersStock,
			const std::vector<size_t>& systemStockIndex,
			const Options_t& options = Options_t::Options())

    : subsystemDynamicsPtrStock_(NUM_SUBSYSTEMS),
      subsystemDerivativesPtrStock_(NUM_SUBSYSTEMS),
      subsystemCostFunctionsPtrStock_(NUM_SUBSYSTEMS),
      subsystemSimulatorsStockPtr_(NUM_SUBSYSTEMS),
      nominalControllersStock_(initialControllersStock),
      nominalTimeTrajectoriesStock_(NUM_SUBSYSTEMS),
      nominalStateTrajectoriesStock_(NUM_SUBSYSTEMS),
      nominalInputTrajectoriesStock_(NUM_SUBSYSTEMS),
      nominalOutputTrajectoriesStock_(NUM_SUBSYSTEMS),
      sensitivityTimeTrajectoryStock_(NUM_SUBSYSTEMS),
      nablaOutputTrajectoryStock_(NUM_SUBSYSTEMS),
      nablaInputTrajectoryStock_(NUM_SUBSYSTEMS),
      nablaqTrajectoryStock_(NUM_SUBSYSTEMS),
      nablaQvTrajectoryStock_(NUM_SUBSYSTEMS),
      nablaRvTrajectoryStock_(NUM_SUBSYSTEMS),
      nablaEvTrajectoryStock_(NUM_SUBSYSTEMS),
      AmTrajectoryStock_(NUM_SUBSYSTEMS),
      BmTrajectoryStock_(NUM_SUBSYSTEMS),
      nc1TrajectoriesStock_(NUM_SUBSYSTEMS),
      EvTrajectoryStock_(NUM_SUBSYSTEMS),
      CmTrajectoryStock_(NUM_SUBSYSTEMS),
      DmTrajectoryStock_(NUM_SUBSYSTEMS),
      qTrajectoryStock_(NUM_SUBSYSTEMS),
      QvTrajectoryStock_(NUM_SUBSYSTEMS),
      QmTrajectoryStock_(NUM_SUBSYSTEMS),
      RvTrajectoryStock_(NUM_SUBSYSTEMS),
      RmTrajectoryStock_(NUM_SUBSYSTEMS),
      PmTrajectoryStock_(NUM_SUBSYSTEMS),
      RmInverseTrajectoryStock_(NUM_SUBSYSTEMS),
      AmConstrainedTrajectoryStock_(NUM_SUBSYSTEMS),
      BmConstrainedTrajectoryStock_(NUM_SUBSYSTEMS),
      QmConstrainedTrajectoryStock_(NUM_SUBSYSTEMS),
      QvConstrainedTrajectoryStock_(NUM_SUBSYSTEMS),
      RvConstrainedTrajectoryStock_(NUM_SUBSYSTEMS),
      PmConstrainedTrajectoryStock_(NUM_SUBSYSTEMS),
      DmDagerTrajectoryStock_(NUM_SUBSYSTEMS),
      RmConstraintProjectionTrajectoryStock_(NUM_SUBSYSTEMS),
      SsTimeTrajectoryStock_(NUM_SUBSYSTEMS),
      sTrajectoryStock_(NUM_SUBSYSTEMS),
      SvTrajectoryStock_(NUM_SUBSYSTEMS),
      SveTrajectoryStock_(NUM_SUBSYSTEMS),
      SmTrajectoryStock_(NUM_SUBSYSTEMS),
      nablasTrajectoryStock_(NUM_SUBSYSTEMS),
      nablaSvTrajectoryStock_(NUM_SUBSYSTEMS),
      nablaSevTrajectoryStock_(NUM_SUBSYSTEMS),
      nablaSmTrajectoryStock_(NUM_SUBSYSTEMS),
      switchingTimes_(NUM_SUBSYSTEMS+1),
      nominalRolloutIsUpdated_(false),
      options_(options)
	{

		if (subsystemDynamicsPtr.size() != subsystemDerivativesPtr.size())
			throw std::runtime_error("Number of subsystem derivaties is not equal to the number of subsystems.");
		if (subsystemDynamicsPtr.size() != subsystemCostFunctionsPtr.size())
			throw std::runtime_error("Number of cost functions is not equal to the number of subsystems.");
		if (subsystemDynamicsPtr.size()-1 < *std::max_element(systemStockIndex.begin(), systemStockIndex.end()))
			throw std::runtime_error("systemStockIndex points to non-existing subsystem");
		if (initialControllersStock.size() != NUM_SUBSYSTEMS)
			throw std::runtime_error("initialControllersStock has less controllers then the number of subsystems");
		if (systemStockIndex.size() != NUM_SUBSYSTEMS)
			throw std::runtime_error("systemStockIndex has less elements then the number of subsystems");

		for (int i=0; i<NUM_SUBSYSTEMS; i++) {

			subsystemDynamicsPtrStock_[i] = subsystemDynamicsPtr[systemStockIndex[i]]->clone();
			subsystemDerivativesPtrStock_[i] = subsystemDerivativesPtr[systemStockIndex[i]]->clone();
			subsystemCostFunctionsPtrStock_[i] = subsystemCostFunctionsPtr[systemStockIndex[i]]->clone();

			subsystemSimulatorsStockPtr_[i] = std::make_shared<ODE45<STATE_DIM> >(subsystemDynamicsPtrStock_[i]);
		}
	}

	~GSLQP() {}

	void rollout(const state_vector_t& initState,
			const std::vector<controller_t>& controllersStock,
			std::vector<scalar_array_t>& timeTrajectoriesStock,
			std::vector<state_vector_array_t>& stateTrajectoriesStock,
			std::vector<control_vector_array_t>& inputTrajectoriesStock,
			std::vector<output_vector_array_t>& outputTrajectoriesStock,
			std::vector<std::vector<size_t> >& nc1TrajectoriesStock,
			std::vector<constraint1_vector_array_t>& EvTrajectoryStock);

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
			const std::vector<output_vector_array_t>& stateTrajectoriesStock,
			const std::vector<control_vector_array_t>& inputTrajectoriesStock,
			scalar_t& totalCost);

	void getRolloutSensitivity2SwitchingTime(std::vector<scalar_array_t>& sensitivityTimeTrajectoriesStock,
		std::vector<nabla_output_matrix_array_t>& sensitivityOutputTrajectoriesStock,
		std::vector<nabla_input_matrix_array_t>& sensitivityInputTrajectoriesStock);

	void getController(std::vector<controller_t>& controllersStock);

	void getValueFuntion(const scalar_t& time, const output_vector_t& output, scalar_t& valueFuntion);

	void getCostFuntionDerivative(const output_vector_t& initOutput, Eigen::Matrix<double,NUM_SUBSYSTEMS-1,1>& costFuntionDerivative);

	void getNominalTrajectories(std::vector<scalar_array_t>& nominalTimeTrajectoriesStock,
			std::vector<state_vector_array_t>& nominalStateTrajectoriesStock,
			std::vector<control_vector_array_t>& nominalInputTrajectoriesStock,
			std::vector<output_vector_array_t>& nominalOutputTrajectoriesStock = std::vector<output_vector_array_t>());

	void run(const state_vector_t& initState, const std::vector<scalar_t>& switchingTimes);


protected:
	void solveSequentialRiccatiEquations(const scalar_t& learningRate);

	void solveFullSequentialRiccatiEquations(const scalar_t& learningRate);

	void approximateOptimalControlProblem();

	void calculatecontroller(std::vector<controller_t>& controllersStock,
			std::vector<control_vector_array_t>& feedForwardControlStock,
			std::vector<control_vector_array_t>& feedForwardConstraintInputStock);

	void lineSearch(const std::vector<control_vector_array_t>& feedForwardControlStock,
			const std::vector<control_vector_array_t>& feedForwardConstraintInputStock,
			scalar_t& learningRateStar);

	void transformeLocalValueFuntion2Global();

	void transformeLocalValueFuntionDerivative2Global();

	void rolloutSensitivity2SwitchingTime(bool nablaSvUpdated);

	void inputIncrementSensitivity2SwitchingTime(std::vector<scalar_array_t>& nablaUffTimeTrajectoryStock,
			std::vector<nabla_input_matrix_array_t>& nablaUffTrajectoryStock);

	template <typename Derived>
	bool makePSD(Eigen::MatrixBase<Derived>& squareMatrix);


private:
	std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > > subsystemDynamicsPtrStock_;
	std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > > subsystemDerivativesPtrStock_;
	std::vector<std::shared_ptr<CostFunctionBase<OUTPUT_DIM, INPUT_DIM> > > subsystemCostFunctionsPtrStock_;

	std::vector<std::shared_ptr<ODE45<STATE_DIM> > > subsystemSimulatorsStockPtr_;

	scalar_t nominalTotalCost_;
	std::vector<controller_t> nominalControllersStock_;
	std::vector<scalar_array_t> nominalTimeTrajectoriesStock_;
	std::vector<state_vector_array_t>   nominalStateTrajectoriesStock_;
	std::vector<control_vector_array_t> nominalInputTrajectoriesStock_;
	std::vector<output_vector_array_t>  nominalOutputTrajectoriesStock_;

	std::vector<scalar_array_t> sensitivityTimeTrajectoryStock_;
	std::vector<nabla_output_matrix_array_t> nablaOutputTrajectoryStock_;
	std::vector<nabla_input_matrix_array_t> nablaInputTrajectoryStock_;
	//
	std::vector<nabla_scalar_rowvector_array_t> nablaqTrajectoryStock_;
	std::vector<nabla_output_matrix_array_t> nablaQvTrajectoryStock_;
	std::vector<nabla_input_matrix_array_t> nablaRvTrajectoryStock_;
	std::vector<nabla_constraint1_matrix_array_t> nablaEvTrajectoryStock_;
	nabla_scalar_rowvector_t nablaqFinal_;
	nabla_output_matrix_t nablaQvFinal_;

	std::vector<state_matrix_array_t>        AmTrajectoryStock_;
	std::vector<control_gain_matrix_array_t> BmTrajectoryStock_;

	std::vector<std::vector<size_t> >       nc1TrajectoriesStock_;  // nc1: Number of the active constraints
	std::vector<constraint1_vector_array_t> EvTrajectoryStock_;
	std::vector<constraint1_state_matrix_array_t>   CmTrajectoryStock_;
	std::vector<constraint1_control_matrix_array_t> DmTrajectoryStock_;

	eigen_scalar_t  qFinal_;
	output_vector_t QvFinal_;
	state_matrix_t  QmFinal_;
	std::vector<eigen_scalar_array_t> qTrajectoryStock_;
	std::vector<output_vector_array_t> QvTrajectoryStock_;
	std::vector<state_matrix_array_t> QmTrajectoryStock_;
	std::vector<control_vector_array_t> RvTrajectoryStock_;
	std::vector<control_matrix_array_t> RmTrajectoryStock_;
	std::vector<control_feedback_array_t> PmTrajectoryStock_;

	std::vector<control_matrix_array_t>      RmInverseTrajectoryStock_;
	std::vector<state_matrix_array_t>        AmConstrainedTrajectoryStock_;
	std::vector<control_gain_matrix_array_t> BmConstrainedTrajectoryStock_;
	std::vector<state_matrix_array_t>        QmConstrainedTrajectoryStock_;
	std::vector<output_vector_array_t>       QvConstrainedTrajectoryStock_;
	std::vector<control_vector_array_t>      RvConstrainedTrajectoryStock_;
	std::vector<control_feedback_array_t>    PmConstrainedTrajectoryStock_;
	std::vector<control_constraint1_matrix_array_t> DmDagerTrajectoryStock_;
	std::vector<constraint1_matrix_array_t>  RmConstraintProjectionTrajectoryStock_;


	std::vector<scalar_array_t> 	  SsTimeTrajectoryStock_;
	std::vector<eigen_scalar_array_t> sTrajectoryStock_;
	std::vector<output_vector_array_t> SvTrajectoryStock_;
	std::vector<output_vector_array_t> SveTrajectoryStock_;
	std::vector<state_matrix_array_t> SmTrajectoryStock_;
	std::vector<nabla_s_array_t>  nablasTrajectoryStock_;
	std::vector<nabla_Sv_array_t> nablaSvTrajectoryStock_;
	std::vector<nabla_Sv_array_t> nablaSevTrajectoryStock_;
	std::vector<nabla_Sm_array_t> nablaSmTrajectoryStock_;

	scalar_array_t switchingTimes_;
	state_vector_t initState_;

	bool nominalRolloutIsUpdated_;

	Options_t options_;
};

#include "implementation/GSLQP.h"

#endif /* GSLQP_H_ */
