/*
 * GSLQP.h
 *
 *  Created on: Dec 18, 2015
 *      Author: farbod
 */

#ifndef GSLQP_OCS2_H_
#define GSLQP_OCS2_H_

#include <vector>
#include <array>
#include <algorithm>
#include <cstddef>
#include <Eigen/Dense>
#include <Eigen/StdVector>

#include "Dimensions.h"

#include "dynamics/ControlledSystemBase.h"
#include "dynamics/DerivativesBase.h"
#include "costs/CostFunctionBaseOCS2.h"

#include "integration/Integrator.h"
#include "misc/LinearInterpolation.h"

#include "GSLQ/SequentialRiccatiEquations.h"
#include "GSLQ/SequentialErrorEquation.h"
#include "GSLQ/FullSequentialRiccatiEquations.h"
#include "GSLQ/RolloutSensitivityEquations.h"
#include "GSLQ/SLQP.h"
#include "GSLQ/SolveBVP.h"


namespace ocs2{

template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
class GSLQP
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	typedef SequentialErrorEquation<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS> ErrorEquation_t;
	typedef FullSequentialRiccatiEquations<OUTPUT_DIM, INPUT_DIM, NUM_SUBSYSTEMS> FullRiccatiEquations_t;

	typedef Dimensions<STATE_DIM, INPUT_DIM, OUTPUT_DIM> DIMENSIONS;
	typedef typename DIMENSIONS::template LinearFunction_t<Eigen::Dynamic> lagrange_t;
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
			const std::vector<std::shared_ptr<CostFunctionBaseOCS2<OUTPUT_DIM, INPUT_DIM> > >& subsystemCostFunctionsPtr,
			const std::vector<controller_t>& initialControllersStock,
			const std::vector<size_t>& systemStockIndex,
			const Options_t& options = Options_t::Options())

    : slqp_(subsystemDynamicsPtr, subsystemDerivativesPtr, subsystemCostFunctionsPtr, initialControllersStock, systemStockIndex, options),
      sensitivityTimeTrajectoryStock_(NUM_SUBSYSTEMS),
      nablaOutputTrajectoryStock_(NUM_SUBSYSTEMS),
      nablaInputTrajectoryStock_(NUM_SUBSYSTEMS),
      nablaqTrajectoryStock_(NUM_SUBSYSTEMS),
      nablaQvTrajectoryStock_(NUM_SUBSYSTEMS),
      nablaRvTrajectoryStock_(NUM_SUBSYSTEMS),
      nablaEvTrajectoryStock_(NUM_SUBSYSTEMS),
      nablasTrajectoryStock_(NUM_SUBSYSTEMS),
      nablaSvTrajectoryStock_(NUM_SUBSYSTEMS),
      nablaSevTrajectoryStock_(NUM_SUBSYSTEMS),
      nablaSmTrajectoryStock_(NUM_SUBSYSTEMS),
      switchingTimes_(NUM_SUBSYSTEMS+1),
      options_(options)
	{}

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

	void calculateCostFunction(const std::vector<scalar_array_t>& timeTrajectoriesStock,
			const std::vector<output_vector_array_t>& stateTrajectoriesStock,
			const std::vector<control_vector_array_t>& inputTrajectoriesStock,
			scalar_t& totalCost);

	void calculateMeritFunction(const std::vector<scalar_array_t>& timeTrajectoriesStock,
			const std::vector<std::vector<size_t> >& nc1TrajectoriesStock,
			const std::vector<constraint1_vector_array_t>& EvTrajectoryStock,
			const std::vector<std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd> > >&  lagrangeTrajectoriesStock,
			const scalar_t& totalCost,
			scalar_t& meritFuntionValue,
			scalar_t& constraintISE);

	double calculateConstraintISE(const std::vector<scalar_array_t>& timeTrajectoriesStock,
			const std::vector<std::vector<size_t>>& nc1TrajectoriesStock,
			const std::vector<constraint1_vector_array_t>& EvTrajectoriesStock,
			scalar_t& constraintISE);

	void getRolloutSensitivity2SwitchingTime(std::vector<scalar_array_t>& sensitivityTimeTrajectoriesStock,
		std::vector<nabla_output_matrix_array_t>& sensitivityOutputTrajectoriesStock,
		std::vector<nabla_input_matrix_array_t>& sensitivityInputTrajectoriesStock);

	void getController(std::vector<controller_t>& controllersStock);

	void getValueFuntion(const scalar_t& time, const output_vector_t& output, scalar_t& valueFuntion);

	void getCostFuntionDerivative(const output_vector_t& initOutput, Eigen::Matrix<double,NUM_SUBSYSTEMS-1,1>& costFuntionDerivative);

	void getNominalTrajectories(std::vector<scalar_array_t>& nominalTimeTrajectoriesStock,
			std::vector<state_vector_array_t>& nominalStateTrajectoriesStock,
			std::vector<control_vector_array_t>& nominalInputTrajectoriesStock,
			std::vector<output_vector_array_t>& nominalOutputTrajectoriesStock);

	void run(const state_vector_t& initState, const std::vector<scalar_t>& switchingTimes);


protected:
	void solveFullSequentialRiccatiEquations(const scalar_t& learningRate);

	void transformLocalValueFuntionDerivative2Global();

	void rolloutSensitivity2SwitchingTime(bool nablaSvUpdated);

	void inputIncrementSensitivity2SwitchingTime(std::vector<scalar_array_t>& nablaUffTimeTrajectoryStock,
			std::vector<nabla_input_matrix_array_t>& nablaUffTrajectoryStock);

private:
	SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS> slqp_;


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

	std::vector<nabla_s_array_t>  nablasTrajectoryStock_;
	std::vector<nabla_Sv_array_t> nablaSvTrajectoryStock_;
	std::vector<nabla_Sv_array_t> nablaSevTrajectoryStock_;
	std::vector<nabla_Sm_array_t> nablaSmTrajectoryStock_;

	scalar_array_t switchingTimes_;
	state_vector_t initState_;
	Options_t options_;
};

} // namespace ocs2

#include "implementation/GSLQP.h"

#endif /* GSLQP_H_ */
