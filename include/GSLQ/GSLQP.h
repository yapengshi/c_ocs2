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
#include "GSLQ/FullSequentialRiccatiEquations.h"
#include "GSLQ/RolloutSensitivityEquations.h"

template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
class GSLQP
{
public:
	typedef SequentialRiccatiEquations<STATE_DIM, INPUT_DIM, NUM_Subsystems> RiccatiEquations;
	typedef FullSequentialRiccatiEquations<STATE_DIM, INPUT_DIM, NUM_Subsystems> FullRiccatiEquations;

	typedef Dimensions<STATE_DIM, INPUT_DIM> DIMENSIONS;
	typedef typename DIMENSIONS::controller_t controller_t;
	typedef typename DIMENSIONS::Options Options_t;
	typedef typename DIMENSIONS::scalar_t 		scalar_t;
	typedef typename DIMENSIONS::scalar_array_t scalar_array_t;
	typedef typename DIMENSIONS::eigen_scalar_t       eigen_scalar_t;
	typedef typename DIMENSIONS::eigen_scalar_array_t eigen_scalar_array_t;
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

	typedef RolloutSensitivityEquations<STATE_DIM, INPUT_DIM, NUM_Subsystems> RolloutSensitivityEquations_t;
	typedef typename RolloutSensitivityEquations_t::nabla_state_matrix_t       nabla_state_matrix_t;
	typedef typename RolloutSensitivityEquations_t::nabla_state_matrix_array_t nabla_state_matrix_array_t;
	typedef typename RolloutSensitivityEquations_t::nabla_input_matrix_t       nabla_input_matrix_t;
	typedef typename RolloutSensitivityEquations_t::nabla_input_matrix_array_t nabla_input_matrix_array_t;
	typedef typename RolloutSensitivityEquations_t::nabla_scalar_rowvector_t       nabla_scalar_rowvector_t;
	typedef typename RolloutSensitivityEquations_t::nabla_scalar_rowvector_array_t nabla_scalar_rowvector_array_t;

	typedef std::array<state_matrix_t, NUM_Subsystems-1> nabla_Sm_t;
	typedef std::array<state_vector_t, NUM_Subsystems-1> nabla_Sv_t;
	typedef std::array<eigen_scalar_t, NUM_Subsystems-1> nabla_s_t;
	typedef std::vector<nabla_Sm_t> nabla_Sm_array_t;
	typedef std::vector<nabla_Sv_t> nabla_Sv_array_t;
    typedef std::vector<nabla_s_t>  nabla_s_array_t;


	GSLQP(const std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM> > >& subsystemDynamicsPtr,
			const std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM> > >& subsystemDerivativesPtr,
			const std::vector<std::shared_ptr<CostFunctionBase<STATE_DIM, INPUT_DIM> > >& subsystemCostFunctionsPtr,
			const std::vector<controller_t>& initialControllersStock,
			const std::vector<size_t>& systemStockIndex,
			const Options_t& options = Options_t::Options())

    : subsystemDynamicsPtrStock_(NUM_Subsystems),
      subsystemDerivativesPtrStock_(NUM_Subsystems),
      subsystemCostFunctionsPtrStock_(NUM_Subsystems),
      subsystemSimulatorsStockPtr_(NUM_Subsystems),
      nominalControllersStock_(initialControllersStock),
      nominalTimeTrajectoriesStock_(NUM_Subsystems),
      nominalStateTrajectoriesStock_(NUM_Subsystems),
      nominalInputTrajectoriesStock_(NUM_Subsystems),
      sensitivityTimeTrajectoryStock_(NUM_Subsystems),
      nablaStateTrajectoryStock_(NUM_Subsystems),
      nablaInputTrajectoryStock_(NUM_Subsystems),
      nablaqTrajectoryStock_(NUM_Subsystems),
      nablaQvTrajectoryStock_(NUM_Subsystems),
      nablaRvTrajectoryStock_(NUM_Subsystems),
      AmTrajectoryStock_(NUM_Subsystems),
      BmTrajectoryStock_(NUM_Subsystems),
      qTrajectoryStock_(NUM_Subsystems),
      QvTrajectoryStock_(NUM_Subsystems),
      QmTrajectoryStock_(NUM_Subsystems),
      RvTrajectoryStock_(NUM_Subsystems),
      RmTrajectoryStock_(NUM_Subsystems),
      PmTrajectoryStock_(NUM_Subsystems),
      SsTimeTrajectoryStock_(NUM_Subsystems),
      sTrajectoryStock_(NUM_Subsystems),
      SvTrajectoryStock_(NUM_Subsystems),
      SmTrajectoryStock_(NUM_Subsystems),
      nablasTrajectoryStock_(NUM_Subsystems),
      nablaSvTrajectoryStock_(NUM_Subsystems),
      nablaSmTrajectoryStock_(NUM_Subsystems),
      switchingTimes_(NUM_Subsystems+1),
      nominalRolloutIsUpdated_(false),
      options_(options)
	{

		if (subsystemDynamicsPtr.size() != subsystemDerivativesPtr.size())
			throw std::runtime_error("Number of subsystem derivaties is not equal to the number of subsystems.");
		if (subsystemDynamicsPtr.size() != subsystemCostFunctionsPtr.size())
			throw std::runtime_error("Number of cost functions is not equal to the number of subsystems.");
		if (subsystemDynamicsPtr.size()-1 < *std::max_element(systemStockIndex.begin(), systemStockIndex.end()))
			throw std::runtime_error("systemStockIndex points to non-existing subsystem");
		if (initialControllersStock.size() != NUM_Subsystems)
			throw std::runtime_error("initialControllersStock has less controllers then the number of subsystems");
		if (systemStockIndex.size() != NUM_Subsystems)
			throw std::runtime_error("systemStockIndex has less elements then the number of subsystems");

		for (int i=0; i<NUM_Subsystems; i++) {

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
			std::vector<control_vector_array_t>& inputTrajectoriesStock);

	void rolloutCost(const std::vector<scalar_array_t>& timeTrajectoriesStock,
			const std::vector<state_vector_array_t>& stateTrajectoriesStock,
			const std::vector<control_vector_array_t>& inputTrajectoriesStock,
			scalar_t& totalCost);

	void getRolloutSensitivity2SwitchingTime(std::vector<scalar_array_t>& sensitivityTimeTrajectoriesStock,
		std::vector<nabla_state_matrix_array_t>& sensitivityStateTrajectoriesStock,
		std::vector<nabla_input_matrix_array_t>& sensitivityInputTrajectoriesStock);

	void getController(std::vector<controller_t>& controllersStock) { controllersStock = nominalControllersStock_;}

	void getValueFuntion(const scalar_t& time, const state_vector_t& state, scalar_t& valueFuntion);

	void getCostFuntionDerivative(const state_vector_t& initState, Eigen::Matrix<double,NUM_Subsystems-1,1>& costFuntionDerivative);

	void getNominalTrajectories(std::vector<scalar_array_t>& nominalTimeTrajectoriesStock,
			std::vector<state_vector_array_t>& nominalStateTrajectoriesStock,
			std::vector<control_vector_array_t>& nominalInputTrajectoriesStock);

	void run(const state_vector_t& initState, const std::vector<scalar_t>& switchingTimes);


protected:
	void SolveSequentialRiccatiEquations(const scalar_t& learningRate);

	void SolveFullSequentialRiccatiEquations(const scalar_t& learningRate);

	void approximateOptimalControlProblem();

	void calculatecontroller(scalar_t& learningRate);

	void lineSearch(const std::vector<controller_t>& controllersStock, const std::vector<control_vector_array_t>& deltaUffStock,
			scalar_t& learningRateStar);

	void transformeLocalValueFuntion2Global();

	void transformeLocalValueFuntionDerivative2Global();

	void rolloutSensitivity2SwitchingTime();

	template <typename Derived>
	bool makePSD(Eigen::MatrixBase<Derived>& squareMatrix);


private:
	std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM> > > subsystemDynamicsPtrStock_;
	std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM> > > subsystemDerivativesPtrStock_;
	std::vector<std::shared_ptr<CostFunctionBase<STATE_DIM, INPUT_DIM> > > subsystemCostFunctionsPtrStock_;

	std::vector<std::shared_ptr<ODE45<STATE_DIM> > > subsystemSimulatorsStockPtr_;

	scalar_t nominalTotalCost_;
	std::vector<controller_t> nominalControllersStock_;
	std::vector<scalar_array_t> nominalTimeTrajectoriesStock_;
	std::vector<state_vector_array_t> nominalStateTrajectoriesStock_;
	std::vector<control_vector_array_t> nominalInputTrajectoriesStock_;

	std::vector<scalar_array_t> sensitivityTimeTrajectoryStock_;
	std::vector<nabla_state_matrix_array_t> nablaStateTrajectoryStock_;
	std::vector<nabla_input_matrix_array_t> nablaInputTrajectoryStock_;
	//
	std::vector<nabla_scalar_rowvector_array_t> nablaqTrajectoryStock_;
	std::vector<nabla_state_matrix_array_t> nablaQvTrajectoryStock_;
	std::vector<nabla_input_matrix_array_t> nablaRvTrajectoryStock_;

	std::vector<state_matrix_array_t>        AmTrajectoryStock_;
	std::vector<control_gain_matrix_array_t> BmTrajectoryStock_;

	eigen_scalar_t qFinal_;
	state_vector_t QvFinal_;
	state_matrix_t QmFinal_;
	std::vector<eigen_scalar_array_t> qTrajectoryStock_;
	std::vector<state_vector_array_t> QvTrajectoryStock_;
	std::vector<state_matrix_array_t> QmTrajectoryStock_;
	std::vector<control_vector_array_t> RvTrajectoryStock_;
	std::vector<control_matrix_array_t> RmTrajectoryStock_;
	std::vector<control_feedback_array_t> PmTrajectoryStock_;

	std::vector<scalar_array_t> 	  SsTimeTrajectoryStock_;
	std::vector<eigen_scalar_array_t> sTrajectoryStock_;
	std::vector<state_vector_array_t> SvTrajectoryStock_;
	std::vector<state_matrix_array_t> SmTrajectoryStock_;
	std::vector<nabla_s_array_t>  nablasTrajectoryStock_;
	std::vector<nabla_Sv_array_t> nablaSvTrajectoryStock_;
	std::vector<nabla_Sm_array_t> nablaSmTrajectoryStock_;

	scalar_array_t switchingTimes_;
	state_vector_t initState_;

	bool nominalRolloutIsUpdated_;

	Options_t options_;
};

#include "implementation/GSLQP.h"

#endif /* GSLQP_H_ */
