/*
 * SLQP.h
 *
 *  Created on: Jun 16, 2016
 *      Author: farbod
 */

#ifndef SLQP_H_
#define SLQP_H_

#include <vector>
#include <array>
#include <algorithm>
#include <cstddef>
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
#include "GSLQ/GSLQP.h"


template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
class SLQP
{
public:
	typedef SequentialRiccatiEquations<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS> RiccatiEquations_t;
	typedef SequentialErrorEquation<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS> ErrorEquation_t;

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


	SLQP(const std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > >& subsystemDynamicsPtr,
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
      nominalcostateTrajectoriesStock_(NUM_SUBSYSTEMS),
      nominalLagrangeTrajectoriesStock_(NUM_SUBSYSTEMS),
      lagrangeControllerStock_(NUM_SUBSYSTEMS),
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
      QmConstrainedTrajectoryStock_(NUM_SUBSYSTEMS),
      QvConstrainedTrajectoryStock_(NUM_SUBSYSTEMS),
      RmConstrainedTrajectoryStock_(NUM_SUBSYSTEMS),
      DmDagerTrajectoryStock_(NUM_SUBSYSTEMS),
      EvProjectedTrajectoryStock_(NUM_SUBSYSTEMS),
      CmProjectedTrajectoryStock_(NUM_SUBSYSTEMS),
      DmProjectedTrajectoryStock_(NUM_SUBSYSTEMS),
      SsTimeTrajectoryStock_(NUM_SUBSYSTEMS),
      sTrajectoryStock_(NUM_SUBSYSTEMS),
      SvTrajectoryStock_(NUM_SUBSYSTEMS),
      SveTrajectoryStock_(NUM_SUBSYSTEMS),
      SmTrajectoryStock_(NUM_SUBSYSTEMS),
      switchingTimes_(NUM_SUBSYSTEMS+1),
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

	~SLQP() {}

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

	void getController(std::vector<controller_t>& controllersStock);

	void getValueFuntion(const scalar_t& time, const output_vector_t& output, scalar_t& valueFuntion);

	void getCostFuntion(const output_vector_t& initOutput, scalar_t& costFunction);

	void getNominalTrajectories(std::vector<scalar_array_t>& nominalTimeTrajectoriesStock,
			std::vector<state_vector_array_t>& nominalStateTrajectoriesStock,
			std::vector<control_vector_array_t>& nominalInputTrajectoriesStock,
			std::vector<output_vector_array_t>& nominalOutputTrajectoriesStock);

	void run(const state_vector_t& initState, const std::vector<scalar_t>& switchingTimes);


protected:
	void solveSequentialRiccatiEquations(const scalar_t& learningRate);

	void approximateOptimalControlProblem();

	void calculateControllerAndLagrangian(std::vector<controller_t>& controllersStock,
			std::vector<lagrange_t>& lagrangeMultiplierFunctionsStock,
			std::vector<control_vector_array_t>& feedForwardConstraintInputStock,
			bool firstCall=true);

	void calculateRolloutLagrangeMultiplier(const std::vector<scalar_array_t>& timeTrajectoriesStock,
			const std::vector<output_vector_array_t>& outputTrajectoriesStock,
			const std::vector<lagrange_t>& lagrangeMultiplierFunctionsStock,
			std::vector<std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd> > >&  lagrangeTrajectoriesStock);

	void calculateRolloutCostate(const std::vector<scalar_array_t>& timeTrajectoriesStock,
			const std::vector<output_vector_array_t>& outputTrajectoriesStock,
			std::vector<output_vector_array_t>& costateTrajectoriesStock);

	void lineSearch(const std::vector<control_vector_array_t>& feedForwardConstraintInputStock,
			scalar_t& learningRateStar,
			scalar_t maxLearningRateStar=1.0);

	void transformLocalValueFuntion2Global();

	template <typename Derived>
	bool makePSD(Eigen::MatrixBase<Derived>& squareMatrix);


private:
	std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > > subsystemDynamicsPtrStock_;
	std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > > subsystemDerivativesPtrStock_;
	std::vector<std::shared_ptr<CostFunctionBase<OUTPUT_DIM, INPUT_DIM> > > subsystemCostFunctionsPtrStock_;

	std::vector<std::shared_ptr<ODE45<STATE_DIM> > > subsystemSimulatorsStockPtr_;

	scalar_t nominalTotalCost_;
	scalar_t nominalTotalMerit_;
	scalar_t nominalConstraint1ISE_;
	std::vector<controller_t> nominalControllersStock_;
	std::vector<scalar_array_t> nominalTimeTrajectoriesStock_;
	std::vector<state_vector_array_t>   nominalStateTrajectoriesStock_;
	std::vector<control_vector_array_t> nominalInputTrajectoriesStock_;
	std::vector<output_vector_array_t>  nominalOutputTrajectoriesStock_;
	std::vector<output_vector_array_t>  nominalcostateTrajectoriesStock_;
	std::vector<std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd> > >  nominalLagrangeTrajectoriesStock_;

	std::vector<lagrange_t> lagrangeControllerStock_;

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

	std::vector<control_matrix_array_t> RmInverseTrajectoryStock_;
	std::vector<state_matrix_array_t>   AmConstrainedTrajectoryStock_;
	std::vector<state_matrix_array_t>   QmConstrainedTrajectoryStock_;
	std::vector<output_vector_array_t>  QvConstrainedTrajectoryStock_;
	std::vector<control_matrix_array_t> RmConstrainedTrajectoryStock_;
	std::vector<control_constraint1_matrix_array_t> DmDagerTrajectoryStock_;
	std::vector<control_vector_array_t>   EvProjectedTrajectoryStock_;  // DmDager * Ev
	std::vector<control_feedback_array_t> CmProjectedTrajectoryStock_;  // DmDager * Cm
	std::vector<control_matrix_array_t>   DmProjectedTrajectoryStock_;  // DmDager * Dm


	std::vector<scalar_array_t> 	  SsTimeTrajectoryStock_;
	std::vector<eigen_scalar_array_t> sTrajectoryStock_;
	std::vector<output_vector_array_t> SvTrajectoryStock_;
	std::vector<output_vector_array_t> SveTrajectoryStock_;
	std::vector<state_matrix_array_t> SmTrajectoryStock_;

	scalar_array_t switchingTimes_;
	state_vector_t initState_;
	size_t iteration_;
	Options_t options_;

public:
	template <size_t GSLQP_STATE_DIM, size_t GSLQP_INPUT_DIM, size_t GSLQP_OUTPUT_DIM, size_t GSLQP_NUM_SUBSYSTEMS>
	friend class GSLQP;
};

#include "implementation/SLQP.h"


#endif /* SLQP_H_ */
