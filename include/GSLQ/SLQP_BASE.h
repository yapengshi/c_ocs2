/*
 * SLQP_BASE.h
 *
 *  Created on: August 14, 2016
 *      Author: mgiftthaler
 */

#ifndef SLQP_BASE_OCS2_H_
#define SLQP_BASE_OCS2_H_

#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
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

#include <chrono>

namespace ocs2{

template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
class SLQP_BASE
{
public:
	typedef std::shared_ptr<SLQP_BASE<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS> > Ptr;
	typedef SequentialRiccatiEquations<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS> RiccatiEquations_t;
	typedef SequentialErrorEquation<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS> ErrorEquation_t;

	typedef Dimensions<STATE_DIM, INPUT_DIM, OUTPUT_DIM> DIMENSIONS;

	typedef typename DIMENSIONS::template LinearFunction_t<Eigen::Dynamic> lagrange_t;

	typedef typename DIMENSIONS::controller_t controller_t;
	typedef typename DIMENSIONS::Options Options_t;
	typedef typename DIMENSIONS::scalar_t scalar_t;
	typedef typename DIMENSIONS::scalar_array_t scalar_array_t;
	typedef typename DIMENSIONS::eigen_scalar_t eigen_scalar_t;
	typedef typename DIMENSIONS::eigen_scalar_array_t eigen_scalar_array_t;
	typedef typename DIMENSIONS::state_vector_t state_vector_t;
	typedef typename DIMENSIONS::state_vector_array_t state_vector_array_t;
	typedef typename DIMENSIONS::output_vector_t output_vector_t;
	typedef typename DIMENSIONS::output_vector_array_t output_vector_array_t;
	typedef typename DIMENSIONS::control_vector_t control_vector_t;
	typedef typename DIMENSIONS::control_vector_array_t control_vector_array_t;
	typedef typename DIMENSIONS::control_feedback_t control_feedback_t;
	typedef typename DIMENSIONS::control_feedback_array_t control_feedback_array_t;
	typedef typename DIMENSIONS::state_matrix_t state_matrix_t;
	typedef typename DIMENSIONS::state_matrix_array_t state_matrix_array_t;
	typedef typename DIMENSIONS::control_matrix_t control_matrix_t;
	typedef typename DIMENSIONS::control_matrix_array_t control_matrix_array_t;
	typedef typename DIMENSIONS::control_gain_matrix_t control_gain_matrix_t;
	typedef typename DIMENSIONS::control_gain_matrix_array_t control_gain_matrix_array_t;
	typedef typename DIMENSIONS::constraint1_vector_t constraint1_vector_t;
	typedef typename DIMENSIONS::constraint1_vector_array_t constraint1_vector_array_t;
	typedef typename DIMENSIONS::constraint1_state_matrix_t constraint1_state_matrix_t;
	typedef typename DIMENSIONS::constraint1_state_matrix_array_t constraint1_state_matrix_array_t;
	typedef typename DIMENSIONS::constraint1_control_matrix_t constraint1_control_matrix_t;
	typedef typename DIMENSIONS::constraint1_control_matrix_array_t constraint1_control_matrix_array_t;
	typedef typename DIMENSIONS::control_constraint1_matrix_t control_constraint1_matrix_t;
	typedef typename DIMENSIONS::control_constraint1_matrix_array_t control_constraint1_matrix_array_t;
	typedef typename DIMENSIONS::constraint2_vector_t       constraint2_vector_t;
	typedef typename DIMENSIONS::constraint2_vector_array_t constraint2_vector_array_t;
	typedef typename DIMENSIONS::constraint2_state_matrix_t       constraint2_state_matrix_t;
	typedef typename DIMENSIONS::constraint2_state_matrix_array_t constraint2_state_matrix_array_t;



	SLQP_BASE(const Options_t& options = Options_t::Options())
    : nominalTimeTrajectoriesStock_(NUM_SUBSYSTEMS),
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
      nc2TrajectoriesStock_(NUM_SUBSYSTEMS),
      HvTrajectoryStock_(NUM_SUBSYSTEMS),
      FmTrajectoryStock_(NUM_SUBSYSTEMS),
      nc2FinalStock_(NUM_SUBSYSTEMS),
      HvFinalStock_(NUM_SUBSYSTEMS),
      FmFinalStock_(NUM_SUBSYSTEMS),
      qFinalStock_(NUM_SUBSYSTEMS),
      QvFinalStock_(NUM_SUBSYSTEMS),
      QmFinalStock_(NUM_SUBSYSTEMS),
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
      iteration_(0),
      options_(options)
	{}


	virtual ~SLQP_BASE() {}

	virtual void rollout(const state_vector_t& initState,
			const std::vector<controller_t>& controllersStock,
			std::vector<scalar_array_t>& timeTrajectoriesStock,
			std::vector<state_vector_array_t>& stateTrajectoriesStock,
			std::vector<control_vector_array_t>& inputTrajectoriesStock,
			std::vector<output_vector_array_t>& outputTrajectoriesStock,
			std::vector<std::vector<size_t> >& nc1TrajectoriesStock,
			std::vector<constraint1_vector_array_t>& EvTrajectoryStock,
			std::vector<std::vector<size_t> >& nc2TrajectoriesStock,
			std::vector<constraint2_vector_array_t>& HvTrajectoryStock,
			std::vector<size_t>& nc2FinalStock,
			std::vector<constraint2_vector_t>& HvFinalStock) = 0;

	virtual void rollout(const state_vector_t& initState,
			const std::vector<controller_t>& controllersStock,
			std::vector<scalar_array_t>& timeTrajectoriesStock,
			std::vector<state_vector_array_t>& stateTrajectoriesStock,
			std::vector<control_vector_array_t>& inputTrajectoriesStock,
			std::vector<output_vector_array_t>& outputTrajectoriesStock) = 0;

	virtual void rollout(const state_vector_t& initState,
			const std::vector<controller_t>& controllersStock,
			std::vector<scalar_array_t>& timeTrajectoriesStock,
			std::vector<state_vector_array_t>& stateTrajectoriesStock,
			std::vector<control_vector_array_t>& inputTrajectoriesStock) = 0;

	virtual void rollout(const state_vector_t& initState,
			const std::vector<controller_t>& controllersStock,
			const double& stoppingTime,
			state_vector_t& stateVectorWhereStopped,
			control_vector_t& controlInputWhereStopped,
			output_vector_t& outputWhereStopped,
			size_t& numSubsystemWhereStopped) = 0;

	virtual void calculateCostFunction(const std::vector<scalar_array_t>& timeTrajectoriesStock,
			const std::vector<output_vector_array_t>& stateTrajectoriesStock,
			const std::vector<control_vector_array_t>& inputTrajectoriesStock,
			scalar_t& totalCost) = 0;

	virtual void calculateCostFunction(const std::vector<scalar_array_t>& timeTrajectoriesStock,
			const std::vector<output_vector_array_t>& stateTrajectoriesStock,
			const std::vector<control_vector_array_t>& inputTrajectoriesStock,
			const std::vector<std::vector<size_t> >& nc2TrajectoriesStock,
			const std::vector<constraint2_vector_array_t>& HvTrajectoryStock,
			const std::vector<size_t>& nc2FinalStock,
			const std::vector<constraint2_vector_t>& HvFinalStock,
			scalar_t& totalCost) = 0;

	virtual void calculateMeritFunction(const std::vector<scalar_array_t>& timeTrajectoriesStock,
			const std::vector<std::vector<size_t> >& nc1TrajectoriesStock,
			const std::vector<constraint1_vector_array_t>& EvTrajectoryStock,
			const std::vector<std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd> > >&  lagrangeTrajectoriesStock,
			const scalar_t& totalCost,
			scalar_t& meritFuntionValue,
			scalar_t& constraintISE) = 0;

	virtual double calculateConstraintISE(const std::vector<scalar_array_t>& timeTrajectoriesStock,
			const std::vector<std::vector<size_t>>& nc1TrajectoriesStock,
			const std::vector<constraint1_vector_array_t>& EvTrajectoriesStock,
			scalar_t& constraintISE) = 0;

	virtual void getController(std::vector<controller_t>& controllersStock) = 0;

	virtual void setController(const std::vector<controller_t>& controllersStock) = 0;

	virtual void getValueFuntion(const scalar_t& time, const output_vector_t& output, scalar_t& valueFuntion) = 0;

	virtual void getCostFuntion(scalar_t& costFunction, scalar_t& constraintISE) = 0;

	virtual void getNominalTrajectories(std::vector<scalar_array_t>& nominalTimeTrajectoriesStock,
			std::vector<state_vector_array_t>& nominalStateTrajectoriesStock,
			std::vector<control_vector_array_t>& nominalInputTrajectoriesStock,
			std::vector<output_vector_array_t>& nominalOutputTrajectoriesStock) = 0;

	virtual void run(const state_vector_t& initState, const std::vector<scalar_t>& switchingTimes,
			const std::vector<controller_t>& initialControllersStock=std::vector<controller_t>()) = 0;

	virtual std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM>>>& getSubsystemDynamicsPtrStock() = 0;

	virtual void setNewCostReferenceState(const output_vector_t& newReference) = 0;

	Options_t& options() {return options_;}

	void getIterationsLog(eigen_scalar_array_t& iterationCost, eigen_scalar_array_t& iterationISE1) const {
		iterationCost = iterationCost_;
		iterationISE1 = iterationISE1_;
	}

	void getSwitchingTimes(scalar_array_t& switchingTimes){switchingTimes = switchingTimes_;};

	size_t getNumIterations() const {return iteration_;}

protected:

	void solveSequentialRiccatiEquations(const scalar_t& learningRate);

	template <typename Derived>
	bool makePSD(Eigen::MatrixBase<Derived>& squareMatrix);

	scalar_t nominalTotalCost_;
	scalar_t nominalTotalMerit_;
	scalar_t nominalConstraint1ISE_;

	std::vector<scalar_array_t> 		nominalTimeTrajectoriesStock_;
	std::vector<state_vector_array_t>   nominalStateTrajectoriesStock_;
	std::vector<control_vector_array_t> nominalInputTrajectoriesStock_;
	std::vector<output_vector_array_t>  nominalOutputTrajectoriesStock_;
	std::vector<output_vector_array_t>  nominalcostateTrajectoriesStock_;
	std::vector<std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd> > >  nominalLagrangeTrajectoriesStock_;

	std::vector<lagrange_t> lagrangeControllerStock_;

	std::vector<state_matrix_array_t>        AmTrajectoryStock_;
	std::vector<control_gain_matrix_array_t> BmTrajectoryStock_;

	std::vector<std::vector<size_t> >       nc1TrajectoriesStock_;  	// nc1: Number of the Type-1  active constraints
	std::vector<constraint1_vector_array_t> EvTrajectoryStock_;
	std::vector<constraint1_state_matrix_array_t>   CmTrajectoryStock_;
	std::vector<constraint1_control_matrix_array_t> DmTrajectoryStock_;

	std::vector<std::vector<size_t> > 		nc2TrajectoriesStock_;  // nc2: Number of the Type-2 active constraints
	std::vector<constraint2_vector_array_t> HvTrajectoryStock_;
	std::vector<constraint2_state_matrix_array_t> FmTrajectoryStock_;
	std::vector<size_t> 					nc2FinalStock_;
	std::vector<constraint2_vector_t> 		HvFinalStock_;
	std::vector<constraint2_state_matrix_t> FmFinalStock_;

	std::vector<eigen_scalar_t>  qFinalStock_;
	std::vector<output_vector_t> QvFinalStock_;
	std::vector<state_matrix_t>  QmFinalStock_;

	std::vector<eigen_scalar_array_t> 	qTrajectoryStock_;
	std::vector<output_vector_array_t> 	QvTrajectoryStock_;
	std::vector<state_matrix_array_t> 	QmTrajectoryStock_;
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

	eigen_scalar_array_t iterationCost_;
	eigen_scalar_array_t iterationISE1_;

	Options_t options_;

public:
	template <size_t GSLQP_STATE_DIM, size_t GSLQP_INPUT_DIM, size_t GSLQP_OUTPUT_DIM, size_t GSLQP_NUM_SUBSYSTEMS>
	friend class GSLQP;
};

} // namespace ocs2

#include <GSLQ/implementation/SLQP_BASE.h>

#endif /* SLQP_H_ */
