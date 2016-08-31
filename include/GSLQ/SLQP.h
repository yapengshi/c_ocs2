/*
 * SLQP.h
 *
 *  Created on: Jun 16, 2016
 *      Author: farbod
 */

#ifndef SLQP_OCS2_H_
#define SLQP_OCS2_H_


#include "Dimensions.h"

#include "GSLQ/SLQP_BASE.h"


namespace ocs2{

template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
class SLQP : public SLQP_BASE<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	typedef SLQP_BASE<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS> BASE;

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
	typedef typename DIMENSIONS::constraint1_matrix_t constraint1_matrix_t;
	typedef typename DIMENSIONS::constraint1_matrix_array_t constraint1_matrix_array_t;
	typedef typename DIMENSIONS::constraint1_state_matrix_t constraint1_state_matrix_t;
	typedef typename DIMENSIONS::constraint1_state_matrix_array_t constraint1_state_matrix_array_t;
	typedef typename DIMENSIONS::constraint1_control_matrix_t constraint1_control_matrix_t;
	typedef typename DIMENSIONS::constraint1_control_matrix_array_t constraint1_control_matrix_array_t;
	typedef typename DIMENSIONS::control_constraint1_matrix_t control_constraint1_matrix_t;
	typedef typename DIMENSIONS::control_constraint1_matrix_array_t control_constraint1_matrix_array_t;


	SLQP(const std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > >& subsystemDynamicsPtr,
			const std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > >& subsystemDerivativesPtr,
			const std::vector<std::shared_ptr<CostFunctionBaseOCS2<OUTPUT_DIM, INPUT_DIM> > >& subsystemCostFunctionsPtr,
			const std::vector<controller_t>& initialControllersStock,
			const std::vector<size_t>& systemStockIndex,
			const Options_t& options = Options_t::Options())
    :
      subsystemDynamicsPtrStock_(NUM_SUBSYSTEMS),
      subsystemDerivativesPtrStock_(NUM_SUBSYSTEMS),
      subsystemCostFunctionsPtrStock_(NUM_SUBSYSTEMS),
      subsystemSimulatorsStockPtr_(NUM_SUBSYSTEMS),
      nominalControllersStock_(initialControllersStock),
      options_(options)
	{
		Eigen::initParallel();

		if (subsystemDynamicsPtr.size() != subsystemDerivativesPtr.size())
			throw std::runtime_error("Number of subsystem derivatives is not equal to the number of subsystems.");
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

	void rollout(const state_vector_t& initState,
			const std::vector<controller_t>& controllersStock,
			const double& stoppingTime,
			state_vector_t& stateVectorWhereStopped,
			control_vector_t& controlInputWhereStopped,
			output_vector_t& outputWhereStopped,
			size_t& numSubsystemWhereStopped) override;

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

	void setController(const std::vector<controller_t>& controllersStock) override;

	void getValueFuntion(const scalar_t& time, const output_vector_t& output, scalar_t& valueFuntion); // fixme: getValueFunction

	void getCostFuntion(const output_vector_t& initOutput, scalar_t& costFunction, scalar_t& constriantCostFunction);

	void getNominalTrajectories(std::vector<scalar_array_t>& nominalTimeTrajectoriesStock,
			std::vector<state_vector_array_t>& nominalStateTrajectoriesStock,
			std::vector<control_vector_array_t>& nominalInputTrajectoriesStock,
			std::vector<output_vector_array_t>& nominalOutputTrajectoriesStock);

	void run(const state_vector_t& initState, const std::vector<scalar_t>& switchingTimes);

	virtual std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM>>>& getSubsystemDynamicsPtrStock() override{
		return subsystemDynamicsPtrStock_;
	}

	void setNewCostReferenceState(const output_vector_t& newReference){	// todo: move to implementation
		for (size_t i = 0; i<subsystemCostFunctionsPtrStock_.size(); i++)
			subsystemCostFunctionsPtrStock_[i]->updateReferenceState(newReference);
	}


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

	template <typename Derived>
	bool makePSD(Eigen::MatrixBase<Derived>& squareMatrix);


private:
	std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM>>> subsystemDynamicsPtrStock_;

	std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM, OUTPUT_DIM> > > subsystemDerivativesPtrStock_;
	std::vector<std::shared_ptr<CostFunctionBaseOCS2<OUTPUT_DIM, INPUT_DIM> > > subsystemCostFunctionsPtrStock_;

	std::vector<std::shared_ptr<ODE45<STATE_DIM> > > subsystemSimulatorsStockPtr_;

	std::vector<controller_t> nominalControllersStock_;

	Options_t options_;

public:
	template <size_t GSLQP_STATE_DIM, size_t GSLQP_INPUT_DIM, size_t GSLQP_OUTPUT_DIM, size_t GSLQP_NUM_SUBSYSTEMS>
	friend class GSLQP;
};

} // namespace ocs2

#include "implementation/SLQP.h"


#endif /* SLQP_H_ */
