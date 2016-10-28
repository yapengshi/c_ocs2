/*
 * Implementation of SLQP.h
 *
 *  Created on: Jan 5, 2016
 *      Author: farbod, markus
 */


namespace ocs2{

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * Forward integrate the system dynamics with given controller:
 * 		inputs:
 * 			+ initState: initial state at time BASE::switchingTimes_[0]
 * 			+ controllersStock: controller for each subsystem
 * 		outputs:
 * 			+ timeTrajectoriesStock:  rollout simulated time steps
 * 			+ stateTrajectoriesStock: rollout states
 * 			+ inputTrajectoriesStock: rollout control inputs
 * 			+ (optional) outputTrajectoriesStock: rollout outputs
 * 			+ (optional) nc1TrajectoriesStock: number of active constraints at each time step
 * 			+ (optional) EvTrajectoryStock: value of the constraint (if the rollout is constrained the value is
 * 											always zero otherwise it is nonzero)
 */

/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::rollout(const state_vector_t& initState,
		const std::vector<controller_t>& controllersStock,
		std::vector<scalar_array_t>& timeTrajectoriesStock,
		std::vector<state_vector_array_t>& stateTrajectoriesStock,
		std::vector<control_vector_array_t>& inputTrajectoriesStock,
		std::vector<output_vector_array_t>& outputTrajectoriesStock)  {

	if (controllersStock.size() != NUM_SUBSYSTEMS)
		throw std::runtime_error("controllersStock has less controllers then the number of subsystems");

	timeTrajectoriesStock.resize(NUM_SUBSYSTEMS);
	stateTrajectoriesStock.resize(NUM_SUBSYSTEMS);
	inputTrajectoriesStock.resize(NUM_SUBSYSTEMS);
	outputTrajectoriesStock.resize(NUM_SUBSYSTEMS);

	state_vector_t x0 = initState;
	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		timeTrajectoriesStock[i].clear();
		stateTrajectoriesStock[i].clear();

		// max number of steps of integration
		size_t maxNumSteps = BASE::options_.maxNumStepsPerSecond_ * std::max( 1.0, BASE::switchingTimes_[i+1]-BASE::switchingTimes_[i]);

		// initialize subsystem i
		subsystemDynamicsPtrStock_[i]->initializeModel(BASE::switchingTimes_, x0, i, "GSLPQ");
		// set controller for subsystem i
		subsystemDynamicsPtrStock_[i]->setController(controllersStock[i]);

		// simulate subsystem i
		subsystemSimulatorsStockPtr_[i]->integrate(x0, BASE::switchingTimes_[i], BASE::switchingTimes_[i+1],
				stateTrajectoriesStock[i], timeTrajectoriesStock[i],
				1e-3, BASE::options_.AbsTolODE_, BASE::options_.RelTolODE_, maxNumSteps);

		if (stateTrajectoriesStock[i].back() != stateTrajectoriesStock[i].back())
			throw std::runtime_error("System became unstable during the SLQP rollout.");

		// compute control trajectory for subsystem i
		inputTrajectoriesStock[i].resize(timeTrajectoriesStock[i].size());
		outputTrajectoriesStock[i].resize(timeTrajectoriesStock[i].size());
		for (int k=0; k<timeTrajectoriesStock[i].size(); k++)
		{
			subsystemDynamicsPtrStock_[i]->computeOutput(timeTrajectoriesStock[i][k], stateTrajectoriesStock[i][k], outputTrajectoriesStock[i][k]);
			subsystemDynamicsPtrStock_[i]->computeInput(timeTrajectoriesStock[i][k], outputTrajectoriesStock[i][k], inputTrajectoriesStock[i][k]);
		}

		// reset the initial state
		x0 = stateTrajectoriesStock[i].back();
	}
}


template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::rollout(
		const state_vector_t& initState,
		const std::vector<controller_t>& controllersStock,
		const double& stoppingTime_in,
		state_vector_t& stateVectorWhereStopped,
		control_vector_t& controlInputWhereStopped,
		output_vector_t& outputWhereStopped,
		size_t& numSubsystemWhereStopped){

	double stoppingTime = stoppingTime_in;

	if(stoppingTime > BASE::switchingTimes_.back()){
		std::cerr << "WARNING: final integration time cannot be bigger than last switching time. Truncating stopping time." << std::endl;
		stoppingTime = BASE::switchingTimes_.back();
	}

	if (controllersStock.size() != NUM_SUBSYSTEMS)
		throw std::runtime_error("controllersStock has less controllers then the number of subsystems");

	std::vector<scalar_array_t> timeTrajectoriesStock;
	std::vector<state_vector_array_t> stateTrajectoriesStock;
	timeTrajectoriesStock.resize(NUM_SUBSYSTEMS);
	stateTrajectoriesStock.resize(NUM_SUBSYSTEMS);

	state_vector_t x0 = initState;
	int i = 0;
	while (i<NUM_SUBSYSTEMS) {

		timeTrajectoriesStock[i].clear();
		stateTrajectoriesStock[i].clear();

		size_t maxNumSteps = BASE::options_.maxNumStepsPerSecond_*(BASE::switchingTimes_[i+1]-BASE::switchingTimes_[i]);
		maxNumSteps = ((1000>maxNumSteps) ? 1000 : maxNumSteps);

		// initialize subsystem i
		subsystemDynamicsPtrStock_[i]->initializeModel(BASE::switchingTimes_, x0, i, "GSLPQ");
		// set controller for subsystem i
		subsystemDynamicsPtrStock_[i]->setController(controllersStock[i]);

		// determine correct stopping time for this subsystem
		double t_stop_thisSubsystem;
		bool stopAfterThisSubsystem = false;
		if(stoppingTime > BASE::switchingTimes_[i+1])
			t_stop_thisSubsystem = BASE::switchingTimes_[i+1];
		else{
			t_stop_thisSubsystem = stoppingTime;
			stopAfterThisSubsystem = true;
		}

		// simulate subsystem i
		subsystemSimulatorsStockPtr_[i]->integrate(
				x0, BASE::switchingTimes_[i], t_stop_thisSubsystem,
				stateTrajectoriesStock[i], timeTrajectoriesStock[i],
				1e-3, BASE::options_.AbsTolODE_, BASE::options_.RelTolODE_, maxNumSteps);

		if (stateTrajectoriesStock[i].back() != stateTrajectoriesStock[i].back())
			throw std::runtime_error("System became unstable during the SLQP rollout.");

		if(stopAfterThisSubsystem == true){
			numSubsystemWhereStopped = i;
			stateVectorWhereStopped = stateTrajectoriesStock[i].back();
			subsystemDynamicsPtrStock_[i]->computeOutput(timeTrajectoriesStock[i].back(), stateVectorWhereStopped, outputWhereStopped);
			subsystemDynamicsPtrStock_[i]->computeInput(timeTrajectoriesStock[i].back(), outputWhereStopped, controlInputWhereStopped);
			break;
		}

		// reset the initial state
		x0 = stateTrajectoriesStock[i].back();

		i++;
	}
}


/******************************************************************************************************/
// rolling out, where output trajectories are of no interest.
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::rollout(
		const state_vector_t& initState,
		const std::vector<controller_t>& controllersStock,
		std::vector<scalar_array_t>& timeTrajectoriesStock,
		std::vector<state_vector_array_t>& stateTrajectoriesStock,
		std::vector<control_vector_array_t>& inputTrajectoriesStock)  {

	std::vector<output_vector_array_t> outputTrajectoriesStock;
	rollout(initState, controllersStock, timeTrajectoriesStock,
			stateTrajectoriesStock, inputTrajectoriesStock, outputTrajectoriesStock);
}

/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::rollout(
		const state_vector_t& initState,
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
		std::vector<constraint2_vector_t>& HvFinalStock) {

	// STEP1 : perform normal rollout
	rollout(initState, controllersStock,
			timeTrajectoriesStock, stateTrajectoriesStock, inputTrajectoriesStock, outputTrajectoriesStock);


	// STEP2 : calculate constraint violations

	// constraint type 1 computations which consists of number of active constraints at each time point
	// and the value of the constraint (if the rollout is constrained the value is always zero otherwise
	// it is nonzero)
	nc1TrajectoriesStock.resize(NUM_SUBSYSTEMS);
	EvTrajectoryStock.resize(NUM_SUBSYSTEMS);

	// constraint type 2 computations which consists of number of active constraints at each time point
	// and the value of the constraint
	nc2TrajectoriesStock.resize(NUM_SUBSYSTEMS);
	HvTrajectoryStock.resize(NUM_SUBSYSTEMS);
	nc2FinalStock.resize(NUM_SUBSYSTEMS);
	HvFinalStock.resize(NUM_SUBSYSTEMS);


	for (int i=0; i<NUM_SUBSYSTEMS; i++)
	{
		size_t N = timeTrajectoriesStock[i].size();
		nc1TrajectoriesStock[i].resize(N);
		EvTrajectoryStock[i].resize(N);
		nc2TrajectoriesStock[i].resize(N);
		HvTrajectoryStock[i].resize(N);

		// compute constraint1 trajectory for subsystem i
		for (int k=0; k<N; k++) {

			// constraint 1 type
			subsystemDynamicsPtrStock_[i]->computeConstriant1(timeTrajectoriesStock[i][k],
					stateTrajectoriesStock[i][k], inputTrajectoriesStock[i][k],
					nc1TrajectoriesStock[i][k], EvTrajectoryStock[i][k]);

			if (nc1TrajectoriesStock[i][k] > INPUT_DIM)
				throw std::runtime_error("Number of active type-1 constraints should be less-equal to the number of input dimension.");

			// constraint type 2
			subsystemDynamicsPtrStock_[i]->computeConstriant2(timeTrajectoriesStock[i][k],
					stateTrajectoriesStock[i][k],
					nc2TrajectoriesStock[i][k], HvTrajectoryStock[i][k]);

			if (nc2TrajectoriesStock[i][k] > INPUT_DIM)
				throw std::runtime_error("Number of active type-2 constraints should be less-equal to the number of input dimension.");

		}  // end of k loop

		subsystemDynamicsPtrStock_[i]->computeFinalConstriant2(timeTrajectoriesStock[i].back(), stateTrajectoriesStock[i].back(),
				nc2FinalStock[i], HvFinalStock[i]);
		if (nc2FinalStock[i] > INPUT_DIM)
			throw std::runtime_error("Number of active type-2 constraints at final time should be less-equal to the number of input dimension.");

	}  // end of i loop

}


/*****************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateCostFunction(
		const std::vector<scalar_array_t>& timeTrajectoriesStock,
		const std::vector<output_vector_array_t>& outputTrajectoriesStock,
		const std::vector<control_vector_array_t>& inputTrajectoriesStock,
		const std::vector<std::vector<size_t> >& nc2TrajectoriesStock,
		const std::vector<constraint2_vector_array_t>& HvTrajectoryStock,
		const std::vector<size_t>& nc2FinalStock,
		const std::vector<constraint2_vector_t>& HvFinalStock,
		scalar_t& totalCost) {

	calculateCostFunction(timeTrajectoriesStock, outputTrajectoriesStock, inputTrajectoriesStock, totalCost);
	double stateConstraintPenalty = BASE::options_.stateConstraintPenaltyCoeff_ * pow(BASE::options_.stateConstraintPenaltyBase_, BASE::iteration_);

	for (int i=0; i<NUM_SUBSYSTEMS; i++) {
		// integrates constraint type 2
		for (int k=0; k<timeTrajectoriesStock[i].size()-1; k++) {
			size_t nc2 = nc2TrajectoriesStock[i][k];
			if (nc2 > 0) {
				double dt = timeTrajectoriesStock[i][k+1]-timeTrajectoriesStock[i][k];
				totalCost += 0.5 * dt * stateConstraintPenalty * HvTrajectoryStock[i][k].head(nc2).squaredNorm();
			}
		}  // end of k loop

		// final constraint type 2
		size_t nc2Final = nc2FinalStock[i];
		if (nc2Final>0)
			totalCost += 0.5 * stateConstraintPenalty * HvFinalStock[i].head(nc2Final).squaredNorm();

	}  // end of i loop
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * compute the cost for a given rollout
 * 		inputs:
 * 			+ timeTrajectoriesStock:  rollout simulated time steps
 * 			+ outputTrajectoriesStock: rollout outputs
 * 			+ inputTrajectoriesStock: rollout control inputs
 *
 * 		outputs:
 * 			+ totalCost: the total cost of the trajectory
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateCostFunction(
		const std::vector<scalar_array_t>& timeTrajectoriesStock,
		const std::vector<output_vector_array_t>& outputTrajectoriesStock,
		const std::vector<control_vector_array_t>& inputTrajectoriesStock,
		scalar_t& totalCost)  {

	totalCost = 0.0;
	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		// integrates the intermediate cost using the trapezoidal approximation method
		scalar_t currentIntermediateCost;
		scalar_t nextIntermediateCost;
		for (int k=0; k<timeTrajectoriesStock[i].size()-1; k++) {

			if (k==0) {
				subsystemCostFunctionsPtrStock_[i]->setCurrentStateAndControl(timeTrajectoriesStock[i][k], outputTrajectoriesStock[i][k], inputTrajectoriesStock[i][k]);
				subsystemCostFunctionsPtrStock_[i]->evaluate(currentIntermediateCost);
			} else {
				currentIntermediateCost = nextIntermediateCost;
			}

			// feed next state and control to cost function
			subsystemCostFunctionsPtrStock_[i]->setCurrentStateAndControl(timeTrajectoriesStock[i][k+1], outputTrajectoriesStock[i][k+1], inputTrajectoriesStock[i][k+1]);
			// evaluate intermediate cost for next time step
			subsystemCostFunctionsPtrStock_[i]->evaluate(nextIntermediateCost);

			totalCost += 0.5*(currentIntermediateCost+nextIntermediateCost)*(timeTrajectoriesStock[i][k+1]-timeTrajectoriesStock[i][k]);
		}  // end of k loop

		// terminal cost
		if (i==NUM_SUBSYSTEMS-1)  {
			scalar_t finalCost;
			subsystemCostFunctionsPtrStock_[i]->setCurrentStateAndControl(timeTrajectoriesStock[i].back(), outputTrajectoriesStock[i].back(), inputTrajectoriesStock[i].back());
			subsystemCostFunctionsPtrStock_[i]->terminalCost(finalCost);
			totalCost += finalCost;
		}

	}  // end of i loop
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * compute the merit function for given rollout
 * 		inputs:
 * 			+ timeTrajectoriesStock: simulation time trajectory
 * 			+ nc1TrajectoriesStock: rollout's number of active constraints in each time step
 * 			+ EvTrajectoryStock: rollout's constraints value
 * 			+ lagrangeTrajectoriesStock: constraint Lagrange multiplier for the given rollout
 * 			+ totalCost: the total cost of the trajectory
 *
 * 		outputs:
 * 			+ meritFuntionValue: the merit function value
 * 			+ constraintISE: integral of Square Error (ISE)
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateMeritFunction(
		const std::vector<scalar_array_t>& timeTrajectoriesStock,
		const std::vector<std::vector<size_t> >& nc1TrajectoriesStock,
		const std::vector<constraint1_vector_array_t>& EvTrajectoryStock,
		const std::vector<std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd> > >&  lagrangeTrajectoriesStock,
		const scalar_t& totalCost,
		scalar_t& meritFunctionValue,
		scalar_t& constraintISE)  {

	// add cost function
	meritFunctionValue = totalCost;

	// add the L2 penalty for constraint violation
	calculateConstraintISE(timeTrajectoriesStock, nc1TrajectoriesStock, EvTrajectoryStock, constraintISE);
	double pho = BASE::iteration_/(BASE::options_.maxIterationGSLQP_-1) * BASE::options_.meritFunctionRho_;	// TODO FIXME - leads ot segfault if maxIterationGSLQP = 1

	meritFunctionValue += 0.5*pho*constraintISE;


	// add the the lagrangian term for the constraint
	scalar_t currentIntermediateMerit;
	scalar_t nextIntermediateMerit;

	for (int i=0; i<NUM_SUBSYSTEMS; i++)
	{
		// integrates the intermediate merit using the trapezoidal approximation method
		currentIntermediateMerit = 0.0;
		nextIntermediateMerit = 0.0;
		for (int k=0; k<timeTrajectoriesStock[i].size()-1; k++)
		{
			if (k==0)
				currentIntermediateMerit = EvTrajectoryStock[i][k].head(nc1TrajectoriesStock[i][k]).transpose() * lagrangeTrajectoriesStock[i][k];
			else
				currentIntermediateMerit = nextIntermediateMerit;

			nextIntermediateMerit = EvTrajectoryStock[i][k+1].head(nc1TrajectoriesStock[i][k+1]).transpose() * lagrangeTrajectoriesStock[i][k+1];

			meritFunctionValue += 0.5*(currentIntermediateMerit+nextIntermediateMerit)*(timeTrajectoriesStock[i][k+1]-timeTrajectoriesStock[i][k]);
		}  // end of k loop
	}  // end of i loop

}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * Constraint's Integral of Squared Error (ISE)
 * 		Inputs:
 * 			+ timeTrajectoriesStock: simulation time trajectory
 * 			+ nc1TrajectoriesStock: rollout's number of active constraints in each time step
 * 			+ EvTrajectoriesStock: rollout's constraints value
 * 		Output:
 * 			+ constraintISE: integral of Square Error (ISE)
 * 		Return:
 * 			+ maximum constraint norm
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
double SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateConstraintISE(
		const std::vector<scalar_array_t>& timeTrajectoriesStock,
		const std::vector<std::vector<size_t>>& nc1TrajectoriesStock,
		const std::vector<constraint1_vector_array_t>& EvTrajectoriesStock,
		scalar_t& constraintISE)  {

	constraintISE = 0.0;

	double maxConstraintNorm = 0.0;
	scalar_t currentSquaredNormError;
	scalar_t nextSquaredNormError;

	for (size_t i=0; i<NUM_SUBSYSTEMS; i++)  {

		currentSquaredNormError = 0.0;
		nextSquaredNormError = 0.0;

		for (size_t k=0; k<timeTrajectoriesStock[i].size()-1; k++)  {

			if (k==0) {
				const size_t& nc1 = nc1TrajectoriesStock[i][0];
				if (nc1>0)
					currentSquaredNormError = EvTrajectoriesStock[i][0].head(nc1).squaredNorm();
				else
					currentSquaredNormError = 0.0;
			} else
				currentSquaredNormError = nextSquaredNormError;

			maxConstraintNorm = ((maxConstraintNorm<currentSquaredNormError)? currentSquaredNormError: maxConstraintNorm);

			const size_t& nc1 = nc1TrajectoriesStock[i][k+1];
			if (nc1>0)
				nextSquaredNormError = EvTrajectoriesStock[i][k+1].head(nc1).squaredNorm();
			else
				nextSquaredNormError = 0.0;

			constraintISE += 0.5 * (currentSquaredNormError+nextSquaredNormError) * (timeTrajectoriesStock[i][k+1]-timeTrajectoriesStock[i][k]);

		}  // end of k loop
	}  // end of i loop

	return sqrt(maxConstraintNorm);
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * approximates the nonlinear problem as a linear-quadratic problem around the nominal
 * state and control trajectories. This method updates the following variables:
 *
 * 		+ linearized system model
 * 		+ dxdt = Am(t)x(t) + Bm(t)u(t)
 * 		+ s.t. Cm(t)x(t) + Dm(t)t(t) + Ev(t) = 0
 * 		+ BASE::AmTrajectoryStock_: Am matrix
 * 		+ BASE::BmTrajectoryStock_: Bm matrix
 * 		+ BASE::CmTrajectoryStock_: Cm matrix
 * 		+ BASE::DmTrajectoryStock_: Dm matrix
 * 		+ BASE::EvTrajectoryStock_: Ev vector
 *
 * 		+ quadratized intermediate cost function
 * 		+ intermediate cost: q(t) + 0.5 y(t)Qm(t)y(t) + y(t)'Qv(t) + u(t)'Pm(t)y(t) + 0.5u(t)'Rm(t)u(t) + u(t)'Rv(t)
 * 		+ BASE::qTrajectoryStock_:  q
 * 		+ BASE::QvTrajectoryStock_: Qv vector
 * 		+ BASE::QmTrajectoryStock_: Qm matrix
 * 		+ BASE::PmTrajectoryStock_: Pm matrix
 * 		+ BASE::RvTrajectoryStock_: Rv vector
 * 		+ BASE::RmTrajectoryStock_: Rm matrix
 * 		+ BASE::RmInverseTrajectoryStock_: inverse of Rm matrix
 *
 * 		+ as well as the constrained coefficients of
 * 			linearized system model
 * 			quadratized intermediate cost function
 * 			quadratized final cost
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::approximateOptimalControlProblem()  {

	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		// initialize subsystem i dynamics derivatives
		subsystemDerivativesPtrStock_[i]->initializeModel(BASE::switchingTimes_, BASE::nominalStateTrajectoriesStock_[i].front(), i, "GSLPQ");

		int N = BASE::nominalTimeTrajectoriesStock_[i].size();

		BASE::AmTrajectoryStock_[i].resize(N);
		BASE::BmTrajectoryStock_[i].resize(N);
		BASE::CmTrajectoryStock_[i].resize(N);
		BASE::DmTrajectoryStock_[i].resize(N);
		BASE::FmTrajectoryStock_[i].resize(N);

		BASE::qTrajectoryStock_[i].resize(N);
		BASE::QvTrajectoryStock_[i].resize(N);
		BASE::QmTrajectoryStock_[i].resize(N);
		BASE::RvTrajectoryStock_[i].resize(N);
		BASE::RmTrajectoryStock_[i].resize(N);
		BASE::RmInverseTrajectoryStock_[i].resize(N);
		BASE::PmTrajectoryStock_[i].resize(N);

		for (int k=0; k<N; k++)
		{
			// LINEARIZE SYSTEM DYNAMICS AND CONSTRAINTS
			subsystemDerivativesPtrStock_[i]->setCurrentStateAndControl(BASE::nominalTimeTrajectoriesStock_[i][k],
					BASE::nominalStateTrajectoriesStock_[i][k], BASE::nominalInputTrajectoriesStock_[i][k], BASE::nominalOutputTrajectoriesStock_[i][k]);
			subsystemDerivativesPtrStock_[i]->getDerivativeState(BASE::AmTrajectoryStock_[i][k]);
			subsystemDerivativesPtrStock_[i]->getDerivativesControl(BASE::BmTrajectoryStock_[i][k]);

			// if constraint type 1 is active
			if (BASE::nc1TrajectoriesStock_[i][k] > 0) {
				subsystemDerivativesPtrStock_[i]->getConstraint1DerivativesState(BASE::CmTrajectoryStock_[i][k]);
				subsystemDerivativesPtrStock_[i]->getConstraint1DerivativesControl(BASE::DmTrajectoryStock_[i][k]);
			}

			// if constraint type 2 is active
			if (BASE::nc2TrajectoriesStock_[i][k] > 0) {
				subsystemDerivativesPtrStock_[i]->getConstraint2DerivativesState(BASE::FmTrajectoryStock_[i][k]);
			}

			// QUADRATIC APPROXIMATION TO THE COST FUNCTION
			subsystemCostFunctionsPtrStock_[i]->setCurrentStateAndControl(BASE::nominalTimeTrajectoriesStock_[i][k],
					BASE::nominalOutputTrajectoriesStock_[i][k], BASE::nominalInputTrajectoriesStock_[i][k]);
			subsystemCostFunctionsPtrStock_[i]->evaluate(BASE::qTrajectoryStock_[i][k](0));
			subsystemCostFunctionsPtrStock_[i]->stateDerivative(BASE::QvTrajectoryStock_[i][k]);
			subsystemCostFunctionsPtrStock_[i]->stateSecondDerivative(BASE::QmTrajectoryStock_[i][k]);
			subsystemCostFunctionsPtrStock_[i]->controlDerivative(BASE::RvTrajectoryStock_[i][k]);
			subsystemCostFunctionsPtrStock_[i]->controlSecondDerivative(BASE::RmTrajectoryStock_[i][k]);
			BASE::RmInverseTrajectoryStock_[i][k] = BASE::RmTrajectoryStock_[i][k].inverse();
			subsystemCostFunctionsPtrStock_[i]->stateControlDerivative(BASE::PmTrajectoryStock_[i][k]);

		}


		if (i==NUM_SUBSYSTEMS-1)  {
			subsystemCostFunctionsPtrStock_[i]->terminalCost(BASE::qFinalStock_[i](0));
			subsystemCostFunctionsPtrStock_[i]->terminalCostStateDerivative(BASE::QvFinalStock_[i]);
			subsystemCostFunctionsPtrStock_[i]->terminalCostStateSecondDerivative(BASE::QmFinalStock_[i]);
			// making sure that Qm remains PSD
			this->makePSD(BASE::QmFinalStock_[i]);
		}
		else {
			BASE::qFinalStock_[i].setZero();
			BASE::QvFinalStock_[i].setZero();
			BASE::QmFinalStock_[i].setZero();
		}
	}


	// constraint type-2 coefficients
	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		int N = BASE::nominalTimeTrajectoriesStock_[i].size();

		// constraint type-2 coefficients
		for (int k=0; k<N; k++) {

			double stateConstraintPenalty = BASE::options_.stateConstraintPenaltyCoeff_ * pow(BASE::options_.stateConstraintPenaltyBase_, BASE::iteration_);

			// constrained type-2 intermediate coefficients
			size_t nc2 = BASE::nc2TrajectoriesStock_[i][k];

			if (nc2 > 0) {

				//subsystemDerivativesPtrStock_[i]->getConstraint2DerivativesState(FmTrajectoryStock_[i][k]);

				BASE::qTrajectoryStock_[i][k]  += 0.5 * stateConstraintPenalty * BASE::HvTrajectoryStock_[i][k].head(nc2).transpose() * BASE::HvTrajectoryStock_[i][k].head(nc2);
				BASE::QvTrajectoryStock_[i][k] += stateConstraintPenalty * BASE::FmTrajectoryStock_[i][k].topRows(nc2).transpose() * BASE::HvTrajectoryStock_[i][k].head(nc2);
				BASE::QmTrajectoryStock_[i][k] += stateConstraintPenalty * BASE::FmTrajectoryStock_[i][k].topRows(nc2).transpose() * BASE::FmTrajectoryStock_[i][k].topRows(nc2);
			}
		}  // end of k loop

		// constrained type-2 final coefficients
		if (BASE::nc2FinalStock_[i] > 0) {
			size_t nc2 = BASE::nc2FinalStock_[i];

			subsystemDerivativesPtrStock_[i]->getFinalConstraint2DerivativesState(BASE::FmFinalStock_[i]);

			double stateConstraintPenalty = BASE::options_.stateConstraintPenaltyCoeff_ * pow(BASE::options_.stateConstraintPenaltyBase_, BASE::iteration_);

			BASE::qFinalStock_[i]  += 0.5 * stateConstraintPenalty * BASE::HvFinalStock_[i].head(nc2).transpose() * BASE::HvFinalStock_[i].head(nc2);
			BASE::QvFinalStock_[i] += stateConstraintPenalty * BASE::FmFinalStock_[i].topRows(nc2).transpose() * BASE::HvFinalStock_[i].head(nc2);
			BASE::QmFinalStock_[i] += stateConstraintPenalty * BASE::FmFinalStock_[i].topRows(nc2).transpose() * BASE::FmFinalStock_[i].topRows(nc2);
		}
	} // end of i loop


	// constraint type 1 coefficients
	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		int N = BASE::nominalTimeTrajectoriesStock_[i].size();

		BASE::RmConstrainedTrajectoryStock_[i].resize(N);
		BASE::DmDagerTrajectoryStock_[i].resize(N);

		BASE::AmConstrainedTrajectoryStock_[i].resize(N);
		BASE::QmConstrainedTrajectoryStock_[i].resize(N);
		BASE::QvConstrainedTrajectoryStock_[i].resize(N);

		BASE::EvProjectedTrajectoryStock_[i].resize(N);
		BASE::CmProjectedTrajectoryStock_[i].resize(N);
		BASE::DmProjectedTrajectoryStock_[i].resize(N);

		for (int k=0; k<N; k++) {
			size_t nc1 = BASE::nc1TrajectoriesStock_[i][k];

			if (nc1 == 0) {
				BASE::DmDagerTrajectoryStock_[i][k].setZero();
				BASE::EvProjectedTrajectoryStock_[i][k].setZero();
				BASE::CmProjectedTrajectoryStock_[i][k].setZero();
				BASE::DmProjectedTrajectoryStock_[i][k].setZero();

				BASE::AmConstrainedTrajectoryStock_[i][k] = BASE::AmTrajectoryStock_[i][k];
				BASE::QmConstrainedTrajectoryStock_[i][k] = BASE::QmTrajectoryStock_[i][k];
				BASE::QvConstrainedTrajectoryStock_[i][k] = BASE::QvTrajectoryStock_[i][k];
				BASE::RmConstrainedTrajectoryStock_[i][k] = BASE::RmTrajectoryStock_[i][k];

			} else {
				Eigen::MatrixXd Cm = BASE::CmTrajectoryStock_[i][k].topRows(nc1);
				Eigen::MatrixXd Dm = BASE::DmTrajectoryStock_[i][k].topRows(nc1);
				Eigen::MatrixXd Ev = BASE::EvTrajectoryStock_[i][k].head(nc1);

				Eigen::MatrixXd RmProjected = ( Dm*BASE::RmInverseTrajectoryStock_[i][k]*Dm.transpose() ).inverse();
				Eigen::MatrixXd DmDager = BASE::RmInverseTrajectoryStock_[i][k] * Dm.transpose() * RmProjected;

				BASE::DmDagerTrajectoryStock_[i][k].leftCols(nc1) = DmDager;
				BASE::EvProjectedTrajectoryStock_[i][k] = DmDager * Ev;
				BASE::CmProjectedTrajectoryStock_[i][k] = DmDager * Cm;
				BASE::DmProjectedTrajectoryStock_[i][k] = DmDager * Dm;

				control_matrix_t DmNullSpaceProjection = control_matrix_t::Identity() - BASE::DmProjectedTrajectoryStock_[i][k];
				state_matrix_t   PmTransDmDagerCm = BASE::PmTrajectoryStock_[i][k].transpose()*BASE::CmProjectedTrajectoryStock_[i][k];

				BASE::AmConstrainedTrajectoryStock_[i][k] = BASE::AmTrajectoryStock_[i][k] - BASE::BmTrajectoryStock_[i][k]*BASE::CmProjectedTrajectoryStock_[i][k];
				BASE::QmConstrainedTrajectoryStock_[i][k] = BASE::QmTrajectoryStock_[i][k] + Cm.transpose()*RmProjected*Cm - PmTransDmDagerCm - PmTransDmDagerCm.transpose();
				BASE::QvConstrainedTrajectoryStock_[i][k] = BASE::QvTrajectoryStock_[i][k] - BASE::CmProjectedTrajectoryStock_[i][k].transpose()*BASE::RvTrajectoryStock_[i][k];
				BASE::RmConstrainedTrajectoryStock_[i][k] = DmNullSpaceProjection.transpose() * BASE::RmTrajectoryStock_[i][k] * DmNullSpaceProjection;
			}

			// making sure that constrained Qm is PSD
			this->makePSD(BASE::QmConstrainedTrajectoryStock_[i][k]);

		}  // end of k loop
	}  // end of i loop
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * calculates the controller and linear function approximation of the type-1 constraint Lagrangian:
 * 		This method uses the following variables:
 * 			+ constrained, linearized model
 * 			+ constrained, quadratized cost
 *
 * 		The method outputs:
 * 			+ controllersStock: the controller that stabilizes the system around the new nominal trajectory and
 * 								improves the constraints as well as the increment to the feedforward control input.
 * 			+ lagrangeMultiplierFunctionsStock: the linear function approximation of the type-1 constraint Lagrangian.
 *			+ feedForwardConstraintInputStock: the increment to the feedforward control input due to constraint type-1
 *			+ firstCall: True if it is the first time that this method is called. This is required for the
 *			             lagrangeMultiplierFunctionsStock's feedforward part.
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateControllerAndLagrangian(
		std::vector<controller_t>& controllersStock,
		std::vector<lagrange_t>& lagrangeMultiplierFunctionsStock,
		std::vector<control_vector_array_t>& feedForwardConstraintInputStock,
		bool firstCall /*=true*/) {

	// functions for controller and lagrane multiplier
	LinearInterpolation<output_vector_t,Eigen::aligned_allocator<output_vector_t> >   nominalOutputFunc;
	LinearInterpolation<control_vector_t,Eigen::aligned_allocator<control_vector_t> > nominalInputFunc;

	LinearInterpolation<control_gain_matrix_t,Eigen::aligned_allocator<control_gain_matrix_t> > BmFunc;
	LinearInterpolation<control_feedback_t,Eigen::aligned_allocator<control_feedback_t> > PmFunc;
	LinearInterpolation<control_matrix_t,Eigen::aligned_allocator<control_matrix_t> >     RmInverseFunc;
	LinearInterpolation<control_vector_t,Eigen::aligned_allocator<control_vector_t> >     RvFunc;
	LinearInterpolation<control_vector_t,Eigen::aligned_allocator<control_vector_t> >     EvProjectedFunc;
	LinearInterpolation<control_feedback_t,Eigen::aligned_allocator<control_feedback_t> > CmProjectedFunc;
	LinearInterpolation<control_matrix_t,Eigen::aligned_allocator<control_matrix_t> >     DmProjectedFunc;

	// functions for lagrane multiplier only
	LinearInterpolation<Eigen::VectorXd,Eigen::aligned_allocator<Eigen::VectorXd> > nominalLagrangeMultiplierFunc;

	LinearInterpolation<control_constraint1_matrix_t,Eigen::aligned_allocator<control_constraint1_matrix_t> > DmDagerFunc;
	LinearInterpolation<control_matrix_t,Eigen::aligned_allocator<control_matrix_t> > RmFunc;


	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		// functions for controller and lagrane multiplier

		nominalOutputFunc.setTimeStamp( &(BASE::nominalTimeTrajectoriesStock_[i]) );
		nominalOutputFunc.setData( &(BASE::nominalOutputTrajectoriesStock_[i]) );

		nominalInputFunc.setTimeStamp( &(BASE::nominalTimeTrajectoriesStock_[i]) );
		nominalInputFunc.setData( &(BASE::nominalInputTrajectoriesStock_[i]) );

		BmFunc.setTimeStamp( &(BASE::nominalTimeTrajectoriesStock_[i]) );
		BmFunc.setData( &(BASE::BmTrajectoryStock_[i]) );

		PmFunc.setTimeStamp( &(BASE::nominalTimeTrajectoriesStock_[i]) );
		PmFunc.setData( &(BASE::PmTrajectoryStock_[i]) );

		RmInverseFunc.setTimeStamp( &(BASE::nominalTimeTrajectoriesStock_[i]) );
		RmInverseFunc.setData( &(BASE::RmInverseTrajectoryStock_[i]) );

		RvFunc.setTimeStamp( &(BASE::nominalTimeTrajectoriesStock_[i]) );
		RvFunc.setData( &(BASE::RvTrajectoryStock_[i]) );

		EvProjectedFunc.setTimeStamp( &(BASE::nominalTimeTrajectoriesStock_[i]) );
		EvProjectedFunc.setData( &(BASE::EvProjectedTrajectoryStock_[i]) );

		CmProjectedFunc.setTimeStamp( &(BASE::nominalTimeTrajectoriesStock_[i]) );
		CmProjectedFunc.setData( &(BASE::CmProjectedTrajectoryStock_[i]) );

		DmProjectedFunc.setTimeStamp( &(BASE::nominalTimeTrajectoriesStock_[i]) );
		DmProjectedFunc.setData( &(BASE::DmProjectedTrajectoryStock_[i]) );

		// functions for lagrange multiplier only

		RmFunc.setTimeStamp( &(BASE::nominalTimeTrajectoriesStock_[i]) );
		RmFunc.setData( &(BASE::RmTrajectoryStock_[i]) );

		DmDagerFunc.setTimeStamp( &(BASE::nominalTimeTrajectoriesStock_[i]) );
		DmDagerFunc.setData( &(BASE::DmDagerTrajectoryStock_[i]) );

		if (firstCall==false) {
			nominalLagrangeMultiplierFunc.setTimeStamp( &(BASE::nominalTimeTrajectoriesStock_[i]) );
			nominalLagrangeMultiplierFunc.setData( &(BASE::nominalLagrangeTrajectoriesStock_[i]) );
		}

		int N = BASE::SsTimeTrajectoryStock_[i].size();

		controllersStock[i].time_ = BASE::SsTimeTrajectoryStock_[i];
		controllersStock[i].k_.resize(N);
		controllersStock[i].uff_.resize(N);
		controllersStock[i].deltaUff_.resize(N);

		feedForwardConstraintInputStock[i].resize(N);

		lagrangeMultiplierFunctionsStock[i].time_ = BASE::SsTimeTrajectoryStock_[i];
		lagrangeMultiplierFunctionsStock[i].k_.resize(N);
		lagrangeMultiplierFunctionsStock[i].uff_.resize(N);
		lagrangeMultiplierFunctionsStock[i].deltaUff_.resize(N);

		for (int k=0; k<N; k++) {

			const double& time = BASE::SsTimeTrajectoryStock_[i][k];
			size_t greatestLessTimeStampIndex;

			output_vector_t nominalOutput;
			nominalOutputFunc.interpolate(time, nominalOutput);
			greatestLessTimeStampIndex = nominalOutputFunc.getGreatestLessTimeStampIndex();
			control_vector_t nominalInput;
			nominalInputFunc.interpolate(time, nominalInput, greatestLessTimeStampIndex);

			control_gain_matrix_t Bm;
			BmFunc.interpolate(time, Bm, greatestLessTimeStampIndex);
			control_feedback_t Pm;
			PmFunc.interpolate(time, Pm, greatestLessTimeStampIndex);
			control_vector_t Rv;
			RvFunc.interpolate(time, Rv, greatestLessTimeStampIndex);
			control_matrix_t RmInverse;
			RmInverseFunc.interpolate(time, RmInverse, greatestLessTimeStampIndex);
			control_vector_t EvProjected;
			EvProjectedFunc.interpolate(time, EvProjected, greatestLessTimeStampIndex);
			control_feedback_t CmProjected;
			CmProjectedFunc.interpolate(time, CmProjected, greatestLessTimeStampIndex);
			control_matrix_t DmProjected;
			DmProjectedFunc.interpolate(time, DmProjected, greatestLessTimeStampIndex);

			control_feedback_t Lm  = RmInverse * (Pm + Bm.transpose()*BASE::SmTrajectoryStock_[i][k]);
			control_vector_t   Lv  = RmInverse * (Rv + Bm.transpose()*BASE::SvTrajectoryStock_[i][k]);
			control_vector_t   Lve = RmInverse * (Bm.transpose()*BASE::SveTrajectoryStock_[i][k]);

			control_matrix_t DmNullProjection = control_matrix_t::Identity()-DmProjected;
			controllersStock[i].k_[k]   = -DmNullProjection*Lm - CmProjected;
			controllersStock[i].uff_[k] = nominalInput - controllersStock[i].k_[k]*nominalOutput;
			controllersStock[i].deltaUff_[k]      = -DmNullProjection*Lv;
			feedForwardConstraintInputStock[i][k] = -DmNullProjection*Lve - EvProjected;

			// checking the numerical stability of the controller parameters
			try {
				if (controllersStock[i].k_[k] != controllersStock[i].k_[k])
					throw std::runtime_error("Feedback gains are unstable.");
				if (controllersStock[i].deltaUff_[k] != controllersStock[i].deltaUff_[k])
					throw std::runtime_error("feedForwardControl is unstable.");
				if (feedForwardConstraintInputStock[i][k] != feedForwardConstraintInputStock[i][k])
					throw std::runtime_error("feedForwardConstraintInput is unstable.");
			}
			catch(const std::exception& error)  {
				std::cerr << "what(): " << error.what() << " at time " << controllersStock[i].time_[k] << " [sec]." << std::endl;
			}


			// lagrane multiplier calculation

			const size_t& nc1 = BASE::nc1TrajectoriesStock_[i][greatestLessTimeStampIndex];

			control_constraint1_matrix_t DmDager;
			DmDagerFunc.interpolate(time, DmDager, greatestLessTimeStampIndex);
			control_matrix_t Rm;
			RmFunc.interpolate(time, Rm, greatestLessTimeStampIndex);

			Eigen::MatrixXd DmDagerTransRm = DmDager.leftCols(nc1).transpose() * Rm;

			lagrangeMultiplierFunctionsStock[i].k_[k]   = DmDagerTransRm * (CmProjected - Lm);
			lagrangeMultiplierFunctionsStock[i].uff_[k] = -lagrangeMultiplierFunctionsStock[i].k_[k]*nominalOutput;
			Eigen::VectorXd localVff = DmDagerTransRm * (EvProjected-Lv-Lve);

			if (firstCall==true || BASE::options_.lineSearchByMeritFuntion_==false) {
				lagrangeMultiplierFunctionsStock[i].uff_[k] += localVff;
				lagrangeMultiplierFunctionsStock[i].deltaUff_[k] = Eigen::VectorXd::Zero(nc1);

			} else {
				Eigen::VectorXd nominalLagrangeMultiplier;
				nominalLagrangeMultiplierFunc.interpolate(time, nominalLagrangeMultiplier, greatestLessTimeStampIndex);
				lagrangeMultiplierFunctionsStock[i].uff_[k] += nominalLagrangeMultiplier;

				lagrangeMultiplierFunctionsStock[i].deltaUff_[k] = localVff - nominalLagrangeMultiplier;
			}

			// checking the numerical stability of the controller parameters
			try {
				if (lagrangeMultiplierFunctionsStock[i].k_[k] != lagrangeMultiplierFunctionsStock[i].k_[k])
					throw std::runtime_error("Feedback lagrangeMultiplier is unstable.");
				if (lagrangeMultiplierFunctionsStock[i].deltaUff_[k] != lagrangeMultiplierFunctionsStock[i].deltaUff_[k])
					throw std::runtime_error("Feedforward lagrangeMultiplier is unstable.");
			}
			catch(const std::exception& error)  {
				std::cerr << "what(): " << error.what() << " at time " << lagrangeMultiplierFunctionsStock[i].time_[k] << " [sec]." << std::endl;
			}


		}  // end of k loop
	}  // end of i loop

}



/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * compute the Lagrage multiplier over the given rollout
 * 		inputs:
 * 			+ timeTrajectoriesStock: rollout simulated time steps
 * 			+ outputTrajectoriesStock: rollout outputs
 * 			+ lagrangeMultiplierFunctionsStock: the coefficients of the linear function for lagrangeMultiplier
 *
 * 		outputs:
 * 			+ lagrangeTrajectoriesStock: lagrangeMultiplier value over the given trajectory
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateRolloutLagrangeMultiplier(
		const std::vector<scalar_array_t>& timeTrajectoriesStock,
		const std::vector<output_vector_array_t>& outputTrajectoriesStock,
		const std::vector<lagrange_t>& lagrangeMultiplierFunctionsStock,
		std::vector<std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd> > >&  lagrangeTrajectoriesStock)  {

	typedef Eigen::Matrix<double, Eigen::Dynamic, 1> constraint_vector_t;
	typedef Eigen::Matrix<double, Eigen::Dynamic, OUTPUT_DIM> constraint_matrix_t;


	LinearInterpolation<constraint_vector_t, Eigen::aligned_allocator<constraint_vector_t> > vffFunc;
	LinearInterpolation<constraint_matrix_t, Eigen::aligned_allocator<constraint_matrix_t> > vfbFunc;

	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		vffFunc.setTimeStamp(&lagrangeMultiplierFunctionsStock[i].time_);
		vffFunc.setData(&lagrangeMultiplierFunctionsStock[i].uff_);

		vfbFunc.setTimeStamp(&lagrangeMultiplierFunctionsStock[i].time_);
		vfbFunc.setData(&lagrangeMultiplierFunctionsStock[i].k_);

		lagrangeTrajectoriesStock[i].resize(timeTrajectoriesStock[i].size());

		for (int k=0; k<timeTrajectoriesStock[i].size(); k++) {

			constraint_vector_t vff;
			vffFunc.interpolate(timeTrajectoriesStock[i][k], vff);
			size_t greatestLessTimeIndex = vffFunc.getGreatestLessTimeStampIndex();

			constraint_matrix_t vfb;
			vfbFunc.interpolate(timeTrajectoriesStock[i][k], vfb, greatestLessTimeIndex);

			lagrangeTrajectoriesStock[i][k] = vff + vfb*outputTrajectoriesStock[i][k];

		}  // end of k loop
	}  // end of i loop
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * compute the co-state over the given rollout
 * 		inputs:
 * 			+ timeTrajectoriesStock: rollout simulated time steps
 * 			+ outputTrajectoriesStock: rollout outputs
 *
 * 		outputs:
 * 			+ costateTrajectoriesStock: co-state vector for the given trajectory
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateRolloutCostate(
		const std::vector<scalar_array_t>& timeTrajectoriesStock,
		const std::vector<output_vector_array_t>& outputTrajectoriesStock,
		std::vector<output_vector_array_t>& costateTrajectoriesStock)  {


	LinearInterpolation<state_matrix_t,Eigen::aligned_allocator<state_matrix_t> >   SmFunc;
	LinearInterpolation<output_vector_t,Eigen::aligned_allocator<output_vector_t> > SvFunc;
	LinearInterpolation<output_vector_t,Eigen::aligned_allocator<output_vector_t> > nominalOutputFunc;

	costateTrajectoriesStock.resize(NUM_SUBSYSTEMS);

	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		SmFunc.setTimeStamp(&BASE::SsTimeTrajectoryStock_[i]);
		SmFunc.setData(&BASE::SmTrajectoryStock_[i]);
		SvFunc.setTimeStamp(&BASE::SsTimeTrajectoryStock_[i]);
		SvFunc.setData(&BASE::SvTrajectoryStock_[i]);
		nominalOutputFunc.setTimeStamp(&BASE::nominalTimeTrajectoriesStock_[i]);
		nominalOutputFunc.setData(&BASE::nominalOutputTrajectoriesStock_[i]);

		size_t N = timeTrajectoriesStock[i].size();
		costateTrajectoriesStock[i].resize(N);

		for (int k=0; k<N; k++) {

			const double& t = timeTrajectoriesStock[i][k];

			state_matrix_t Sm;
			SmFunc.interpolate(t, Sm);
			size_t greatestLessTimeStampIndex = SmFunc.getGreatestLessTimeStampIndex();
			output_vector_t Sv;
			SvFunc.interpolate(t, Sv, greatestLessTimeStampIndex);

			output_vector_t nominalOutput;
			nominalOutputFunc.interpolate(t, nominalOutput);

			costateTrajectoriesStock[i][k] = Sv + Sm*(outputTrajectoriesStock[i][k]-nominalOutput);

		}  // end of k loop
	}  // end of i loop
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * line search on the feedforwrd parts of the controller and lagrange multipliers.
 * Based on the option flag lineSearchByMeritFuntion_ it uses two different approaches for line search:
 * 		+ lineSearchByMeritFuntion_=TRUE: it uses the merit function to choose the best stepSize for the
 * 		feedforward elements of controller and lagrangeMultiplierFuntion
 * 		ineSearchByMeritFuntion_=FALSE: the constraint correction term is added by a user defined stepSize.
 * 		The line search uses the pure cost function for choosing the best stepSize.
 *
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::lineSearch(
		const std::vector<control_vector_array_t>& feedForwardConstraintInputStock,
		scalar_t& learningRateStar,
		scalar_t maxLearningRateStar/*=1.0*/)  {

	// display
	if (BASE::options_.dispayGSLQP_)  {
		// less-equal operator for eigen vectors
		auto eigenVectorLessEqual = [] (const control_vector_t& u1, const control_vector_t& u2){ return u1.norm() < u2.norm(); };

		std::vector<control_vector_t> maxDeltaUffStock(NUM_SUBSYSTEMS);
		std::vector<control_vector_t> maxDeltaUffeStock(NUM_SUBSYSTEMS);
		for (size_t i=0; i<NUM_SUBSYSTEMS; i++)  {
			maxDeltaUffStock[i]  = *std::max_element(nominalControllersStock_[i].deltaUff_.begin(), nominalControllersStock_[i].deltaUff_.template end(), eigenVectorLessEqual);
			maxDeltaUffeStock[i] = *std::max_element(feedForwardConstraintInputStock[i].begin(), feedForwardConstraintInputStock[i].end(), eigenVectorLessEqual);
		}
		control_vector_t maxDeltaUff  = *std::max_element(maxDeltaUffStock.begin(), maxDeltaUffStock.end(), eigenVectorLessEqual);
		control_vector_t maxDeltaUffe = *std::max_element(maxDeltaUffeStock.begin(), maxDeltaUffeStock.end(), eigenVectorLessEqual);

		std::cerr << "max delta_uff norm: " << maxDeltaUff.norm()  << std::endl;
		std::cerr << "max uff_error norm: " << maxDeltaUffe.norm() << std::endl;
	}


	for (int i=0; i<NUM_SUBSYSTEMS; i++)
		for (int k=0; k<nominalControllersStock_[i].time_.size(); k++)
			if (BASE::options_.lineSearchByMeritFuntion_==true)
				nominalControllersStock_[i].deltaUff_[k] += feedForwardConstraintInputStock[i][k];
			else{
				nominalControllersStock_[i].uff_[k] += BASE::options_.constraintStepSize_*feedForwardConstraintInputStock[i][k];
			}

	// perform one rollout while the input correction for the type-1 constraint is considered.
	rollout(BASE::initState_, nominalControllersStock_, BASE::nominalTimeTrajectoriesStock_,
			BASE::nominalStateTrajectoriesStock_, BASE::nominalInputTrajectoriesStock_, BASE::nominalOutputTrajectoriesStock_,
			BASE::nc1TrajectoriesStock_, BASE::EvTrajectoryStock_, BASE::nc2TrajectoriesStock_,
			BASE::HvTrajectoryStock_, BASE::nc2FinalStock_, BASE::HvFinalStock_);
	calculateCostFunction(BASE::nominalTimeTrajectoriesStock_, BASE::nominalOutputTrajectoriesStock_, BASE::nominalInputTrajectoriesStock_,
			BASE::nc2TrajectoriesStock_, BASE::HvTrajectoryStock_, BASE::nc2FinalStock_, BASE::HvFinalStock_,
			BASE::nominalTotalCost_);

	// calculate the merit function
	if (BASE::options_.lineSearchByMeritFuntion_==true) {
		// calculate the lagrange multiplier with learningRate zero
		calculateRolloutLagrangeMultiplier(BASE::nominalTimeTrajectoriesStock_, BASE::nominalOutputTrajectoriesStock_, BASE::lagrangeControllerStock_,
				BASE::nominalLagrangeTrajectoriesStock_);
		// calculate the merit function
		calculateMeritFunction(BASE::nominalTimeTrajectoriesStock_, BASE::nc1TrajectoriesStock_, BASE::EvTrajectoryStock_, BASE::nominalLagrangeTrajectoriesStock_, BASE::nominalTotalCost_,
				BASE::nominalTotalMerit_, BASE::nominalConstraint1ISE_);
	} else {
		BASE::nominalTotalMerit_ = BASE::nominalTotalCost_;
		calculateConstraintISE(BASE::nominalTimeTrajectoriesStock_, BASE::nc1TrajectoriesStock_, BASE::EvTrajectoryStock_, BASE::nominalConstraint1ISE_);
	}

	// display
	if (BASE::options_.dispayGSLQP_)  {
		std::cerr << "\t learningRate 0.0 \t cost: " << BASE::nominalTotalCost_ << " \t merit: " << BASE::nominalTotalMerit_ <<
				" \t constraint ISE: " << BASE::nominalConstraint1ISE_ << std::endl;
		if (std::accumulate(BASE::nc2FinalStock_.begin(), BASE::nc2FinalStock_.end(), 0) > 0) {
			std::cerr << "\t final constraint type-2:   ";
			for(size_t i=0; i<NUM_SUBSYSTEMS; i++) std::cerr << "[" << i  << "]: " << BASE::HvFinalStock_[i].head(BASE::nc2FinalStock_[i]).transpose() << ",  ";
			std::cerr << std::endl;
		}

	}
	scalar_t learningRate = maxLearningRateStar;
	const std::vector<controller_t> controllersStock = nominalControllersStock_;
	const std::vector<lagrange_t> lagrangeMultiplierFunctionsStock = BASE::lagrangeControllerStock_;

	// local search forward simulation's variables
	scalar_t lsTotalCost;
	scalar_t lsTotalMerit;
	scalar_t lsConstraint1ISE;
	std::vector<controller_t>           lsControllersStock(NUM_SUBSYSTEMS);
	std::vector<lagrange_t>             lsLagrangeControllersStock(NUM_SUBSYSTEMS);
	std::vector<scalar_array_t>         lsTimeTrajectoriesStock(NUM_SUBSYSTEMS);
	std::vector<state_vector_array_t>   lsStateTrajectoriesStock(NUM_SUBSYSTEMS);
	std::vector<control_vector_array_t> lsInputTrajectoriesStock(NUM_SUBSYSTEMS);
	std::vector<output_vector_array_t>  lsOutputTrajectoriesStock(NUM_SUBSYSTEMS);
	std::vector<std::vector<size_t> >   lsNc1TrajectoriesStock(NUM_SUBSYSTEMS);
	std::vector<constraint1_vector_array_t> lsEvTrajectoryStock(NUM_SUBSYSTEMS);
	std::vector<std::vector<size_t> >   lsNc2TrajectoriesStock(NUM_SUBSYSTEMS);
	std::vector<constraint2_vector_array_t> lsHvTrajectoryStock(NUM_SUBSYSTEMS);
	std::vector<size_t>               	lsNc2FinalStock(NUM_SUBSYSTEMS);
	std::vector<constraint2_vector_t> 	lsHvFinalStock(NUM_SUBSYSTEMS);
	std::vector<std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd> > >  lsLagrangeTrajectoriesStock(NUM_SUBSYSTEMS);

	while (learningRate >= BASE::options_.minLearningRateGSLQP_)  {
		// modifying uff by the local increments
		lsControllersStock = controllersStock;
		for (int i=0; i<NUM_SUBSYSTEMS; i++)
			for (int k=0; k<lsControllersStock[i].time_.size(); k++)
				lsControllersStock[i].uff_[k] += learningRate*lsControllersStock[i].deltaUff_[k];

		// modifying vff by the local increments
		lsLagrangeControllersStock = lagrangeMultiplierFunctionsStock;
		for (int i=0; i<NUM_SUBSYSTEMS; i++)
			for (int k=0; k<lsLagrangeControllersStock[i].time_.size(); k++)
				lsLagrangeControllersStock[i].uff_[k] += learningRate*lsLagrangeControllersStock[i].deltaUff_[k];

		// perform rollout
		try {
			rollout(BASE::initState_, lsControllersStock, lsTimeTrajectoriesStock, lsStateTrajectoriesStock, lsInputTrajectoriesStock, lsOutputTrajectoriesStock,
					lsNc1TrajectoriesStock, lsEvTrajectoryStock,  lsNc2TrajectoriesStock, lsHvTrajectoryStock, lsNc2FinalStock, lsHvFinalStock);
			// calculate rollout cost
			calculateCostFunction(lsTimeTrajectoriesStock, lsOutputTrajectoriesStock, lsInputTrajectoriesStock, lsNc2TrajectoriesStock, lsHvTrajectoryStock, lsNc2FinalStock, lsHvFinalStock, lsTotalCost);

			// calculate the merit function
			if (BASE::options_.lineSearchByMeritFuntion_==true) {
				// calculate the lagrange multiplier
				calculateRolloutLagrangeMultiplier(lsTimeTrajectoriesStock, lsOutputTrajectoriesStock, lsLagrangeControllersStock,
						lsLagrangeTrajectoriesStock);
				// calculate the merit function value
				calculateMeritFunction(lsTimeTrajectoriesStock, lsNc1TrajectoriesStock, lsEvTrajectoryStock, lsLagrangeTrajectoriesStock, lsTotalCost,
						lsTotalMerit, lsConstraint1ISE);
			} else {
				lsTotalMerit = lsTotalCost;
				calculateConstraintISE(lsTimeTrajectoriesStock, lsNc1TrajectoriesStock, lsEvTrajectoryStock, lsConstraint1ISE);
			}

			// display
			if (BASE::options_.dispayGSLQP_)  std::cerr << "\t learningRate " << learningRate << " \t cost: " << lsTotalCost << " \t merit: " << lsTotalMerit <<
					" \t constraint ISE: " << lsConstraint1ISE << std::endl;
		}
		catch(const std::exception& error)
		{
			std::cerr << "\t rollout with learningRate " << learningRate << " is terminated due to the slow simulation!" << std::endl;
			lsTotalCost  = std::numeric_limits<scalar_t>::max();
			lsTotalMerit = std::numeric_limits<scalar_t>::max();
		}

		// break condition 1: it exits with largest learningRate that its cost is smaller than nominal cost.
		if (lsTotalMerit < BASE::nominalTotalMerit_*(1-1e-3*learningRate))
			break;  // exit while loop
		else
			learningRate = 0.5*learningRate;

	}  // end of while


	if (learningRate >= BASE::options_.minLearningRateGSLQP_)  {
		BASE::nominalTotalCost_      = lsTotalCost;
		BASE::nominalTotalMerit_     = lsTotalMerit;
		BASE::nominalConstraint1ISE_ = lsConstraint1ISE;
		learningRateStar = learningRate;

		BASE::nc2FinalStock_.swap(lsNc2FinalStock);
		BASE::HvFinalStock_.swap(lsHvFinalStock);

		for (size_t i = 0; i<NUM_SUBSYSTEMS; i++)	// swapping where possible for efficiency
		{
			nominalControllersStock_[i].swap(lsControllersStock[i]);
			BASE::nominalTimeTrajectoriesStock_[i].swap(lsTimeTrajectoriesStock[i]);
			BASE::nominalStateTrajectoriesStock_[i].swap(lsStateTrajectoriesStock[i]);
			BASE::nominalInputTrajectoriesStock_[i].swap(lsInputTrajectoriesStock[i]);
			BASE::nominalOutputTrajectoriesStock_[i].swap(lsOutputTrajectoriesStock[i]);
			BASE::nc1TrajectoriesStock_[i].swap(lsNc1TrajectoriesStock[i]);
			BASE::EvTrajectoryStock_[i].swap(lsEvTrajectoryStock[i]);
			BASE::nc2TrajectoriesStock_[i].swap(lsNc2TrajectoriesStock[i]);
			BASE::HvTrajectoryStock_[i].swap(lsHvTrajectoryStock[i]);
			BASE::lagrangeControllerStock_[i].swap(lsLagrangeControllersStock[i]);;
			BASE::nominalLagrangeTrajectoriesStock_[i].swap(lsLagrangeTrajectoriesStock[i]);
		}

	} else // since the open loop input is not change, the nominal trajectories will be unchanged
		learningRateStar = 0.0;

	// clear the feedforward increments
	for (size_t i=0; i<NUM_SUBSYSTEMS; i++) {
		nominalControllersStock_[i].deltaUff_.clear();
		BASE::lagrangeControllerStock_[i].deltaUff_.clear();
	}

	// display
	if (BASE::options_.dispayGSLQP_)  std::cerr << "The chosen learningRate is: " << learningRateStar << std::endl;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * get the calculated optimal controller structure
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getController(std::vector<controller_t>& controllersStock) {

	controllersStock = nominalControllersStock_;
}


template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::setController(const std::vector<controller_t>& controllersStock){
	nominalControllersStock_ = controllersStock;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * calculate the value function for the given time and output vector
 * 		inputs
 * 			+ time: inquiry time
 * 			+ output: inquiry output
 *
 * 		output:
 * 			+ valueFuntion: value function at the inquiry time and output
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getValueFuntion(const scalar_t& time, const output_vector_t& output, scalar_t& valueFuntion)  {

	int activeSubsystem = -1;
	for (int i=0; i<NUM_SUBSYSTEMS; i++)  {
		activeSubsystem = i;
		if (BASE::switchingTimes_[i]<=time && time<BASE::switchingTimes_[i+1])
			break;
	}

	state_matrix_t Sm;
	LinearInterpolation<state_matrix_t,Eigen::aligned_allocator<state_matrix_t> > SmFunc(
			&BASE::SsTimeTrajectoryStock_[activeSubsystem], &BASE::SmTrajectoryStock_[activeSubsystem]);
	SmFunc.interpolate(time, Sm);
	size_t greatestLessTimeStampIndex = SmFunc.getGreatestLessTimeStampIndex();

	output_vector_t Sv;
	LinearInterpolation<output_vector_t,Eigen::aligned_allocator<output_vector_t> > SvFunc(
			&BASE::SsTimeTrajectoryStock_[activeSubsystem], &BASE::SvTrajectoryStock_[activeSubsystem]);
	SvFunc.interpolate(time, Sv, greatestLessTimeStampIndex);

	eigen_scalar_t s;
	LinearInterpolation<eigen_scalar_t,Eigen::aligned_allocator<eigen_scalar_t> > sFunc(
			&BASE::SsTimeTrajectoryStock_[activeSubsystem], &BASE::sTrajectoryStock_[activeSubsystem]);
	sFunc.interpolate(time, s, greatestLessTimeStampIndex);

	output_vector_t xNominal;
	LinearInterpolation<output_vector_t,Eigen::aligned_allocator<output_vector_t> > xNominalFunc(&BASE::nominalTimeTrajectoriesStock_[activeSubsystem], &BASE::nominalOutputTrajectoriesStock_[activeSubsystem]);
	xNominalFunc.interpolate(time, xNominal);

	output_vector_t deltaX = output-xNominal;
	valueFuntion = (s + deltaX.transpose()*Sv + 0.5*deltaX.transpose()*Sm*deltaX).eval()(0);
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * calculate the cost function at the initial time
 * 		inputs
 * 			+ initOutput: initial output
 *
 * 		output:
 * 			+ cost function value
 * 			+ cost function value plus the constraint ISE multiplied by pho
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getCostFuntion(scalar_t& costFunction, scalar_t& constriantISE) {
	costFunction = BASE::nominalTotalCost_;
	constriantISE = BASE::nominalConstraint1ISE_;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * get the optimal state, output, and input trajectories
 * 		output
 * 			+ BASE::nominalTimeTrajectoriesStock_: optimal time trajectory
 * 			+ BASE::nominalStateTrajectoriesStock_: optimal state trajectory
 * 			+ BASE::nominalInputTrajectoriesStock_: optimal control input trajectory
 * 			+ BASE::nominalOutputTrajectoriesStock_: optimal output trajectory
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getNominalTrajectories(std::vector<scalar_array_t>& nominalTimeTrajectoriesStock,
		std::vector<state_vector_array_t>& nominalStateTrajectoriesStock,
		std::vector<control_vector_array_t>& nominalInputTrajectoriesStock,
		std::vector<output_vector_array_t>& nominalOutputTrajectoriesStock)   {

	nominalTimeTrajectoriesStock   = BASE::nominalTimeTrajectoriesStock_;
	nominalStateTrajectoriesStock  = BASE::nominalStateTrajectoriesStock_;
	nominalInputTrajectoriesStock  = BASE::nominalInputTrajectoriesStock_;
	nominalOutputTrajectoriesStock = BASE::nominalOutputTrajectoriesStock_;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * solve the SLQP Riccati differential equations:
 * 		input:
 * 			+ learningRate: the feeadforward learningRate
 *
 * 		uses:
 * 			+ linearized dynamics
 * 			+ quadratized cost
 *
 * 		modifies:
 * 			V(t,y) = y^T*Sm*y + y^T*(Sv+Sve) + s
 * 			+ BASE::SsTimeTrajectoryStock_: time stamp
 * 			+ BASE::SmTrajectoryStock_: Sm matrix
 * 			+ BASE::SvTrajectoryStock_: Sv vector
 * 			+ BASE::SveTrajectoryStock_: Sve vector
 * 			+ BASE::sTrajectoryStock_: s scalar
 */
//template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
//void SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::solveSequentialRiccatiEquations(const scalar_t& learningRate)  {
//
//	LinearInterpolation<state_matrix_t, Eigen::aligned_allocator<state_matrix_t> > SmFunc;
//
//	// final value for the last Riccati equations
//	typename RiccatiEquations_t::s_vector_t allSsFinal;
//	allSsFinal.setZero();
//	// final value for the last error equation
//	output_vector_t SveFinal = output_vector_t::Zero();
//
//	for (int i=NUM_SUBSYSTEMS-1; i>=0; i--) {
//
//		// final cost of the subsystem is added to the following subsystem solution
//		typename RiccatiEquations_t::s_vector_t allCostFinal;
//		RiccatiEquations_t::convert2Vector(BASE::QmFinalStock_[i], BASE::QvFinalStock_[i], BASE::qFinalStock_[i], allCostFinal);
//		allSsFinal += allCostFinal;
//
//		// set data for Riccati equations
//		std::shared_ptr<RiccatiEquations_t> riccatiEquationsPtr( new RiccatiEquations_t() );
//		riccatiEquationsPtr->setData(learningRate, i, BASE::switchingTimes_[i], BASE::switchingTimes_[i+1],
//				&BASE::nominalTimeTrajectoriesStock_[i],
//				&BASE::AmConstrainedTrajectoryStock_[i], &BASE::BmTrajectoryStock_[i],
//				&BASE::qTrajectoryStock_[i], &BASE::QvConstrainedTrajectoryStock_[i], &BASE::QmConstrainedTrajectoryStock_[i],
//				&BASE::RvTrajectoryStock_[i], &BASE::RmInverseTrajectoryStock_[i], &BASE::RmConstrainedTrajectoryStock_[i], &BASE::PmTrajectoryStock_[i]);
//
//		// max number of steps of integration
//		size_t maxNumSteps = BASE::options_.maxNumStepsPerSecond_ * std::max( 1.0, BASE::switchingTimes_[i+1]-BASE::switchingTimes_[i] );
//		// integrating the Riccati equations
//		ODE45<RiccatiEquations_t::S_DIM_> ode45(riccatiEquationsPtr);
//		std::vector<double> normalizedTimeTrajectory;
//		std::vector<typename RiccatiEquations_t::s_vector_t, Eigen::aligned_allocator<typename RiccatiEquations_t::s_vector_t> > allSsTrajectory;
//		ode45.integrate(allSsFinal, i, i+1, allSsTrajectory, normalizedTimeTrajectory,
//				1e-5, BASE::options_.AbsTolODE_, BASE::options_.RelTolODE_,  maxNumSteps);
//
//		// denormalizing time and constructing 'Sm', 'Sv', and 's'
//		int N = normalizedTimeTrajectory.size();
//		BASE::SsTimeTrajectoryStock_[i].resize(N);
//		BASE::SmTrajectoryStock_[i].resize(N);
//		BASE::SvTrajectoryStock_[i].resize(N);
//		BASE::sTrajectoryStock_[i].resize(N);
//		for (int k=0; k<N; k++) {
//			RiccatiEquations_t::convert2Matrix(allSsTrajectory[N-1-k], BASE::SmTrajectoryStock_[i][k], BASE::SvTrajectoryStock_[i][k], BASE::sTrajectoryStock_[i][k]);
//			BASE::SsTimeTrajectoryStock_[i][k] = (BASE::switchingTimes_[i]-BASE::switchingTimes_[i+1])*(normalizedTimeTrajectory[N-1-k]-i) + BASE::switchingTimes_[i+1];
//		}  // end of k loop
//
//		// testing the numerical stability of the Riccati equations
//		for (int k=N-1; k>=0; k--) {
//			try {
//				if (BASE::SmTrajectoryStock_[i][k] != BASE::SmTrajectoryStock_[i][k])  throw std::runtime_error("Sm is unstable.");
//				if (BASE::SvTrajectoryStock_[i][k] != BASE::SvTrajectoryStock_[i][k])  throw std::runtime_error("Sv is unstable.");
//				if (BASE::sTrajectoryStock_[i][k] != BASE::sTrajectoryStock_[i][k])    throw std::runtime_error("s is unstable.");
//			}
//			catch(const std::exception& error)
//			{
//				std::cerr << "what(): " << error.what() << " at time " << BASE::SsTimeTrajectoryStock_[i][k] << " [sec]." << std::endl;
//				for (int kp=k; kp<k+10; kp++)  {
//					if (kp >= N) continue;
//					std::cerr << "Sm[" << BASE::SsTimeTrajectoryStock_[i][kp] << "]:\n"<< BASE::SmTrajectoryStock_[i][kp].transpose() << std::endl;
//					std::cerr << "Sv[" << BASE::SsTimeTrajectoryStock_[i][kp] << "]:\t"<< BASE::SvTrajectoryStock_[i][kp].transpose() << std::endl;
//					std::cerr << "s[" << BASE::SsTimeTrajectoryStock_[i][kp] << "]: \t"<< BASE::sTrajectoryStock_[i][kp].transpose() << std::endl;
//				}
//				exit(1);
//			}
//
//		}  // end of k loop
//
//		// set the final value for next Riccati equation
//		allSsFinal = allSsTrajectory.back();
//
//		/*
//		 * Type_1 constraints error correction compensation
//		 */
//
//		// Skip calculation of the error correction term Sve if the constrained simulation is used for forward simulation
//		if (BASE::options_.simulationIsConstrained_) {
//			BASE::SveTrajectoryStock_[i].resize(N);
//			for (int k=0; k<N; k++)
//				BASE::SveTrajectoryStock_[i][k].setZero();
//			continue;
//		}
//
//		// Calculating the coefficients of the error equation
//		SmFunc.setTimeStamp( &(BASE::SsTimeTrajectoryStock_[i]) );
//		SmFunc.setData( &(BASE::SmTrajectoryStock_[i]) );
//
//		output_vector_array_t GvTrajectory(BASE::nominalTimeTrajectoriesStock_[i].size());
//		state_matrix_array_t  GmTrajectory(BASE::nominalTimeTrajectoriesStock_[i].size());
//
//		for (int k=0; k<BASE::nominalTimeTrajectoriesStock_[i].size(); k++) {
//			state_matrix_t Sm;
//			SmFunc.interpolate(BASE::nominalTimeTrajectoriesStock_[i][k], Sm);
//
//			control_feedback_t Lm = BASE::RmInverseTrajectoryStock_[i][k]*(BASE::PmTrajectoryStock_[i][k]+BASE::BmTrajectoryStock_[i][k].transpose()*Sm);
//
//			GmTrajectory[k] = BASE::AmConstrainedTrajectoryStock_[i][k] -
//					BASE::BmTrajectoryStock_[i][k]*BASE::RmInverseTrajectoryStock_[i][k]*BASE::RmConstrainedTrajectoryStock_[i][k]*Lm;
//
//			GvTrajectory[k] = (BASE::CmProjectedTrajectoryStock_[i][k]-Lm).transpose()*
//					BASE::RmTrajectoryStock_[i][k]*BASE::EvProjectedTrajectoryStock_[i][k];
//
//		}  // end of k loop
//
//		// set data for error equations
//		std::shared_ptr<ErrorEquation_t> errorEquationPtr( new ErrorEquation_t() );
//		errorEquationPtr->setData(i, BASE::switchingTimes_[i], BASE::switchingTimes_[i+1],
//				&BASE::nominalTimeTrajectoriesStock_[i], &GvTrajectory, &GmTrajectory);
//
//		// integrating the Riccati equations
//		ODE45<OUTPUT_DIM> errorOde45(errorEquationPtr);
//		output_vector_array_t SveTrajectory;
//		errorOde45.integrate(SveFinal, normalizedTimeTrajectory, SveTrajectory, 1e-3, BASE::options_.AbsTolODE_, BASE::options_.RelTolODE_);
//
//		// reset the final value for next Riccati equation
//		SveFinal = SveTrajectory.back();
//
//		BASE::SveTrajectoryStock_[i].resize(N);
//		for (int k=0; k<N; k++) {
//			BASE::SveTrajectoryStock_[i][k] = SveTrajectory[N-1-k];
//
//			// testing the numerical stability of the Riccati error equation
//			try {
//				if (BASE::SveTrajectoryStock_[i][k] != BASE::SveTrajectoryStock_[i][k])  throw std::runtime_error("Sve is unstable");
//			}
//			catch(const std::exception& error) 	{
//				std::cerr << "what(): " << error.what() << " at time " << BASE::SsTimeTrajectoryStock_[i][k] << " [sec]." << std::endl;
//				for (int kp=k; kp<N; kp++)   std::cerr << "Sve[" << BASE::SsTimeTrajectoryStock_[i][kp] << "]:\t"<< BASE::SveTrajectoryStock_[i][kp].transpose() << std::endl;
//				exit(1);
//			}
//		}
//
//	}  // end of i loop
//}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * run the SLQP algorithm for a given state and switching times
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::run(const state_vector_t& initState,
		const std::vector<scalar_t>& switchingTimes, const std::vector<controller_t>& initialControllersStock)  {

	BASE::iterationCost_.clear();
	BASE::iterationISE1_.clear();

	if (switchingTimes.size() != NUM_SUBSYSTEMS+1)
		throw std::runtime_error("Number of switching times should be one plus the number of subsystems.");

	if (initialControllersStock.empty()==false) {
		if (initialControllersStock.size() != NUM_SUBSYSTEMS)  throw std::runtime_error("initialControllersStock has less controllers than the number of subsystems");
		nominalControllersStock_ = initialControllersStock;
	} else
		if (nominalControllersStock_.empty()==true)
			throw std::runtime_error("initialControllersStock should be provided since nominalControllersStock is empty.");


	BASE::switchingTimes_ = switchingTimes;
	BASE::initState_ = initState;

	// display
	if (BASE::options_.dispayGSLQP_)
	{
		std::cerr << "\n#### SLQP solver starts with switching times [" << switchingTimes[0];
		for (size_t i=1; i<=NUM_SUBSYSTEMS; i++)   std::cerr << ", " << switchingTimes[i];
		std::cerr << "] ..." << std::endl << std::endl;
	}

	BASE::iteration_ = 0;
	double relCost;
	double learningRateStar;
	double relConstraint1ISE;
	bool isConstraint1Satisfied  = false;
	bool isCostFunctionConverged = false;
	bool isOptimizationConverged = false;
	bool nominalLagrangeMultiplierUpdated = false;

	// initial controller rollout
	rollout(BASE::initState_, nominalControllersStock_,
			BASE::nominalTimeTrajectoriesStock_, BASE::nominalStateTrajectoriesStock_, BASE::nominalInputTrajectoriesStock_, BASE::nominalOutputTrajectoriesStock_,
			BASE::nc1TrajectoriesStock_, BASE::EvTrajectoryStock_, BASE::nc2TrajectoriesStock_, BASE::HvTrajectoryStock_, BASE::nc2FinalStock_, BASE::HvFinalStock_);

	// initial controller cost
	calculateCostFunction(BASE::nominalTimeTrajectoriesStock_, BASE::nominalOutputTrajectoriesStock_, BASE::nominalInputTrajectoriesStock_,
			BASE::nc2TrajectoriesStock_, BASE::HvTrajectoryStock_, BASE::nc2FinalStock_, BASE::HvFinalStock_,
			BASE::nominalTotalCost_);

	// initial controller merit
	BASE::nominalTotalMerit_ = BASE::nominalTotalCost_;

	// initial controller constraint type-1 ISE
	calculateConstraintISE(BASE::nominalTimeTrajectoriesStock_, BASE::nc1TrajectoriesStock_, BASE::EvTrajectoryStock_, BASE::nominalConstraint1ISE_);

	// display
	if (BASE::options_.dispayGSLQP_)  std::cerr << "\n#### Initial controller: \n cost: " << BASE::nominalTotalCost_ << " \t constraint ISE: " << BASE::nominalConstraint1ISE_ << std::endl;

	BASE::iterationCost_.push_back( (Eigen::VectorXd(1) << BASE::nominalTotalCost_).finished() );
	BASE::iterationISE1_.push_back( (Eigen::VectorXd(1) << BASE::nominalConstraint1ISE_).finished() );

	// SLQP main loop
	while (BASE::iteration_<BASE::options_.maxIterationGSLQP_ && isOptimizationConverged==false)  {

		double costCashed = BASE::nominalTotalCost_;
		double constraint1ISECashed = BASE::nominalConstraint1ISE_;

		// display
		if (BASE::options_.dispayGSLQP_)  std::cerr << "\n#### Iteration " <<  BASE::iteration_+1 << std::endl;

		// linearizing the dynamics and quadratizing the cost function along nominal trajectories
		approximateOptimalControlProblem();

		// solve Riccati equations
		this->solveSequentialRiccatiEquations(1.0 /*nominal learningRate*/);

		// calculate controller and lagrange multiplier
		std::vector<control_vector_array_t> feedForwardConstraintInputStock(NUM_SUBSYSTEMS);
		calculateControllerAndLagrangian(nominalControllersStock_, BASE::lagrangeControllerStock_,
				feedForwardConstraintInputStock, ~nominalLagrangeMultiplierUpdated);

		nominalLagrangeMultiplierUpdated = true;

		// finding the optimal learningRate
		lineSearch(feedForwardConstraintInputStock, learningRateStar, BASE::options_.maxLearningRateGSLQP_);

		// calculates type-1 constraint ISE and maximum norm
		double constraint1MaxNorm = calculateConstraintISE(BASE::nominalTimeTrajectoriesStock_, BASE::nc1TrajectoriesStock_, BASE::EvTrajectoryStock_, BASE::nominalConstraint1ISE_);

		// calculates type-2 constraint ISE and maximum norm
		double nominalConstraint2ISE, constraint2MaxNorm;
		constraint2MaxNorm = calculateConstraintISE(BASE::nominalTimeTrajectoriesStock_, BASE::nc2TrajectoriesStock_, BASE::HvTrajectoryStock_, nominalConstraint2ISE);

		// loop variables
		BASE::iteration_++;
		relCost = fabs(BASE::nominalTotalCost_-costCashed);
		relConstraint1ISE = fabs(BASE::nominalConstraint1ISE_-constraint1ISECashed);
		isConstraint1Satisfied  = BASE::nominalConstraint1ISE_<=BASE::options_.minAbsConstraint1ISE_ || relConstraint1ISE<=BASE::options_.minRelConstraint1ISE_;
		isCostFunctionConverged = learningRateStar==0 || relCost<=BASE::options_.minRelCostGSLQP_;
		isOptimizationConverged = isCostFunctionConverged==true && isConstraint1Satisfied==true;

		// display
		if (BASE::options_.dispayGSLQP_)  {
			std::cerr << "optimization cost:         " << BASE::nominalTotalCost_ << std::endl;
			std::cerr << "constraint type-1 ISE:     " << BASE::nominalConstraint1ISE_ << std::endl;
			std::cerr << "constraint type-1 MaxNorm: " << constraint1MaxNorm << std::endl;
			std::cerr << "constraint type-2 ISE:     " << nominalConstraint2ISE << std::endl;
			std::cerr << "constraint type-2 MaxNorm: " << constraint2MaxNorm << std::endl;
			//			if (std::accumulate(BASE::nc2FinalStock_.begin(), BASE::nc2FinalStock_.end(), 0) > 0) {
			std::cerr << "final constraint type-2: 	 ";
			for(size_t i=0; i<NUM_SUBSYSTEMS; i++) std::cerr << "[" << i  << "]: " << BASE::HvFinalStock_[i].head(BASE::nc2FinalStock_[i]).transpose() << ",  ";
			std::cerr << std::endl;
			//			}
		}
		if(BASE::options_.displayShortSummary_){
			std::cout << "#### Iter " << BASE::iteration_-1 << ".   opt. cost: " << BASE::nominalTotalCost_ << ".    constraint ISE: " << BASE::nominalConstraint1ISE_ <<
					".   constr MaxNorm: " << constraint1MaxNorm << std::endl;
		}
	}  // end of while loop


	// linearizing the dynamics and quadratizing the cost function along nominal trajectories
	approximateOptimalControlProblem();

	// solve Riccati equations with learningRate zero
	this->solveSequentialRiccatiEquations(0.0 /*learningRate*/);

	// calculate the nominal co-state
	calculateRolloutCostate(BASE::nominalTimeTrajectoriesStock_, BASE::nominalOutputTrajectoriesStock_,
			BASE::nominalcostateTrajectoriesStock_);


	// display
	if (BASE::options_.dispayGSLQP_ )  {
		std::cerr << "\n+++++++++++++++++++++++++++++++++++" << std::endl;
		std::cerr <<   "+++++++ SLQP solver is ended ++++++" << std::endl;
		std::cerr <<   "+++++++++++++++++++++++++++++++++++" << std::endl;
		if (isOptimizationConverged) {
			if (learningRateStar==0)
				std::cerr << "SLQP successfully terminates as learningRate reduced to zero." << std::endl;
			else
				std::cerr << "SLQP successfully terminates as cost relative change (relCost=" << relCost <<") reached to the minimum value." << std::endl;

			if (BASE::nominalConstraint1ISE_<=BASE::options_.minAbsConstraint1ISE_)
				std::cerr << "Type-1 constraint absolute ISE (absConstraint1ISE=" << BASE::nominalConstraint1ISE_ << ") reached to the minimum value." << std::endl;
			else
				std::cerr << "Type-1 constraint relative ISE (relConstraint1ISE=" << relConstraint1ISE << ") reached to the minimum value." << std::endl;
		} else
			std::cerr << "Maximum number of iterations has reached." << std::endl;
	}

}

} // namespace ocs2

