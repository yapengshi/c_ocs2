/*
 * Implementation of SLQP_MP.h
 *
 *  Created on: Jun 20, 2016
 *      Author: markus, farbod
 */


namespace ocs2{


template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::~SLQP_MP()
{
	workersActive_ = false;
	workerTask_ = SHUTDOWN;

	workerWakeUpCondition_.notify_all();

	if(mp_options_.debugPrintMP_)
	{
		printString("Shutting down workers");
	}

	for (size_t i=0; i<workerThreads_.size(); i++)
	{
		workerThreads_[i].join();
	}

	if(mp_options_.debugPrintMP_)
	{
		printString("All workers shut down");
	}
}


template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::rollout(
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
		std::vector<constraint2_vector_t>& HvFinalStock){

	rollout(mp_options_.nThreads_, initState,
			controllersStock, timeTrajectoriesStock, stateTrajectoriesStock,
			inputTrajectoriesStock, outputTrajectoriesStock,
			nc1TrajectoriesStock, EvTrajectoryStock,
			nc2TrajectoriesStock, HvTrajectoryStock, nc2FinalStock, HvFinalStock);
}

template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::rollout(
		const state_vector_t& initState,
		const std::vector<controller_t>& controllersStock,
		std::vector<scalar_array_t>& timeTrajectoriesStock,
		std::vector<state_vector_array_t>& stateTrajectoriesStock,
		std::vector<control_vector_array_t>& inputTrajectoriesStock,
		std::vector<output_vector_array_t>& outputTrajectoriesStock){

	rollout(mp_options_.nThreads_, initState,
			controllersStock, timeTrajectoriesStock, stateTrajectoriesStock,
			inputTrajectoriesStock, outputTrajectoriesStock);
}

template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::rollout(
		const state_vector_t& initState,
		const std::vector<controller_t>& controllersStock,
		std::vector<scalar_array_t>& timeTrajectoriesStock,
		std::vector<state_vector_array_t>& stateTrajectoriesStock,
		std::vector<control_vector_array_t>& inputTrajectoriesStock){

	rollout(mp_options_.nThreads_, initState,
			controllersStock, timeTrajectoriesStock, stateTrajectoriesStock, inputTrajectoriesStock);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * Forward integrate the system dynamics with given controller:
 * 		inputs:
 * 			+ threadId: 		identifier of chosen thread
 * 			+ initState: 		initial state at time BASE::switchingTimes_[0]
 * 			+ controller_local: controller for each subsystem
 * 		outputs:
 * 			+ t_local: 			rollout simulated time steps
 * 			+ x_local: 			rollout states
 * 			+ u_local: 			rollout control inputs
 * 			+ (optional) y_local: rollout outputs
 * 			+ (optional) nc1TrajectoriesStock: 	number of active constraints at each time step
 * 			+ (optional) EvTrajectoryStock: value of the constraint (if the rollout is constrained the value is
 * 											always zero otherwise it is nonzero)
 */

/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::rollout(
		const size_t& threadId,
		const state_vector_t& initState,
		const std::vector<controller_t>& controller_local,
		std::vector<scalar_array_t>& t_local,
		std::vector<state_vector_array_t>& x_local,
		std::vector<control_vector_array_t>& u_local,
		std::vector<output_vector_array_t>& y_local){

	if (controller_local.size() != NUM_SUBSYSTEMS)
		throw std::runtime_error("controller_local has less controllers then the number of subsystems");

	t_local.resize(NUM_SUBSYSTEMS);
	x_local.resize(NUM_SUBSYSTEMS);
	u_local.resize(NUM_SUBSYSTEMS);
	y_local.resize(NUM_SUBSYSTEMS);

	state_vector_t x0 = initState;
	for (int i=0; i<NUM_SUBSYSTEMS; i++)
	{
		t_local[i].clear();
		x_local[i].clear();

		// max number of steps of integration
		size_t maxNumSteps = BASE::options_.maxNumStepsPerSecond_ * std::max( 1.0, BASE::switchingTimes_[i+1]-BASE::switchingTimes_[i]);

		// initialize subsystem i
		dynamics_[threadId][i]->initializeModel(BASE::switchingTimes_, x0, i, "GSLPQ");

		// set controller for subsystem i
		dynamics_[threadId][i]->setController(controller_local[i]);

		// simulate subsystem i
		integratorsODE45_[threadId][i]->integrate(
				x0, BASE::switchingTimes_[i], BASE::switchingTimes_[i+1],
				x_local[i], t_local[i],
				1e-3, BASE::options_.AbsTolODE_, BASE::options_.RelTolODE_, maxNumSteps);

		if (x_local[i].back() != x_local[i].back())
			throw std::runtime_error("System became unstable during the SLQP_MP rollout.");


		// compute control trajectory for subsystem i
		u_local[i].resize(t_local[i].size());
		y_local[i].resize(t_local[i].size());
		for (int k=0; k<t_local[i].size(); k++)
		{
			dynamics_[threadId][i]->computeOutput(t_local[i][k], x_local[i][k], y_local[i][k]);
			dynamics_[threadId][i]->computeInput(t_local[i][k], y_local[i][k], u_local[i][k]);
		}

		// reset the initial state
		x0 = x_local[i].back();
	}
}



template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::rollout(
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
		dynamics_[mp_options_.nThreads_][i]->initializeModel(BASE::switchingTimes_, x0, i, "GSLPQ");

		// set controller for subsystem i
		dynamics_[mp_options_.nThreads_][i]->setController(controllersStock[i]);

		// determine correct stopping time for this subsystem
		double t_stop_thisSubsystem;
		bool stopAfterThisSubsystem = false;
		if(stoppingTime > BASE::switchingTimes_[i+1])
			t_stop_thisSubsystem = BASE::switchingTimes_[i+1];
		else{
			t_stop_thisSubsystem = stoppingTime;
			stopAfterThisSubsystem = true;
		}

		// integrate subsystem i
		integratorsODE45_[mp_options_.nThreads_][i]->integrate(
				x0, BASE::switchingTimes_[i], t_stop_thisSubsystem,
				stateTrajectoriesStock[i], timeTrajectoriesStock[i],
				1e-3, BASE::options_.AbsTolODE_, BASE::options_.RelTolODE_, maxNumSteps);

		if (stateTrajectoriesStock[i].back() != stateTrajectoriesStock[i].back())
			throw std::runtime_error("System became unstable during the SLQP rollout.");

		if(stopAfterThisSubsystem == true){
			numSubsystemWhereStopped = i;
			stateVectorWhereStopped = stateTrajectoriesStock[i].back();
			dynamics_[mp_options_.nThreads_][i]->computeOutput(timeTrajectoriesStock[i].back(), stateVectorWhereStopped, outputWhereStopped);
			dynamics_[mp_options_.nThreads_][i]->computeInput(timeTrajectoriesStock[i].back(), outputWhereStopped, controlInputWhereStopped);
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
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::rollout(
		const size_t& threadId,
		const state_vector_t& initState,
		const std::vector<controller_t>& controller_local,
		std::vector<scalar_array_t>& t_local,
		std::vector<state_vector_array_t>& x_local,
		std::vector<control_vector_array_t>& u_local)  {

	std::vector<output_vector_array_t> y_local;
	rollout(threadId, initState, controller_local, t_local, x_local, u_local, y_local);
}


/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::rollout(
		const size_t& threadId,
		const state_vector_t& initState,
		const std::vector<controller_t>& controller_local,
		std::vector<scalar_array_t>& t_local,
		std::vector<state_vector_array_t>& x_local,
		std::vector<control_vector_array_t>& u_local,
		std::vector<output_vector_array_t>& y_local,
		std::vector<std::vector<size_t> >& nc1TrajectoriesStock,
		std::vector<constraint1_vector_array_t>& EvTrajectoryStock,
		std::vector<std::vector<size_t> >& nc2TrajectoriesStock,
		std::vector<constraint2_vector_array_t>& HvTrajectoryStock,
		std::vector<size_t>& nc2FinalStock,
		std::vector<constraint2_vector_t>& HvFinalStock) {

	// STEP1 : do a rollout
	rollout(threadId, initState, controller_local, t_local, x_local, u_local, y_local);


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
		size_t N = t_local[i].size();
		nc1TrajectoriesStock[i].resize(N);
		EvTrajectoryStock[i].resize(N);
		nc2TrajectoriesStock[i].resize(N);
		HvTrajectoryStock[i].resize(N);


		// compute constraint1 trajectory for subsystem i
		for (int k=0; k<N; k++)
		{
			// constraint 1 type
			dynamics_[threadId][i]->computeConstriant1(
					t_local[i][k], x_local[i][k], u_local[i][k],
					nc1TrajectoriesStock[i][k], EvTrajectoryStock[i][k]);

			if (nc1TrajectoriesStock[i][k] > INPUT_DIM)
				throw std::runtime_error("Number of active type-1 constraints should be less-equal to the number of input dimension.");


			// constraint type 2
			dynamics_[threadId][i]->computeConstriant2(t_local[i][k],
					x_local[i][k],
					nc2TrajectoriesStock[i][k], HvTrajectoryStock[i][k]);

			if (nc2TrajectoriesStock[i][k] > INPUT_DIM)
				throw std::runtime_error("Number of active type-2 constraints should be less-equal to the number of input dimension.");

		}  // end of k loop

		dynamics_[threadId][i]->computeFinalConstriant2(t_local[i].back(), x_local[i].back(),
				nc2FinalStock[i], HvFinalStock[i]);
		if (nc2FinalStock[i] > INPUT_DIM)
			throw std::runtime_error("Number of active type-2 constraints at final time should be less-equal to the number of input dimension.");

	}  // end of i loop

}


template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateCostFunction(
		const std::vector<scalar_array_t>& timeTrajectoriesStock,
		const std::vector<output_vector_array_t>& stateTrajectoriesStock,
		const std::vector<control_vector_array_t>& inputTrajectoriesStock,
		scalar_t& totalCost) {

	calculateCostFunction(timeTrajectoriesStock, stateTrajectoriesStock, inputTrajectoriesStock, totalCost, mp_options_.nThreads_);
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
 *			+ threadId: working thread, defaults to the thread with lowest id, thus this is the default thread for single-core cost computation
 *				(allows to let method be called from the outside)
 * 		outputs:
 * 			+ totalCost: the total cost of the trajectory
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateCostFunction(
		const std::vector<scalar_array_t>& timeTrajectoriesStock,
		const std::vector<output_vector_array_t>& outputTrajectoriesStock,
		const std::vector<control_vector_array_t>& inputTrajectoriesStock,
		scalar_t& totalCost,
		size_t threadId)  {

	totalCost = 0.0;
	for (size_t i=0; i<NUM_SUBSYSTEMS; i++)
	{

		// integrates the intermediate cost using the trapezoidal approximation method
		scalar_t currentIntermediateCost;
		scalar_t nextIntermediateCost;
		for (int k=0; k<timeTrajectoriesStock[i].size()-1; k++) {

			if (k==0) {
				costFunctions_[threadId][i]->setCurrentStateAndControl(timeTrajectoriesStock[i][k], outputTrajectoriesStock[i][k], inputTrajectoriesStock[i][k]);
				costFunctions_[threadId][i]->evaluate(currentIntermediateCost);
			} else
			{
				currentIntermediateCost = nextIntermediateCost;
			}

			// feed next state and control to cost function
			costFunctions_[threadId][i]->setCurrentStateAndControl(timeTrajectoriesStock[i][k+1], outputTrajectoriesStock[i][k+1], inputTrajectoriesStock[i][k+1]);
			// evaluate intermediate cost for next time step
			costFunctions_[threadId][i]->evaluate(nextIntermediateCost);

			totalCost += 0.5*(currentIntermediateCost+nextIntermediateCost)*(timeTrajectoriesStock[i][k+1]-timeTrajectoriesStock[i][k]);
		}  // end of k loop

		// terminal cost
		if (i==NUM_SUBSYSTEMS-1)
		{
			scalar_t finalCost;
			costFunctions_[threadId][i]->setCurrentStateAndControl(timeTrajectoriesStock[i].back(), outputTrajectoriesStock[i].back(), inputTrajectoriesStock[i].back());
			costFunctions_[threadId][i]->terminalCost(finalCost);
			totalCost += finalCost;
		}

	}  // end of i loop
}


template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateCostFunction(const std::vector<scalar_array_t>& timeTrajectoriesStock,
		const std::vector<output_vector_array_t>& stateTrajectoriesStock,
		const std::vector<control_vector_array_t>& inputTrajectoriesStock,
		const std::vector<std::vector<size_t> >& nc2TrajectoriesStock,
		const std::vector<constraint2_vector_array_t>& HvTrajectoryStock,
		const std::vector<size_t>& nc2FinalStock,
		const std::vector<constraint2_vector_t>& HvFinalStock,
		scalar_t& totalCost){

	calculateCostFunction(timeTrajectoriesStock,
			stateTrajectoriesStock,
			inputTrajectoriesStock,
			nc2TrajectoriesStock,
			HvTrajectoryStock,
			nc2FinalStock,
			HvFinalStock,
			totalCost,
			mp_options_.nThreads_);
}

/*****************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateCostFunction(
		const std::vector<scalar_array_t>& timeTrajectoriesStock,
		const std::vector<output_vector_array_t>& outputTrajectoriesStock,
		const std::vector<control_vector_array_t>& inputTrajectoriesStock,
		const std::vector<std::vector<size_t> >& nc2TrajectoriesStock,
		const std::vector<constraint2_vector_array_t>& HvTrajectoryStock,
		const std::vector<size_t>& nc2FinalStock,
		const std::vector<constraint2_vector_t>& HvFinalStock,
		scalar_t& totalCost,
		size_t threadId) {

	calculateCostFunction(timeTrajectoriesStock, outputTrajectoriesStock, inputTrajectoriesStock, totalCost, threadId);
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
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateMeritFunction(
		const std::vector<scalar_array_t>& timeTrajectoriesStock,
		const std::vector<std::vector<size_t> >& nc1TrajectoriesStock,
		const std::vector<constraint1_vector_array_t>& EvTrajectoryStock,
		const std::vector<std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd> > >&  lagrangeTrajectoriesStock,
		const scalar_t& totalCost,
		scalar_t& meritFuntionValue,
		scalar_t& constraintISE)  {

	// add cost function
	meritFuntionValue = totalCost;

	// add the L2 penalty for constraint violation
	calculateConstraintISE(timeTrajectoriesStock, nc1TrajectoriesStock, EvTrajectoryStock, constraintISE);
	double pho = BASE::iteration_/(BASE::options_.maxIterationGSLQP_-1) * BASE::options_.meritFunctionRho_;
	meritFuntionValue += 0.5*pho*constraintISE;

	// add the the lagrangian term for the constraint
	scalar_t currentIntermediateMerit;
	scalar_t nextIntermediateMerit;
	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		// integrates the intermediate merit using the trapezoidal approximation method
		currentIntermediateMerit = 0.0;
		nextIntermediateMerit = 0.0;
		for (int k=0; k<timeTrajectoriesStock[i].size()-1; k++) {

			if (k==0)
				currentIntermediateMerit = EvTrajectoryStock[i][k].head(nc1TrajectoriesStock[i][k]).transpose() * lagrangeTrajectoriesStock[i][k];
			else
				currentIntermediateMerit = nextIntermediateMerit;

			nextIntermediateMerit = EvTrajectoryStock[i][k+1].head(nc1TrajectoriesStock[i][k+1]).transpose() * lagrangeTrajectoriesStock[i][k+1];

			meritFuntionValue += 0.5*(currentIntermediateMerit+nextIntermediateMerit)*(timeTrajectoriesStock[i][k+1]-timeTrajectoriesStock[i][k]);
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
double SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateConstraintISE(
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
 * 		+ quadratized final cost in the last subsystem: qFinal(t) + 0.5 y(t)QmFinal(t)y(t) + y(t)'QvFinal(t)
 * 		+ BASE::qFinal_: qFinal
 * 		+ BASE::qFinal_: QvFinal vector
 * 		+ BASE::qFinal_: QmFinal matrix
 *
 * 		+ as well as the constrained coefficients of
 * 			linearized system model
 * 			quadratized intermediate cost function
 * 			quadratized final cost
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::approximateOptimalControlProblem()
{
	for (int i=0; i<NUM_SUBSYSTEMS; i++)
	{
		subsystemProcessed_ = i;

		if(mp_options_.debugPrintMP_)
			printString("[MP] Starting approximation of subsystem " + std::to_string(i) + " out of " + std::to_string( (size_t) NUM_SUBSYSTEMS-1));

		approximateSubsystemLQ(i);

		if(mp_options_.debugPrintMP_)
			printString("[MP] ended approximation of subsystem " + std::to_string(i));


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
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateControllerAndLagrangian(
		std::vector<controller_t>& controllersStock,
		std::vector<lagrange_t>& lagrangeMultiplierFunctionsStock,
		std::vector<control_vector_array_t>& feedForwardConstraintInputStock,
		bool firstCall /*=true*/) {
	// this needs to be implemented!

	std::cout << "protected calculateControllerAndLagrangian(...) for calls from outside needs to be implemented" << std::endl;
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
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateRolloutLagrangeMultiplier(
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
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateRolloutCostate(
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
 * line search on the feedforward parts of the controller and lagrange multipliers.
 * Based on the option flag lineSearchByMeritFuntion_ it uses two different approaches for line search:
 * 		+ lineSearchByMeritFuntion_=TRUE: it uses the merit function to choose the best stepSize for the
 * 		feedforward elements of controller and lagrangeMultiplierFunction
 * 		+ lineSearchByMeritFuntion_=FALSE: the constraint correction term is added by a user defined stepSize.
 * 		The line search uses the pure cost function for choosing the best stepSize.
 *
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::lineSearch() {

	learningRateStar_ = 0.0;	// default learning rate as zero

	// display
	if (BASE::options_.dispayGSLQP_)
	{
		// less-equal operator for eigen vectors
		auto eigenVectorLessEqual = [] (const control_vector_t& u1, const control_vector_t& u2){ return u1.norm() < u2.norm(); };

		std::vector<control_vector_t> maxDeltaUffStock(NUM_SUBSYSTEMS);
		std::vector<control_vector_t> maxDeltaUffeStock(NUM_SUBSYSTEMS);
		for (size_t i=0; i<NUM_SUBSYSTEMS; i++)
		{
			maxDeltaUffStock[i]  = *std::max_element(controllers_[mp_options_.nThreads_][i].deltaUff_.begin(), controllers_[mp_options_.nThreads_][i].deltaUff_.template end(), eigenVectorLessEqual);
			maxDeltaUffeStock[i] = *std::max_element(feedForwardConstraintInputStock_[i].begin(), feedForwardConstraintInputStock_[i].end(), eigenVectorLessEqual);
		}
		control_vector_t maxDeltaUff  = *std::max_element(maxDeltaUffStock.begin(), maxDeltaUffStock.end(), eigenVectorLessEqual);
		control_vector_t maxDeltaUffe = *std::max_element(maxDeltaUffeStock.begin(), maxDeltaUffeStock.end(), eigenVectorLessEqual);

		std::cerr << "max delta_uff norm: " << maxDeltaUff.norm()  << std::endl;
		std::cerr << "max uff_error norm: " << maxDeltaUffe.norm() << std::endl;
	}


	for (int i=0; i<NUM_SUBSYSTEMS; i++)
		for (int k=0; k<controllers_[mp_options_.nThreads_][i].time_.size(); k++)
			if (BASE::options_.lineSearchByMeritFuntion_==true)
				controllers_[mp_options_.nThreads_][i].deltaUff_[k] += feedForwardConstraintInputStock_[i][k];
			else
				controllers_[mp_options_.nThreads_][i].uff_[k] += BASE::options_.constraintStepSize_*feedForwardConstraintInputStock_[i][k];


	// perform one rollout while the input correction for the type-1 constraint is considered.
	rollout(mp_options_.nThreads_, BASE::initState_, controllers_[mp_options_.nThreads_], BASE::nominalTimeTrajectoriesStock_,
			BASE::nominalStateTrajectoriesStock_, BASE::nominalInputTrajectoriesStock_, BASE::nominalOutputTrajectoriesStock_,
			BASE::nc1TrajectoriesStock_, BASE::EvTrajectoryStock_,BASE::nc2TrajectoriesStock_,
			BASE::HvTrajectoryStock_, BASE::nc2FinalStock_, BASE::HvFinalStock_);

	calculateCostFunction( BASE::nominalTimeTrajectoriesStock_, BASE::nominalOutputTrajectoriesStock_, BASE::nominalInputTrajectoriesStock_,
			BASE::nc2TrajectoriesStock_, BASE::HvTrajectoryStock_, BASE::nc2FinalStock_, BASE::HvFinalStock_, BASE::nominalTotalCost_, mp_options_.nThreads_);


	// calculate the merit function
	if (BASE::options_.lineSearchByMeritFuntion_==true)
	{
		// calculate the lagrange multiplier with learningRate zero
		calculateRolloutLagrangeMultiplier(BASE::nominalTimeTrajectoriesStock_, BASE::nominalOutputTrajectoriesStock_, BASE::lagrangeControllerStock_,
				BASE::nominalLagrangeTrajectoriesStock_);
		// calculate the merit function
		calculateMeritFunction(BASE::nominalTimeTrajectoriesStock_, BASE::nc1TrajectoriesStock_, BASE::EvTrajectoryStock_, BASE::nominalLagrangeTrajectoriesStock_, BASE::nominalTotalCost_,
				BASE::nominalTotalMerit_, BASE::nominalConstraint1ISE_);
	}
	else
	{
		BASE::nominalTotalMerit_ = BASE::nominalTotalCost_;
		calculateConstraintISE(BASE::nominalTimeTrajectoriesStock_, BASE::nc1TrajectoriesStock_, BASE::EvTrajectoryStock_, BASE::nominalConstraint1ISE_);
	}

	lowestTotalMerit_ = BASE::nominalTotalMerit_;
	lowestTotalCost_  = BASE::nominalTotalCost_;
	lowestConstraint1ISE_ = BASE::nominalConstraint1ISE_;


	// display
	if (BASE::options_.dispayGSLQP_)  {std::cerr << "\t learningRate 0.0 \t cost: " << BASE::nominalTotalCost_ << " \t merit: " << BASE::nominalTotalMerit_ <<
		" \t constraint ISE: " << BASE::nominalConstraint1ISE_ << std::endl;}


	initLScontrollersStock_ = controllers_[mp_options_.nThreads_];		// this will serve to init the workers
	initLSlagrangeMultiplierFunctionsStock_ = BASE::lagrangeControllerStock_;

	subsystemProcessed_ = 0; // not required for linesearch, but assign to not let it dangle around
	alphaProcessed_.clear();
	alphaTaken_ = 0;
	alphaBestFound_ = false;
	lsWorkerCompleted_ = 0;

	size_t maxNumOfLineSearches =  (int) (log(BASE::options_.minLearningRateGSLQP_/BASE::options_.maxLearningRateGSLQP_) / log(BASE::options_.lineSearchContractionRate_)) +1;
	alphaExpMax_ = maxNumOfLineSearches;
	alphaExpBest_ = maxNumOfLineSearches;
	alphaProcessed_.resize(maxNumOfLineSearches, 0);

	if(mp_options_.debugPrintMP_)
		printString("[MP]: calculated maximum number of line searches " + std::to_string(alphaExpMax_));

	if(mp_options_.debugPrintMP_)
		printString("[MP] Waking up workers for line search ");

	workerTask_ = LINE_SEARCH;

	std::unique_lock<std::mutex> lock (workerWakeUpMutex_);
	workerWakeUpCondition_.notify_all();
	lock.unlock();

	if(mp_options_.debugPrintMP_)
		printString("[MP] Will sleep now until we have results ");


	std::unique_lock<std::mutex> waitLock(alphaBestFoundMutex_);
	while(lsWorkerCompleted_.load() < mp_options_.nThreads_)
		alphaBestFoundCondition_.wait(waitLock);

	waitLock.unlock();

	workerTask_ = IDLE;

	if(mp_options_.debugPrintMP_)
		printString("[MP]: Woke up again, should have results now.");



	// clear the feedforward increments
	for (int j=0; j<mp_options_.nThreads_+1; j++){
		for (size_t i=0; i<NUM_SUBSYSTEMS; i++)
		{
			controllers_[mp_options_.nThreads_][i].deltaUff_.clear();
			BASE::lagrangeControllerStock_[i].deltaUff_.clear();
		}
	}

	// reset integrator events
	killIntegrationEventHandler_->resetEvent();	// reset all integrations

	// display
	if (BASE::options_.dispayGSLQP_)
		printString("The chosen learningRate is: " + std::to_string(learningRateStar_));

	BASE::nominalTotalMerit_ = lowestTotalMerit_;
	BASE::nominalTotalCost_ =	lowestTotalCost_;
	BASE::nominalConstraint1ISE_ = lowestConstraint1ISE_;
}



/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * get the calculated optimal controller structure
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getController(std::vector<controller_t>& controllersStock) {

	controllersStock = controllers_[mp_options_.nThreads_];
}


template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::setController(const std::vector<controller_t>& controllersStock){
	controllers_[mp_options_.nThreads_] = controllersStock;
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
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getValueFuntion(
		const scalar_t& time, const output_vector_t& output, scalar_t& valueFuntion)  {

	int activeSubsystem = -1;

	for (int i=0; i<NUM_SUBSYSTEMS; i++)
	{
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
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getCostFuntion(
		scalar_t& costFunction, scalar_t& constriantISE)  {

	costFunction = BASE::nominalTotalCost_;
	constriantISE = BASE::nominalConstraint1ISE_;}


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
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getNominalTrajectories(
		std::vector<scalar_array_t>& nominalTimeTrajectoriesStock,
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
//void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::solveSequentialRiccatiEquations(const scalar_t& learningRate)  {
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
//		std::shared_ptr<RiccatiEquations_t> riccatiEquationsPtr( new RiccatiEquations_t() );		riccatiEquationsPtr->setData(learningRate, i, BASE::switchingTimes_[i], BASE::switchingTimes_[i+1],
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
//				1e-3, BASE::options_.AbsTolODE_, BASE::options_.RelTolODE_,  maxNumSteps);
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
//		std::shared_ptr<ErrorEquation_t> errorEquationPtr( new ErrorEquation_t());
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
 * make the given square matrix psd
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
template <typename Derived>
bool SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::makePSD(Eigen::MatrixBase<Derived>& squareMatrix) {

	if (squareMatrix.rows() != squareMatrix.cols())
		throw std::runtime_error("Not a square matrix: makePSD() method is for square matrix.");

	Eigen::SelfAdjointEigenSolver<Derived> eig(squareMatrix);
	Eigen::VectorXd lambda = eig.eigenvalues();

	bool hasNegativeEigenValue = false;
	for (size_t j=0; j<lambda.size() ; j++)
		if (lambda(j) < 0.0) {
			hasNegativeEigenValue = true;
			lambda(j) = 0.0;
		}

	if (hasNegativeEigenValue)
		squareMatrix = eig.eigenvectors() * lambda.asDiagonal() * eig.eigenvectors().inverse();
	//	else
	//		squareMatrix = 0.5*(squareMatrix+squareMatrix.transpose()).eval();

	return hasNegativeEigenValue;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * run the SLQP algorithm for a given state and switching times
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::run(const state_vector_t& initState, const std::vector<scalar_t>& switchingTimes,
		const std::vector<controller_t>& initialControllersStock) {

	BASE::iterationCost_.clear();
	BASE::iterationISE1_.clear();

	if (switchingTimes.size() != NUM_SUBSYSTEMS+1)
		throw std::runtime_error("Number of switching times should be one plus the number of subsystems.");


	if (initialControllersStock.empty()==false) {
		if (initialControllersStock.size() != NUM_SUBSYSTEMS)
			throw std::runtime_error("initialControllersStock has less controllers than the number of subsystems");
		for(size_t i = 0; i<mp_options_.nThreads_+1; i++)
			controllers_[i] = initialControllersStock;
	}
	else
		if (controllers_[mp_options_.nThreads_].empty()==true)
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
	double relConstraint1ISE;
	bool isConstraint1Satisfied  = false;
	bool isCostFunctionConverged = false;
	bool isOptimizationConverged = false;
	nominalLagrangeMultiplierUpdated_ = false;

	// initial controller rollout
	rollout(mp_options_.nThreads_,
			BASE::initState_,
			controllers_[mp_options_.nThreads_],
			BASE::nominalTimeTrajectoriesStock_, BASE::nominalStateTrajectoriesStock_, BASE::nominalInputTrajectoriesStock_, BASE::nominalOutputTrajectoriesStock_,
			BASE::nc1TrajectoriesStock_, BASE::EvTrajectoryStock_,
			BASE::nc2TrajectoriesStock_, BASE::HvTrajectoryStock_, BASE::nc2FinalStock_, BASE::HvFinalStock_);

	// initial controller cost
	calculateCostFunction(BASE::nominalTimeTrajectoriesStock_, BASE::nominalOutputTrajectoriesStock_, BASE::nominalInputTrajectoriesStock_,
			BASE::nc2TrajectoriesStock_, BASE::HvTrajectoryStock_, BASE::nc2FinalStock_, BASE::HvFinalStock_,
			BASE::nominalTotalCost_, mp_options_.nThreads_);

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
		auto start = std::chrono::high_resolution_clock::now();

		Eigen::setNbThreads(1); // disable Eigen multi-threading

		approximateOptimalControlProblem();

		Eigen::setNbThreads(0); // restore default Eigen thread number

		auto end = std::chrono::high_resolution_clock::now();
		if(mp_options_.debugPrintMP_){
			std::chrono::duration<double, std::milli> diff = end - start;
			printString("iteration " + std::to_string(BASE::iteration_) + ": mp LQ approximation took " + std::to_string(diff.count()) + "ms");
		}

		// solve Riccati equations
		auto start2 = std::chrono::high_resolution_clock::now();

		this->solveSequentialRiccatiEquations(1.0 /*nominal learningRate*/);

		auto end2 = std::chrono::high_resolution_clock::now();
		if(mp_options_.debugPrintMP_){
			std::chrono::duration<double, std::milli> diff2 = end2 - start2;
			printString("iteration " + std::to_string(BASE::iteration_) + ": mp solve riccati took " + std::to_string(diff2.count()) + "ms");
		}

		auto start3 = std::chrono::high_resolution_clock::now();
		Eigen::setNbThreads(1); // disable Eigen multi-threading

		calculateControllerAndLagrangian();

		Eigen::setNbThreads(0); // restore default Eigen thread number
		auto end3 = std::chrono::high_resolution_clock::now();

		if(mp_options_.debugPrintMP_){
			std::chrono::duration<double, std::milli> diff3 = end3 - start3;
			printString("iteration " + std::to_string(BASE::iteration_) + ": mp calc controller and lagrangian took " + std::to_string(diff3.count()) + "ms");
		}

		nominalLagrangeMultiplierUpdated_ = true;

		// finding the optimal learningRate
		auto start4 = std::chrono::high_resolution_clock::now();
		Eigen::setNbThreads(1); // disable Eigen multi-threading

		lineSearch();

		Eigen::setNbThreads(0); // restore default Eigen thread number
		auto end4 = std::chrono::high_resolution_clock::now();

		if(mp_options_.debugPrintMP_){
			std::chrono::duration<double, std::milli> diff4 = end4 - start4;
			printString("iteration " + std::to_string(BASE::iteration_) + ": mp line search took " + std::to_string(diff4.count()) + "ms");
		}

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
		isCostFunctionConverged = learningRateStar_==0 || relCost<=BASE::options_.minRelCostGSLQP_;
		isOptimizationConverged = isCostFunctionConverged==true && isConstraint1Satisfied==true;

		// display
		if (BASE::options_.dispayGSLQP_)  {
			std::cerr << "optimization cost:         " << BASE::nominalTotalCost_ << std::endl;
			std::cerr << "constraint type-1 ISE:     " << BASE::nominalConstraint1ISE_ << std::endl;
			std::cerr << "constraint type-1 MaxNorm: " << constraint1MaxNorm << std::endl;
			std::cerr << "constraint type-2 ISE:     " << nominalConstraint2ISE << std::endl;
			std::cerr << "constraint type-2 MaxNorm: " << constraint2MaxNorm << std::endl;
			if (std::accumulate(BASE::nc2FinalStock_.begin(), BASE::nc2FinalStock_.end(), 0) > 0) {
				std::cerr << "final constraint type-2: 	 ";
				for(size_t i=0; i<NUM_SUBSYSTEMS; i++) std::cerr << "[" << i  << "]: " << BASE::HvFinalStock_[i].head(BASE::nc2FinalStock_[i]).transpose() << ",  ";
				std::cerr << std::endl;
			}
		}
		if(BASE::options_.displayShortSummary_){
			printString("#### Iter " + std::to_string(BASE::iteration_-1) + ".   opt. cost: " + std::to_string(BASE::nominalTotalCost_)+ ".    constraint ISE: "
					+std::to_string(BASE::nominalConstraint1ISE_) +".   constr MaxNorm: " + std::to_string(constraint1MaxNorm));
		}
	}  // end of while loop

	// linearizing the dynamics and quadratizing the cost function along nominal trajectories
	approximateOptimalControlProblem();

	// solve Riccati equations with learningRate zero
	this->solveSequentialRiccatiEquations(0.0 /*learningRate*/);

	// calculate the nominal co-state
	calculateRolloutCostate(BASE::nominalTimeTrajectoriesStock_, BASE::nominalOutputTrajectoriesStock_, BASE::nominalcostateTrajectoriesStock_);

	// display
	if (BASE::options_.dispayGSLQP_)  {
		std::cout << "\n+++++++++++++++++++++++++++++++++++" << std::endl;
		std::cout <<   "++++++ SLQP solver has ended ++++++" << std::endl;
		std::cout <<   "+++++++++++++++++++++++++++++++++++" << std::endl;
		if (isOptimizationConverged) {
			if (learningRateStar_==0)
				printString("SLQP successfully terminates as learningRate reduced to zero.");
			else
				printString("SLQP successfully terminates as cost relative change (relCost=" + std::to_string(relCost) + ") reached to the minimum value.");

			if (BASE::nominalConstraint1ISE_<=BASE::options_.minAbsConstraint1ISE_)
				printString("Type-1 constraint absolute ISE (absConstraint1ISE=" + std::to_string(BASE::nominalConstraint1ISE_) + ") reached to the minimum value.");
			else
				printString("Type-1 constraint relative ISE (relConstraint1ISE=" + std::to_string(relConstraint1ISE) + ") reached to the minimum value.");
		} else
			printString("Maximum number of iterations is reached.");
	}

}


template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::launchWorkerThreads()
{
	workersActive_ = true;
	workerTask_ = IDLE;

	for (size_t i=0; i < mp_options_.nThreads_; i++)
	{
		workerThreads_.push_back(std::thread(&SLQP_MP::threadWork, this, i));
	}
}


template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::threadWork(size_t threadId)
{
	if(mp_options_.debugPrintMP_)
		printString("[Thread " + std::to_string(threadId) + "]: launched");

	// local variables
	size_t uniqueProcessID = 0;
	size_t subsystemProcessed_local = 0;
	size_t iteration_local = BASE::iteration_;
	int workerTask_local = IDLE;

	while(workersActive_)
	{
		subsystemProcessed_local = subsystemProcessed_.load();
		workerTask_local = workerTask_.load();
		iteration_local = BASE::iteration_;

		// display
		if(mp_options_.debugPrintMP_){
			printString("[Thread " + std::to_string(threadId) + "]: previous procId: " + std::to_string(uniqueProcessID) +
					", current procId: " +std::to_string(generateUniqueProcessID(iteration_local, (int) workerTask_local, (int) subsystemProcessed_local)));
		}

		/* We want to put the worker to sleep if
		 * - the workerTask_ is IDLE
		 * - or we are finished both workerTask_ is not yet reset, thus the process ID is still the same
		 * */
		if ( workerTask_local == IDLE || uniqueProcessID == generateUniqueProcessID(iteration_local, (int) workerTask_local, (int) subsystemProcessed_local))
		{
			if(mp_options_.debugPrintMP_)
				printString("[Thread " + std::to_string(threadId) + "]: going to sleep !");

			// sleep until the state is not IDLE any more and we have a different process ID than before
			std::unique_lock<std::mutex> waitLock(workerWakeUpMutex_);
			while(workerTask_ == IDLE ||  (uniqueProcessID == generateUniqueProcessID(BASE::iteration_, (int)workerTask_.load(), (int) subsystemProcessed_.load()))){
				workerWakeUpCondition_.wait(waitLock);
			}
			waitLock.unlock();

			subsystemProcessed_local = subsystemProcessed_.load();
			workerTask_local = workerTask_.load();
			iteration_local = BASE::iteration_;

			if(mp_options_.debugPrintMP_)
				printString("[Thread " + std::to_string(threadId) + "]: woke up !");
		}

		if (!workersActive_)
			break;

		switch(workerTask_local)
		{
		case APPROXIMATE_LQ:
		{
			if(mp_options_.debugPrintMP_)
				printString("[Thread " + std::to_string(threadId) + "]: now busy with APPROXIMATE_LQ on subsystem " + std::to_string(subsystemProcessed_local));

			approximateSubsystemLQWorker(threadId, subsystemProcessed_local);
			uniqueProcessID = generateUniqueProcessID (iteration_local, APPROXIMATE_LQ, subsystemProcessed_local);

			break;
		}
		case CALCULATE_CONTROLLER_AND_LAGRANGIAN:
		{
			if(mp_options_.debugPrintMP_)
				printString("[Thread " + std::to_string(threadId) + "]: now busy with CALCULATE_CONTROLLER_AND_LAGRANGIAN !");

			calculateControllerAndLagrangianWorker(threadId, subsystemProcessed_local);
			uniqueProcessID = generateUniqueProcessID (iteration_local, CALCULATE_CONTROLLER_AND_LAGRANGIAN, subsystemProcessed_local);

			break;
		}
		case LINE_SEARCH:
		{
			if(mp_options_.debugPrintMP_)
				printString("[Thread " + std::to_string(threadId) + "]: now busy with LINE_SEARCH !");

			lineSearchWorker(threadId);
			uniqueProcessID = generateUniqueProcessID (iteration_local, LINE_SEARCH, subsystemProcessed_local);
			break;
		}
		case SHUTDOWN:
		{
			if(mp_options_.debugPrintMP_)
				printString("[Thread "+ std::to_string(threadId) +"]: now shutting down!");
			return;
		}
		}

		if(mp_options_.debugPrintMP_)
			printString("[Thread " + std::to_string(threadId) +"]: done with job. Will wait for next now!");
	}
}


template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::approximateSubsystemLQ(const size_t i)
{
	size_t N =   BASE::nominalTimeTrajectoriesStock_[i].size();

	// initialize subsystem i dynamics derivatives
	for(size_t j = 0; j< mp_options_.nThreads_+1; j++)
	{
		assert( BASE::nominalTimeTrajectoriesStock_[i].size() == BASE::nominalStateTrajectoriesStock_[i].size());
		linearizedSystems_[j][i]->initializeModel(BASE::switchingTimes_, BASE::nominalStateTrajectoriesStock_[i].front(), i, "GSLPQ");
	}

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

	// for constraints
	BASE::RmConstrainedTrajectoryStock_[i].resize(N);
	BASE::DmDagerTrajectoryStock_[i].resize(N);
	BASE::AmConstrainedTrajectoryStock_[i].resize(N);
	BASE::QmConstrainedTrajectoryStock_[i].resize(N);
	BASE::QvConstrainedTrajectoryStock_[i].resize(N);
	BASE::EvProjectedTrajectoryStock_[i].resize(N);
	BASE::CmProjectedTrajectoryStock_[i].resize(N);
	BASE::DmProjectedTrajectoryStock_[i].resize(N);


	kTaken_approx_[i] = 0;
	kCompleted_approx_[i]= 0;
	KMax_subsystem_approx_[i] = N;

	if(mp_options_.debugPrintMP_)
		printString("[MP]: Waking up workers to do linearisation for subsystem " + std::to_string(i));

	workerTask_ = APPROXIMATE_LQ;
	std::unique_lock<std::mutex> lock (workerWakeUpMutex_);
	workerWakeUpCondition_.notify_all();
	lock.unlock();

	if(mp_options_.debugPrintMP_)
		printString("[MP]: Will wait now until workers have linearized dynamics of subsystem " + std::to_string(i));

	std::unique_lock<std::mutex> waitLock(kCompletedMutex_);
	while(kCompleted_approx_[i].load() < KMax_subsystem_approx_[i]){
		kCompletedCondition_.wait(waitLock);
	}

	waitLock.unlock();
	workerTask_ = IDLE;


	if (i==NUM_SUBSYSTEMS-1) // if last subsystem, set terminal cost
	{
		if(mp_options_.debugPrintMP_)
			printString("[MP]: Approximating terminal cost with single thread, subsystem  " + std::to_string(i));

		costFunctions_[mp_options_.nThreads_][i]->setCurrentStateAndControl(BASE::nominalTimeTrajectoriesStock_[i].back(),
				BASE::nominalOutputTrajectoriesStock_[i].back(), BASE::nominalInputTrajectoriesStock_[i].back());

		costFunctions_[mp_options_.nThreads_][i]->terminalCost(BASE::qFinalStock_[i](0));
		costFunctions_[mp_options_.nThreads_][i]->terminalCostStateDerivative(BASE::QvFinalStock_[i]);
		costFunctions_[mp_options_.nThreads_][i]->terminalCostStateSecondDerivative(BASE::QmFinalStock_[i]);

		// making sure that Qm remains PSD
		makePSD(BASE::QmFinalStock_[i]);
	}
	else {
		BASE::qFinalStock_[i].setZero();
		BASE::QvFinalStock_[i].setZero();
		BASE::QmFinalStock_[i].setZero();
	}

	// constrained type-2 final coefficients
	if (BASE::nc2FinalStock_[i] > 0) {
		size_t nc2 = BASE::nc2FinalStock_[i];

		linearizedSystems_[mp_options_.nThreads_][i]->setCurrentStateAndControl(
				BASE::nominalTimeTrajectoriesStock_[i].back(),
				BASE::nominalStateTrajectoriesStock_[i].back(),
				BASE::nominalInputTrajectoriesStock_[i].back(),
				BASE::nominalOutputTrajectoriesStock_[i].back());

		linearizedSystems_[mp_options_.nThreads_][i]->getFinalConstraint2DerivativesState(BASE::FmFinalStock_[i]);

		double stateConstraintPenalty = BASE::options_.stateConstraintPenaltyCoeff_ * pow(BASE::options_.stateConstraintPenaltyBase_, BASE::iteration_);

		BASE::qFinalStock_[i]  += 0.5 * stateConstraintPenalty * BASE::HvFinalStock_[i].head(nc2).transpose() * BASE::HvFinalStock_[i].head(nc2);
		BASE::QvFinalStock_[i] += stateConstraintPenalty * BASE::FmFinalStock_[i].topRows(nc2).transpose() * BASE::HvFinalStock_[i].head(nc2);
		BASE::QmFinalStock_[i] += stateConstraintPenalty * BASE::FmFinalStock_[i].topRows(nc2).transpose() * BASE::FmFinalStock_[i].topRows(nc2);
	}

}



template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateControllerAndLagrangian() {

	for (int i=0; i<NUM_SUBSYSTEMS; i++)
	{
		subsystemProcessed_ =  i;

		kTaken_ctrl_[i] = 0;
		kCompleted_ctrl_[i] = 0;
		KMax_subsystem_ctrl_[i]  = BASE::SsTimeTrajectoryStock_[i].size(); // number of elements in the trajectory of this subsystem

		// initialize interpolators
		for(size_t n = 0; n< mp_options_.nThreads_+1; n++)
		{
			// functions for controller and lagrange-multiplier
			nominalOutputFunc_[n].setTimeStamp( &(BASE::nominalTimeTrajectoriesStock_[i]) );
			nominalOutputFunc_[n].setData( &(BASE::nominalOutputTrajectoriesStock_[i]) );

			nominalInputFunc_[n].setTimeStamp( &(BASE::nominalTimeTrajectoriesStock_[i]) );
			nominalInputFunc_[n].setData( &(BASE::nominalInputTrajectoriesStock_[i]) );

			BmFunc_[n].setTimeStamp( &(BASE::nominalTimeTrajectoriesStock_[i]) );
			BmFunc_[n].setData( &(BASE::BmTrajectoryStock_[i]) );

			PmFunc_[n].setTimeStamp( &(BASE::nominalTimeTrajectoriesStock_[i]) );
			PmFunc_[n].setData( &(BASE::PmTrajectoryStock_[i]) );

			RmInverseFunc_[n].setTimeStamp( &(BASE::nominalTimeTrajectoriesStock_[i]) );
			RmInverseFunc_[n].setData( &(BASE::RmInverseTrajectoryStock_[i]) );

			RvFunc_[n].setTimeStamp( &(BASE::nominalTimeTrajectoriesStock_[i]) );
			RvFunc_[n].setData( &(BASE::RvTrajectoryStock_[i]) );

			EvProjectedFunc_[n].setTimeStamp( &(BASE::nominalTimeTrajectoriesStock_[i]) );
			EvProjectedFunc_[n].setData( &(BASE::EvProjectedTrajectoryStock_[i]) );

			CmProjectedFunc_[n].setTimeStamp( &(BASE::nominalTimeTrajectoriesStock_[i]) );
			CmProjectedFunc_[n].setData( &(BASE::CmProjectedTrajectoryStock_[i]) );

			DmProjectedFunc_[n].setTimeStamp( &(BASE::nominalTimeTrajectoriesStock_[i]) );
			DmProjectedFunc_[n].setData( &(BASE::DmProjectedTrajectoryStock_[i]) );

			// functions for lagrange multiplier only
			RmFunc_[n].setTimeStamp( &(BASE::nominalTimeTrajectoriesStock_[i]) );
			RmFunc_[n].setData( &(BASE::RmTrajectoryStock_[i]) );

			DmDagerFunc_[n].setTimeStamp( &(BASE::nominalTimeTrajectoriesStock_[i]) );
			DmDagerFunc_[n].setData( &(BASE::DmDagerTrajectoryStock_[i]) );

			if (nominalLagrangeMultiplierUpdated_ == true) {
				nominalLagrangeMultiplierFunc_[n].setTimeStamp( &(BASE::nominalTimeTrajectoriesStock_[i]) );
				nominalLagrangeMultiplierFunc_[n].setData( &(BASE::nominalLagrangeTrajectoriesStock_[i]) );
			}
		}

		controllers_[mp_options_.nThreads_][i].time_ = BASE::SsTimeTrajectoryStock_[i];
		controllers_[mp_options_.nThreads_][i].k_.resize(KMax_subsystem_ctrl_[i]);
		controllers_[mp_options_.nThreads_][i].uff_.resize(KMax_subsystem_ctrl_[i]);
		controllers_[mp_options_.nThreads_][i].deltaUff_.resize(KMax_subsystem_ctrl_[i]);

		feedForwardConstraintInputStock_[i].resize(KMax_subsystem_ctrl_[i]);

		BASE::lagrangeControllerStock_[i].time_ = BASE::SsTimeTrajectoryStock_[i];
		BASE::lagrangeControllerStock_[i].k_.resize(KMax_subsystem_ctrl_[i]);
		BASE::lagrangeControllerStock_[i].uff_.resize(KMax_subsystem_ctrl_[i]);
		BASE::lagrangeControllerStock_[i].deltaUff_.resize(KMax_subsystem_ctrl_[i]);


		if(mp_options_.debugPrintMP_)
			printString("[MP]: Waking up workers to calc. controller for subsystem " + std::to_string(i));

		workerTask_ = CALCULATE_CONTROLLER_AND_LAGRANGIAN;
		std::unique_lock<std::mutex> lock (workerWakeUpMutex_);
		workerWakeUpCondition_.notify_all();
		lock.unlock();

		if(mp_options_.debugPrintMP_)
			printString("[MP]: Will wait now controllers have been calculated for subsystem " + std::to_string(i));


		std::unique_lock<std::mutex> waitLock(kCompletedMutex_);

		while(kCompleted_ctrl_[i].load() < KMax_subsystem_ctrl_[i] ){
			kCompletedCondition_.wait(waitLock);
		}

		waitLock.unlock();
		workerTask_ = IDLE;

		if(mp_options_.debugPrintMP_)
			printString("[MP]: Back to main thread, workers should now have designed controllers for subsystem " + std::to_string(i));

	}  // end of i loop

}

template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
size_t SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::approximateSubsystemLQWorker(size_t threadId, size_t subsystemProcessed)
{

	size_t k = 0;
	size_t kCompleted_local = 0;

	while(true)
	{
		k = kTaken_approx_[subsystemProcessed]++;

		if(k < KMax_subsystem_approx_[subsystemProcessed]){

			if(mp_options_.debugPrintMP_){
				if (k%10 == 0) {
					printString("[Thread " + std::to_string(threadId) + "], subsystem " + std::to_string(subsystemProcessed)
					+ ":Start approximating system LQ on index k = " + std::to_string(k) + " out of " + std::to_string(KMax_subsystem_approx_[subsystemProcessed]-1));
				}
			}

			executeApproximateSubsystemLQ(threadId, k, subsystemProcessed);
			kCompleted_local = ++kCompleted_approx_[subsystemProcessed];
		}

		if (k >= KMax_subsystem_approx_[subsystemProcessed]-1) // if all k's are already covered, notify and return
		{
			if(kCompleted_local >=KMax_subsystem_approx_[subsystemProcessed])
			{
				if(mp_options_.debugPrintMP_){
					printString("[Thread " + std::to_string(threadId) + "], subsystem "
							+ std::to_string(subsystemProcessed) + ", k " + std::to_string(k)
					+ ", kCompleted_local " + std::to_string(kCompleted_local)
					+ ", KMax_subsystem_approx_ " + std::to_string(KMax_subsystem_approx_[subsystemProcessed])
					+ ": leaving approximateSubsystemLQWorker AND NOTIFYING ");
				}
				std::unique_lock<std::mutex> lock (kCompletedMutex_);
				kCompletedCondition_.notify_all();
				lock.unlock();
			}
			else{
				if(mp_options_.debugPrintMP_){
					printString("[Thread " + std::to_string(threadId) + "], subsystem "
							+ std::to_string(subsystemProcessed) + ", k " + std::to_string(k) + ", kCompleted_local " + std::to_string(kCompleted_local)
					+ ", KMax_subsystem_approx_ " + std::to_string(KMax_subsystem_approx_[subsystemProcessed])
					+ ": leaving approximateSubsystemLQWorker but NOT notifying ");
				}
			}

			return subsystemProcessed;
		}
	}

	return subsystemProcessed;
}


template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
size_t SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateControllerAndLagrangianWorker(size_t threadId, size_t subsystemProcessed)
{

	size_t k = 0;
	size_t kCompleted_local = 0;

	while(true)
	{
		k = kTaken_ctrl_[subsystemProcessed]++;

		if(k < KMax_subsystem_ctrl_[subsystemProcessed]){

			if(mp_options_.debugPrintMP_){
				if (k%10 == 0) {
					printString("[Thread " + std::to_string(threadId) + "]: Start calculating controller on index k = " + std::to_string(k) +
							" out of " + std::to_string(KMax_subsystem_ctrl_[subsystemProcessed]-1));
				}
			}

			executeCalculateControllerAndLagrangian(threadId, k, subsystemProcessed);
			kCompleted_local = ++kCompleted_ctrl_[subsystemProcessed];
		}


		if (k >= KMax_subsystem_ctrl_[subsystemProcessed]-1)	// if all k's are already covered, notify and return
		{
			if(kCompleted_local>=KMax_subsystem_ctrl_[subsystemProcessed])
			{
				if(mp_options_.debugPrintMP_)
					printString("[Thread " + std::to_string(threadId) + "], subsystem " + std::to_string(subsystemProcessed) + ": leaving calculateControllerAndLagrangianWorker() AND NOTIFYING ");

				std::unique_lock<std::mutex> lock (kCompletedMutex_);
				kCompletedCondition_.notify_all();
				lock.unlock();
			}else
			{
				if(mp_options_.debugPrintMP_)
					printString("[Thread " + std::to_string(threadId) + "], subsystem " + std::to_string(subsystemProcessed) + ": leaving calculateControllerAndLagrangianWorker() but NOT notifying ");				}

			return subsystemProcessed;
		}
	}

	return subsystemProcessed;
}


template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
size_t SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::executeApproximateSubsystemLQ(size_t threadId, size_t k, size_t subsystemProcessed)
{
	const size_t i = subsystemProcessed;

	// LINEARIZE SYSTEM DYNAMICS AND CONSTRAINTS
	linearizedSystems_[threadId][i]->setCurrentStateAndControl(
			BASE::nominalTimeTrajectoriesStock_[i][k],
			BASE::nominalStateTrajectoriesStock_[i][k],
			BASE::nominalInputTrajectoriesStock_[i][k],
			BASE::nominalOutputTrajectoriesStock_[i][k]);

	linearizedSystems_[threadId][i]->getDerivativeState(BASE::AmTrajectoryStock_[i][k]);
	linearizedSystems_[threadId][i]->getDerivativesControl(BASE::BmTrajectoryStock_[i][k]);

	// if constraint type 1 is active
	if (BASE::nc1TrajectoriesStock_[i][k] > 0)
	{
		linearizedSystems_[threadId][i]->getConstraint1DerivativesState(BASE::CmTrajectoryStock_[i][k]);
		linearizedSystems_[threadId][i]->getConstraint1DerivativesControl(BASE::DmTrajectoryStock_[i][k]);
	}

	// if constraint type 2 is active
	if (BASE::nc2TrajectoriesStock_[i][k] > 0) {
		linearizedSystems_[threadId][i]->getConstraint2DerivativesState(BASE::FmTrajectoryStock_[i][k]);
	}

	// QUADRATIC APPROXIMATION TO THE COST FUNCTION
	costFunctions_[threadId][i]->setCurrentStateAndControl(
			BASE::nominalTimeTrajectoriesStock_[i][k],
			BASE::nominalOutputTrajectoriesStock_[i][k],
			BASE::nominalInputTrajectoriesStock_[i][k]);
	costFunctions_[threadId][i]->evaluate(BASE::qTrajectoryStock_[i][k](0));
	costFunctions_[threadId][i]->stateDerivative(BASE::QvTrajectoryStock_[i][k]);
	costFunctions_[threadId][i]->stateSecondDerivative(BASE::QmTrajectoryStock_[i][k]);
	costFunctions_[threadId][i]->controlDerivative(BASE::RvTrajectoryStock_[i][k]);
	costFunctions_[threadId][i]->controlSecondDerivative(BASE::RmTrajectoryStock_[i][k]);
	BASE::RmInverseTrajectoryStock_[i][k] = BASE::RmTrajectoryStock_[i][k].inverse();
	costFunctions_[threadId][i]->stateControlDerivative(BASE::PmTrajectoryStock_[i][k]);


	// constraint type 2 coefficients
	double stateConstraintPenalty = BASE::options_.stateConstraintPenaltyCoeff_ * pow(BASE::options_.stateConstraintPenaltyBase_, BASE::iteration_);
	size_t nc2 = BASE::nc2TrajectoriesStock_[i][k];

	if (nc2 > 0) {
		//				subsystemDerivativesPtrStock_[i]->getConstraint2DerivativesState(FmTrajectoryStock_[i][k]);
		BASE::qTrajectoryStock_[i][k]  += 0.5 * stateConstraintPenalty * BASE::HvTrajectoryStock_[i][k].head(nc2).transpose() * BASE::HvTrajectoryStock_[i][k].head(nc2);
		BASE::QvTrajectoryStock_[i][k] += stateConstraintPenalty * BASE::FmTrajectoryStock_[i][k].topRows(nc2).transpose() * BASE::HvTrajectoryStock_[i][k].head(nc2);
		BASE::QmTrajectoryStock_[i][k] += stateConstraintPenalty * BASE::FmTrajectoryStock_[i][k].topRows(nc2).transpose() * BASE::FmTrajectoryStock_[i][k].topRows(nc2);
	}

	// constraint type 1 coefficients
	size_t nc1 = BASE::nc1TrajectoriesStock_[i][k];

	if (nc1 == 0)
	{
		BASE::DmDagerTrajectoryStock_[i][k].setZero();
		BASE::EvProjectedTrajectoryStock_[i][k].setZero();
		BASE::CmProjectedTrajectoryStock_[i][k].setZero();
		BASE::DmProjectedTrajectoryStock_[i][k].setZero();

		BASE::AmConstrainedTrajectoryStock_[i][k] = BASE::AmTrajectoryStock_[i][k];
		BASE::QmConstrainedTrajectoryStock_[i][k] = BASE::QmTrajectoryStock_[i][k];
		BASE::QvConstrainedTrajectoryStock_[i][k] = BASE::QvTrajectoryStock_[i][k];
		BASE::RmConstrainedTrajectoryStock_[i][k] = BASE::RmTrajectoryStock_[i][k];
	}
	else
	{
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
	makePSD(BASE::QmConstrainedTrajectoryStock_[i][k]);

	return i;
}


template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
size_t SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::executeCalculateControllerAndLagrangian(size_t threadId, size_t k, size_t subsystemProcessed)
{

	const size_t i = subsystemProcessed;

	const double time = BASE::SsTimeTrajectoryStock_[i][k];
	size_t greatestLessTimeStampIndex;

	// local variables
	output_vector_t nominalOutput;
	control_vector_t nominalInput;
	control_gain_matrix_t Bm;
	control_feedback_t Pm;
	control_vector_t Rv;
	control_matrix_t RmInverse;
	control_vector_t EvProjected;
	control_feedback_t CmProjected;
	control_matrix_t DmProjected;
	control_constraint1_matrix_t DmDager;
	control_matrix_t Rm;

	nominalOutputFunc_[threadId].interpolate(time, nominalOutput);
	greatestLessTimeStampIndex = nominalOutputFunc_[threadId].getGreatestLessTimeStampIndex();

	nominalInputFunc_[threadId].interpolate(time, nominalInput, greatestLessTimeStampIndex);

	// interpolate
	BmFunc_[threadId].interpolate(time, Bm, greatestLessTimeStampIndex);
	PmFunc_[threadId].interpolate(time, Pm, greatestLessTimeStampIndex);
	RvFunc_[threadId].interpolate(time, Rv, greatestLessTimeStampIndex);
	RmInverseFunc_[threadId].interpolate(time, RmInverse, greatestLessTimeStampIndex);
	EvProjectedFunc_[threadId].interpolate(time, EvProjected, greatestLessTimeStampIndex);
	CmProjectedFunc_[threadId].interpolate(time, CmProjected, greatestLessTimeStampIndex);
	DmProjectedFunc_[threadId].interpolate(time, DmProjected, greatestLessTimeStampIndex);

	control_feedback_t Lm  = RmInverse * (Pm + Bm.transpose()*BASE::SmTrajectoryStock_[i][k]);
	control_vector_t   Lv  = RmInverse * (Rv + Bm.transpose()*BASE::SvTrajectoryStock_[i][k]);
	control_vector_t   Lve = RmInverse * (Bm.transpose()*BASE::SveTrajectoryStock_[i][k]);

	control_matrix_t DmNullProjection = control_matrix_t::Identity()-DmProjected;
	controllers_[mp_options_.nThreads_][i].k_[k]   = -DmNullProjection*Lm - CmProjected;
	controllers_[mp_options_.nThreads_][i].uff_[k] = nominalInput - controllers_[mp_options_.nThreads_][i].k_[k]*nominalOutput;
	controllers_[mp_options_.nThreads_][i].deltaUff_[k]      = -DmNullProjection*Lv;
	feedForwardConstraintInputStock_[i][k] = -DmNullProjection*Lve - EvProjected;


	// checking the numerical stability of the controller parameters
	try {
		if (controllers_[mp_options_.nThreads_][i].k_[k] != controllers_[mp_options_.nThreads_][i].k_[k])
			throw std::runtime_error("Feedback gains are unstable.");
		if (controllers_[mp_options_.nThreads_][i].deltaUff_[k] != controllers_[mp_options_.nThreads_][i].deltaUff_[k])
			throw std::runtime_error("feedForwardControl is unstable.");
		if (feedForwardConstraintInputStock_[i][k] != feedForwardConstraintInputStock_[i][k])
			throw std::runtime_error("feedForwardConstraintInput is unstable.");
	}
	catch(const std::exception& error)  {
		std::cerr << "what(): " << error.what() << " at time " << controllers_[mp_options_.nThreads_][i].time_[k] << " [sec]." << std::endl;
	}


	// lagrange multiplier calculation
	const size_t& nc1 = BASE::nc1TrajectoriesStock_[i][greatestLessTimeStampIndex];
	DmDagerFunc_[threadId].interpolate(time, DmDager, greatestLessTimeStampIndex);
	RmFunc_[threadId].interpolate(time, Rm, greatestLessTimeStampIndex);
	Eigen::MatrixXd DmDagerTransRm = DmDager.leftCols(nc1).transpose() * Rm;


	BASE::lagrangeControllerStock_[i].k_[k]   = DmDagerTransRm * (CmProjected - Lm);
	BASE::lagrangeControllerStock_[i].uff_[k] = -BASE::lagrangeControllerStock_[i].k_[k]*nominalOutput;
	Eigen::VectorXd localVff = DmDagerTransRm * (EvProjected-Lv-Lve);

	if (nominalLagrangeMultiplierUpdated_ == false || BASE::options_.lineSearchByMeritFuntion_==false)
	{
		BASE::lagrangeControllerStock_[i].uff_[k] += localVff;
		BASE::lagrangeControllerStock_[i].deltaUff_[k] = Eigen::VectorXd::Zero(nc1);
	}
	else
	{
		Eigen::VectorXd nominalLagrangeMultiplier;

		nominalLagrangeMultiplierFunc_[threadId].interpolate(time, nominalLagrangeMultiplier, greatestLessTimeStampIndex);

		BASE::lagrangeControllerStock_[i].uff_[k] += nominalLagrangeMultiplier;

		BASE::lagrangeControllerStock_[i].deltaUff_[k] = localVff - nominalLagrangeMultiplier;
	}

	// checking the numerical stability of the controller parameters
	try {
		if (BASE::lagrangeControllerStock_[i].k_[k] != BASE::lagrangeControllerStock_[i].k_[k])
			throw std::runtime_error("Feedback lagrangeMultiplier is unstable.");
		if (BASE::lagrangeControllerStock_[i].deltaUff_[k] != BASE::lagrangeControllerStock_[i].deltaUff_[k])
			throw std::runtime_error("Feedforward lagrangeMultiplier is unstable.");
	}
	catch(const std::exception& error)  {
		std::cerr << "what(): " << error.what() << " at time " << BASE::lagrangeControllerStock_[i].time_[k] << " [sec]." << std::endl;
	}

	return i;
}



template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::lineSearchWorker(size_t threadId)
{
	if(mp_options_.debugPrintMP_)
		printString("[Thread " + std::to_string(threadId) + "]: Starting lineSearchWorker. ");

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
	std::vector<size_t>               lsNc2FinalStock(NUM_SUBSYSTEMS);
	std::vector<constraint2_vector_t> lsHvFinalStock(NUM_SUBSYSTEMS);
	std::vector<std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd> >> lsLagrangeTrajectoriesStock(NUM_SUBSYSTEMS);


	while(true)
	{
		size_t alphaExp = alphaTaken_++;

		scalar_t learningRate  = BASE::options_.maxLearningRateGSLQP_ * std::pow(BASE::options_.lineSearchContractionRate_, alphaExp);


		if (learningRate < BASE::options_.minLearningRateGSLQP_ || alphaBestFound_.load() == true)
		{
			if(mp_options_.debugPrintMP_)
			{
				if (alphaBestFound_.load() == true)
					printString("[Thread " + std::to_string(threadId) + "]: Leaving lineSearchWorker because best alpha is found OR no improvement for any alpha");
				else
					printString("[Thread "+ std::to_string(threadId) +"]: Leaving lineSearchWorker because learningRate < BASE::options_.minLearningRateGSLQP_");
			}

			break;
		}

		if(mp_options_.debugPrintMP_)
			printString("[Thread " + std::to_string(threadId) + "]: Trying learningRate " + std::to_string(learningRate));


		lsControllersStock 			= initLScontrollersStock_;
		lsLagrangeControllersStock 	= initLSlagrangeMultiplierFunctionsStock_;


		executeLineSearch(
				threadId,
				learningRate,
				lsTotalCost,
				lsTotalMerit,
				lsConstraint1ISE,
				lsControllersStock,
				lsLagrangeControllersStock,
				lsTimeTrajectoriesStock,
				lsStateTrajectoriesStock,
				lsInputTrajectoriesStock,
				lsOutputTrajectoriesStock,
				lsNc1TrajectoriesStock,
				lsEvTrajectoryStock,
				lsNc2TrajectoriesStock,
				lsHvTrajectoryStock,
				lsNc2FinalStock,
				lsHvFinalStock,
				lsLagrangeTrajectoriesStock);


		// make sure we do not alter an existing result
		if (alphaBestFound_.load() == true)
		{
			if(mp_options_.debugPrintMP_)
				printString("[Thread " + std::to_string(threadId) + "]: Leaving lineSearchWorker because best alpha already found by another thread.");

			break;
		}


		lineSearchResultMutex_.lock();

		bool updatePolicy = false;
		if(mp_options_.lsStepsizeGreedy_ == true)
		{
			// act stepsize greedy, merit should be better than in last iteration but learning rate should be as high as possible
			if(lsTotalMerit <  (BASE::nominalTotalMerit_ * (1-1e-3*learningRate)) && learningRate > learningRateStar_)
			{
				updatePolicy = true;
				if(mp_options_.debugPrintMP_){
					printString("[LS, Thread " + std::to_string(threadId) + "]: stepsize-greedy mode : better stepsize and merit found: " + std::to_string(lsTotalMerit)
					+ " at learningRate: " + std::to_string(learningRate));
				}
			}
			else{
				if(mp_options_.debugPrintMP_){
					printString("[LS, Thread " + std::to_string(threadId) + "]: stepsize-greedy mode : no better combination found, merit " + std::to_string(lsTotalMerit)
					+ " at learningRate: " + std::to_string(learningRate));
				}
			}
		}
		else // line search acts merit greedy, minimize merit as much as possible
		{
			if(lsTotalMerit <  (lowestTotalMerit_ * (1-1e-3*learningRate)))
			{
				updatePolicy = true;
				if(mp_options_.debugPrintMP_){
					printString("[LS, Thread " + std::to_string(threadId) + "]: merit-greedy mode : better merit found: " + std::to_string(lsTotalMerit)
					+ " at learningRate: " + std::to_string(learningRate));
				}
			}
			else{
				if(mp_options_.debugPrintMP_){
					printString("[LS, Thread " + std::to_string(threadId) + "]: merit-greedy mode : no better merit found, merit " + std::to_string(lsTotalMerit)
					+ " at learningRate: " + std::to_string(learningRate) + ". Best merit was " + std::to_string(lowestTotalMerit_));
				}
			}
		}


		if (updatePolicy == true)
		{
			alphaExpBest_ 	  	= alphaExp;
			lowestTotalMerit_ 	= lsTotalMerit;
			lowestTotalCost_ 	= lsTotalCost;
			lowestConstraint1ISE_ = lsConstraint1ISE;
			learningRateStar_ 	= learningRate;

			for (size_t i = 0; i<NUM_SUBSYSTEMS; i++)	// swapping where possible for improved efficiency
			{
				controllers_[mp_options_.nThreads_][i].swap(lsControllersStock[i]);
				BASE::nominalTimeTrajectoriesStock_[i].swap(lsTimeTrajectoriesStock[i]);
				BASE::nominalStateTrajectoriesStock_[i].swap(lsStateTrajectoriesStock[i]);
				BASE::nominalInputTrajectoriesStock_[i].swap(lsInputTrajectoriesStock[i]);
				BASE::nominalOutputTrajectoriesStock_[i].swap(lsOutputTrajectoriesStock[i]);
				BASE::nc1TrajectoriesStock_[i].swap(lsNc1TrajectoriesStock[i]);
				BASE::EvTrajectoryStock_[i].swap(lsEvTrajectoryStock[i]);
				BASE::nc2TrajectoriesStock_[i].swap(lsNc2TrajectoriesStock[i]);
				BASE::HvTrajectoryStock_[i].swap(lsHvTrajectoryStock[i]);
				BASE::nc2FinalStock_[i] = lsNc2FinalStock[i];
				BASE::HvFinalStock_[i] = lsHvFinalStock[i];
				BASE::lagrangeControllerStock_[i].swap(lsLagrangeControllersStock[i]);;
				BASE::nominalLagrangeTrajectoriesStock_[i].swap(lsLagrangeTrajectoriesStock[i]);
			}
		}

		alphaProcessed_[alphaExp] = 1;

		// we now check if all alphas prior to the best have been processed, this also covers the case that there is no better alpha
		bool allPreviousAlphasProcessed = true;
		for (size_t i=0; i<alphaExpBest_; i++)
		{
			if (alphaProcessed_[i] != 1)
			{
				allPreviousAlphasProcessed = false;
				break;
			}
		}
		if (allPreviousAlphasProcessed)
		{
			alphaBestFound_ = true;
			killIntegrationEventHandler_->setEvent();	// kill all integrators
			if (BASE::options_.dispayGSLQP_) {
				printString("\t LS: terminate other rollouts with different alphas. alpha_best found or terminating without improvement. ");
			}
		}

		lineSearchResultMutex_.unlock();

	}

	lsWorkerCompleted_++;

	if(mp_options_.debugPrintMP_)
		printString("[Thread " + std::to_string(threadId) + "]: Leaving lineSearchWorker ");

	if (lsWorkerCompleted_.load() >= mp_options_.nThreads_)
	{
		std::unique_lock<std::mutex> lock (alphaBestFoundMutex_);
		alphaBestFoundCondition_.notify_all();
		lock.unlock();

		if(mp_options_.debugPrintMP_)
			printString("NOTIFYING by LS WORKER since all workers are now done ");
	}
}


template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::executeLineSearch(
		size_t threadId,
		double learningRate,
		scalar_t& lsTotalCost,
		scalar_t& lsTotalMerit,
		scalar_t& lsConstraint1ISE,
		std::vector<controller_t>& lsControllersStock,
		std::vector<lagrange_t>& lsLagrangeControllersStock,
		std::vector<scalar_array_t>& lsTimeTrajectoriesStock,
		std::vector<state_vector_array_t>& lsStateTrajectoriesStock,
		std::vector<control_vector_array_t>& lsInputTrajectoriesStock,
		std::vector<output_vector_array_t>& lsOutputTrajectoriesStock,
		std::vector<std::vector<size_t> >& lsNc1TrajectoriesStock,
		std::vector<constraint1_vector_array_t>& lsEvTrajectoryStock,
		std::vector<std::vector<size_t> >&   lsNc2TrajectoriesStock,
		std::vector<constraint2_vector_array_t>& lsHvTrajectoryStock,
		std::vector<size_t>&               	lsNc2FinalStock,
		std::vector<constraint2_vector_t>& 	lsHvFinalStock,
		std::vector<std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd> >>& lsLagrangeTrajectoriesStock
){


	// modifying uff by local increments
	for (int i=0; i<NUM_SUBSYSTEMS; i++)
		for (int k=0; k<lsControllersStock[i].time_.size(); k++)
			lsControllersStock[i].uff_[k] += learningRate * lsControllersStock[i].deltaUff_[k];


	// modifying vff by the local increments
	for (int i=0; i<NUM_SUBSYSTEMS; i++)
		for (int k=0; k < lsLagrangeControllersStock[i].time_.size(); k++)
			lsLagrangeControllersStock[i].uff_[k] += learningRate * lsLagrangeControllersStock[i].deltaUff_[k];


	try {
		rollout(threadId, BASE::initState_, lsControllersStock, lsTimeTrajectoriesStock,
				lsStateTrajectoriesStock, lsInputTrajectoriesStock, lsOutputTrajectoriesStock,
				lsNc1TrajectoriesStock, lsEvTrajectoryStock,
				lsNc2TrajectoriesStock, lsHvTrajectoryStock, lsNc2FinalStock, lsHvFinalStock);

		// calculate rollout cost
		calculateCostFunction(lsTimeTrajectoriesStock, lsOutputTrajectoriesStock, lsInputTrajectoriesStock,
				lsNc2TrajectoriesStock, lsHvTrajectoryStock, lsNc2FinalStock, lsHvFinalStock, lsTotalCost, threadId);

		// calculate the merit function
		if (BASE::options_.lineSearchByMeritFuntion_==true)
		{
			// calculate the lagrange multiplier
			calculateRolloutLagrangeMultiplier(lsTimeTrajectoriesStock, lsOutputTrajectoriesStock, lsLagrangeControllersStock, lsLagrangeTrajectoriesStock);
			// calculate the merit function value
			calculateMeritFunction(lsTimeTrajectoriesStock, lsNc1TrajectoriesStock, lsEvTrajectoryStock, lsLagrangeTrajectoriesStock, lsTotalCost, lsTotalMerit, lsConstraint1ISE);
		} else
		{
			lsTotalMerit = lsTotalCost;
			calculateConstraintISE(lsTimeTrajectoriesStock, lsNc1TrajectoriesStock, lsEvTrajectoryStock, lsConstraint1ISE);
		}

		// display
		if (BASE::options_.dispayGSLQP_){
			printString("\t [Thread" + std::to_string(threadId) + "] - learningRate " + std::to_string(learningRate) + " \t cost: " + std::to_string(lsTotalCost) +" \t merit: " + std::to_string(lsTotalMerit) +
					" \t constraint ISE: " + std::to_string(lsConstraint1ISE));
			if (std::accumulate(lsNc2FinalStock.begin(), lsNc2FinalStock.end(), 0) > 0) {
				std::cerr << "\t final constraint type-2:   ";
				for(size_t i=0; i<NUM_SUBSYSTEMS; i++) std::cerr << "[" << i  << "]: " << lsHvFinalStock[i].head(lsNc2FinalStock[i]).transpose() << ",  ";
				std::cerr << std::endl;
			}

		}
	}
	catch(const std::exception& error)
	{
		if(mp_options_.debugPrintMP_)
			std::cout << error.what() << std::endl;

		lsTotalMerit = std::numeric_limits<scalar_t>::max();
		lsTotalCost  = std::numeric_limits<scalar_t>::max();
	}
}



template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::setNewCostReferenceState(const output_vector_t& newReference){	// todo: move to implementation
	// for all threads + 1
	for (size_t i=0; i<mp_options_.nThreads_+1; i++)
	{
		// .. initialize all subsystems, etc.
		for(size_t j = 0; j<NUM_SUBSYSTEMS; j++)
		{
			costFunctions_[i][j]->updateReferenceState(newReference);
		}
	}
}



} // namespace ocs2
