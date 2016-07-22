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
		std::cout<<"Shutting down workers"<<std::endl;
	}

	for (size_t i=0; i<workerThreads_.size(); i++)
	{
		workerThreads_[i].join();
	}

	if(mp_options_.debugPrintMP_)
	{
		std::cout<<"All workers shut down"<<std::endl;
	}
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * Forward integrate the system dynamics with given controller:
 * 		inputs:
 * 			+ threadId: 		identifier of chosen thread
 * 			+ initState: 		initial state at time switchingTimes_[0]
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
		const size_t threadId,
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

		size_t maxNumSteps = options_.maxNumStepsPerSecond_*(switchingTimes_[i+1]-switchingTimes_[i]);
		maxNumSteps = ((1000>maxNumSteps) ? 1000 : maxNumSteps);

		// initialize subsystem i
		dynamics_[threadId][i]->initializeModel(switchingTimes_, x0, i, "GSLPQ");

		// set controller for subsystem i
		dynamics_[threadId][i]->setController(controller_local[i]);

		// simulate subsystem i
		integratorsODE45_[threadId][i]->integrate(
				x0, switchingTimes_[i], switchingTimes_[i+1],
				x_local[i], t_local[i],
				1e-3, options_.AbsTolODE_, options_.RelTolODE_, maxNumSteps);

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


/******************************************************************************************************/
// rolling out, where output trajectories are of no interest.
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::rollout(
		const size_t threadId,
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
		const size_t threadId,
		const state_vector_t& initState,
		const std::vector<controller_t>& controller_local,
		std::vector<scalar_array_t>& t_local,
		std::vector<state_vector_array_t>& x_local,
		std::vector<control_vector_array_t>& u_local,
		std::vector<output_vector_array_t>& y_local,
		std::vector<std::vector<size_t> >& nc1TrajectoriesStock,
		std::vector<constraint1_vector_array_t>& EvTrajectoryStock)  {

	// STEP1 : perform normal rollout
	rollout(threadId, initState, controller_local, t_local, x_local, u_local, y_local);


	// STEP2 : calculate constraint violations

	// constraint type 1 computations which consists of number of active constraints at each time point
	// and the value of the constraint (if the rollout is constrained the value is always zero otherwise
	// it is nonzero)
	nc1TrajectoriesStock.resize(NUM_SUBSYSTEMS);
	EvTrajectoryStock.resize(NUM_SUBSYSTEMS);

	for (int i=0; i<NUM_SUBSYSTEMS; i++)
	{
		size_t N = t_local[i].size();
		nc1TrajectoriesStock[i].resize(N);
		EvTrajectoryStock[i].resize(N);

		// compute constraint1 trajectory for subsystem i
		for (int k=0; k<N; k++)
		{
			dynamics_[threadId][i]->computeConstriant1(
					t_local[i][k], x_local[i][k], u_local[i][k],
					nc1TrajectoriesStock[i][k], EvTrajectoryStock[i][k]);

			if (nc1TrajectoriesStock[i][k] > INPUT_DIM)
				throw std::runtime_error("Number of active type-1 constraints should be less-equal to the number of input dimension.");
		}  // end of k loop
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
 *			+ threadId: working thread, defaults to the thread with lowest id, thus this is the default thread for single-core cost computation
 *				(allows to let method be called from the outside)
 * 		outputs:
 * 			+ totalCost: the total cost of the trajectory
 */
// TODO: shall we parallelize this method as well? Or does this not make sense because it only gets called from already parallel instances, e.g. in lineSearch?
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
	double pho = iteration_/(options_.maxIterationGSLQP_-1) * options_.meritFunctionRho_;
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
 * 		+ AmTrajectoryStock_: Am matrix
 * 		+ BmTrajectoryStock_: Bm matrix
 * 		+ CmTrajectoryStock_: Cm matrix
 * 		+ DmTrajectoryStock_: Dm matrix
 * 		+ EvTrajectoryStock_: Ev vector
 *
 * 		+ quadratized intermediate cost function
 * 		+ intermediate cost: q(t) + 0.5 y(t)Qm(t)y(t) + y(t)'Qv(t) + u(t)'Pm(t)y(t) + 0.5u(t)'Rm(t)u(t) + u(t)'Rv(t)
 * 		+ qTrajectoryStock_:  q
 * 		+ QvTrajectoryStock_: Qv vector
 * 		+ QmTrajectoryStock_: Qm matrix
 * 		+ PmTrajectoryStock_: Pm matrix
 * 		+ RvTrajectoryStock_: Rv vector
 * 		+ RmTrajectoryStock_: Rm matrix
 * 		+ RmInverseTrajectoryStock_: inverse of Rm matrix
 *
 * 		+ quadratized final cost in the last subsystem: qFinal(t) + 0.5 y(t)QmFinal(t)y(t) + y(t)'QvFinal(t)
 * 		+ qFinal_: qFinal
 * 		+ qFinal_: QvFinal vector
 * 		+ qFinal_: QmFinal matrix
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
		{
			std::cout << "[MP] Starting approximation of subsystem "<< i << " out of " << (size_t) NUM_SUBSYSTEMS-1 << std::endl;
		}
		auto start = std::chrono::steady_clock::now();

		approximateSubsystemLQ(i);

		auto end = std::chrono::steady_clock::now();
		auto diff = end - start;

		if(mp_options_.debugPrintMP_)
		{
			std::cout << "[MP] ended approximation of subsystem "<< i << std::endl;
			std::cout << "[MP] approximation took "<<std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;
		}
	}

	// todo: continue parallelization here
	// constrained coefficients
	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		int N = nominalTimeTrajectoriesStock_[i].size();

		RmConstrainedTrajectoryStock_[i].resize(N);
		DmDagerTrajectoryStock_[i].resize(N);

		AmConstrainedTrajectoryStock_[i].resize(N);
		QmConstrainedTrajectoryStock_[i].resize(N);
		QvConstrainedTrajectoryStock_[i].resize(N);

		EvProjectedTrajectoryStock_[i].resize(N);
		CmProjectedTrajectoryStock_[i].resize(N);
		DmProjectedTrajectoryStock_[i].resize(N);

		for (int k=0; k<N; k++) {
			size_t nc1 = nc1TrajectoriesStock_[i][k];

			if (nc1 == 0)
			{
				DmDagerTrajectoryStock_[i][k].setZero();
				EvProjectedTrajectoryStock_[i][k].setZero();
				CmProjectedTrajectoryStock_[i][k].setZero();
				DmProjectedTrajectoryStock_[i][k].setZero();

				AmConstrainedTrajectoryStock_[i][k] = AmTrajectoryStock_[i][k];
				QmConstrainedTrajectoryStock_[i][k] = QmTrajectoryStock_[i][k];
				QvConstrainedTrajectoryStock_[i][k] = QvTrajectoryStock_[i][k];
				RmConstrainedTrajectoryStock_[i][k] = RmTrajectoryStock_[i][k];

			}
			else
			{
				Eigen::MatrixXd Cm = CmTrajectoryStock_[i][k].topRows(nc1);
				Eigen::MatrixXd Dm = DmTrajectoryStock_[i][k].topRows(nc1);
				Eigen::MatrixXd Ev = EvTrajectoryStock_[i][k].head(nc1);

				Eigen::MatrixXd RmProjected = ( Dm*RmInverseTrajectoryStock_[i][k]*Dm.transpose() ).inverse();
				Eigen::MatrixXd DmDager = RmInverseTrajectoryStock_[i][k] * Dm.transpose() * RmProjected;

				DmDagerTrajectoryStock_[i][k].leftCols(nc1) = DmDager;
				EvProjectedTrajectoryStock_[i][k] = DmDager * Ev;
				CmProjectedTrajectoryStock_[i][k] = DmDager * Cm;
				DmProjectedTrajectoryStock_[i][k] = DmDager * Dm;

				control_matrix_t DmNullSpaceProjection = control_matrix_t::Identity() - DmProjectedTrajectoryStock_[i][k];
				state_matrix_t   PmTransDmDagerCm = PmTrajectoryStock_[i][k].transpose()*CmProjectedTrajectoryStock_[i][k];

				AmConstrainedTrajectoryStock_[i][k] = AmTrajectoryStock_[i][k] - BmTrajectoryStock_[i][k]*CmProjectedTrajectoryStock_[i][k];
				QmConstrainedTrajectoryStock_[i][k] = QmTrajectoryStock_[i][k] + Cm.transpose()*RmProjected*Cm - PmTransDmDagerCm - PmTransDmDagerCm.transpose();
				QvConstrainedTrajectoryStock_[i][k] = QvTrajectoryStock_[i][k] - CmProjectedTrajectoryStock_[i][k].transpose()*RvTrajectoryStock_[i][k];
				RmConstrainedTrajectoryStock_[i][k] = DmNullSpaceProjection.transpose() * RmTrajectoryStock_[i][k] * DmNullSpaceProjection;
			}

			// making sure that constrained Qm is PSD
			makePSD(QmConstrainedTrajectoryStock_[i][k]);

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

		SmFunc.setTimeStamp(&SsTimeTrajectoryStock_[i]);
		SmFunc.setData(&SmTrajectoryStock_[i]);
		SvFunc.setTimeStamp(&SsTimeTrajectoryStock_[i]);
		SvFunc.setData(&SvTrajectoryStock_[i]);
		nominalOutputFunc.setTimeStamp(&nominalTimeTrajectoriesStock_[i]);
		nominalOutputFunc.setData(&nominalOutputTrajectoriesStock_[i]);

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
	if (options_.dispayGSLQP_)
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
			if (options_.lineSearchByMeritFuntion_==true)
				controllers_[mp_options_.nThreads_][i].deltaUff_[k] += feedForwardConstraintInputStock_[i][k];
			else
				controllers_[mp_options_.nThreads_][i].uff_[k] += options_.constraintStepSize_*feedForwardConstraintInputStock_[i][k];


	// perform one rollout while the input correction for the type-1 constraint is considered.
	rollout(mp_options_.nThreads_, initState_, controllers_[mp_options_.nThreads_], nominalTimeTrajectoriesStock_,
			nominalStateTrajectoriesStock_, nominalInputTrajectoriesStock_, nominalOutputTrajectoriesStock_,
			nc1TrajectoriesStock_, EvTrajectoryStock_);

	calculateCostFunction( nominalTimeTrajectoriesStock_, nominalOutputTrajectoriesStock_, nominalInputTrajectoriesStock_, nominalTotalCost_, mp_options_.nThreads_);


	// calculate the merit function
	if (options_.lineSearchByMeritFuntion_==true)
	{
		// calculate the lagrange multiplier with learningRate zero
		calculateRolloutLagrangeMultiplier(nominalTimeTrajectoriesStock_, nominalOutputTrajectoriesStock_, lagrangeControllerStock_,
				nominalLagrangeTrajectoriesStock_);
		// calculate the merit function
		calculateMeritFunction(nominalTimeTrajectoriesStock_, nc1TrajectoriesStock_, EvTrajectoryStock_, nominalLagrangeTrajectoriesStock_, nominalTotalCost_,
				nominalTotalMerit_, nominalConstraint1ISE_);
	}
	else
	{
		nominalTotalMerit_ = nominalTotalCost_;
		calculateConstraintISE(nominalTimeTrajectoriesStock_, nc1TrajectoriesStock_, EvTrajectoryStock_, nominalConstraint1ISE_);
	}

	lowestTotalMerit_ = nominalTotalMerit_;


	// display
	if (options_.dispayGSLQP_)  {std::cerr << "\t learningRate 0.0 \t cost: " << nominalTotalCost_ << " \t merit: " << nominalTotalMerit_ <<
			" \t constraint ISE: " << nominalConstraint1ISE_ << std::endl;}


	initLScontrollersStock_ = controllers_[mp_options_.nThreads_];		// this will serve to init the workers
	initLSlagrangeMultiplierFunctionsStock_ = lagrangeControllerStock_;

	subsystemProcessed_ = 0; // not required for linesearch, but assign to not let it dangle around
	alphaProcessed_.clear();
	alphaTaken_ = 0;
	alphaBestFound_ = false;
	lsWorkerCompleted_ = 0;

	size_t maxNumOfLineSearches =  (int) (log(options_.minLearningRateGSLQP_/options_.maxLearningRateGSLQP_) / log(options_.lineSearchContractionRate_)) +1;
	alphaExpMax_ = maxNumOfLineSearches;
	alphaExpBest_ = maxNumOfLineSearches;
	alphaProcessed_.resize(maxNumOfLineSearches, 0);

	if(mp_options_.debugPrintMP_)
		{std::cout<<"[MP]: calculated maximum number of line searches "<< alphaExpMax_ <<std::endl;}

	if(mp_options_.debugPrintMP_)
		{std::cout << "[MP] Waking up workers for line search " << std::endl;}

	workerTask_ = LINE_SEARCH;
	workerWakeUpCondition_.notify_all();

	if(mp_options_.debugPrintMP_)
		{std::cout << "[MP] Will sleep now until we have results " << std::endl;}


	std::unique_lock<std::mutex> waitLock(alphaBestFoundMutex_);
	alphaBestFoundCondition_.wait(waitLock, [this]{return ( /*alphaBestFound_.load() &&*/ (lsWorkerCompleted_.load() >= mp_options_.nThreads_));});
	waitLock.unlock();

	workerTask_ = IDLE;

	if(mp_options_.debugPrintMP_) {std::cout << "[MP]: Woke up again, should have results now." << std::endl;}



	// clear the feedforward increments
	for (int j=0; j<mp_options_.nThreads_+1; j++){
		for (size_t i=0; i<NUM_SUBSYSTEMS; i++)
		{
			controllers_[mp_options_.nThreads_][i].deltaUff_.clear();
			lagrangeControllerStock_[i].deltaUff_.clear();
		}
	}

	// display
	if (options_.dispayGSLQP_)  {std::cerr << "The chosen learningRate is: " << learningRateStar_ << std::endl;}
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
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getValueFuntion(const scalar_t& time, const output_vector_t& output, scalar_t& valueFuntion)  {

	int activeSubsystem = -1;

	for (int i=0; i<NUM_SUBSYSTEMS; i++)
	{
		activeSubsystem = i;
		if (switchingTimes_[i]<=time && time<switchingTimes_[i+1])
			break;
	}

	state_matrix_t Sm;
	LinearInterpolation<state_matrix_t,Eigen::aligned_allocator<state_matrix_t> > SmFunc(
			&SsTimeTrajectoryStock_[activeSubsystem], &SmTrajectoryStock_[activeSubsystem]);
	SmFunc.interpolate(time, Sm);
	size_t greatestLessTimeStampIndex = SmFunc.getGreatestLessTimeStampIndex();

	output_vector_t Sv;
	LinearInterpolation<output_vector_t,Eigen::aligned_allocator<output_vector_t> > SvFunc(
			&SsTimeTrajectoryStock_[activeSubsystem], &SvTrajectoryStock_[activeSubsystem]);
	SvFunc.interpolate(time, Sv, greatestLessTimeStampIndex);

	eigen_scalar_t s;
	LinearInterpolation<eigen_scalar_t,Eigen::aligned_allocator<eigen_scalar_t> > sFunc(
			&SsTimeTrajectoryStock_[activeSubsystem], &sTrajectoryStock_[activeSubsystem]);
	sFunc.interpolate(time, s, greatestLessTimeStampIndex);

	output_vector_t xNominal;
	LinearInterpolation<output_vector_t,Eigen::aligned_allocator<output_vector_t> > xNominalFunc(&nominalTimeTrajectoriesStock_[activeSubsystem], &nominalOutputTrajectoriesStock_[activeSubsystem]);
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
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getCostFuntion(const output_vector_t& initOutput, scalar_t& costFunction, scalar_t& constriantCostFunction)  {

	const state_matrix_t&  Sm = SmTrajectoryStock_[0][0];
	const output_vector_t& Sv = SvTrajectoryStock_[0][0];
	const eigen_scalar_t&  s  = sTrajectoryStock_[0][0];

	output_vector_t deltaX = initOutput-nominalOutputTrajectoriesStock_[0][0];

	costFunction = (s + deltaX.transpose()*Sv + 0.5*deltaX.transpose()*Sm*deltaX).eval()(0);


	double pho = iteration_/(options_.maxIterationGSLQP_-1) * options_.meritFunctionRho_;
	constriantCostFunction = costFunction + pho*nominalConstraint1ISE_;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * get the optimal state, output, and input trajectories
 * 		output
 * 			+ nominalTimeTrajectoriesStock_: optimal time trajectory
 * 			+ nominalStateTrajectoriesStock_: optimal state trajectory
 * 			+ nominalInputTrajectoriesStock_: optimal control input trajectory
 * 			+ nominalOutputTrajectoriesStock_: optimal output trajectory
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getNominalTrajectories(std::vector<scalar_array_t>& nominalTimeTrajectoriesStock,
		std::vector<state_vector_array_t>& nominalStateTrajectoriesStock,
		std::vector<control_vector_array_t>& nominalInputTrajectoriesStock,
		std::vector<output_vector_array_t>& nominalOutputTrajectoriesStock)   {

	nominalTimeTrajectoriesStock   = nominalTimeTrajectoriesStock_;
	nominalStateTrajectoriesStock  = nominalStateTrajectoriesStock_;
	nominalInputTrajectoriesStock  = nominalInputTrajectoriesStock_;
	nominalOutputTrajectoriesStock = nominalOutputTrajectoriesStock_;
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
 * 			+ SsTimeTrajectoryStock_: time stamp
 * 			+ SmTrajectoryStock_: Sm matrix
 * 			+ SvTrajectoryStock_: Sv vector
 * 			+ SveTrajectoryStock_: Sve vector
 * 			+ sTrajectoryStock_: s scalar
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::solveSequentialRiccatiEquations(const scalar_t& learningRate)  {


	LinearInterpolation<state_matrix_t, Eigen::aligned_allocator<state_matrix_t> > SmFunc;

	// final value for the last Riccati equations
	typename RiccatiEquations_t::s_vector_t allSsFinal;
	RiccatiEquations_t::convert2Vector(QmFinal_, QvFinal_, qFinal_, allSsFinal);

	// final value for the last error equation
	output_vector_t SveFinal = output_vector_t::Zero();

	for (int i=NUM_SUBSYSTEMS-1; i>=0; i--) {

		// set data for Riccati equations
		auto riccatiEquationsPtr = std::make_shared<RiccatiEquations_t>();
		riccatiEquationsPtr->setData(learningRate, i, switchingTimes_[i], switchingTimes_[i+1],
				&nominalTimeTrajectoriesStock_[i],
				&AmConstrainedTrajectoryStock_[i], &BmTrajectoryStock_[i],
				&qTrajectoryStock_[i], &QvConstrainedTrajectoryStock_[i], &QmConstrainedTrajectoryStock_[i],
				&RvTrajectoryStock_[i], &RmInverseTrajectoryStock_[i], &RmConstrainedTrajectoryStock_[i], &PmTrajectoryStock_[i]);

		// integrating the Riccati equations
		ODE45<RiccatiEquations_t::S_DIM_> ode45(riccatiEquationsPtr);
		std::vector<double> normalizedTimeTrajectory;
		std::vector<typename RiccatiEquations_t::s_vector_t, Eigen::aligned_allocator<typename RiccatiEquations_t::s_vector_t> > allSsTrajectory;
		ode45.integrate(allSsFinal, i, i+1, allSsTrajectory, normalizedTimeTrajectory,
				1e-3, options_.AbsTolODE_, options_.RelTolODE_);

		// denormalizing time and constructing 'Sm', 'Sv', and 's'
		int N = normalizedTimeTrajectory.size();
		SsTimeTrajectoryStock_[i].resize(N);
		SmTrajectoryStock_[i].resize(N);
		SvTrajectoryStock_[i].resize(N);
		sTrajectoryStock_[i].resize(N);
		for (int k=0; k<N; k++) {
			RiccatiEquations_t::convert2Matrix(allSsTrajectory[N-1-k], SmTrajectoryStock_[i][k], SvTrajectoryStock_[i][k], sTrajectoryStock_[i][k]);
			SsTimeTrajectoryStock_[i][k] = (switchingTimes_[i]-switchingTimes_[i+1])*(normalizedTimeTrajectory[N-1-k]-i) + switchingTimes_[i+1];
		}  // end of k loop

		// testing the numerical stability of the Riccati equations
		for (int k=N-1; k>=0; k--) {
			try {
				if (SmTrajectoryStock_[i][k] != SmTrajectoryStock_[i][k])  throw std::runtime_error("Sm is unstable.");
				if (SvTrajectoryStock_[i][k] != SvTrajectoryStock_[i][k])  throw std::runtime_error("Sv is unstable.");
				if (sTrajectoryStock_[i][k] != sTrajectoryStock_[i][k])    throw std::runtime_error("s is unstable.");
			}
			catch(const std::exception& error)
			{
				std::cerr << "what(): " << error.what() << " at time " << SsTimeTrajectoryStock_[i][k] << " [sec]." << std::endl;
				for (int kp=k; kp<k+10; kp++)  {
					if (kp >= N) continue;
					std::cerr << "Sm[" << SsTimeTrajectoryStock_[i][kp] << "]:\n"<< SmTrajectoryStock_[i][kp].transpose() << std::endl;
					std::cerr << "Sv[" << SsTimeTrajectoryStock_[i][kp] << "]:\t"<< SvTrajectoryStock_[i][kp].transpose() << std::endl;
					std::cerr << "s[" << SsTimeTrajectoryStock_[i][kp] << "]: \t"<< sTrajectoryStock_[i][kp].transpose() << std::endl;
				}
				exit(1);
			}

		}  // end of k loop

		// set the final value for next Riccati equation
		allSsFinal = allSsTrajectory.back();

		/*
		 * Type_1 constraints error correction compensation
		 */

		// Skip calculation of the error correction term Sve if the constrained simulation is used for forward simulation
		if (options_.simulationIsConstrained_) {
			SveTrajectoryStock_[i].resize(N);
			for (int k=0; k<N; k++)
				SveTrajectoryStock_[i][k].setZero();
			continue;
		}

		// Calculating the coefficients of the error equation
		SmFunc.setTimeStamp( &(SsTimeTrajectoryStock_[i]) );
		SmFunc.setData( &(SmTrajectoryStock_[i]) );

		output_vector_array_t GvTrajectory(nominalTimeTrajectoriesStock_[i].size());
		state_matrix_array_t  GmTrajectory(nominalTimeTrajectoriesStock_[i].size());

		for (int k=0; k<nominalTimeTrajectoriesStock_[i].size(); k++) {
			state_matrix_t Sm;
			SmFunc.interpolate(nominalTimeTrajectoriesStock_[i][k], Sm);

			control_feedback_t Lm = RmInverseTrajectoryStock_[i][k]*(PmTrajectoryStock_[i][k]+BmTrajectoryStock_[i][k].transpose()*Sm);

			GmTrajectory[k] = AmConstrainedTrajectoryStock_[i][k] -
					BmTrajectoryStock_[i][k]*RmInverseTrajectoryStock_[i][k]*RmConstrainedTrajectoryStock_[i][k]*Lm;

			GvTrajectory[k] = (CmProjectedTrajectoryStock_[i][k]-DmProjectedTrajectoryStock_[i][k]*Lm).transpose()*
					RmTrajectoryStock_[i][k]*EvProjectedTrajectoryStock_[i][k];

		}  // end of k loop

		// set data for error equations
		auto errorEquationPtr = std::make_shared<ErrorEquation_t>();
		errorEquationPtr->setData(i, switchingTimes_[i], switchingTimes_[i+1],
				&nominalTimeTrajectoriesStock_[i], &GvTrajectory, &GmTrajectory);

		// integrating the Riccati equations
		ODE45<OUTPUT_DIM> errorOde45(errorEquationPtr);
		output_vector_array_t SveTrajectory;
		errorOde45.integrate(SveFinal, normalizedTimeTrajectory, SveTrajectory, 1e-3, options_.AbsTolODE_, options_.RelTolODE_);

		// reset the final value for next Riccati equation
		SveFinal = SveTrajectory.back();

		SveTrajectoryStock_[i].resize(N);
		for (int k=0; k<N; k++) {
			SveTrajectoryStock_[i][k] = SveTrajectory[N-1-k];

			// testing the numerical stability of the Riccati error equation
			try {
				if (SveTrajectoryStock_[i][k] != SveTrajectoryStock_[i][k])  throw std::runtime_error("Sve is unstable");
			}
			catch(const std::exception& error) 	{
				std::cerr << "what(): " << error.what() << " at time " << SsTimeTrajectoryStock_[i][k] << " [sec]." << std::endl;
				for (int kp=k; kp<N; kp++)   std::cerr << "Sve[" << SsTimeTrajectoryStock_[i][kp] << "]:\t"<< SveTrajectoryStock_[i][kp].transpose() << std::endl;
				exit(1);
			}
		}

	}  // end of i loop
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * make the given square matrix psd
 */
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
template <typename Derived>
bool SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::makePSD(Eigen::MatrixBase<Derived>& squareMatrix) {

	if (squareMatrix.rows() != squareMatrix.cols())  throw std::runtime_error("Not a square matrix: makePSD() method is for square matrix.");

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
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::run(const state_vector_t& initState, const std::vector<scalar_t>& switchingTimes) {

	if (switchingTimes.size() != NUM_SUBSYSTEMS+1)
		throw std::runtime_error("Number of switching times should be one plus the number of subsystems.");

	switchingTimes_ = switchingTimes;
	initState_ = initState;

	// display
	if (options_.dispayGSLQP_)
	{
		std::cerr << "\n#### SLQP solver starts with switching times [" << switchingTimes[0];
		for (size_t i=1; i<=NUM_SUBSYSTEMS; i++)   std::cerr << ", " << switchingTimes[i];
		std::cerr << "] ..." << std::endl << std::endl;
	}

	iteration_ = 0;
	double relCost;
	double relConstraint1ISE;
	bool isConstraint1Satisfied  = false;
	bool isCostFunctionConverged = false;
	bool isOptimizationConverged = false;
	nominalLagrangeMultiplierUpdated_ = false;

	// initial controller rollout
	rollout(mp_options_.nThreads_,
			initState_,
			controllers_[mp_options_.nThreads_],
			nominalTimeTrajectoriesStock_, nominalStateTrajectoriesStock_, nominalInputTrajectoriesStock_, nominalOutputTrajectoriesStock_,
			nc1TrajectoriesStock_, EvTrajectoryStock_);

	// initial controller cost
	calculateCostFunction(nominalTimeTrajectoriesStock_, nominalOutputTrajectoriesStock_, nominalInputTrajectoriesStock_, nominalTotalCost_, mp_options_.nThreads_); // todo: parallelize ?

	// initial controller merit
	nominalTotalMerit_ = nominalTotalCost_;

	// initial controller constraint type-1 ISE
	calculateConstraintISE(nominalTimeTrajectoriesStock_, nc1TrajectoriesStock_, EvTrajectoryStock_, nominalConstraint1ISE_);

	// display
	if (options_.dispayGSLQP_)  std::cerr << "\n#### Initial controller: \n cost: " << nominalTotalCost_ << " \t constraint ISE: " << nominalConstraint1ISE_ << std::endl;

	// SLQP main loop
	while (iteration_<options_.maxIterationGSLQP_ && isOptimizationConverged==false)  {

		double costCashed = nominalTotalCost_;
		double constraint1ISECashed = nominalConstraint1ISE_;

		// display
		if (options_.dispayGSLQP_)  std::cerr << "\n#### Iteration " <<  iteration_ << std::endl;

		// linearizing the dynamics and quadratizing the cost function along nominal trajectories
		approximateOptimalControlProblem(); // todo parallelize, working here

		// solve Riccati equations
		solveSequentialRiccatiEquations(1.0 /*nominal learningRate*/); //todo do not parallize

		calculateControllerAndLagrangian();
		nominalLagrangeMultiplierUpdated_ = true;

		// finding the optimal learningRate
		lineSearch();


		// calculates type-1 constraint ISE and maximum norm
		double constraint1MaxNorm = calculateConstraintISE(nominalTimeTrajectoriesStock_, nc1TrajectoriesStock_, EvTrajectoryStock_, nominalConstraint1ISE_);

		// loop variables
		iteration_++;
		relCost = fabs(nominalTotalCost_-costCashed);
		relConstraint1ISE = fabs(nominalConstraint1ISE_-constraint1ISECashed);
		isConstraint1Satisfied  = nominalConstraint1ISE_<=options_.minAbsConstraint1ISE_ || relConstraint1ISE<=options_.minRelConstraint1ISE_;
		isCostFunctionConverged = learningRateStar_==0 || relCost<=options_.minRelCostGSLQP_;
		isOptimizationConverged = isCostFunctionConverged==true && isConstraint1Satisfied==true;

		// display
		if (options_.dispayGSLQP_)  {
			std::cerr << "optimization cost:  " << nominalTotalCost_ << std::endl;
			std::cerr << "constraint ISE:     " << nominalConstraint1ISE_ << std::endl;
			std::cerr << "constraint MaxNorm: " << constraint1MaxNorm << std::endl;
		}
	}  // end of while loop

	// linearizing the dynamics and quadratizing the cost function along nominal trajectories
	approximateOptimalControlProblem();

	// solve Riccati equations with learningRate zero
	solveSequentialRiccatiEquations(0.0 /*learningRate*/);

	// calculate the nominal co-state
	calculateRolloutCostate(nominalTimeTrajectoriesStock_, nominalOutputTrajectoriesStock_, nominalcostateTrajectoriesStock_);

	// display
	if (options_.dispayGSLQP_ )  {
		std::cout << "\n+++++++++++++++++++++++++++++++++++" << std::endl;
		std::cout <<   "+++++++ SLQP solver is ended ++++++" << std::endl;
		std::cout <<   "+++++++++++++++++++++++++++++++++++" << std::endl;
		if (isOptimizationConverged) {
			if (learningRateStar_==0)
				std::cerr << "SLQP successfully terminates as learningRate reduced to zero." << std::endl;
			else
				std::cerr << "SLQP successfully terminates as cost relative change (relCost=" << relCost <<") reached to the minimum value." << std::endl;

			if (nominalConstraint1ISE_<=options_.minAbsConstraint1ISE_)
				std::cerr << "Type-1 constraint absolute ISE (absConstraint1ISE=" << nominalConstraint1ISE_ << ") reached to the minimum value." << std::endl;
			else
				std::cerr << "Type-1 constraint relative ISE (relConstraint1ISE=" << relConstraint1ISE << ") reached to the minimum value." << std::endl;
		} else
			std::cerr << "Maximum number of iterations has reached." << std::endl;
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
	if(mp_options_.debugPrintMP_)	std::cout<<"[Thread "<<threadId<<"]: launched"<<std::endl;

	size_t uniqueProcessID = 0;

	while(workersActive_)
	{
		if(mp_options_.debugPrintMP_){
			std::cout<<"[Thread " << threadId << "]: previous procId: " <<  uniqueProcessID << std::endl;
			std::cout<<"[Thread " << threadId << "]: current procId : " <<  generateUniqueProcessID(iteration_, workerTask_, subsystemProcessed_) << std::endl;
		}


		/* We want to put the worker to sleep if
		 * - the workerTask_ is IDLE
		 * - or we are finished both workerTask_ is not yet reset, thus the process ID is still the same
		 * */
		if ( workerTask_ == IDLE || uniqueProcessID == generateUniqueProcessID(iteration_, workerTask_, subsystemProcessed_))
		{
			if(mp_options_.debugPrintMP_){
				std::string output;	output = "[Thread " + std::to_string(threadId) + "]: going to sleep !";
				std::cout << output << std::endl;
			}

			std::unique_lock<std::mutex> waitLock(workerWakeUpMutex_);

			// sleep until the state is not IDLE any more and we have a different process ID than before

			workerWakeUpCondition_.wait(waitLock, [this , uniqueProcessID]{
				return (workerTask_ != IDLE &&  (uniqueProcessID != generateUniqueProcessID(iteration_, workerTask_, subsystemProcessed_) ) );
			});

			waitLock.unlock();
		}

		if(mp_options_.debugPrintMP_){
			std::string output;	output = "[Thread " + std::to_string(threadId) + "]: woke up !";
			std::cout << output << std::endl;		}

		if (!workersActive_)
			break;

		switch(workerTask_)
		{
		case APPROXIMATE_LQ:
		{
			if(mp_options_.debugPrintMP_){
				std::string output;	output = "[Thread " + std::to_string(threadId) + "]: now busy with APPROXIMATE_LQ !";
				std::cout << output << std::endl;
			}

			size_t currentIteration = iteration_;
			size_t currentlyProcessedSubsys = approximateSubsystemLQWorker(threadId);
			uniqueProcessID = generateUniqueProcessID (currentIteration, APPROXIMATE_LQ, currentlyProcessedSubsys);

			break;
		}
		case CALCULATE_CONTROLLER_AND_LAGRANGIAN:
		{
			if(mp_options_.debugPrintMP_){
				std::string output;	output = "[Thread " + std::to_string(threadId) + "]: now busy with CALCULATE_CONTROLLER_AND_LAGRANGIAN !";
				std::cout << output << std::endl;
			}

			size_t currentIteration = iteration_;
			size_t currentlyProcessedSubsys = calculateControllerAndLagrangianWorker(threadId);
			uniqueProcessID = generateUniqueProcessID (currentIteration, CALCULATE_CONTROLLER_AND_LAGRANGIAN, currentlyProcessedSubsys);
			break;
		}
		case LINE_SEARCH:
		{
			if(mp_options_.debugPrintMP_) {
				std::string output;	output = "[Thread " + std::to_string(threadId) + "]: now busy with LINE_SEARCH !";
				std::cout << output << std::endl;}

			size_t currentIteration = iteration_;
			lineSearchWorker(threadId);
			uniqueProcessID = generateUniqueProcessID (currentIteration, LINE_SEARCH, subsystemProcessed_);
			break;
		}
		case SHUTDOWN:
		{
			if(mp_options_.debugPrintMP_){
				std::cout<<"[Thread "<<threadId<<"]: now shutting down!"<<std::endl;}
			return;
			break;
		}
		}

		if(mp_options_.debugPrintMP_){
			std::cout<<"[Thread "<<threadId<<"]: done with job. Will wait for next now!"<<std::endl;}
	}
}


template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::approximateSubsystemLQ(const size_t i)
{
	kTaken_ = 0;
	kCompleted_= 0;
	KMax_subsystem_approx_  = nominalTimeTrajectoriesStock_[i].size(); // number of elements in the trajectory of this subsystem
	// todo: no atomic variable

	// initialize subsystem i dynamics derivatives
	for(size_t j = 0; j< mp_options_.nThreads_; j++)
	{
		assert( nominalTimeTrajectoriesStock_[i].size() == nominalStateTrajectoriesStock_[i].size());
		linearizedSystems_[j][i]->initializeModel(switchingTimes_, nominalStateTrajectoriesStock_[i].front(), i, "GSLPQ");
	}

	AmTrajectoryStock_[i].resize(KMax_subsystem_approx_);
	BmTrajectoryStock_[i].resize(KMax_subsystem_approx_);
	CmTrajectoryStock_[i].resize(KMax_subsystem_approx_);
	DmTrajectoryStock_[i].resize(KMax_subsystem_approx_);

	qTrajectoryStock_[i].resize(KMax_subsystem_approx_);
	QvTrajectoryStock_[i].resize(KMax_subsystem_approx_);
	QmTrajectoryStock_[i].resize(KMax_subsystem_approx_);
	RvTrajectoryStock_[i].resize(KMax_subsystem_approx_);
	RmTrajectoryStock_[i].resize(KMax_subsystem_approx_);
	RmInverseTrajectoryStock_[i].resize(KMax_subsystem_approx_);
	PmTrajectoryStock_[i].resize(KMax_subsystem_approx_);


	if(mp_options_.debugPrintMP_){
		std::string output;	output = "[MP]: Waking up workers to do linearisation for subsystem " + std::to_string(i);
		std::cout << output << std::endl;
	}

	workerTask_ = APPROXIMATE_LQ;
	workerWakeUpCondition_.notify_all();

	if(mp_options_.debugPrintMP_){
		std::string output;	output = "[MP]: Will wait now until workers have linearized dynamics of subsystem " + std::to_string(i);
		std::cout << output << std::endl;
	}

	std::unique_lock<std::mutex> waitLock(kCompletedMutex_);
	// with timeout 20 000 milliseconds

	if(kCompletedCondition_.wait_for(waitLock, std::chrono::milliseconds(20000),[this]{return (kCompleted_.load() >= KMax_subsystem_approx_) ;}))
	{
		if(mp_options_.debugPrintMP_){
			std::string output;	output = "[MP]: Back to main thread, workers should now have linearized dynamics of subsystem " + std::to_string(i);
			std::cout << output << std::endl;
		}
	}
	else
	{
		throw(std::runtime_error("[MP FATAL]: Back to main thread, because TIMEOUT occured. Designing controller should not take longer than 20 sec."));
	}

	waitLock.unlock();

	workerTask_ = IDLE;


	if (i==NUM_SUBSYSTEMS-1) // if last subsystem, set terminal cost
	{
		if(mp_options_.debugPrintMP_)	std::cout<<"[MP]: Approximating terminal cost with single thread, subsystem  "<< i << std::endl;

		costFunctions_[mp_options_.nThreads_][i]->setCurrentStateAndControl(nominalTimeTrajectoriesStock_[i].back(), nominalOutputTrajectoriesStock_[i].back(), nominalInputTrajectoriesStock_[i].back());

		costFunctions_[mp_options_.nThreads_][i]->terminalCost(qFinal_(0));
		costFunctions_[mp_options_.nThreads_][i]->terminalCostStateDerivative(QvFinal_);
		costFunctions_[mp_options_.nThreads_][i]->terminalCostStateSecondDerivative(QmFinal_);

		// making sure that Qm remains PSD
		makePSD(QmFinal_);
	}
}



template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateControllerAndLagrangian() {

	for (int i=0; i<NUM_SUBSYSTEMS; i++)
	{
		subsystemProcessed_ =  i;

		kTaken_ = 0;
		kCompleted_ = 0;
		KMax_subsystem_ctrl_  = SsTimeTrajectoryStock_[i].size();; // number of elements in the trajectory of this subsystem

		// initialize interpolators
		for(size_t n = 0; n< mp_options_.nThreads_+1; n++)
		{
			// functions for controller and lagrange-multiplier
			nominalOutputFunc_[n].setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
			nominalOutputFunc_[n].setData( &(nominalOutputTrajectoriesStock_[i]) );

			nominalInputFunc_[n].setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
			nominalInputFunc_[n].setData( &(nominalInputTrajectoriesStock_[i]) );

			BmFunc_[n].setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
			BmFunc_[n].setData( &(BmTrajectoryStock_[i]) );

			PmFunc_[n].setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
			PmFunc_[n].setData( &(PmTrajectoryStock_[i]) );

			RmInverseFunc_[n].setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
			RmInverseFunc_[n].setData( &(RmInverseTrajectoryStock_[i]) );

			RvFunc_[n].setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
			RvFunc_[n].setData( &(RvTrajectoryStock_[i]) );

			EvProjectedFunc_[n].setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
			EvProjectedFunc_[n].setData( &(EvProjectedTrajectoryStock_[i]) );

			CmProjectedFunc_[n].setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
			CmProjectedFunc_[n].setData( &(CmProjectedTrajectoryStock_[i]) );

			DmProjectedFunc_[n].setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
			DmProjectedFunc_[n].setData( &(DmProjectedTrajectoryStock_[i]) );

			// functions for lagrange multiplier only
			RmFunc_[n].setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
			RmFunc_[n].setData( &(RmTrajectoryStock_[i]) );

			DmDagerFunc_[n].setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
			DmDagerFunc_[n].setData( &(DmDagerTrajectoryStock_[i]) );

			if (nominalLagrangeMultiplierUpdated_ == true) {
				nominalLagrangeMultiplierFunc_[n].setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
				nominalLagrangeMultiplierFunc_[n].setData( &(nominalLagrangeTrajectoriesStock_[i]) );
			}
		}

		controllers_[mp_options_.nThreads_][i].time_ = SsTimeTrajectoryStock_[i];
		controllers_[mp_options_.nThreads_][i].k_.resize(KMax_subsystem_ctrl_);
		controllers_[mp_options_.nThreads_][i].uff_.resize(KMax_subsystem_ctrl_);
		controllers_[mp_options_.nThreads_][i].deltaUff_.resize(KMax_subsystem_ctrl_);

		feedForwardConstraintInputStock_[i].resize(KMax_subsystem_ctrl_);

		lagrangeControllerStock_[i].time_ = SsTimeTrajectoryStock_[i];
		lagrangeControllerStock_[i].k_.resize(KMax_subsystem_ctrl_);
		lagrangeControllerStock_[i].uff_.resize(KMax_subsystem_ctrl_);
		lagrangeControllerStock_[i].deltaUff_.resize(KMax_subsystem_ctrl_);


		if(mp_options_.debugPrintMP_)
		{
			std::cout<<"[MP]: Waking up workers to calc. controller for subsystem "<< i << std::endl;
		}

		workerTask_ = CALCULATE_CONTROLLER_AND_LAGRANGIAN;
		workerWakeUpCondition_.notify_all();

		if(mp_options_.debugPrintMP_)
		{
			std::cout<<"[MP]: Will wait now until workers have calculated controller for subsystem " << i <<std::endl;
		}

		std::unique_lock<std::mutex> waitLock(kCompletedMutex_);
		// .. with timeout for 20 000 milliseconds
		if (kCompletedCondition_.wait_for(waitLock,  std::chrono::milliseconds(20000), [this]{return kCompleted_.load() >= KMax_subsystem_ctrl_ ;})) // problem: this does sometimes not get notified
		{
			if(mp_options_.debugPrintMP_) std::cout<<"[MP]: Back to main thread, workers should now have designed controllers for subsystem " << i <<std::endl;
		}
		else
			throw(std::runtime_error("[MP FATAL]: Back to main thread, because TIMEOUT occured. LQ Approximating system should not take longer than 20 sec"));

		waitLock.unlock();

		workerTask_ = IDLE;


	}  // end of i loop

}

template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
size_t SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::approximateSubsystemLQWorker(size_t threadId)
{
	size_t subsystem = subsystemProcessed_;

	while(true)
	{
		size_t k = kTaken_++;


//		if (k >= KMax_subsystem_approx_) // if all k's are already covered, notify and return
//		{
//			if(kCompleted_.load() >=KMax_subsystem_approx_)
//			{
//				kCompletedCondition_.notify_all();
//				std::string test;	test = "thread " + std::to_string(threadId) + " leaving AND NOTIFYING "; // todo: fix everywhere
//				std::cout << test << std::endl;
//			}
//			else
//				std::cout << "thread " << threadId << " leaving but not notifying " << std::endl;
//
//			return subsystem; // todo: may be wrong value here
//		}

		if(k < KMax_subsystem_approx_){

			if(mp_options_.debugPrintMP_){
	//			if (k%10 == 0) {
					std::string output;	output = "[Thread " + std::to_string(threadId) + "]: Start approximating system LQ on index k = " + std::to_string(k) +
					" out of " + std::to_string(KMax_subsystem_approx_-1);
					std::cout << output << std::endl;
	//			}
			}

			subsystem = executeApproximateSubsystemLQ(threadId, k);
			kCompleted_++;
		}

		if (k >= KMax_subsystem_approx_-1) // if all k's are already covered, notify and return
		{
			if(kCompleted_.load() >=KMax_subsystem_approx_)
			{
				kCompletedCondition_.notify_all();
				if(mp_options_.debugPrintMP_){
					std::string output;	output = "[Thread " + std::to_string(threadId) + "]: leaving AND NOTIFYING "; // todo: fix everywhere
					std::cout << output << std::endl;
				}
			}
			else{
				{
					if(mp_options_.debugPrintMP_){
					std::string output;	output = "[Thread " + std::to_string(threadId) + "]: leaving but NOT notifying "; // todo: fix everywhere
					std::cout << output << std::endl;
					}
				}
			}

			return subsystem; // todo: may be wrong value here
		}
	}

	return subsystem;
}


template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
size_t SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateControllerAndLagrangianWorker(size_t threadId)
{

	size_t subsystem = subsystemProcessed_;

	while(true)
	{
		size_t k = kTaken_++;

		if(k < KMax_subsystem_ctrl_){

			if(mp_options_.debugPrintMP_){
				//			if (k%10 == 0) {
				std::string output;	output = "[Thread " + std::to_string(threadId) + "]: Start calculating controller on index k = " + std::to_string(k) +
				" out of " + std::to_string(KMax_subsystem_ctrl_-1);
				std::cout << output << std::endl;
				//			}
			}

			subsystem = executeCalculateControllerAndLagrangian(threadId, k);
			kCompleted_++;
		}


		if (k >= KMax_subsystem_ctrl_-1)	// if all k's are already covered, notify and return
		{
			if(kCompleted_.load()>=KMax_subsystem_ctrl_){
				kCompletedCondition_.notify_all();
				if(mp_options_.debugPrintMP_){
					std::string output;	output = "[Thread " + std::to_string(threadId) + "]: leaving AND NOTIFYING ";
					std::cout << output << std::endl;
				}
			}else
			{
				{
					if(mp_options_.debugPrintMP_){
					std::string output;	output = "[Thread " + std::to_string(threadId) + "]: leaving but NOT notifying ";
					std::cout << output << std::endl;
					}
				}
			}

			return subsystem;
		}
	}

	return subsystem;
}


template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
size_t SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::executeApproximateSubsystemLQ(size_t threadId, size_t k)
{
	const size_t i = subsystemProcessed_;

	// LINEARIZE SYSTEM DYNAMICS AND CONSTRAINTS
	linearizedSystems_[threadId][i]->setCurrentStateAndControl(
			nominalTimeTrajectoriesStock_[i][k],
			nominalStateTrajectoriesStock_[i][k],
			nominalInputTrajectoriesStock_[i][k],
			nominalOutputTrajectoriesStock_[i][k]);

	linearizedSystems_[threadId][i]->getDerivativeState(AmTrajectoryStock_[i][k]);
	linearizedSystems_[threadId][i]->getDerivativesControl(BmTrajectoryStock_[i][k]);

	// if constraint type 1 is active
	if (nc1TrajectoriesStock_[i][k] > 0)
	{
		linearizedSystems_[threadId][i]->getConstraint1DerivativesState(CmTrajectoryStock_[i][k]);
		linearizedSystems_[threadId][i]->getConstraint1DerivativesControl(DmTrajectoryStock_[i][k]);
	}

	// QUADRATIC APPROXIMATION TO THE COST FUNCTION
	costFunctions_[threadId][i]->setCurrentStateAndControl(
			nominalTimeTrajectoriesStock_[i][k],
			nominalOutputTrajectoriesStock_[i][k],
			nominalInputTrajectoriesStock_[i][k]);
	costFunctions_[threadId][i]->evaluate(qTrajectoryStock_[i][k](0));
	costFunctions_[threadId][i]->stateDerivative(QvTrajectoryStock_[i][k]);
	costFunctions_[threadId][i]->stateSecondDerivative(QmTrajectoryStock_[i][k]);
	costFunctions_[threadId][i]->controlDerivative(RvTrajectoryStock_[i][k]);
	costFunctions_[threadId][i]->controlSecondDerivative(RmTrajectoryStock_[i][k]);
	RmInverseTrajectoryStock_[i][k] = RmTrajectoryStock_[i][k].inverse();
	costFunctions_[threadId][i]->stateControlDerivative(PmTrajectoryStock_[i][k]);

	return i;
}


template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
size_t SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::executeCalculateControllerAndLagrangian(size_t threadId, size_t k)
{

	const size_t i = subsystemProcessed_;

	const double time = SsTimeTrajectoryStock_[i][k];
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

	control_feedback_t Lm  = RmInverse * (Pm + Bm.transpose()*SmTrajectoryStock_[i][k]);
	control_vector_t   Lv  = RmInverse * (Rv + Bm.transpose()*SvTrajectoryStock_[i][k]);
	control_vector_t   Lve = RmInverse * (Bm.transpose()*SveTrajectoryStock_[i][k]);

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
	const size_t& nc1 = nc1TrajectoriesStock_[i][greatestLessTimeStampIndex];
	DmDagerFunc_[threadId].interpolate(time, DmDager, greatestLessTimeStampIndex);
	RmFunc_[threadId].interpolate(time, Rm, greatestLessTimeStampIndex);
	Eigen::MatrixXd DmDagerTransRm = DmDager.leftCols(nc1).transpose() * Rm;


	lagrangeControllerStock_[i].k_[k]   = DmDagerTransRm * (CmProjected - Lm);
	lagrangeControllerStock_[i].uff_[k] = -lagrangeControllerStock_[i].k_[k]*nominalOutput;
	Eigen::VectorXd localVff = DmDagerTransRm * (EvProjected-Lv-Lve);

	if (nominalLagrangeMultiplierUpdated_ == false || options_.lineSearchByMeritFuntion_==false)
	{
		lagrangeControllerStock_[i].uff_[k] += localVff;
		lagrangeControllerStock_[i].deltaUff_[k] = Eigen::VectorXd::Zero(nc1);
	}
	else
	{
		Eigen::VectorXd nominalLagrangeMultiplier;

		nominalLagrangeMultiplierFunc_[threadId].interpolate(time, nominalLagrangeMultiplier, greatestLessTimeStampIndex);

		lagrangeControllerStock_[i].uff_[k] += nominalLagrangeMultiplier;

		lagrangeControllerStock_[i].deltaUff_[k] = localVff - nominalLagrangeMultiplier;
	}

	// checking the numerical stability of the controller parameters
	try {
		if (lagrangeControllerStock_[i].k_[k] != lagrangeControllerStock_[i].k_[k])
			throw std::runtime_error("Feedback lagrangeMultiplier is unstable.");
		if (lagrangeControllerStock_[i].deltaUff_[k] != lagrangeControllerStock_[i].deltaUff_[k])
			throw std::runtime_error("Feedforward lagrangeMultiplier is unstable.");
	}
	catch(const std::exception& error)  {
		std::cerr << "what(): " << error.what() << " at time " << lagrangeControllerStock_[i].time_[k] << " [sec]." << std::endl;
	}

	return i;
}



template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP_MP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::lineSearchWorker(size_t threadId)
{
	if(mp_options_.debugPrintMP_)
		{std::cout<<"[Thread "<<threadId<<"]: Starting lineSearchWorker " <<std::endl;}


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
	std::vector<std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd> >> lsLagrangeTrajectoriesStock(NUM_SUBSYSTEMS);


	while(true)
	{
		size_t alphaExp = alphaTaken_++;

		scalar_t learningRate  = options_.maxLearningRateGSLQP_ * std::pow(options_.lineSearchContractionRate_, alphaExp);


		if (learningRate < options_.minLearningRateGSLQP_ || alphaBestFound_.load() == true)
		{
			if(mp_options_.debugPrintMP_)
			{
				if (alphaBestFound_.load() == true)
					std::cout<<"[Thread "<<threadId<<"]: Leaving lineSearchWorker because best alpha is found ..." <<std::endl;
				else
					std::cout<<"[Thread "<<threadId<<"]: Leaving lineSearchWorker because learningRate < options_.minLearningRateGSLQP_" <<std::endl;
			}

			break;
		}

		if(mp_options_.debugPrintMP_){
			std::cout<<"[Thread "<<threadId<<"]: Trying learningRate " << learningRate <<std::endl;}


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
				lsLagrangeTrajectoriesStock);

		// wait for the thread with a lower alphaExp index to finish first. Why?
		// observations tell us it is better to use a bigger stepsize!
		while(std::accumulate(alphaProcessed_.begin(), std::next(alphaProcessed_.begin(), alphaExp), 0) < alphaExp && alphaBestFound_.load() == false)
		{
			std::chrono::milliseconds dura( static_cast<int>(5) );	// Sleep 5ms until we check again. TODO: which time should we select?
			std::this_thread::sleep_for( dura );
		}

		// make sure we do not alter an existing result
		if (alphaBestFound_.load() == true)
		{
			if(mp_options_.debugPrintMP_){
				std::cout<<"[Thread "<<threadId<<"]: Leaving lineSearchWorker because best alpha already found by another thread." <<std::endl;}

			break;
		}


		lineSearchResultMutex_.lock();
		if (lsTotalMerit <  (nominalTotalMerit_ * (1-1e-3*learningRate)))	// that is a common criterion in optimization
		{
			if(mp_options_.debugPrintMP_){
				std::cout<<"[LineSearch, Thread "<<threadId<<"]: Lower cost found: "<<lsTotalMerit<<" at learningRate: "<<learningRate<<std::endl;}

			alphaExpBest_ 	  = alphaExp;
			lowestTotalMerit_ = lsTotalMerit;

			nominalTotalCost_      				= lsTotalCost;
			nominalTotalMerit_     				= lsTotalMerit;
			nominalConstraint1ISE_ 				= lsConstraint1ISE;
			controllers_[mp_options_.nThreads_] = lsControllersStock;
			nominalTimeTrajectoriesStock_  		= lsTimeTrajectoriesStock;
			nominalStateTrajectoriesStock_  	= lsStateTrajectoriesStock;		// todo: swap to save time
			nominalInputTrajectoriesStock_  	= lsInputTrajectoriesStock;
			nominalOutputTrajectoriesStock_ 	= lsOutputTrajectoriesStock;
			nc1TrajectoriesStock_ 				= lsNc1TrajectoriesStock;
			EvTrajectoryStock_ 					= lsEvTrajectoryStock;
			learningRateStar_ 					= learningRate;
			lagrangeControllerStock_ 			= lsLagrangeControllersStock;
			nominalLagrangeTrajectoriesStock_ 	= lsLagrangeTrajectoriesStock;

		}
		else
		{
			if(mp_options_.debugPrintMP_)
				std::cout<<"[LineSearch, Thread "<<threadId<<"]: No lower cost found, cost "<<lsTotalMerit<<" at learningRate "<<learningRate<<" . Best cost was "<<lowestTotalMerit_ <<std::endl;
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
		}

		lineSearchResultMutex_.unlock();

	}

	lsWorkerCompleted_++; // todo hack


	if (lsWorkerCompleted_.load() >= mp_options_.nThreads_)
	{
		// only the very last thread leaving line search notifies. TODO: improve, that might be time consuming
		alphaBestFoundCondition_.notify_all();

		if(mp_options_.debugPrintMP_)
			std::cout << "NOTIFYING by LS WORKER since all workers are now done " << std::endl;
	}

	if(mp_options_.debugPrintMP_)
		std::cout<<"[Thread "<<threadId<<"]: Leaving lineSearchWorker " <<std::endl;
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
		std::vector<std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd> >>& lsLagrangeTrajectoriesStock){

	// modifying uff by local increments
	for (int i=0; i<NUM_SUBSYSTEMS; i++)
		for (int k=0; k<lsControllersStock[i].time_.size(); k++)
			lsControllersStock[i].uff_[k] += learningRate * lsControllersStock[i].deltaUff_[k];


	// modifying vff by the local increments
	for (int i=0; i<NUM_SUBSYSTEMS; i++)
		for (int k=0; k < lsLagrangeControllersStock[i].time_.size(); k++)
			lsLagrangeControllersStock[i].uff_[k] += learningRate * lsLagrangeControllersStock[i].deltaUff_[k];


	try {
		rollout(threadId, initState_, lsControllersStock, lsTimeTrajectoriesStock,
				lsStateTrajectoriesStock, lsInputTrajectoriesStock, lsOutputTrajectoriesStock,
				lsNc1TrajectoriesStock, lsEvTrajectoryStock);

		// calculate rollout cost
		calculateCostFunction(lsTimeTrajectoriesStock, lsOutputTrajectoriesStock, lsInputTrajectoriesStock, lsTotalCost, threadId);

		// calculate the merit function
		if (options_.lineSearchByMeritFuntion_==true)
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
		if (options_.dispayGSLQP_)  std::cerr << "\t learningRate " << learningRate << " \t cost: " << lsTotalCost << " \t merit: " << lsTotalMerit <<
				" \t constraint ISE: " << lsConstraint1ISE << std::endl;
	}
	catch(const std::exception& error)
	{
		std::cerr << "\t rollout with learningRate " << learningRate << " is terminated due to the slow simulation!" << std::endl;
		lsTotalMerit = std::numeric_limits<scalar_t>::max();
		lsTotalCost  = std::numeric_limits<scalar_t>::max();
	}
}

} // namespace ocs2
