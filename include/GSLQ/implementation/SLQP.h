/*
 * Implementation of SLQP.h
 *
 *  Created on: Jan 5, 2016
 *      Author: farbod
 */


namespace ocs2{

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/*
 * Forward integrate the system dynamics with given controller:
 * 		inputs:
 * 			+ initState: initial state at time switchingTimes_[0]
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

		size_t maxNumSteps = options_.maxNumStepsPerSecond_*(switchingTimes_[i+1]-switchingTimes_[i]);
		maxNumSteps = ((1000>maxNumSteps) ? 1000 : maxNumSteps);

		// initialize subsystem i
		subsystemDynamicsPtrStock_[i]->initializeModel(switchingTimes_, x0, i, "GSLPQ");
		// set controller for subsystem i
		subsystemDynamicsPtrStock_[i]->setController(controllersStock[i]);
		// simulate subsystem i
		subsystemSimulatorsStockPtr_[i]->integrate(x0, switchingTimes_[i], switchingTimes_[i+1],
				stateTrajectoriesStock[i], timeTrajectoriesStock[i],
				1e-3, options_.AbsTolODE_, options_.RelTolODE_, maxNumSteps);

		if (stateTrajectoriesStock[i].back() != stateTrajectoriesStock[i].back())
				throw std::runtime_error("System became unstable during the SLQP rollout.");

		// compute control trajectory for subsystem i
		inputTrajectoriesStock[i].resize(timeTrajectoriesStock[i].size());
		outputTrajectoriesStock[i].resize(timeTrajectoriesStock[i].size());
		for (int k=0; k<timeTrajectoriesStock[i].size(); k++) {
			subsystemDynamicsPtrStock_[i]->computeOutput(timeTrajectoriesStock[i][k], stateTrajectoriesStock[i][k], outputTrajectoriesStock[i][k]);
			subsystemDynamicsPtrStock_[i]->computeInput(timeTrajectoriesStock[i][k], outputTrajectoriesStock[i][k], inputTrajectoriesStock[i][k]);
		}

		// reset the initial state
		x0 = stateTrajectoriesStock[i].back();
	}
}


/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::rollout(const state_vector_t& initState,
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
void SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::rollout(const state_vector_t& initState,
		const std::vector<controller_t>& controllersStock,
		std::vector<scalar_array_t>& timeTrajectoriesStock,
		std::vector<state_vector_array_t>& stateTrajectoriesStock,
		std::vector<control_vector_array_t>& inputTrajectoriesStock,
		std::vector<output_vector_array_t>& outputTrajectoriesStock,
		std::vector<std::vector<size_t> >& nc1TrajectoriesStock,
		std::vector<constraint1_vector_array_t>& EvTrajectoryStock)  {

	rollout(initState, controllersStock,
			timeTrajectoriesStock, stateTrajectoriesStock, inputTrajectoriesStock, outputTrajectoriesStock);

	// constraint type 1 computations which consists of number of active constraints at each time point
	// and the value of the constraint (if the rollout is constrained the value is always zero otherwise
	// it is nonzero)
	nc1TrajectoriesStock.resize(NUM_SUBSYSTEMS);
	EvTrajectoryStock.resize(NUM_SUBSYSTEMS);

	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		size_t N = timeTrajectoriesStock[i].size();
		nc1TrajectoriesStock[i].resize(N);
		EvTrajectoryStock[i].resize(N);

		// compute constraint1 trajectory for subsystem i
		for (int k=0; k<N; k++) {
			subsystemDynamicsPtrStock_[i]->computeConstriant1(timeTrajectoriesStock[i][k],
					stateTrajectoriesStock[i][k], inputTrajectoriesStock[i][k],
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
		scalar_t& meritFuntionValue,
		scalar_t& constraintISE)  {

	// add cost function
	meritFuntionValue = totalCost;

	// add the L2 penalty for constraint violation
	calculateConstraintISE(timeTrajectoriesStock, nc1TrajectoriesStock, EvTrajectoryStock, constraintISE);
	double pho = iteration_/(options_.maxIterationGSLQP_-1) * options_.meritFunctionRho_;
	meritFuntionValue += 0.5*pho*constraintISE;

	// add the the lagrangian term for the constraint
	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		// integrates the intermediate merit using the trapezoidal approximation method
		scalar_t currentIntermediateMerit;
		scalar_t nextIntermediateMerit;
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
double SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculateConstraintISE(
		const std::vector<scalar_array_t>& timeTrajectoriesStock,
		const std::vector<std::vector<size_t>>& nc1TrajectoriesStock,
		const std::vector<constraint1_vector_array_t>& EvTrajectoriesStock,
		scalar_t& constraintISE)  {

	constraintISE = 0.0;
	double maxConstraintNorm = 0.0;

	for (size_t i=0; i<NUM_SUBSYSTEMS; i++)  {

		scalar_t currentSquaredNormError;
		scalar_t nextSquaredNormError;

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
void SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::approximateOptimalControlProblem()  {

	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		// initialize subsystem i dynamics derivatives
		subsystemDerivativesPtrStock_[i]->initializeModel(switchingTimes_, nominalStateTrajectoriesStock_[i].front(), i, "GSLPQ");

		int N = nominalTimeTrajectoriesStock_[i].size();

		AmTrajectoryStock_[i].resize(N);
		BmTrajectoryStock_[i].resize(N);
		CmTrajectoryStock_[i].resize(N);
		DmTrajectoryStock_[i].resize(N);

		qTrajectoryStock_[i].resize(N);
		QvTrajectoryStock_[i].resize(N);
		QmTrajectoryStock_[i].resize(N);
		RvTrajectoryStock_[i].resize(N);
		RmTrajectoryStock_[i].resize(N);
		RmInverseTrajectoryStock_[i].resize(N);
		PmTrajectoryStock_[i].resize(N);

		for (int k=0; k<N; k++) {

			subsystemDerivativesPtrStock_[i]->setCurrentStateAndControl(nominalTimeTrajectoriesStock_[i][k],
					nominalStateTrajectoriesStock_[i][k], nominalInputTrajectoriesStock_[i][k], nominalOutputTrajectoriesStock_[i][k]);
			subsystemDerivativesPtrStock_[i]->getDerivativeState(AmTrajectoryStock_[i][k]);
			subsystemDerivativesPtrStock_[i]->getDerivativesControl(BmTrajectoryStock_[i][k]);
			// if constraint type 1 is active
			if (nc1TrajectoriesStock_[i][k] > 0) {
				subsystemDerivativesPtrStock_[i]->getConstraint1DerivativesState(CmTrajectoryStock_[i][k]);
				subsystemDerivativesPtrStock_[i]->getConstraint1DerivativesControl(DmTrajectoryStock_[i][k]);
			}

			subsystemCostFunctionsPtrStock_[i]->setCurrentStateAndControl(nominalTimeTrajectoriesStock_[i][k],
					nominalOutputTrajectoriesStock_[i][k], nominalInputTrajectoriesStock_[i][k]);
			subsystemCostFunctionsPtrStock_[i]->evaluate(qTrajectoryStock_[i][k](0));
			subsystemCostFunctionsPtrStock_[i]->stateDerivative(QvTrajectoryStock_[i][k]);
			subsystemCostFunctionsPtrStock_[i]->stateSecondDerivative(QmTrajectoryStock_[i][k]);
			subsystemCostFunctionsPtrStock_[i]->controlDerivative(RvTrajectoryStock_[i][k]);
			subsystemCostFunctionsPtrStock_[i]->controlSecondDerivative(RmTrajectoryStock_[i][k]);
			RmInverseTrajectoryStock_[i][k] = RmTrajectoryStock_[i][k].inverse();
			subsystemCostFunctionsPtrStock_[i]->stateControlDerivative(PmTrajectoryStock_[i][k]);

		}

		if (i==NUM_SUBSYSTEMS-1)  {
			subsystemCostFunctionsPtrStock_[i]->terminalCost(qFinal_(0));
			subsystemCostFunctionsPtrStock_[i]->terminalCostStateDerivative(QvFinal_);
			subsystemCostFunctionsPtrStock_[i]->terminalCostStateSecondDerivative(QmFinal_);
			// making sure that Qm remains PSD
			makePSD(QmFinal_);
		}
	}

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

			if (nc1 == 0) {
				DmDagerTrajectoryStock_[i][k].setZero();
				EvProjectedTrajectoryStock_[i][k].setZero();
				CmProjectedTrajectoryStock_[i][k].setZero();
				DmProjectedTrajectoryStock_[i][k].setZero();

				AmConstrainedTrajectoryStock_[i][k] = AmTrajectoryStock_[i][k];
				QmConstrainedTrajectoryStock_[i][k] = QmTrajectoryStock_[i][k];
				QvConstrainedTrajectoryStock_[i][k] = QvTrajectoryStock_[i][k];
				RmConstrainedTrajectoryStock_[i][k] = RmTrajectoryStock_[i][k];

			} else {
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

		nominalOutputFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		nominalOutputFunc.setData( &(nominalOutputTrajectoriesStock_[i]) );

		nominalInputFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		nominalInputFunc.setData( &(nominalInputTrajectoriesStock_[i]) );

		BmFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		BmFunc.setData( &(BmTrajectoryStock_[i]) );

		PmFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		PmFunc.setData( &(PmTrajectoryStock_[i]) );

		RmInverseFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		RmInverseFunc.setData( &(RmInverseTrajectoryStock_[i]) );

		RvFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		RvFunc.setData( &(RvTrajectoryStock_[i]) );

		EvProjectedFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		EvProjectedFunc.setData( &(EvProjectedTrajectoryStock_[i]) );

		CmProjectedFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		CmProjectedFunc.setData( &(CmProjectedTrajectoryStock_[i]) );

		DmProjectedFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		DmProjectedFunc.setData( &(DmProjectedTrajectoryStock_[i]) );

		// functions for lagrane multiplier only

		RmFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		RmFunc.setData( &(RmTrajectoryStock_[i]) );

		DmDagerFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		DmDagerFunc.setData( &(DmDagerTrajectoryStock_[i]) );

		if (firstCall==false) {
			nominalLagrangeMultiplierFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
			nominalLagrangeMultiplierFunc.setData( &(nominalLagrangeTrajectoriesStock_[i]) );
		}

		int N = SsTimeTrajectoryStock_[i].size();

		controllersStock[i].time_ = SsTimeTrajectoryStock_[i];
		controllersStock[i].k_.resize(N);
		controllersStock[i].uff_.resize(N);
		controllersStock[i].deltaUff_.resize(N);

		feedForwardConstraintInputStock[i].resize(N);

		lagrangeMultiplierFunctionsStock[i].time_ = SsTimeTrajectoryStock_[i];
		lagrangeMultiplierFunctionsStock[i].k_.resize(N);
		lagrangeMultiplierFunctionsStock[i].uff_.resize(N);
		lagrangeMultiplierFunctionsStock[i].deltaUff_.resize(N);

		for (int k=0; k<N; k++) {

			const double& time = SsTimeTrajectoryStock_[i][k];
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

			control_feedback_t Lm  = RmInverse * (Pm + Bm.transpose()*SmTrajectoryStock_[i][k]);
			control_vector_t   Lv  = RmInverse * (Rv + Bm.transpose()*SvTrajectoryStock_[i][k]);
			control_vector_t   Lve = RmInverse * (Bm.transpose()*SveTrajectoryStock_[i][k]);

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

			const size_t& nc1 = nc1TrajectoriesStock_[i][greatestLessTimeStampIndex];

			control_constraint1_matrix_t DmDager;
			DmDagerFunc.interpolate(time, DmDager, greatestLessTimeStampIndex);
			control_matrix_t Rm;
			RmFunc.interpolate(time, Rm, greatestLessTimeStampIndex);

			Eigen::MatrixXd DmDagerTransRm = DmDager.leftCols(nc1).transpose() * Rm;

			lagrangeMultiplierFunctionsStock[i].k_[k]   = DmDagerTransRm * (CmProjected - Lm);
			lagrangeMultiplierFunctionsStock[i].uff_[k] = -lagrangeMultiplierFunctionsStock[i].k_[k]*nominalOutput;
			Eigen::VectorXd localVff = DmDagerTransRm * (EvProjected-Lv-Lve);

			if (firstCall==true || options_.lineSearchByMeritFuntion_==false) {
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
	if (options_.dispayGSLQP_)  {
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
			if (options_.lineSearchByMeritFuntion_==true)
				nominalControllersStock_[i].deltaUff_[k] += feedForwardConstraintInputStock[i][k];
			else
				nominalControllersStock_[i].uff_[k] += options_.constraintStepSize_*feedForwardConstraintInputStock[i][k];

	// perform one rollout while the input correction for the type-1 constraint is considered.
	rollout(initState_, nominalControllersStock_, nominalTimeTrajectoriesStock_,
			nominalStateTrajectoriesStock_, nominalInputTrajectoriesStock_, nominalOutputTrajectoriesStock_,
			nc1TrajectoriesStock_, EvTrajectoryStock_);
	calculateCostFunction(nominalTimeTrajectoriesStock_, nominalOutputTrajectoriesStock_, nominalInputTrajectoriesStock_,
			nominalTotalCost_);

	// calculate the merit function
	if (options_.lineSearchByMeritFuntion_==true) {
		// calculate the lagrange multiplier with learningRate zero
		calculateRolloutLagrangeMultiplier(nominalTimeTrajectoriesStock_, nominalOutputTrajectoriesStock_, lagrangeControllerStock_,
				nominalLagrangeTrajectoriesStock_);
		// calculate the merit function
		calculateMeritFunction(nominalTimeTrajectoriesStock_, nc1TrajectoriesStock_, EvTrajectoryStock_, nominalLagrangeTrajectoriesStock_, nominalTotalCost_,
				nominalTotalMerit_, nominalConstraint1ISE_);
	} else {
		nominalTotalMerit_ = nominalTotalCost_;
		calculateConstraintISE(nominalTimeTrajectoriesStock_, nc1TrajectoriesStock_, EvTrajectoryStock_, nominalConstraint1ISE_);
	}

	// display
	if (options_.dispayGSLQP_)  std::cerr << "\t learningRate 0.0 \t cost: " << nominalTotalCost_ << " \t merit: " << nominalTotalMerit_ <<
			" \t constraint ISE: " << nominalConstraint1ISE_ << std::endl;

	scalar_t learningRate = maxLearningRateStar;
	const std::vector<controller_t> controllersStock = nominalControllersStock_;
	const std::vector<lagrange_t> lagrangeMultiplierFunctionsStock = lagrangeControllerStock_;

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
	std::vector<std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd> > >  lsLagrangeTrajectoriesStock(NUM_SUBSYSTEMS);

	while (learningRate >= options_.minLearningRateGSLQP_)  {
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
			rollout(initState_, lsControllersStock, lsTimeTrajectoriesStock, lsStateTrajectoriesStock, lsInputTrajectoriesStock, lsOutputTrajectoriesStock,
					lsNc1TrajectoriesStock, lsEvTrajectoryStock);
			// calculate rollout cost
			calculateCostFunction(lsTimeTrajectoriesStock, lsOutputTrajectoriesStock, lsInputTrajectoriesStock, lsTotalCost);

			// calculate the merit function
			if (options_.lineSearchByMeritFuntion_==true) {
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
			if (options_.dispayGSLQP_)  std::cerr << "\t learningRate " << learningRate << " \t cost: " << lsTotalCost << " \t merit: " << lsTotalMerit <<
					" \t constraint ISE: " << lsConstraint1ISE << std::endl;
		}
		catch(const std::exception& error)
		{
			std::cerr << "\t rollout with learningRate " << learningRate << " is terminated due to the slow simulation!" << std::endl;
			lsTotalCost  = std::numeric_limits<scalar_t>::max();
			lsTotalMerit = std::numeric_limits<scalar_t>::max();
		}

		// break condition 1: it exits with largest learningRate that its cost is smaller than nominal cost.
		if (lsTotalMerit < nominalTotalMerit_*(1-1e-3*learningRate))
			break;  // exit while loop
		else
			learningRate = 0.5*learningRate;

	}  // end of while


	if (learningRate >= options_.minLearningRateGSLQP_)  {
		nominalTotalCost_      = lsTotalCost;
		nominalTotalMerit_     = lsTotalMerit;
		nominalConstraint1ISE_ = lsConstraint1ISE;
		nominalControllersStock_ = lsControllersStock;
		nominalTimeTrajectoriesStock_  = lsTimeTrajectoriesStock;
		nominalStateTrajectoriesStock_  = lsStateTrajectoriesStock;
		nominalInputTrajectoriesStock_  = lsInputTrajectoriesStock;
		nominalOutputTrajectoriesStock_ = lsOutputTrajectoriesStock;
		nc1TrajectoriesStock_ = lsNc1TrajectoriesStock;
		EvTrajectoryStock_ = lsEvTrajectoryStock;
		learningRateStar = learningRate;
		lagrangeControllerStock_ = lsLagrangeControllersStock;
		nominalLagrangeTrajectoriesStock_ = lsLagrangeTrajectoriesStock;

	} else // since the open loop input is not change, the nominal trajectories will be unchanged
		learningRateStar = 0.0;

	// clear the feedforward increments
	for (size_t i=0; i<NUM_SUBSYSTEMS; i++) {
		nominalControllersStock_[i].deltaUff_.clear();
		lagrangeControllerStock_[i].deltaUff_.clear();
	}

	// display
	if (options_.dispayGSLQP_)  std::cerr << "The chosen learningRate is: " << learningRateStar << std::endl;
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
		if (switchingTimes_[i]<=time && time<switchingTimes_[i+1])
			break;
	}

	state_matrix_t Sm;
	LinearInterpolation<state_matrix_t,Eigen::aligned_allocator<state_matrix_t> > SmFunc(&SsTimeTrajectoryStock_[activeSubsystem], &SmTrajectoryStock_[activeSubsystem]);
	SmFunc.interpolate(time, Sm);
	size_t greatestLessTimeStampIndex = SmFunc.getGreatestLessTimeStampIndex();
	output_vector_t Sv;
	LinearInterpolation<output_vector_t,Eigen::aligned_allocator<output_vector_t> > SvFunc(&SsTimeTrajectoryStock_[activeSubsystem], &SvTrajectoryStock_[activeSubsystem]);
	SvFunc.interpolate(time, Sv, greatestLessTimeStampIndex);
	eigen_scalar_t s;
	LinearInterpolation<eigen_scalar_t,Eigen::aligned_allocator<eigen_scalar_t> > sFunc(&SsTimeTrajectoryStock_[activeSubsystem], &sTrajectoryStock_[activeSubsystem]);
	sFunc.interpolate(time, s, greatestLessTimeStampIndex);

	output_vector_t xNomilnal;
	LinearInterpolation<output_vector_t,Eigen::aligned_allocator<output_vector_t> > xNominalFunc(&nominalTimeTrajectoriesStock_[activeSubsystem], &nominalOutputTrajectoriesStock_[activeSubsystem]);
	xNominalFunc.interpolate(time, xNomilnal);

	output_vector_t deltaX = output-xNomilnal;
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
void SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getCostFuntion(const output_vector_t& initOutput, scalar_t& costFunction, scalar_t& constriantCostFunction)  {

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
void SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getNominalTrajectories(std::vector<scalar_array_t>& nominalTimeTrajectoriesStock,
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
void SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::solveSequentialRiccatiEquations(const scalar_t& learningRate)  {

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
bool SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::makePSD(Eigen::MatrixBase<Derived>& squareMatrix) {

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
void SLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::run(const state_vector_t& initState, const std::vector<scalar_t>& switchingTimes)  {

	if (switchingTimes.size() != NUM_SUBSYSTEMS+1)
		throw std::runtime_error("Number of switching times should be one plus the number of subsystems.");

	switchingTimes_ = switchingTimes;
	initState_ = initState;

	// display
	if (options_.dispayGSLQP_) {
		std::cerr << "\n#### SLQP solver starts with switching times [" << switchingTimes[0];
		for (size_t i=1; i<=NUM_SUBSYSTEMS; i++)   std::cerr << ", " << switchingTimes[i];
		std::cerr << "] ..." << std::endl << std::endl;
	}

	iteration_ = 0;
	double relCost;
	double learningRateStar;
	double relConstraint1ISE;
	bool isConstraint1Satisfied  = false;
	bool isCostFunctionConverged = false;
	bool isOptimizationConverged = false;
	bool nominalLagrangeMultiplierUpdated = false;

	// initial controller rollout
	rollout(initState_, nominalControllersStock_,
			nominalTimeTrajectoriesStock_, nominalStateTrajectoriesStock_, nominalInputTrajectoriesStock_, nominalOutputTrajectoriesStock_,
			nc1TrajectoriesStock_, EvTrajectoryStock_);
	// initial controller cost
	calculateCostFunction(nominalTimeTrajectoriesStock_, nominalOutputTrajectoriesStock_, nominalInputTrajectoriesStock_,
			nominalTotalCost_);
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
		approximateOptimalControlProblem(); // todo parallelize

		// solve Riccati equations
		solveSequentialRiccatiEquations(1.0 /*nominal learningRate*/); //todo do not touch

		// calculate controller and lagrange multiplier
		std::vector<control_vector_array_t> feedForwardConstraintInputStock(NUM_SUBSYSTEMS);
		calculateControllerAndLagrangian(nominalControllersStock_, lagrangeControllerStock_, feedForwardConstraintInputStock, ~nominalLagrangeMultiplierUpdated); // todo parallelize
		nominalLagrangeMultiplierUpdated = true;

		// finding the optimal learningRate
		lineSearch(feedForwardConstraintInputStock, learningRateStar, options_.maxLearningRateGSLQP_); // todo parallelize

		// calculates type-1 constraint ISE and maximum norm
		double constraint1MaxNorm = calculateConstraintISE(nominalTimeTrajectoriesStock_, nc1TrajectoriesStock_, EvTrajectoryStock_, nominalConstraint1ISE_);

		// loop variables
		iteration_++;
		relCost = fabs(nominalTotalCost_-costCashed);
		relConstraint1ISE = fabs(nominalConstraint1ISE_-constraint1ISECashed);
		isConstraint1Satisfied  = nominalConstraint1ISE_<=options_.minAbsConstraint1ISE_ || relConstraint1ISE<=options_.minRelConstraint1ISE_;
		isCostFunctionConverged = learningRateStar==0 || relCost<=options_.minRelCostGSLQP_;
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
	calculateRolloutCostate(nominalTimeTrajectoriesStock_, nominalOutputTrajectoriesStock_,
			nominalcostateTrajectoriesStock_);

//	// transform from local value function to global representation (avoid this function)
//	transformLocalValueFuntion2Global();

	// display
	if (options_.dispayGSLQP_ )  {
		std::cerr << "\n+++++++++++++++++++++++++++++++++++" << std::endl;
		std::cerr <<   "+++++++ SLQP solver is ended ++++++" << std::endl;
		std::cerr <<   "+++++++++++++++++++++++++++++++++++" << std::endl;
		if (isOptimizationConverged) {
			if (learningRateStar==0)
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

} // namespace ocs2

