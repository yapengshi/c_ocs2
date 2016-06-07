/*
 * Implementation of GSLQP.h
 *
 *  Created on: Jan 5, 2016
 *      Author: farbod
 */


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::rollout(const state_vector_t& initState,
		const std::vector<controller_t>& controllersStock,
		std::vector<scalar_array_t>& timeTrajectoriesStock,
		std::vector<state_vector_array_t>& stateTrajectoriesStock,
		std::vector<control_vector_array_t>& inputTrajectoriesStock,
		std::vector<output_vector_array_t>& outputTrajectoriesStock,
		std::vector<std::vector<size_t> >& nc1TrajectoriesStock,
		std::vector<constraint1_vector_array_t>& EvTrajectoryStock)  {

	rollout(initState, controllersStock,
			timeTrajectoriesStock, stateTrajectoriesStock, inputTrajectoriesStock, outputTrajectoriesStock);

	nc1TrajectoriesStock.resize(NUM_SUBSYSTEMS);
	EvTrajectoryStock.resize(NUM_SUBSYSTEMS);

	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		// compute constraint1 trajectory for subsystem i
		size_t N = timeTrajectoriesStock[i].size();
		nc1TrajectoriesStock[i].resize(N);
		EvTrajectoryStock[i].resize(N);

		for (int k=0; k<N; k++) {
			subsystemDynamicsPtrStock_[i]->computeConstriant1(timeTrajectoriesStock[i][k], stateTrajectoriesStock[i][k], inputTrajectoriesStock[i][k],
					nc1TrajectoriesStock[i][k], EvTrajectoryStock[i][k]);
			if (nc1TrajectoriesStock[i][k] > INPUT_DIM)
				throw std::runtime_error("Number of active type-1 constraints should be less-equal to the number of input dimension.");
		}  // end of k loop
	}  // end of i loop

}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::rollout(const state_vector_t& initState,
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

		// initialize subsystem i
		subsystemDynamicsPtrStock_[i]->initializeModel(switchingTimes_[i], x0, switchingTimes_[i+1], "GSLQP");
		// set controller for subsystem i
		subsystemDynamicsPtrStock_[i]->setController(controllersStock[i]);
		// simulate subsystem i
		subsystemSimulatorsStockPtr_[i]->integrate(x0, switchingTimes_[i], switchingTimes_[i+1],
				stateTrajectoriesStock[i], timeTrajectoriesStock[i],
				1e-3, options_.AbsTolODE_, options_.RelTolODE_);

		if (stateTrajectoriesStock[i].back() != stateTrajectoriesStock[i].back())
				throw std::runtime_error("System became unstable during the GSLQP roullouts.");

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
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::rollout(const state_vector_t& initState,
		const std::vector<controller_t>& controllersStock,
		std::vector<scalar_array_t>& timeTrajectoriesStock,
		std::vector<state_vector_array_t>& stateTrajectoriesStock,
		std::vector<control_vector_array_t>& inputTrajectoriesStock)  {

	std::vector<output_vector_array_t> outputTrajectoriesStock;
	rollout(initState, controllersStock, timeTrajectoriesStock,
			stateTrajectoriesStock, inputTrajectoriesStock, outputTrajectoriesStock);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::rolloutCost(const std::vector<scalar_array_t>& timeTrajectoriesStock,
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
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::approximateOptimalControlProblem()  {

	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

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

		// initialize subsystem i
		subsystemDerivativesPtrStock_[i]->initializeModel(nominalTimeTrajectoriesStock_[i].front(),
				nominalStateTrajectoriesStock_[i].front(), nominalTimeTrajectoriesStock_[i].back(), "GSLQP");

		for (int k=0; k<N; k++) {

			subsystemDerivativesPtrStock_[i]->setCurrentStateAndControl(nominalTimeTrajectoriesStock_[i][k],
					nominalStateTrajectoriesStock_[i][k], nominalInputTrajectoriesStock_[i][k], nominalOutputTrajectoriesStock_[i][k]);
			subsystemDerivativesPtrStock_[i]->getDerivativeState(AmTrajectoryStock_[i][k]);
			subsystemDerivativesPtrStock_[i]->getDerivativesControl(BmTrajectoryStock_[i][k]);

			size_t nc1;
			subsystemDerivativesPtrStock_[i]->getConstraint1DerivativesState(nc1, CmTrajectoryStock_[i][k]);
			if (nc1 != nc1TrajectoriesStock_[i][k])
				throw std::runtime_error("First dimension of the Cm matrix is not compatible with number of active type-1 constraints in nominal rollout.");

			subsystemDerivativesPtrStock_[i]->getConstraint1DerivativesControl(nc1, DmTrajectoryStock_[i][k]);
			if (nc1 != nc1TrajectoriesStock_[i][k])
				throw std::runtime_error("First dimension of the Dm matrix is not compatible with number of active type-1 constraints in nominal rollout.");

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
			// making sure that Qm is PSD
			makePSD(QmFinal_);
		}
	}

	// constrained coefficients
	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		int N = nominalTimeTrajectoriesStock_[i].size();

		RmConstraintProjectionTrajectoryStock_[i].resize(N);
		DmDagerTrajectoryStock_[i].resize(N);
		AmConstrainedTrajectoryStock_[i].resize(N);
		BmConstrainedTrajectoryStock_[i].resize(N);
		QmConstrainedTrajectoryStock_[i].resize(N);
		QvConstrainedTrajectoryStock_[i].resize(N);
		RvConstrainedTrajectoryStock_[i].resize(N);
		PmConstrainedTrajectoryStock_[i].resize(N);

		for (int k=0; k<N; k++) {
			size_t nc1 = nc1TrajectoriesStock_[i][k];

			if (nc1 == 0) {

				RmConstraintProjectionTrajectoryStock_[i][k].setZero();
				DmDagerTrajectoryStock_[i][k].setZero();
				AmConstrainedTrajectoryStock_[i][k] = AmTrajectoryStock_[i][k];
				BmConstrainedTrajectoryStock_[i][k] = BmTrajectoryStock_[i][k];
				QmConstrainedTrajectoryStock_[i][k] = QmTrajectoryStock_[i][k];
				QvConstrainedTrajectoryStock_[i][k] = QvTrajectoryStock_[i][k];
				RvConstrainedTrajectoryStock_[i][k] = RvTrajectoryStock_[i][k];
				PmConstrainedTrajectoryStock_[i][k] = PmTrajectoryStock_[i][k];
			} else {

				Eigen::MatrixXd Cm = CmTrajectoryStock_[i][k].topRows(nc1);
				Eigen::MatrixXd Dm = DmTrajectoryStock_[i][k].topRows(nc1);
				Eigen::MatrixXd RmConstraintProjection = ( Dm*RmInverseTrajectoryStock_[i][k]*Dm.transpose() ).inverse();
				Eigen::MatrixXd DmDager = RmInverseTrajectoryStock_[i][k]*Dm.transpose()*RmConstraintProjection;
				control_matrix_t DmNullSpaceProjection = control_matrix_t::Identity() - DmDager*Dm;
				state_matrix_t PmTransDmDagerCm = PmTrajectoryStock_[i][k].transpose()*DmDager*Cm;

				RmConstraintProjectionTrajectoryStock_[i][k].topLeftCorner(nc1,nc1) = RmConstraintProjection;
				DmDagerTrajectoryStock_[i][k].leftCols(nc1) = DmDager;
				AmConstrainedTrajectoryStock_[i][k] = AmTrajectoryStock_[i][k] - BmTrajectoryStock_[i][k]*DmDager*Cm;
				BmConstrainedTrajectoryStock_[i][k] = BmTrajectoryStock_[i][k]*DmNullSpaceProjection;
				QmConstrainedTrajectoryStock_[i][k] = QmTrajectoryStock_[i][k] + Cm.transpose()*RmConstraintProjection*Cm - PmTransDmDagerCm - PmTransDmDagerCm.transpose();
				QvConstrainedTrajectoryStock_[i][k] = QvTrajectoryStock_[i][k] - (DmDager*Cm).transpose()*RvTrajectoryStock_[i][k];
				RvConstrainedTrajectoryStock_[i][k] = DmNullSpaceProjection.transpose()*RvTrajectoryStock_[i][k];
				PmConstrainedTrajectoryStock_[i][k] = DmNullSpaceProjection.transpose()*PmTrajectoryStock_[i][k];
			}

			// making sure that Qm is PSD
			makePSD(QmConstrainedTrajectoryStock_[i][k]);
		}  // end of k loop
	}  // end of i loop
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::calculatecontroller(scalar_t& learningRateStar) {

	double maxUffe = 0.0;

	std::vector<controller_t> controllersStock(NUM_SUBSYSTEMS);

	LinearInterpolation<control_gain_matrix_t,Eigen::aligned_allocator<control_gain_matrix_t> > BmFunc;
	LinearInterpolation<constraint1_state_matrix_t,Eigen::aligned_allocator<constraint1_state_matrix_t> > CmFunc;
	LinearInterpolation<control_constraint1_matrix_t,Eigen::aligned_allocator<control_constraint1_matrix_t> > DmDagerFunc;
	LinearInterpolation<constraint1_vector_t,Eigen::aligned_allocator<constraint1_vector_t> > EvFunc;

	LinearInterpolation<control_vector_t,Eigen::aligned_allocator<control_vector_t> > RvFunc;
	LinearInterpolation<control_matrix_t,Eigen::aligned_allocator<control_matrix_t> > RmInverseFunc;
	LinearInterpolation<control_feedback_t,Eigen::aligned_allocator<control_feedback_t> > PmFunc;

	LinearInterpolation<output_vector_t,Eigen::aligned_allocator<output_vector_t> > nominalOutputFunc;
	LinearInterpolation<control_vector_t,Eigen::aligned_allocator<control_vector_t> > nominalInputFunc;

	std::vector<control_vector_array_t> deltaUffStock(NUM_SUBSYSTEMS);
	std::vector<control_vector_t> maxDeltaUffStock(NUM_SUBSYSTEMS);

	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		BmFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		BmFunc.setData( &(BmConstrainedTrajectoryStock_[i]) );

		CmFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		CmFunc.setData( &(CmTrajectoryStock_[i]) );

		DmDagerFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		DmDagerFunc.setData( &(DmDagerTrajectoryStock_[i]) );

		EvFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		EvFunc.setData( &(EvTrajectoryStock_[i]) );

		RvFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		RvFunc.setData( &(RvConstrainedTrajectoryStock_[i]) );

		RmInverseFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		RmInverseFunc.setData( &(RmInverseTrajectoryStock_[i]) );

		PmFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		PmFunc.setData( &(PmConstrainedTrajectoryStock_[i]) );

		nominalOutputFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		nominalOutputFunc.setData( &(nominalOutputTrajectoriesStock_[i]) );

		nominalInputFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		nominalInputFunc.setData( &(nominalInputTrajectoriesStock_[i]) );

		int N = SsTimeTrajectoryStock_[i].size();
		controllersStock[i].time_ = SsTimeTrajectoryStock_[i];
		controllersStock[i].k_.resize(N);
		controllersStock[i].uff_.resize(N);
		deltaUffStock[i].resize(N);
		for (int k=0; k<N; k++) {

			control_gain_matrix_t Bm;
			BmFunc.interpolate(SsTimeTrajectoryStock_[i][k], Bm);
			size_t greatestLessTimeStampIndex = BmFunc.getGreatestLessTimeStampIndex();

			control_vector_t Rv;
			RvFunc.interpolate(SsTimeTrajectoryStock_[i][k], Rv, greatestLessTimeStampIndex);
			control_matrix_t RmInverse;
			RmInverseFunc.interpolate(SsTimeTrajectoryStock_[i][k], RmInverse, greatestLessTimeStampIndex);
			control_feedback_t Pm;
			PmFunc.interpolate(SsTimeTrajectoryStock_[i][k], Pm, greatestLessTimeStampIndex);

			control_vector_t uffe = -RmInverse*Bm.transpose()*SveTrajectoryStock_[i][k];
			control_feedback_t ke = control_feedback_t::Zero();
			size_t nc1 = nc1TrajectoriesStock_[i][greatestLessTimeStampIndex];
			if (nc1 != 0) {
				constraint1_state_matrix_t Cm;
				CmFunc.interpolate(SsTimeTrajectoryStock_[i][k], Cm, greatestLessTimeStampIndex);
				control_constraint1_matrix_t DmDager;
				DmDagerFunc.interpolate(SsTimeTrajectoryStock_[i][k], DmDager, greatestLessTimeStampIndex);
				constraint1_vector_t Ev;
				EvFunc.interpolate(SsTimeTrajectoryStock_[i][k], Ev, greatestLessTimeStampIndex);

				uffe += -DmDager.leftCols(nc1)*Ev.head(nc1);
				ke = -DmDager.leftCols(nc1)*Cm.topRows(nc1);
			}

			output_vector_t nominalOutput;
			nominalOutputFunc.interpolate(SsTimeTrajectoryStock_[i][k], nominalOutput);
			control_vector_t nominalInput;
			nominalInputFunc.interpolate(SsTimeTrajectoryStock_[i][k], nominalInput);

			controllersStock[i].k_[k]   = ke - RmInverse * (Pm + Bm.transpose()*SmTrajectoryStock_[i][k]);
			controllersStock[i].uff_[k] = nominalInput + uffe - controllersStock[i].k_[k]*nominalOutput;
			deltaUffStock[i][k] = -RmInverse * (Rv + Bm.transpose()*SvTrajectoryStock_[i][k]);

			if (options_.dispayGSLQP_) maxUffe = (maxUffe<uffe.norm()) ? uffe.norm() : maxUffe;

			// testing the numerical stability of the controller gains
			try {
				if (uffe != uffe)  throw std::runtime_error("uffe is unstable.");
				if (ke != ke)      throw std::runtime_error("ke is unstable.");
				if (controllersStock[i].k_[k] != controllersStock[i].k_[k])  throw std::runtime_error("K is unstable.");
				if (deltaUffStock[i][k] != deltaUffStock[i][k])              throw std::runtime_error("deltaUff is unstable.");
			}
			catch(std::exception const& error)
			{
			    std::cerr << "what(): " << error.what() << " at time " << controllersStock[i].time_[k] << " [sec]." << std::endl;
			}
		}  // end of k loop

		// display
		if (options_.dispayGSLQP_)  maxDeltaUffStock[i] = *std::max_element(deltaUffStock[i].begin(), deltaUffStock[i].end(),
				[] (const control_vector_t& u1, const control_vector_t& u2){ return u1.norm() < u2.norm(); });

	}  // end of i loop

	// display
	if (options_.dispayGSLQP_)  {
		control_vector_t maxDeltaUff = *std::max_element(maxDeltaUffStock.begin(), maxDeltaUffStock.end(),
				[] (const control_vector_t& u1, const control_vector_t& u2){ return u1.norm() < u2.norm(); });
		std::cerr << "max delta_uff norm: " << maxDeltaUff.norm() << std::endl;
		std::cerr << "max uff_error norm: " << maxUffe << std::endl;
	}


	if (!options_.simulationIsConstrained_)  {
		nominalRolloutIsUpdated_ = true;
		nominalControllersStock_ = controllersStock;
		rollout(initState_, nominalControllersStock_, nominalTimeTrajectoriesStock_,
				nominalStateTrajectoriesStock_, nominalInputTrajectoriesStock_, nominalOutputTrajectoriesStock_,
				nc1TrajectoriesStock_, EvTrajectoryStock_);
		rolloutCost(nominalTimeTrajectoriesStock_, nominalOutputTrajectoriesStock_, nominalInputTrajectoriesStock_,
				nominalTotalCost_);
		// display
		if (options_.dispayGSLQP_)  std::cerr << "\t learningRate 0.0 cost:\t" << nominalTotalCost_ << std::endl;
	}


	// finding the optimal learningRate
	lineSearch(controllersStock, deltaUffStock, learningRateStar);

	// calculating the nominal controller
	nominalControllersStock_ = controllersStock;
	if (learningRateStar>0) {
		for (int i=0; i<NUM_SUBSYSTEMS; i++)
			for (int k=0; k<SsTimeTrajectoryStock_[i].size(); k++)
				nominalControllersStock_[i].uff_[k] += learningRateStar*deltaUffStock[i][k];
	}

}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::lineSearch(const std::vector<controller_t>& controllersStock,
		const std::vector<control_vector_array_t>& deltaUffStock,
		scalar_t& learningRateStar)  {

	scalar_t learningRate = learningRateStar;

	scalar_t lsTotalCost;
	std::vector<controller_t> lsControllersStock(NUM_SUBSYSTEMS);
	std::vector<scalar_array_t> lsTimeTrajectoriesStock(NUM_SUBSYSTEMS);
	std::vector<state_vector_array_t> lsStateTrajectoriesStock(NUM_SUBSYSTEMS);
	std::vector<control_vector_array_t> lsInputTrajectoriesStock(NUM_SUBSYSTEMS);
	std::vector<output_vector_array_t> lsOutputTrajectoriesStock(NUM_SUBSYSTEMS);
	std::vector<std::vector<size_t> > lsNc1TrajectoriesStock(NUM_SUBSYSTEMS);
	std::vector<constraint1_vector_array_t> lsEvTrajectoryStock(NUM_SUBSYSTEMS);

	while (learningRate > options_.minLearningRateGSLQP_)  {
		// modifying uff by the local increamant
		lsControllersStock = controllersStock;
		for (int i=0; i<NUM_SUBSYSTEMS; i++)
			for (int k=0; k<SsTimeTrajectoryStock_[i].size(); k++)
				lsControllersStock[i].uff_[k] += learningRate*deltaUffStock[i][k];

		// rollout
		rollout(initState_, lsControllersStock, lsTimeTrajectoriesStock, lsStateTrajectoriesStock, lsInputTrajectoriesStock, lsOutputTrajectoriesStock,
				lsNc1TrajectoriesStock, lsEvTrajectoryStock);
		// calculate rollout cost
		rolloutCost(lsTimeTrajectoriesStock, lsOutputTrajectoriesStock, lsInputTrajectoriesStock, lsTotalCost);

		// display
		if (options_.dispayGSLQP_)  std::cerr << "\t learningRate " << learningRate << " cost:\t" << lsTotalCost << std::endl;

		// break condition 1: it exits with largest learningRate that its cost is smaller than nominal cost.
		if (lsTotalCost < nominalTotalCost_*(1-1e-3*learningRate))  {
			nominalRolloutIsUpdated_ = true;
			nominalTotalCost_ = lsTotalCost;
			nominalControllersStock_ = lsControllersStock;
			nominalTimeTrajectoriesStock_  = lsTimeTrajectoriesStock;
			nominalStateTrajectoriesStock_  = lsStateTrajectoriesStock;
			nominalInputTrajectoriesStock_  = lsInputTrajectoriesStock;
			nominalOutputTrajectoriesStock_ = lsOutputTrajectoriesStock;
			nc1TrajectoriesStock_ = lsNc1TrajectoriesStock;
			EvTrajectoryStock_ = lsEvTrajectoryStock;
			learningRateStar = learningRate;
			break;
		} else {
			learningRate = 0.5*learningRate;
		}

	}  // end of while

	if (learningRate <= options_.minLearningRateGSLQP_) {
		nominalRolloutIsUpdated_ = true;  // since the open loop input will not change, the nominal trajectories will be constatnt (no disturbance effect has been assumed)
		learningRateStar = 0.0;
	}

	// display
	if (options_.dispayGSLQP_)  std::cerr << "The chosen learningRate is: " << learningRateStar << std::endl;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::transformeLocalValueFuntion2Global() {

	LinearInterpolation<output_vector_t,Eigen::aligned_allocator<output_vector_t> > nominalOutputFunc;

	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		nominalOutputFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		nominalOutputFunc.setData( &(nominalOutputTrajectoriesStock_[i]) );

		for (int k=0; k<SsTimeTrajectoryStock_[i].size(); k++) {

			output_vector_t nominalOutput;
			nominalOutputFunc.interpolate(SsTimeTrajectoryStock_[i][k], nominalOutput);

			sTrajectoryStock_[i][k]  += - nominalOutput.transpose()*SvTrajectoryStock_[i][k] + 0.5*nominalOutput.transpose()*SmTrajectoryStock_[i][k]*nominalOutput;
			SvTrajectoryStock_[i][k] += - SmTrajectoryStock_[i][k]*nominalOutput;
		}  // end of k loop
	}  // enf of i loop
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::transformeLocalValueFuntionDerivative2Global() {

	LinearInterpolation<output_vector_t,Eigen::aligned_allocator<output_vector_t> > nominalOutputFunc;

	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		nominalOutputFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		nominalOutputFunc.setData( &(nominalOutputTrajectoriesStock_[i]) );

		for (int k=0; k<SsTimeTrajectoryStock_[i].size(); k++) {

			output_vector_t nominalOutput;
			nominalOutputFunc.interpolate(SsTimeTrajectoryStock_[i][k], nominalOutput);

			for (int j=0; j<NUM_SUBSYSTEMS-1; j++)  {

				nablasTrajectoryStock_[i][k][j] += - nominalOutput.transpose()*nablaSvTrajectoryStock_[i][k][j] +
						0.5*nominalOutput.transpose()*nablaSmTrajectoryStock_[i][k][j]*nominalOutput;
				nablaSvTrajectoryStock_[i][k][j]+= - nablaSmTrajectoryStock_[i][k][j]*nominalOutput;
			}  // end of j loop
		}  // end of k loop
	}  // enf of i loop
}



/******************************************************************************************************/ //????????????????????????????????
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getRolloutSensitivity2SwitchingTime(
		std::vector<scalar_array_t>& sensitivityTimeTrajectoriesStock,
		std::vector<nabla_output_matrix_array_t>& sensitivityOutputTrajectoriesStock,
		std::vector<nabla_input_matrix_array_t>& sensitivityInputTrajectoriesStock)  {

	sensitivityTimeTrajectoriesStock   = sensitivityTimeTrajectoryStock_;
	sensitivityOutputTrajectoriesStock = nablaOutputTrajectoryStock_;
	sensitivityInputTrajectoriesStock  = nablaInputTrajectoryStock_;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getController(std::vector<controller_t>& controllersStock) {

	controllersStock = nominalControllersStock_;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getValueFuntion(const scalar_t& time, const output_vector_t& output, scalar_t& valueFuntion)  {

	int activeSubsystem = -1;
	for (int i=0; i<NUM_SUBSYSTEMS; i++)  {
		activeSubsystem = i;
		if (switchingTimes_[i]<=time && time<switchingTimes_[i+1])
			break;
	}

	state_matrix_t Sm;
	LinearInterpolation<state_matrix_t,Eigen::aligned_allocator<state_matrix_t> > SmFunc(&SsTimeTrajectoryStock_[activeSubsystem], &SmTrajectoryStock_[activeSubsystem]);
	SmFunc.interpolate(time, Sm);
	output_vector_t Sv;
	LinearInterpolation<output_vector_t,Eigen::aligned_allocator<output_vector_t> > SvFunc(&SsTimeTrajectoryStock_[activeSubsystem], &SvTrajectoryStock_[activeSubsystem]);
	SvFunc.interpolate(time, Sv);
	eigen_scalar_t s;
	LinearInterpolation<eigen_scalar_t,Eigen::aligned_allocator<eigen_scalar_t> > sFunc(&SsTimeTrajectoryStock_[activeSubsystem], &sTrajectoryStock_[activeSubsystem]);
	sFunc.interpolate(time, s);

	valueFuntion = (s + output.transpose()*Sv + 0.5*output.transpose()*Sm*output).eval()(0);
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getCostFuntionDerivative(const output_vector_t& initOutput,
		Eigen::Matrix<double,NUM_SUBSYSTEMS-1,1>& costFuntionDerivative)  {


	for (int j=0; j<NUM_SUBSYSTEMS-1; j++)  {

		state_matrix_t dSm  = nablaSmTrajectoryStock_[0][0][j];
		output_vector_t dSv = nablaSvTrajectoryStock_[0][0][j];
		eigen_scalar_t ds   = nablasTrajectoryStock_[0][0][j];

		costFuntionDerivative(j) = (ds + initOutput.transpose()*dSv + 0.5*initOutput.transpose()*dSm*initOutput).eval()(0);
	}
}



/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::getNominalTrajectories(std::vector<scalar_array_t>& nominalTimeTrajectoriesStock,
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
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::solveSequentialRiccatiEquations(const scalar_t& learningRate)  {

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
				&AmConstrainedTrajectoryStock_[i], &BmConstrainedTrajectoryStock_[i],
				&qTrajectoryStock_[i], &QvConstrainedTrajectoryStock_[i], &QmConstrainedTrajectoryStock_[i],
				&RvConstrainedTrajectoryStock_[i], &RmInverseTrajectoryStock_[i], &PmConstrainedTrajectoryStock_[i]);

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
			catch(std::exception const& error)
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
		 * Type_1 constriants error correction
		 */

		// Skip calculation of the error correction term Sve if the constrained simulation is used for forward simulation
		if (options_.simulationIsConstrained_) {
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

			GmTrajectory[k] = AmConstrainedTrajectoryStock_[i][k] - BmConstrainedTrajectoryStock_[i][k]*RmInverseTrajectoryStock_[i][k]*(
					PmConstrainedTrajectoryStock_[i][k]+BmConstrainedTrajectoryStock_[i][k].transpose()*Sm);

			size_t nc1 = nc1TrajectoriesStock_[i][k];
			if (nc1 == 0)
				GvTrajectory[k].setZero();
			else {
				Eigen::MatrixXd Cm = CmTrajectoryStock_[i][k].topRows(nc1);
				Eigen::MatrixXd Ev = EvTrajectoryStock_[i][k].head(nc1);
				Eigen::MatrixXd RmConstraintProjection = RmConstraintProjectionTrajectoryStock_[i][k].topLeftCorner(nc1,nc1);
				Eigen::MatrixXd DmDager = DmDagerTrajectoryStock_[i][k].leftCols(nc1);

				GvTrajectory[k] = Cm.transpose()*RmConstraintProjection*Ev -
						( PmTrajectoryStock_[i][k]+BmTrajectoryStock_[i][k].transpose()*Sm ).transpose() * DmDager * Ev;
			}
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
			catch(std::exception const& error) 	{
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
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::solveFullSequentialRiccatiEquations(const scalar_t& learningRate)  {

	// final value for the last Riccati equations
	typename FullRiccatiEquations_t::all_s_vector_t allSsFinal;
	nabla_Sm_t nablaQmFinal;  nablaQmFinal.fill(state_matrix_t::Zero());
	nabla_Sv_t nablaQvFinal;  nablaQvFinal.fill(output_vector_t::Zero());
	nabla_s_t  nablaqFinal;   nablaqFinal.fill(eigen_scalar_t::Zero());
	for (int i=0; i<NUM_SUBSYSTEMS-1; i++) {
		nablaQvFinal[i] = nablaQvFinal_.template block<OUTPUT_DIM,1>(0,i);
		nablaqFinal[i](0) = nablaqFinal_(i);
	}
	FullRiccatiEquations_t::convert2Vector(QmFinal_, QvFinal_, qFinal_, nablaQmFinal, nablaQvFinal, nablaqFinal, allSsFinal);

	for (int i=NUM_SUBSYSTEMS-1; i>=0; i--) {

		// set data for Riccati equations
		auto riccatiEquationsPtr = std::make_shared<FullRiccatiEquations_t>();
		riccatiEquationsPtr->setData(learningRate,
				i, switchingTimes_[i], switchingTimes_[i+1],
				&nominalTimeTrajectoriesStock_[i],
				&AmConstrainedTrajectoryStock_[i], &BmConstrainedTrajectoryStock_[i],
				&qTrajectoryStock_[i], &QvConstrainedTrajectoryStock_[i], &QmConstrainedTrajectoryStock_[i],
				&RvConstrainedTrajectoryStock_[i], &RmInverseTrajectoryStock_[i], &PmConstrainedTrajectoryStock_[i],
				&sensitivityTimeTrajectoryStock_[i], &nablaqTrajectoryStock_[i],
				&nablaQvTrajectoryStock_[i], &nablaRvTrajectoryStock_[i]);

		// integrating the Riccati equations
		ODE45<FullRiccatiEquations_t::S_DIM_*NUM_SUBSYSTEMS> ode45(riccatiEquationsPtr);
		std::vector<double> normalizedTimeTrajectory;
		std::vector<typename FullRiccatiEquations_t::all_s_vector_t, Eigen::aligned_allocator<typename FullRiccatiEquations_t::all_s_vector_t> > allSsTrajectory;
		ode45.integrate(allSsFinal, i, i+1, allSsTrajectory, normalizedTimeTrajectory,
				1e-3, options_.AbsTolODE_, options_.RelTolODE_);

		// denormalizing time and constructing 'Sm', 'Sv', and 's'
		int N = normalizedTimeTrajectory.size();
		SsTimeTrajectoryStock_[i].resize(N);
		SmTrajectoryStock_[i].resize(N);
		SvTrajectoryStock_[i].resize(N);
		sTrajectoryStock_[i].resize(N);
		nablaSmTrajectoryStock_[i].resize(N);
		nablaSvTrajectoryStock_[i].resize(N);
		nablasTrajectoryStock_[i].resize(N);
		for (int k=0; k<normalizedTimeTrajectory.size(); k++) {

			FullRiccatiEquations_t::convert2Matrix(allSsTrajectory[N-1-k],
					SmTrajectoryStock_[i][k], SvTrajectoryStock_[i][k], sTrajectoryStock_[i][k],
					nablaSmTrajectoryStock_[i][k], nablaSvTrajectoryStock_[i][k], nablasTrajectoryStock_[i][k]);
			SsTimeTrajectoryStock_[i][k] = (switchingTimes_[i]-switchingTimes_[i+1])*(normalizedTimeTrajectory[N-1-k]-i) + switchingTimes_[i+1];
		}  // end of k loop

		// reset the final value for next Riccati equation
		allSsFinal = allSsTrajectory.back();

	}  // end of i loop

}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::rolloutSensitivity2SwitchingTime()  {

	auto rolloutSensitivityEquationsPtr = std::make_shared<RolloutSensitivityEquations_t>();

	nabla_output_vector_t nabla_YmInit;
	RolloutSensitivityEquations_t::convert2Vector(nabla_output_matrix_t::Zero(), nabla_YmInit);

	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		// initialize subsystem i
		subsystemDynamicsPtrStock_[i]->initializeModel(nominalTimeTrajectoriesStock_[i].front(),
				nominalStateTrajectoriesStock_[i].front(), nominalTimeTrajectoriesStock_[i].back(), "GSLQP");

		rolloutSensitivityEquationsPtr->setData(i, switchingTimes_, subsystemDynamicsPtrStock_[i], &nominalControllersStock_[i],
				&nominalTimeTrajectoriesStock_[i], &nominalStateTrajectoriesStock_[i], &nominalInputTrajectoriesStock_[i],
				&AmConstrainedTrajectoryStock_[i], &BmConstrainedTrajectoryStock_[i]);

		// integrating
		scalar_array_t normalizedSensitivityTimeTrajectory;
		nabla_output_vector_array_t sensitivityOutputTrajectory;
		ODE45<(NUM_SUBSYSTEMS-1)*OUTPUT_DIM> ode45(rolloutSensitivityEquationsPtr);
		ode45.integrate(nabla_YmInit, i, i+1, sensitivityOutputTrajectory, normalizedSensitivityTimeTrajectory,
				1e-3, options_.AbsTolODE_, options_.RelTolODE_);

		// denormalizing time and constructing SensitivityStateTrajectory and computing control trajectory sensitivity for subsystem i
		int N = normalizedSensitivityTimeTrajectory.size();
		sensitivityTimeTrajectoryStock_[i].resize(N);
		nablaOutputTrajectoryStock_[i].resize(N);
		nablaInputTrajectoryStock_[i].resize(N);
		for (int k=0; k<N; k++) {

			sensitivityTimeTrajectoryStock_[i][k] = switchingTimes_[i] + (switchingTimes_[i+1]-switchingTimes_[i])*(normalizedSensitivityTimeTrajectory[k]-i);
			RolloutSensitivityEquations_t::convert2Matrix(sensitivityOutputTrajectory[k], nablaOutputTrajectoryStock_[i][k]);
			rolloutSensitivityEquationsPtr->computeInputSensitivity(sensitivityTimeTrajectoryStock_[i][k], nablaOutputTrajectoryStock_[i][k],
					nablaInputTrajectoryStock_[i][k]);
		}

		// reset the initial state
		nabla_YmInit = sensitivityOutputTrajectory.back();
	}

	// calculate nabla_q, nabla_Qv, nabla_Rv
	LinearInterpolation<constraint1_state_matrix_t,Eigen::aligned_allocator<constraint1_state_matrix_t> > CmFunc;
	LinearInterpolation<constraint1_control_matrix_t,Eigen::aligned_allocator<constraint1_control_matrix_t> > DmFunc;
	LinearInterpolation<output_vector_t,Eigen::aligned_allocator<output_vector_t> > QvFunc;
	LinearInterpolation<state_matrix_t,Eigen::aligned_allocator<state_matrix_t> > QmFunc;
	LinearInterpolation<control_vector_t,Eigen::aligned_allocator<control_vector_t> > RvFunc;
	LinearInterpolation<control_matrix_t,Eigen::aligned_allocator<control_matrix_t> > RmFunc;
	LinearInterpolation<control_feedback_t,Eigen::aligned_allocator<control_feedback_t> > PmFunc;

	for (int i=0; i<NUM_SUBSYSTEMS; i++) {

		CmFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		CmFunc.setData( &(CmTrajectoryStock_[i]) );
		DmFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		DmFunc.setData( &(DmTrajectoryStock_[i]) );
		QvFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		QvFunc.setData( &(QvTrajectoryStock_[i]) );
		QmFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		QmFunc.setData( &(QmTrajectoryStock_[i]) );
		RvFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		RvFunc.setData( &(RvTrajectoryStock_[i]) );
		RmFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		RmFunc.setData( &(RmTrajectoryStock_[i]) );
		PmFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		PmFunc.setData( &(PmTrajectoryStock_[i]) );

		int N = sensitivityTimeTrajectoryStock_[i].size();
		nablaqTrajectoryStock_[i].resize(N);
		nablaQvTrajectoryStock_[i].resize(N);
		nablaRvTrajectoryStock_[i].resize(N);

		for (int k=0; k<N; k++) {

			constraint1_state_matrix_t Cm;
			CmFunc.interpolate(sensitivityTimeTrajectoryStock_[i][k], Cm);
			constraint1_control_matrix_t Dm;
			DmFunc.interpolate(sensitivityTimeTrajectoryStock_[i][k], Dm);
			output_vector_t Qv;
			QvFunc.interpolate(sensitivityTimeTrajectoryStock_[i][k], Qv);
			state_matrix_t Qm;
			QmFunc.interpolate(sensitivityTimeTrajectoryStock_[i][k], Qm);
			control_vector_t Rv;
			RvFunc.interpolate(sensitivityTimeTrajectoryStock_[i][k], Rv);
			control_matrix_t Rm;
			RmFunc.interpolate(sensitivityTimeTrajectoryStock_[i][k], Rm);
			control_feedback_t Pm;
			PmFunc.interpolate(sensitivityTimeTrajectoryStock_[i][k], Pm);

			nablaqTrajectoryStock_[i][k]  = Qv.transpose()*nablaOutputTrajectoryStock_[i][k] + Rv.transpose()*nablaInputTrajectoryStock_[i][k];
			nablaQvTrajectoryStock_[i][k] = Qm*nablaOutputTrajectoryStock_[i][k] + Pm.transpose()*nablaInputTrajectoryStock_[i][k];
			nablaRvTrajectoryStock_[i][k] = Pm*nablaOutputTrajectoryStock_[i][k] + Rm*nablaInputTrajectoryStock_[i][k];
//			nablaEvTrajectoryStock_[i][k] = Cm*nablaOutputTrajectoryStock_[i][k] + Dm*nablaInputTrajectoryStock_[i][k];
		}

		if (i==NUM_SUBSYSTEMS-1)  {
			output_vector_t Qv = QvTrajectoryStock_[NUM_SUBSYSTEMS-1].back();
			state_matrix_t Qm  = QmTrajectoryStock_[NUM_SUBSYSTEMS-1].back();

			nablaqFinal_  = Qv.transpose()*nablaOutputTrajectoryStock_[NUM_SUBSYSTEMS-1].back();
			nablaQvFinal_ = Qm*nablaOutputTrajectoryStock_[NUM_SUBSYSTEMS-1].back();
		}
	}  // end of i loop

}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
template <typename Derived>
bool GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::makePSD(Eigen::MatrixBase<Derived>& squareMatrix) {

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
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_SUBSYSTEMS>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_SUBSYSTEMS>::run(const state_vector_t& initState, const std::vector<scalar_t>& switchingTimes)  {


	if (switchingTimes.size() != NUM_SUBSYSTEMS+1)
		throw std::runtime_error("Number of switching times should be one plus the number of subsystems.");

	// display
	if (options_.dispayGSLQP_)  std::cerr << "\n#### GSLQP solver starts ..." << std::endl << std::endl;

	switchingTimes_ = switchingTimes;
	initState_ = initState;

	scalar_t learningRateStar = 1.0;  // resetting learningRateStar
	size_t iteration = 0;
	double relCost = 0.0;
	double absConstraint1RMSE = 0.0;//options_.minAbsConstraint1RMSE_ + 1;
	double relConstraint1RMSE = 0.0;//options_.minRelConstraint1RMSE_ + 1;
	bool isConstraint1Satisfied = true;
	bool isCostValueConverged = false;
	bool isOptimizationConverged = false;
	while (iteration<options_.maxIterationGSLQP_ && isOptimizationConverged==false)  {

		double absConstraint1RMSECashed = absConstraint1RMSE;
		double costCashed = nominalTotalCost_;

		// do a rollout if nominalRolloutIsUpdated_ is fale.
		if (nominalRolloutIsUpdated_== false)  {
			rollout(initState_, nominalControllersStock_, nominalTimeTrajectoriesStock_,
					nominalStateTrajectoriesStock_, nominalInputTrajectoriesStock_, nominalOutputTrajectoriesStock_,
					nc1TrajectoriesStock_, EvTrajectoryStock_);
			rolloutCost(nominalTimeTrajectoriesStock_, nominalOutputTrajectoriesStock_, nominalInputTrajectoriesStock_,
					nominalTotalCost_);
			costCashed = nominalTotalCost_;
			// display
			if (options_.dispayGSLQP_ && iteration==0)  std::cerr << "\n#### Initial controller: \ncost: " << nominalTotalCost_ << std::endl;
		}

		// display
		if (options_.dispayGSLQP_)  std::cerr << "\n#### Iteration " <<  iteration << std::endl;

		// linearizing the dynamics and quadratizing the cost function along nominal trajectories
		approximateOptimalControlProblem();

		// solve Riccati equations
		solveSequentialRiccatiEquations(1.0 /*nominal learningRate*/);

		// calculate controller
		nominalRolloutIsUpdated_ = false;
		learningRateStar = 1.0;  // resetting learningRateStar
		calculatecontroller(learningRateStar);  // uses a line search scheme to determine the optimal learningRate

		// calculates type-1 constraints RMSE
		double constraint1SSE = 0.0;
		double errorSquaredNorm;
		size_t totalNumData = 0;
		double constraint1Max = 0.0;
		for (int i=0; i<NUM_SUBSYSTEMS; i++) {
			for (int k=0; k<nominalTimeTrajectoriesStock_[i].size(); k++) {
				size_t nc1 = nc1TrajectoriesStock_[i][k];
				if (nc1 != 0) {
					errorSquaredNorm = EvTrajectoryStock_[i][k].head(nc1).squaredNorm();
					constraint1SSE += errorSquaredNorm;
					totalNumData++;
					if (constraint1Max < errorSquaredNorm) constraint1Max=errorSquaredNorm;
				}
			}  // end of k loop
		}  // end of i loop
		if (totalNumData == 0)
			absConstraint1RMSE = 0.0;
		else
			absConstraint1RMSE = sqrt(constraint1SSE/(double)totalNumData);

		// display
		if (options_.dispayGSLQP_)  std::cerr << "constraint1 RMSE: " << absConstraint1RMSE << std::endl;
		if (options_.dispayGSLQP_)  std::cerr << "constraint1 Max:  " << sqrt(constraint1Max) << std::endl;
		if (options_.dispayGSLQP_)  std::cerr << "cost: " << nominalTotalCost_ << std::endl;

		// loop variables
		iteration++;
		relConstraint1RMSE = fabs(absConstraint1RMSE-absConstraint1RMSECashed);
		isConstraint1Satisfied = absConstraint1RMSE<=options_.minAbsConstraint1RMSE_ || relConstraint1RMSE<=options_.minRelConstraint1RMSE_;
		relCost = fabs(nominalTotalCost_-costCashed);
		isCostValueConverged = learningRateStar==0 || relCost<=options_.minRelCostGSLQP_;
		isOptimizationConverged = (isCostValueConverged==true) && (isConstraint1Satisfied==true);

	}  // end of while loop

	// display
	if (options_.dispayGSLQP_ && isOptimizationConverged) {
		if (learningRateStar==0)  std::cerr << "GSLQP successfully termintes as learningRate reduced to zero." << std::endl;
		else std::cerr << "GSLQP successfully termintes as cost relative change (relCost=" << relCost <<") reached to the min value." << std::endl;

		if (absConstraint1RMSE<=options_.minAbsConstraint1RMSE_)  std::cerr << "constraint1 absolute RMSE (absConstraint1RMSE=" << absConstraint1RMSE <<") reached to the min value." << std::endl;
		else std::cerr << "constraint1 relative RMSE  (relConstraint1RMSE=" << relConstraint1RMSE << ") reached to the min value." << std::endl;
	}
	if (options_.dispayGSLQP_)  std::cerr << "\n#### Final iteration" << std::endl;

	// linearizing the dynamics and quadratizing the cost function along nominal trajectories
	approximateOptimalControlProblem();
	// calculate nominal rollout sensitivity to switching times
	rolloutSensitivity2SwitchingTime();

	// solve Riccati equations
	learningRateStar = 0.0;  // prevents the changes in the nominal trajectories and just update the gains
	solveFullSequentialRiccatiEquations(learningRateStar);
	// calculate controller
	calculatecontroller(learningRateStar);

 	// transform from local value function and local derivatives to global representation
//	transformeLocalValueFuntion2Global();
	transformeLocalValueFuntionDerivative2Global();

	// display
	if (options_.dispayGSLQP_)  std::cerr << "\n#### GSLQP solver ends." << std::endl;
}

