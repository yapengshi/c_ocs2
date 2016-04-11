/*
 * Implementation of GSLQP.h
 *
 *  Created on: Jan 5, 2016
 *      Author: farbod
 */

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::rollout(const state_vector_t& initState,
		const std::vector<controller_t>& controllersStock,
		std::vector<scalar_array_t>& timeTrajectoriesStock,
		std::vector<state_vector_array_t>& stateTrajectoriesStock,
		std::vector<control_vector_array_t>& inputTrajectoriesStock,
		std::vector<output_vector_array_t>& outputTrajectoriesStock)  {

	if (controllersStock.size() != NUM_Subsystems)
		throw std::runtime_error("controllersStock has less controllers then the number of subsystems");

	timeTrajectoriesStock.resize(NUM_Subsystems);
	stateTrajectoriesStock.resize(NUM_Subsystems);
	inputTrajectoriesStock.resize(NUM_Subsystems);
	outputTrajectoriesStock.resize(NUM_Subsystems);

	state_vector_t x0 = initState;
	for (int i=0; i<NUM_Subsystems; i++) {

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
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::rollout(const state_vector_t& initState,
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
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::rolloutCost(const std::vector<scalar_array_t>& timeTrajectoriesStock,
		const std::vector<output_vector_array_t>& outputTrajectoriesStock,
		const std::vector<control_vector_array_t>& inputTrajectoriesStock,
		scalar_t& totalCost)  {

	totalCost = 0.0;
	for (int i=0; i<NUM_Subsystems; i++) {

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
		}

		// terminal cost
		if (i==NUM_Subsystems-1)  {
			scalar_t finalCost;
			subsystemCostFunctionsPtrStock_[i]->setCurrentStateAndControl(timeTrajectoriesStock[i].back(), outputTrajectoriesStock[i].back(), inputTrajectoriesStock[i].back());
			subsystemCostFunctionsPtrStock_[i]->terminalCost(finalCost);
			totalCost += finalCost;
		}
	}

}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::approximateOptimalControlProblem()  {

	for (int i=0; i<NUM_Subsystems; i++) {

		int N = nominalTimeTrajectoriesStock_[i].size();

		AmTrajectoryStock_[i].resize(N);
		BmTrajectoryStock_[i].resize(N);

		qTrajectoryStock_[i].resize(N);
		QvTrajectoryStock_[i].resize(N);
		QmTrajectoryStock_[i].resize(N);
		RvTrajectoryStock_[i].resize(N);
		RmTrajectoryStock_[i].resize(N);
		PmTrajectoryStock_[i].resize(N);

		// initialize subsystem i
		subsystemDerivativesPtrStock_[i]->initializeModel(nominalTimeTrajectoriesStock_[i].front(),
				nominalStateTrajectoriesStock_[i].front(), nominalTimeTrajectoriesStock_[i].back(), "GSLQP");

		for (int k=0; k<N; k++) {

			subsystemDerivativesPtrStock_[i]->setCurrentStateAndControl(nominalTimeTrajectoriesStock_[i][k],
					nominalStateTrajectoriesStock_[i][k], nominalInputTrajectoriesStock_[i][k], nominalOutputTrajectoriesStock_[i][k]);
			subsystemDerivativesPtrStock_[i]->getDerivativeState(AmTrajectoryStock_[i][k]);
			subsystemDerivativesPtrStock_[i]->getDerivativesControl(BmTrajectoryStock_[i][k]);

			subsystemCostFunctionsPtrStock_[i]->setCurrentStateAndControl(nominalTimeTrajectoriesStock_[i][k],
					nominalOutputTrajectoriesStock_[i][k], nominalInputTrajectoriesStock_[i][k]);
			subsystemCostFunctionsPtrStock_[i]->evaluate(qTrajectoryStock_[i][k](0));
			subsystemCostFunctionsPtrStock_[i]->stateDerivative(QvTrajectoryStock_[i][k]);
			subsystemCostFunctionsPtrStock_[i]->stateSecondDerivative(QmTrajectoryStock_[i][k]);
			subsystemCostFunctionsPtrStock_[i]->controlDerivative(RvTrajectoryStock_[i][k]);
			subsystemCostFunctionsPtrStock_[i]->controlSecondDerivative(RmTrajectoryStock_[i][k]);
			subsystemCostFunctionsPtrStock_[i]->stateControlDerivative(PmTrajectoryStock_[i][k]);

			// making sure that Qm is PSD
			makePSD(QmTrajectoryStock_[i][k]);
		}

		if (i==NUM_Subsystems-1)  {
			subsystemCostFunctionsPtrStock_[i]->terminalCost(qFinal_(0));
			subsystemCostFunctionsPtrStock_[i]->terminalCostStateDerivative(QvFinal_);
			subsystemCostFunctionsPtrStock_[i]->terminalCostStateSecondDerivative(QmFinal_);
			// making sure that Qm is PSD
			makePSD(QmFinal_);
		}
	}
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::calculatecontroller(scalar_t& learningRateStar) {

	std::vector<controller_t> controllersStock(NUM_Subsystems);

	LinearInterpolation<state_matrix_t,Eigen::aligned_allocator<state_matrix_t> > AmFunc;
	LinearInterpolation<control_gain_matrix_t,Eigen::aligned_allocator<control_gain_matrix_t> > BmFunc;

	LinearInterpolation<eigen_scalar_t,Eigen::aligned_allocator<eigen_scalar_t> > qFunc;
	LinearInterpolation<output_vector_t,Eigen::aligned_allocator<output_vector_t> > QvFunc;
	LinearInterpolation<state_matrix_t,Eigen::aligned_allocator<state_matrix_t> > QmFunc;
	LinearInterpolation<control_vector_t,Eigen::aligned_allocator<control_vector_t> > RvFunc;
	LinearInterpolation<control_matrix_t,Eigen::aligned_allocator<control_matrix_t> > RmFunc;
	LinearInterpolation<control_feedback_t,Eigen::aligned_allocator<control_feedback_t> > PmFunc;

	LinearInterpolation<output_vector_t,Eigen::aligned_allocator<output_vector_t> > nominalOutputFunc;
	LinearInterpolation<control_vector_t,Eigen::aligned_allocator<control_vector_t> > nominalInputFunc;

	std::vector<control_vector_array_t> deltaUffStock(NUM_Subsystems);
	std::vector<control_vector_t> maxDeltaUffStock(NUM_Subsystems);

	for (int i=0; i<NUM_Subsystems; i++) {

		AmFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		AmFunc.setData( &(AmTrajectoryStock_[i]) );

		BmFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		BmFunc.setData( &(BmTrajectoryStock_[i]) );

		qFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		qFunc.setData( &(qTrajectoryStock_[i]) );

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

			state_matrix_t Am;
			AmFunc.interpolate(SsTimeTrajectoryStock_[i][k], Am);
			control_gain_matrix_t Bm;
			BmFunc.interpolate(SsTimeTrajectoryStock_[i][k], Bm);

			eigen_scalar_t q;
			qFunc.interpolate(SsTimeTrajectoryStock_[i][k], q);
			output_vector_t Qv;
			QvFunc.interpolate(SsTimeTrajectoryStock_[i][k], Qv);
			state_matrix_t Qm;
			QmFunc.interpolate(SsTimeTrajectoryStock_[i][k], Qm);
			control_vector_t Rv;
			RvFunc.interpolate(SsTimeTrajectoryStock_[i][k], Rv);
			control_matrix_t Rm;
			RmFunc.interpolate(SsTimeTrajectoryStock_[i][k], Rm);
			control_matrix_t RmInverse = Rm.inverse();
			control_feedback_t Pm;
			PmFunc.interpolate(SsTimeTrajectoryStock_[i][k], Pm);

			output_vector_t nominalOutput;
			nominalOutputFunc.interpolate(SsTimeTrajectoryStock_[i][k], nominalOutput);
			control_vector_t nominalInput;
			nominalInputFunc.interpolate(SsTimeTrajectoryStock_[i][k], nominalInput);

			controllersStock[i].k_[k]   = -RmInverse * (Pm + Bm.transpose()*SmTrajectoryStock_[i][k]);
			controllersStock[i].uff_[k] = nominalInput - controllersStock[i].k_[k]*nominalOutput;
			deltaUffStock[i][k] = -RmInverse * (Rv + Bm.transpose()*SvTrajectoryStock_[i][k]);
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
	}

	// finding the optimal learningRate
	lineSearch(controllersStock, deltaUffStock, learningRateStar);

	// calculating the nominal controller
	nominalControllersStock_ = controllersStock;
	if (learningRateStar>0) {
		for (int i=0; i<NUM_Subsystems; i++)
			for (int k=0; k<SsTimeTrajectoryStock_[i].size(); k++)
				nominalControllersStock_[i].uff_[k] += learningRateStar*deltaUffStock[i][k];
	}

}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::lineSearch(const std::vector<controller_t>& controllersStock,
		const std::vector<control_vector_array_t>& deltaUffStock,
		scalar_t& learningRateStar)  {

	scalar_t learningRate = learningRateStar;

	scalar_t lsTotalCost;
	std::vector<controller_t> lsControllersStock(NUM_Subsystems);
	std::vector<scalar_array_t> lsTimeTrajectoriesStock(NUM_Subsystems);
	std::vector<state_vector_array_t> lsStateTrajectoriesStock(NUM_Subsystems);
	std::vector<control_vector_array_t> lsInputTrajectoriesStock(NUM_Subsystems);
	std::vector<output_vector_array_t> lsOutputTrajectoriesStock(NUM_Subsystems);

	while (learningRate > options_.minLearningRateGSLQP_)  {
		// modifying uff by the local increamant
		lsControllersStock = controllersStock;
		for (int i=0; i<NUM_Subsystems; i++)
			for (int k=0; k<SsTimeTrajectoryStock_[i].size(); k++)
				lsControllersStock[i].uff_[k] += learningRate*deltaUffStock[i][k];

		// rollout
		rollout(initState_, lsControllersStock, lsTimeTrajectoriesStock, lsStateTrajectoriesStock, lsInputTrajectoriesStock, lsOutputTrajectoriesStock);
		// calculate rollout cost
		rolloutCost(lsTimeTrajectoriesStock, lsOutputTrajectoriesStock, lsInputTrajectoriesStock, lsTotalCost);

		// display
		if (options_.dispayGSLQP_)  std::cerr << "\t learningRate " << learningRate << " cost: " << nominalTotalCost_ << std::endl;

		// break condition 1: it exits with largest learningRate that its cost is smaller than nominal cost.
		if (lsTotalCost < nominalTotalCost_*(1-1e-3*learningRate))  {
			nominalRolloutIsUpdated_ = true;
			nominalTotalCost_ = lsTotalCost;
			nominalControllersStock_ = lsControllersStock;
			nominalTimeTrajectoriesStock_  = lsTimeTrajectoriesStock;
			nominalStateTrajectoriesStock_  = lsStateTrajectoriesStock;
			nominalInputTrajectoriesStock_  = lsInputTrajectoriesStock;
			nominalOutputTrajectoriesStock_ = lsOutputTrajectoriesStock;
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
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::transformeLocalValueFuntion2Global() {

	LinearInterpolation<output_vector_t,Eigen::aligned_allocator<output_vector_t> > nominalOutputFunc;

	for (int i=0; i<NUM_Subsystems; i++) {

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
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::transformeLocalValueFuntionDerivative2Global() {

	LinearInterpolation<output_vector_t,Eigen::aligned_allocator<output_vector_t> > nominalOutputFunc;

	for (int i=0; i<NUM_Subsystems; i++) {

		nominalOutputFunc.setTimeStamp( &(nominalTimeTrajectoriesStock_[i]) );
		nominalOutputFunc.setData( &(nominalOutputTrajectoriesStock_[i]) );

		for (int k=0; k<SsTimeTrajectoryStock_[i].size(); k++) {

			output_vector_t nominalOutput;
			nominalOutputFunc.interpolate(SsTimeTrajectoryStock_[i][k], nominalOutput);

			for (int j=0; j<NUM_Subsystems-1; j++)  {

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
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::getRolloutSensitivity2SwitchingTime(
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
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::getController(std::vector<controller_t>& controllersStock) {

	controllersStock = nominalControllersStock_;
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::getValueFuntion(const scalar_t& time, const output_vector_t& output, scalar_t& valueFuntion)  {

	int activeSubsystem = -1;
	for (int i=0; i<NUM_Subsystems; i++)  {
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
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::getCostFuntionDerivative(const output_vector_t& initOutput,
		Eigen::Matrix<double,NUM_Subsystems-1,1>& costFuntionDerivative)  {


	for (int j=0; j<NUM_Subsystems-1; j++)  {

		state_matrix_t dSm  = nablaSmTrajectoryStock_[0][0][j];
		output_vector_t dSv = nablaSvTrajectoryStock_[0][0][j];
		eigen_scalar_t ds   = nablasTrajectoryStock_[0][0][j];

		costFuntionDerivative(j) = (ds + initOutput.transpose()*dSv + 0.5*initOutput.transpose()*dSm*initOutput).eval()(0);
	}
}



/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::getNominalTrajectories(std::vector<scalar_array_t>& nominalTimeTrajectoriesStock,
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
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::SolveSequentialRiccatiEquations(const scalar_t& learningRate)  {

	// final value for the last Riccati equations
	typename RiccatiEquations_t::s_vector_t allSsFinal;
	RiccatiEquations_t::convert2Vector(QmFinal_, QvFinal_, qFinal_, allSsFinal);

	for (int i=NUM_Subsystems-1; i>=0; i--) {

		// set data for Riccati equations
		auto riccatiEquationsPtr = std::make_shared<RiccatiEquations_t>();
		riccatiEquationsPtr->setData(learningRate,
				i, switchingTimes_[i], switchingTimes_[i+1],
				&nominalTimeTrajectoriesStock_[i],
				&AmTrajectoryStock_[i], &BmTrajectoryStock_[i],
				&qTrajectoryStock_[i], &QvTrajectoryStock_[i], &QmTrajectoryStock_[i],
				&RvTrajectoryStock_[i], &RmTrajectoryStock_[i], &PmTrajectoryStock_[i]);

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
		for (int k=0; k<normalizedTimeTrajectory.size(); k++) {

			RiccatiEquations_t::convert2Matrix(allSsTrajectory[N-1-k], SmTrajectoryStock_[i][k], SvTrajectoryStock_[i][k], sTrajectoryStock_[i][k]);
			SsTimeTrajectoryStock_[i][k] = (switchingTimes_[i]-switchingTimes_[i+1])*(normalizedTimeTrajectory[N-1-k]-i) + switchingTimes_[i+1];
		}

		// reset the final value for next Riccati equation
		allSsFinal = allSsTrajectory.back();
	}

}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::SolveFullSequentialRiccatiEquations(const scalar_t& learningRate)  {

	// final value for the last Riccati equations
	typename FullRiccatiEquations_t::all_s_vector_t allSsFinal;
	nabla_Sm_t nablaQmFinal;  nablaQmFinal.fill(state_matrix_t::Zero());
	nabla_Sv_t nablaQvFinal;  nablaQvFinal.fill(output_vector_t::Zero());
	nabla_s_t  nablaqFinal;   nablaqFinal.fill(eigen_scalar_t::Zero());
	for (int i=0; i<NUM_Subsystems-1; i++) {
		nablaQvFinal[i] = nablaQvFinal_.template block<OUTPUT_DIM,1>(0,i);
		nablaqFinal[i](0) = nablaqFinal_(i);
	}
	FullRiccatiEquations_t::convert2Vector(QmFinal_, QvFinal_, qFinal_, nablaQmFinal, nablaQvFinal, nablaqFinal, allSsFinal);

	for (int i=NUM_Subsystems-1; i>=0; i--) {

		// set data for Riccati equations
		auto riccatiEquationsPtr = std::make_shared<FullRiccatiEquations_t>();
		riccatiEquationsPtr->setData(learningRate,
				i, switchingTimes_[i], switchingTimes_[i+1],
				&nominalTimeTrajectoriesStock_[i],
				&AmTrajectoryStock_[i], &BmTrajectoryStock_[i],
				&qTrajectoryStock_[i], &QvTrajectoryStock_[i], &QmTrajectoryStock_[i],
				&RvTrajectoryStock_[i], &RmTrajectoryStock_[i], &PmTrajectoryStock_[i],
				&sensitivityTimeTrajectoryStock_[i], &nablaqTrajectoryStock_[i],
				&nablaQvTrajectoryStock_[i], &nablaRvTrajectoryStock_[i]);

		// integrating the Riccati equations
		ODE45<FullRiccatiEquations_t::S_DIM_*NUM_Subsystems> ode45(riccatiEquationsPtr);
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
		}

		// reset the final value for next Riccati equation
		allSsFinal = allSsTrajectory.back();
	}

}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::rolloutSensitivity2SwitchingTime()  {

	auto rolloutSensitivityEquationsPtr = std::make_shared<RolloutSensitivityEquations_t>();

	typename RolloutSensitivityEquations_t::nabla_output_vector_t nabla_YmInit;
	RolloutSensitivityEquations_t::convert2Vector(nabla_output_matrix_t::Zero(), nabla_YmInit);

	for (int i=0; i<NUM_Subsystems; i++) {

		// initialize subsystem i
		subsystemDynamicsPtrStock_[i]->initializeModel(nominalTimeTrajectoriesStock_[i].front(),
				nominalStateTrajectoriesStock_[i].front(), nominalTimeTrajectoriesStock_[i].back(), "GSLQP");

		rolloutSensitivityEquationsPtr->setData(i, switchingTimes_, subsystemDynamicsPtrStock_[i], &nominalControllersStock_[i],
				&nominalTimeTrajectoriesStock_[i], &nominalStateTrajectoriesStock_[i], &nominalInputTrajectoriesStock_[i],
				&AmTrajectoryStock_[i], &BmTrajectoryStock_[i]);

		scalar_array_t normalizedSensitivityTimeTrajectory;
		std::vector<typename RolloutSensitivityEquations_t::nabla_output_vector_t,
			Eigen::aligned_allocator<typename RolloutSensitivityEquations_t::nabla_output_vector_t> > sensitivityOutputTrajectory;

		// integrating
		ODE45<(NUM_Subsystems-1)*OUTPUT_DIM> ode45(rolloutSensitivityEquationsPtr);
		ode45.integrate(nabla_YmInit, i, i+1, sensitivityOutputTrajectory, normalizedSensitivityTimeTrajectory,
				1e-3, options_.AbsTolODE_, options_.RelTolODE_);

		// denormalizing time and constructing SensitivityStateTrajectory and computing control trajectory sensitivity for subsystem i
		int N = sensitivityOutputTrajectory.size();
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
	LinearInterpolation<output_vector_t,Eigen::aligned_allocator<output_vector_t> > QvFunc;
	LinearInterpolation<state_matrix_t,Eigen::aligned_allocator<state_matrix_t> > QmFunc;
	LinearInterpolation<control_vector_t,Eigen::aligned_allocator<control_vector_t> > RvFunc;
	LinearInterpolation<control_matrix_t,Eigen::aligned_allocator<control_matrix_t> > RmFunc;
	LinearInterpolation<control_feedback_t,Eigen::aligned_allocator<control_feedback_t> > PmFunc;

	for (int i=0; i<NUM_Subsystems; i++) {

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
		}

		if (i==NUM_Subsystems-1)  {
			output_vector_t Qv = QvTrajectoryStock_[NUM_Subsystems-1].back();
			state_matrix_t Qm  = QmTrajectoryStock_[NUM_Subsystems-1].back();

			nablaqFinal_  = Qv.transpose()*nablaOutputTrajectoryStock_[NUM_Subsystems-1].back();
			nablaQvFinal_ = Qm*nablaOutputTrajectoryStock_[NUM_Subsystems-1].back();
		}
	}

}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
template <typename Derived>
bool GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::makePSD(Eigen::MatrixBase<Derived>& squareMatrix) {

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
template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
void GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>::run(const state_vector_t& initState, const std::vector<scalar_t>& switchingTimes)  {

	// display
	if (options_.dispayGSLQP_)  std::cerr << "\n#### GSLQP solver starts ..." << std::endl << std::endl;

	if (switchingTimes.size() != NUM_Subsystems+1)
		throw std::runtime_error("Number of switching times should be one plus the number of subsystems.");
	switchingTimes_ = switchingTimes;
	initState_ = initState;

	scalar_t learningRateStar = 1.0;  // resetting learningRateStar
	size_t iteration = 0;
	while (iteration<options_.maxIterationGSLQP_ && learningRateStar>0)  {

		// do a rollout if nominalRolloutIsUpdated_ is fale.
		if (nominalRolloutIsUpdated_==false)  {
			rollout(initState_, nominalControllersStock_, nominalTimeTrajectoriesStock_,
					nominalStateTrajectoriesStock_, nominalInputTrajectoriesStock_, nominalOutputTrajectoriesStock_);
			rolloutCost(nominalTimeTrajectoriesStock_, nominalOutputTrajectoriesStock_, nominalInputTrajectoriesStock_,
					nominalTotalCost_);
			// display
			if (options_.dispayGSLQP_ && iteration==0)  std::cerr << "\n#### Initial controller: \ncost: " << nominalTotalCost_ << std::endl;
		}

		// display
		if (options_.dispayGSLQP_)  std::cerr << "\n#### Iteration " <<  iteration << std::endl;

		// linearizing the dynamics and quadratizing the cost funtion along nominal trajectories
		approximateOptimalControlProblem();

		// solve Riccati equations
		SolveSequentialRiccatiEquations(1.0 /*nominal learningRate*/);

		// calculate controller
		nominalRolloutIsUpdated_ = false;
		learningRateStar = 1.0;  // resetting learningRateStar
		calculatecontroller(learningRateStar);

		// display
		if (options_.dispayGSLQP_)  std::cerr << "cost: " << nominalTotalCost_ << std::endl;

		iteration++;

	}  // end of j loop

	// display
	if (options_.dispayGSLQP_)  std::cerr << "\n#### Final iteration" << std::endl;

	// linearizing the dynamics and quadratizing the cost funtion along nominal trajectories
	approximateOptimalControlProblem();
	// calculate nominal rollout sensitivity to switching times
	rolloutSensitivity2SwitchingTime();

	// solve Riccati equations
	learningRateStar = 0.0;  // prevents the changes in the nominal trajectories and just update the gains
	SolveFullSequentialRiccatiEquations(learningRateStar);
	// calculate controller
	calculatecontroller(learningRateStar);

 	// transforme from local value funtion and local derivatives to global representation
	transformeLocalValueFuntion2Global();
	transformeLocalValueFuntionDerivative2Global();

	// display
	if (options_.dispayGSLQP_)  std::cerr << "\n#### GSLQP solver ends." << std::endl;
}

