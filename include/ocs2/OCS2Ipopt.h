/*
 * OCS2Ipopt.h
 *
 *  Created on: Jan 12, 2016
 *      Author: farbod
 */

#ifndef OCS2IPOPT_H_
#define OCS2IPOPT_H_

#include "IpIpoptApplication.hpp"

#include "ocs2/IpoptCostFunntion.h"

using namespace Ipopt;

template <size_t STATE_DIM, size_t INPUT_DIM, size_t OUTPUT_DIM, size_t NUM_Subsystems>
class OCS2Ipopt
{
public:
	typedef GSLQP<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems> GSLQP_t;

	typedef Dimensions<STATE_DIM, INPUT_DIM, OUTPUT_DIM> DIMENSIONS;
	typedef typename DIMENSIONS::Options Options_t;
	typedef typename DIMENSIONS::controller_t controller_t;
	typedef typename DIMENSIONS::scalar_t 		scalar_t;
	typedef typename DIMENSIONS::scalar_array_t scalar_array_t;
	typedef typename DIMENSIONS::eigen_scalar_t       eigen_scalar_t;
	typedef typename DIMENSIONS::eigen_scalar_array_t eigen_scalar_array_t;
	typedef typename DIMENSIONS::state_vector_t 	  state_vector_t;
	typedef typename DIMENSIONS::state_vector_array_t state_vector_array_t;
	typedef typename DIMENSIONS::control_vector_t 		control_vector_t;
	typedef typename DIMENSIONS::control_vector_array_t control_vector_array_t;
	typedef typename DIMENSIONS::output_vector_t 	  output_vector_t;
	typedef typename DIMENSIONS::output_vector_array_t output_vector_array_t;
	typedef typename DIMENSIONS::control_feedback_t 	  control_feedback_t;
	typedef typename DIMENSIONS::control_feedback_array_t control_feedback_array_t;
	typedef typename DIMENSIONS::state_matrix_t 	  state_matrix_t;
	typedef typename DIMENSIONS::state_matrix_array_t state_matrix_array_t;
	typedef typename DIMENSIONS::control_matrix_t 		control_matrix_t;
	typedef typename DIMENSIONS::control_matrix_array_t control_matrix_array_t;
	typedef typename DIMENSIONS::control_gain_matrix_t 		 control_gain_matrix_t;
	typedef typename DIMENSIONS::control_gain_matrix_array_t control_gain_matrix_array_t;

	OCS2Ipopt (const std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM> > >& subsystemDynamicsPtr,
			const std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM> > >& subsystemDerivativesPtr,
			const std::vector<std::shared_ptr<CostFunctionBase<STATE_DIM, INPUT_DIM> > >& subsystemCostFunctionsPtr,
			const state_vector_array_t&   stateOperatingPoints,
			const control_vector_array_t& inputOperatingPoints,
			const std::vector<size_t>& systemStockIndex,
			const scalar_array_t& initSwitchingTimes,
			const state_vector_t& initState,
			const Options_t& options = Options_t::Options() )

	: subsystemDynamicsPtr_(subsystemDynamicsPtr),
	  subsystemDerivativesPtr_(subsystemDerivativesPtr),
	  subsystemCostFunctionsPtr_(subsystemCostFunctionsPtr),
	  stateOperatingPoints_(stateOperatingPoints),
	  inputOperatingPoints_(inputOperatingPoints),
	  systemStockIndex_(systemStockIndex),
	  initSwitchingTimes_(initSwitchingTimes),
	  initState_(initState),
	  options_(options)
	{
		// Create a new instance of your nlp
		//  (use a SmartPtr, not raw)
		ocs2Nlp_ = new IpoptCostFunntion<STATE_DIM, INPUT_DIM, OUTPUT_DIM, NUM_Subsystems>(subsystemDynamicsPtr_, subsystemDerivativesPtr_, subsystemCostFunctionsPtr_,
				stateOperatingPoints_, inputOperatingPoints_, systemStockIndex_, initSwitchingTimes_, initState_, options_);
	}


	virtual ~OCS2Ipopt() {}


	bool run ();


	void getCost(scalar_t& optimizedTotalCost,
			Eigen::Matrix<double,NUM_Subsystems-1,1>& optimizedTotalCostDerivative,
			scalar_t& ipoptOptimizedTotalCost=0)  const;


	void getController(scalar_array_t& optimizedSwitchingTimes,
			std::vector<controller_t>& optimizedControllersStock) const;


	void getTrajectories(std::vector<scalar_array_t>& optimizedTimeTrajectoriesStock,
			std::vector<state_vector_array_t>& optimizedStateTrajectoriesStock,
			std::vector<control_vector_array_t>& optimizedInputTrajectoriesStock,
			std::vector<output_vector_array_t>& optimizedOutputTrajectoriesStock)  const;

	void getTrajectories(std::vector<scalar_array_t>& optimizedTimeTrajectoriesStock,
			std::vector<state_vector_array_t>& optimizedStateTrajectoriesStock,
			std::vector<control_vector_array_t>& optimizedInputTrajectoriesStock)  const;


protected:
	void rolloutSolution();


private:
	std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM> > > subsystemDynamicsPtr_;
	std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM> > > subsystemDerivativesPtr_;
	std::vector<std::shared_ptr<CostFunctionBase<STATE_DIM, INPUT_DIM> > > subsystemCostFunctionsPtr_;

	state_vector_array_t   stateOperatingPoints_;
	control_vector_array_t inputOperatingPoints_;

	std::vector<size_t> systemStockIndex_;

	scalar_array_t initSwitchingTimes_;
	state_vector_t initState_;

	Options_t options_;

	// NLP
	SmartPtr<TNLP> ocs2Nlp_;

	// IPOPT optimized solution
	scalar_t ipoptOptimizedTotalCost_;
	std::vector<controller_t> ipoptOptimizedControllersStock_;
	scalar_array_t ipoptOptimizedSwitchingTimes_;
	// Simulated results
	scalar_t optimizedTotalCost_;
	//
	Eigen::Matrix<double,NUM_Subsystems-1,1> optimizedTotalCostDerivative_;
	//
	std::vector<scalar_array_t> optimizedTimeTrajectoriesStock_;
	std::vector<state_vector_array_t>    optimizedStateTrajectoriesStock_;
	std::vector<control_vector_array_t>  optimizedInputTrajectoriesStock_;
	std::vector<output_vector_array_t>   optimizedOutputTrajectoriesStock_;

};


#include "implementation/OCS2Ipopt.h"


#endif /* OCS2IPOPT_H_ */
