/*
 * IpopotCostFunntion.h
 *
 *  Created on: Jan 12, 2016
 *      Author: farbod
 */

#ifndef IPOPOTCOSTFUNNTION_H_
#define IPOPOTCOSTFUNNTION_H_

#include "IpTNLP.hpp"

#include "GSLQ/GLQP.h"
#include "GSLQ/GSLQP.h"

using namespace Ipopt;

template <size_t STATE_DIM, size_t INPUT_DIM, size_t NUM_Subsystems>
class IpopotCostFunntion : public TNLP
{
public:
	enum {NumParameters_=NUM_Subsystems-1, NumConstraints_=NUM_Subsystems-2};
	typedef GLQP<STATE_DIM, INPUT_DIM, NUM_Subsystems> GLQP_t;
	typedef GSLQP<STATE_DIM, INPUT_DIM, NUM_Subsystems> GSLQP_t;

	typedef Dimensions<STATE_DIM, INPUT_DIM> DIMENSIONS;
	typedef typename DIMENSIONS::controller_t controller_t;
	typedef typename DIMENSIONS::scalar_t 		scalar_t;
	typedef typename DIMENSIONS::scalar_array_t scalar_array_t;
	typedef typename DIMENSIONS::eigen_scalar_t       eigen_scalar_t;
	typedef typename DIMENSIONS::eigen_scalar_array_t eigen_scalar_array_t;
	typedef typename DIMENSIONS::state_vector_t 	  state_vector_t;
	typedef typename DIMENSIONS::state_vector_array_t state_vector_array_t;
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

	IpopotCostFunntion(const std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM> > >& subsystemDynamicsPtr,
			const std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM> > >& subsystemDerivativesPtr,
			const std::vector<std::shared_ptr<CostFunctionBase<STATE_DIM, INPUT_DIM> > >& subsystemCostFunctionsPtr,
			const state_vector_array_t&   stateOperatingPoints,
			const control_vector_array_t& inputOperatingPoints,
			const std::vector<size_t>& systemStockIndex,
			const scalar_array_t& initSwitchingTimes,
			const state_vector_t& initState,
			const typename GSLQP_t::Options& options = GSLQP_t::Options() )
		: subsystemDynamicsPtr_(subsystemDynamicsPtr),
		  subsystemDerivativesPtr_(subsystemDerivativesPtr),
		  subsystemCostFunctionsPtr_(subsystemCostFunctionsPtr),
		  stateOperatingPoints_(stateOperatingPoints),
		  inputOperatingPoints_(inputOperatingPoints),
		  systemStockIndex_(systemStockIndex),
		  initSwitchingTimes_(initSwitchingTimes),
		  initState_(initState),
		  options_(options),
		  optimizedSwitchingTimes_(NUM_Subsystems+1),
		  numFuntionCall_(0)
	{}

	virtual ~IpopotCostFunntion() {}

	/**@name Overloaded from TNLP */
	//@{
	/** Method to return some info about the nlp */
	virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
			Index& nnz_h_lag, IndexStyleEnum& index_style);

	/** Method to return the bounds for my problem */
	virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
			Index m, Number* g_l, Number* g_u);

	/** Method to return the starting point for the algorithm */
	virtual bool get_starting_point(Index n, bool init_x, Number* x,
			bool init_z, Number* z_L, Number* z_U,
			Index m, bool init_lambda,
			Number* lambda);

	/** Method to return the objective value */
	virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

	/** Method to return the gradient of the objective */
	virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

	/** Method to return the constraint residuals */
	virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

	/** Method to return:
	 *   1) The structure of the jacobian (if "values" is NULL)
	 *   2) The values of the jacobian (if "values" is not NULL)
	 */
	virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
			Index m, Index nele_jac, Index* iRow, Index *jCol,
			Number* values);

	/** Method to return:
	 *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
	 *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
	 */
	virtual bool eval_h(Index n, const Number* x, bool new_x,
			Number obj_factor, Index m, const Number* lambda,
			bool new_lambda, Index nele_hess, Index* iRow,
			Index* jCol, Number* values);

	//@}

	/** @name Solution Methods */
	//@{
	/** This method is called when the algorithm is complete so the TNLP can store/write the solution */
	virtual void finalize_solution(SolverReturn status,
			Index n, const Number* x, const Number* z_L, const Number* z_U,
			Index m, const Number* g, const Number* lambda,
			Number obj_value,
			const IpoptData* ip_data,
			IpoptCalculatedQuantities* ip_cq);
	//@}

protected:
	IpopotCostFunntion();
	IpopotCostFunntion(const IpopotCostFunntion&);
	IpopotCostFunntion& operator=(const IpopotCostFunntion&);

	void solveGSLQP(const Number* x);

private:
	std::vector<std::shared_ptr<ControlledSystemBase<STATE_DIM, INPUT_DIM> > > subsystemDynamicsPtr_;
	std::vector<std::shared_ptr<DerivativesBase<STATE_DIM, INPUT_DIM> > > subsystemDerivativesPtr_;
	std::vector<std::shared_ptr<CostFunctionBase<STATE_DIM, INPUT_DIM> > > subsystemCostFunctionsPtr_;

	state_vector_array_t   stateOperatingPoints_;
	control_vector_array_t inputOperatingPoints_;

	std::vector<size_t> systemStockIndex_;

	scalar_array_t initSwitchingTimes_;
	state_vector_t initState_;

	typename GSLQP_t::Options options_;

	scalar_t optimizedTotalCost_;
	scalar_array_t optimizedSwitchingTimes_;

	scalar_t currentTotalCost_;
	Eigen::Matrix<double,NumParameters_,1> currentCostFuntionDerivative_;

	size_t numFuntionCall_;

};

#include "implementation/IpopotCostFunntion.h"

#endif /* IPOPOTCOSTFUNNTION_H_ */
