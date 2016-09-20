/*
 * CeresTest.cpp
 *
 *  Created on: Jan 11, 2016
 *  Updated June 7, 2016
 *      Author: farbod, mgiftthaler
 */

#include <gtest/gtest.h>

#include "ceres/ceres.h"
#include "glog/logging.h"


// A CostFunction implementing analytically derivatives for the
// function f(x) = 10 - x.
class QuadraticCostFunction : public ceres::SizedCostFunction<1 /* number of residuals */, 1 /* size of first parameter */>
{
public:
	virtual ~QuadraticCostFunction() {}
	virtual bool Evaluate(double const* const* parameters,
			double* residuals,
			double** jacobians) const {
		double x = parameters[0][0];
		// f(x) = 10 - x.
		residuals[0] = 10 - x;
		// f'(x) = -1. Since there's only 1 parameter and that parameter
		// has 1 dimension, there is only 1 element to fill in the
		// jacobians.
		//
		// Since the Evaluate function can be called with the jacobians
		// pointer equal to NULL, the Evaluate function must check to see
		// if jacobians need to be computed.
		//
		// For this simple problem it is overkill to check if jacobians[0]
		// is NULL, but in general when writing more complex
		// CostFunctions, it is possible that Ceres may only demand the
		// derivatives w.r.t. a subset of the parameter blocks.
		if (jacobians != NULL && jacobians[0] != NULL) {
			jacobians[0][0] = -1;
		}
		return true;
	}
};


TEST(CeresTest, QuadraticCostFunction)
{
	std::cout << "CERES TEST"  << std::endl;
	std::cout << "========================================" << std::endl;
	std::cout << "========================================" << std::endl;


	// The variable to solve for with its initial value. It will be
	// mutated in place by the solver.
	double x = 0.5;
	const double initial_x = x;

	// Set up the problem.
	ceres::Problem problem;

	// Set up the only cost function (also known as residual).
	ceres::CostFunction* cost_function = new QuadraticCostFunction;
	problem.AddResidualBlock(cost_function, NULL, &x);

	// Run the solver!
	ceres::Solver::Options options;
	options.minimizer_progress_to_stdout = true;
	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);

	std::cout << summary.BriefReport() << "\n";
	std::cout << "x : " << initial_x
			<< " -> " << x << "\n";

	ASSERT_TRUE(summary.IsSolutionUsable());
}


int main(int argc, char** argv)
{
	google::InitGoogleLogging(argv[0]);
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}


