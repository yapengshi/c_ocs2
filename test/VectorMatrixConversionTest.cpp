/*
 * MatrixVectorConversion.cpp
 *
 *  Created on: Oct 19, 2016
 *      Author: mgiftthaler@ethz.ch
 */


#include <iostream>
#include <cstdlib>

#include <gtest/gtest.h>

#include "GSLQ/SequentialRiccatiEquations.h"

using namespace ocs2;

const size_t state_dim = 4;
const size_t input_dim =  3;
const size_t num_subsystem = 1;

TEST(MatrixToVectorConversionTest, MatrixToVectorConversionTest)
{
	typedef SequentialRiccatiEquations<state_dim, input_dim, state_dim, num_subsystem> SRE;
	typedef typename SRE::state_matrix_t state_matrix_t;
	typedef typename SRE::state_vector_t state_vector_t;
	typedef typename SRE::eigen_scalar_t eigen_scalar_t;
	typedef typename SRE::s_vector_t s_vector_t;

	const size_t nTests = 10;

	for(size_t i = 0; i<nTests; i++){

		// convert a random state matrix and vector and scalar to an s_vector_t
		state_matrix_t testMatrix, backtransformedTestMatrix;
		testMatrix.setRandom();

		// make the testMatrix symmetric
		testMatrix = testMatrix.selfadjointView<Eigen::Upper>();

		state_vector_t testVector, backtransformedTestVector; testVector.setRandom();
		eigen_scalar_t testScalar, backtransformedTestScalar; testScalar.setRandom();

		//	std::cout << "starting data: matrix " << std::endl;
		//	std::cout << testMatrix << std::endl;
		//	std::cout << "starting data: vector " << std::endl;
		//	std::cout << testVector.transpose() << std::endl;
		//	std::cout << "starting data: scalar " << std::endl;
		//	std::cout << testScalar << std::endl;

		// do conversion
		s_vector_t stackedVector;
		SRE::convert2Vector(testMatrix, testVector, testScalar, stackedVector);

		//	std::cout << "stackedVector" << std::endl << stackedVector.transpose() << std::endl;

		SRE::convert2Matrix(stackedVector, backtransformedTestMatrix, backtransformedTestVector, backtransformedTestScalar);

		//	std::cout << "backtransformed data: matrix " << std::endl;
		//	std::cout << backtransformedTestMatrix << std::endl;
		//	std::cout << "backtransformed data: vector " << std::endl;
		//	std::cout << backtransformedTestVector.transpose() << std::endl;
		//	std::cout << "backtransformed data: scalar " << std::endl;
		//	std::cout << backtransformedTestScalar << std::endl;

		ASSERT_TRUE(backtransformedTestMatrix.isApprox(testMatrix));
		ASSERT_TRUE(backtransformedTestVector.isApprox(backtransformedTestVector));
		ASSERT_TRUE(backtransformedTestScalar.isApprox(testScalar));
	}
}


int main(int argc, char** argv)
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
