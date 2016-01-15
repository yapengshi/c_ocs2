/*
 * SystemDynamicsDerivativeTest.cpp
 *
 *  Created on: Jan 13, 2016
 *      Author: farbod
 */

#include <iostream>

#include "dynamics/SystemDynamicsLinearizer.h"

#include "test/EXP1.h"


int main(int argc, char* argv[])
{

	typedef Dimensions<2, 1> DIMENSIONS;

	// subsystem dynamics
	std::vector<std::shared_ptr<ControlledSystemBase<2,1> > > subsystemDynamicsPtr {std::make_shared<EXP1_Sys1>(), std::make_shared<EXP1_Sys2>(), std::make_shared<EXP1_Sys3>()};

	// subsystem derivatives
	std::vector<std::shared_ptr<DerivativesBase<2,1> > > subsystemDerivativesPtr {std::make_shared<EXP1_SysDerivative1>(), std::make_shared<EXP1_SysDerivative2>(), std::make_shared<EXP1_SysDerivative3>()};

	// subsystem numerical derivatives
	std::vector<std::shared_ptr<DerivativesBase<2,1> > > subsystemNumericalDerivativesPtr(3);
	for (size_t i=0; i<3; i++)
		subsystemNumericalDerivativesPtr[i] = std::make_shared<SystemDynamicsLinearizer<2,1> >(subsystemDynamicsPtr[i]);


	DIMENSIONS::state_vector_array_t statePoints;
	DIMENSIONS::control_vector_array_t inputPoints;


	for (size_t k=0; k<100; k++) {
		statePoints.push_back(DIMENSIONS::state_vector_t::Random()*10);
		inputPoints.push_back(DIMENSIONS::control_vector_t::Random()*10);

//		std::cout << "x: " << statePoints.back().transpose() << std::endl;
//		std::cout << "u: " << inputPoints.back().transpose() << std::endl;
	}

	std::array<DIMENSIONS::state_matrix_array_t, 3> APoints1;
	std::array<DIMENSIONS::control_gain_matrix_array_t, 3> BPoints1;

	std::array<DIMENSIONS::state_matrix_array_t, 3> APoints2;
	std::array<DIMENSIONS::control_gain_matrix_array_t, 3> BPoints2;

	size_t numErrorsA = 0;
	size_t numErrorsB = 0;

	const double tol = 1e-7;
	for (size_t i=0; i<3; i++)  {

		size_t N = statePoints.size();
		APoints1[i].resize(N);
		BPoints1[i].resize(N);
		APoints2[i].resize(N);
		BPoints2[i].resize(N);

		for (size_t k=0; k<N; k++) {

			subsystemDerivativesPtr[i]->setCurrentStateAndControl(0.0, statePoints[k], inputPoints[k]);
			subsystemDerivativesPtr[i]->getDerivativeState(APoints1[i][k]);
			subsystemDerivativesPtr[i]->getDerivativesControl(BPoints1[i][k]);

			subsystemNumericalDerivativesPtr[i]->setCurrentStateAndControl(0.0, statePoints[k], inputPoints[k]);
			subsystemNumericalDerivativesPtr[i]->getDerivativeState(APoints2[i][k]);
			subsystemNumericalDerivativesPtr[i]->getDerivativesControl(BPoints2[i][k]);

			if ( !(APoints1[i][k]-APoints2[i][k]).isZero(tol) ) {
				numErrorsA++;
				std::cout << "delta_A[" << i << "] " << APoints1[i][k]-APoints2[i][k] << std::endl;
			}

			if ( !(BPoints1[i][k]-BPoints2[i][k]).isZero(tol) ) {
				numErrorsB++;
				std::cout << "delta_B[" << i << "] " << BPoints1[i][k]-BPoints2[i][k] << std::endl;
			}

		}
	}


	if (numErrorsA == 0)  std::cout << "#### State derivative is correct." << std::endl;
	if (numErrorsB == 0)  std::cout << "#### Input derivative is correct." << std::endl;

}
