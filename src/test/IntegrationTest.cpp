/*
 * IntegrationTest.cpp
 *
 *  Created on: Dec 30, 2015
 *      Author: farbod
 */

#include <memory>

#include "integration/Integrator.h"
#include "test/SecondOrderSystem.h"


int main (int argc, char* argv[])
{
	std::shared_ptr<SecondOrderSystem> sys = std::make_shared<SecondOrderSystem>();

	std::vector<double> cntTimeStamp {0, 10};
	std::vector<Eigen::Matrix<double,1,1>, Eigen::aligned_allocator<Eigen::Matrix<double,1,1> > > uff(2, Eigen::Matrix<double,1,1>::Ones());
	std::vector<Eigen::Matrix<double,1,2>, Eigen::aligned_allocator<Eigen::Matrix<double,1,2> > > k(2, Eigen::Matrix<double,1,2>::Zero());

	sys->setController(cntTimeStamp, uff, k);

	ODE45<2> ode45(sys);

	std::vector<double> timeTrajectory;
	std::vector<Eigen::Matrix<double,2,1>, Eigen::aligned_allocator<Eigen::Matrix<double,2,1>> > stateTrajectory;

	Eigen::Matrix<double,2,1> x0;
	x0.setZero();

	bool flag = ode45.integrate(x0, 0.0, 10.0, stateTrajectory, timeTrajectory, 0.001);

	for (size_t i=0; i<timeTrajectory.size(); i++) {
		std::cout << "At time " <<  timeTrajectory[i] << "\t state is: " << stateTrajectory[i].transpose() << std::endl;
	}

}

