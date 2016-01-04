/*
 * IntegrationTest.cpp
 *
 *  Created on: Dec 30, 2015
 *      Author: farbod
 */

#include <memory>
#include <fstream>

#include <cereal/archives/xml.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/Eigen.hpp>

#include "integration/Integrator.h"
#include "dynamics/ControlledSystemBase.h"


class SecondOrderSystem : public ControlledSystemBase<2,1>
{
public:

	SecondOrderSystem() {}
	~SecondOrderSystem() {}

	void computeDerivative(
			const double& t,
			const Eigen::Matrix<double,2,1>& x,
			const Eigen::Matrix<double,1,1>& u,
			Eigen::Matrix<double,2,1>& dxdt) {

		Eigen::Matrix2d A;
		A << -2, -1, 1, 0;

		Eigen::Vector2d B;
		B << 1, 0;

		dxdt = A*x + B*u;
	}

	std::shared_ptr<ControlledSystemBase<2, 1> > clone() const {
		return std::make_shared<SecondOrderSystem>(*this);
	}

private:

};



int main (int argc, char* argv[])
{
	std::shared_ptr<SecondOrderSystem> sys = std::make_shared<SecondOrderSystem>();

	SecondOrderSystem::scalar_array_t cntTimeStamp {0, 10};
	SecondOrderSystem::control_vector_array_t uff(2, Eigen::Matrix<double,1,1>::Ones());
	SecondOrderSystem::control_feedback_array_t k(2, Eigen::Matrix<double,1,2>::Zero());

	SecondOrderSystem::controller_t controller;
	controller.time_ = cntTimeStamp;
	controller.uff_ = uff;
	controller.k_ = k;

	sys->setController(controller);

	//
	std::shared_ptr<ControlledSystemBase<2, 1> > sysClone1 = sys->clone();
	std::cout << "The cloned pointer is unique: " << std::boolalpha << sysClone1.unique() << std::noboolalpha << std::endl << std::endl;

	ODE45<2> ode45(sys);

	std::vector<double> timeTrajectory;
	std::vector<Eigen::Matrix<double,2,1>, Eigen::aligned_allocator<Eigen::Matrix<double,2,1>> > stateTrajectory;
	std::vector<Eigen::Matrix<double,1,1>, Eigen::aligned_allocator<Eigen::Matrix<double,1,1>> > timeEigenTrajectory;

	Eigen::Matrix<double,2,1> x0;
	x0.setZero();

	bool flag = ode45.integrate(x0, 0.0, 10.0, stateTrajectory, timeTrajectory);

	timeEigenTrajectory.resize(timeTrajectory.size());
	for (size_t i=0; i<timeTrajectory.size(); i++) {
		std::cout << "At time " <<  timeTrajectory[i] << "\t state is: " << stateTrajectory[i].transpose() << std::endl;
		timeEigenTrajectory[i](0) = timeTrajectory[i];
	}

	std::string resultDir = "/home/farbod/Programs/ct_ws/src/c_ocs2/cereal/test/integration_test";
	std::string secondOrderStateFile = resultDir + "/secondOrderState.xml";
	std::string secondOrderTimeFile = resultDir + "/secondOrderTime.xml";


	{ // we need these brackets to make sure the archive goes out of scope and flushes
		std::ofstream xmlSecondOrderState(secondOrderStateFile);
		cereal::XMLOutputArchive archive_secondOrderState_xml(xmlSecondOrderState);
		archive_secondOrderState_xml(CEREAL_NVP(stateTrajectory));

		std::ofstream xmlSecondOrderTime(secondOrderTimeFile);
		cereal::XMLOutputArchive archive_secondOrderTime_xml(xmlSecondOrderTime);
		archive_secondOrderTime_xml(CEREAL_NVP(timeEigenTrajectory));
	}

}

