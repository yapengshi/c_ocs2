/*
 * A unit test
 *
 *  Created on: Sept 20, 2016
 *      Author: mgiftthaler
 */

#include <memory>
#include <fstream>

#include <cereal/archives/xml.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/Eigen.hpp>

#include "integration/Integrator.h"
#include "dynamics/ControlledSystemBase.h"

#include <gtest/gtest.h>

#include <PathTweaker.h>

using namespace ocs2;

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



TEST(IntegrationTest, SecondOrderSystem_ODE45)
{
	std::cout << "INTEGRATION TEST"  << std::endl;
	std::cout << "========================================" << std::endl;
	std::cout << "========================================" << std::endl;

	std::shared_ptr<SecondOrderSystem> sys = std::make_shared<SecondOrderSystem>();

	SecondOrderSystem::scalar_array_t cntTimeStamp {0, 10};
	SecondOrderSystem::control_vector_array_t uff(2, Eigen::Matrix<double,1,1>::Ones());
	SecondOrderSystem::control_feedback_array_t k(2, Eigen::Matrix<double,1,2>::Zero());

	SecondOrderSystem::controller_t controller;
	controller.time_ = cntTimeStamp;
	controller.uff_ = uff;
	controller.k_ = k;

	sys->setController(controller);

	std::shared_ptr<ControlledSystemBase<2, 1> > sysClone1 = sys->clone();
	std::cout << "The cloned pointer is unique: " << std::boolalpha << sysClone1.unique() << std::noboolalpha << std::endl << std::endl;

	ODE45<2> ode45(sys);
	ODE45<2> ode45_2(sys);

	std::vector<double> timeTrajectory, timeTrajectory2;
	std::vector<Eigen::Matrix<double,2,1>, Eigen::aligned_allocator<Eigen::Matrix<double,2,1>> > stateTrajectory, stateTrajectory2;
	std::vector<Eigen::Matrix<double,1,1>, Eigen::aligned_allocator<Eigen::Matrix<double,1,1>> > timeEigenTrajectory, timeEigenTrajectory2;

	Eigen::Matrix<double,2,1> x0, x02;
	x0.setZero();
	x02.setZero();

	ode45.integrate(x0, 0.0, 10.0, stateTrajectory, timeTrajectory);	// integrate adaptive
	double dt = 0.05;
	ode45_2.integrate(x02, 0.0, 10.0, dt, stateTrajectory2, timeTrajectory2);	//integrate const


	timeEigenTrajectory.resize(timeTrajectory.size());
	timeEigenTrajectory2.resize(timeTrajectory2.size());

//	std::cout << "adaptive integration" << std::endl;
//	for (size_t i=0; i<timeTrajectory.size(); i++)
//	{
//		std::cout << "At time " <<  timeTrajectory[i] << "\t state is: " << stateTrajectory[i].transpose() << std::endl;
//		timeEigenTrajectory[i](0) = timeTrajectory[i];
//	}

//	std::cout << "const integration" << std::endl;
//	for (size_t i=0; i<timeTrajectory2.size(); i++)
//	{
//		std::cout << "At time " <<  timeTrajectory2[i] << "\t state is: " << stateTrajectory2[i].transpose() << std::endl;
//		timeEigenTrajectory2[i](0) = timeTrajectory2[i];
//	}

	bool resultsGood = true;
	if(fabs(timeTrajectory.back() - 10.0) > 1e-6)
		resultsGood = false;

	if(fabs(stateTrajectory.back()[1] - 1.0) > 1e-3)
		resultsGood = false;

	if(fabs(timeTrajectory2.back() - 10.0) > 1e-6)
		resultsGood = false;

	if(fabs(stateTrajectory2.back()[1] - 1.0) > 1e-3)
		resultsGood = false;


	ASSERT_TRUE(resultsGood);
}


TEST(IntegrationTest, SecondOrderSystem_AdamsBashfort)
{
	std::cout << "INTEGRATION TEST"  << std::endl;
	std::cout << "========================================" << std::endl;
	std::cout << "========================================" << std::endl;

	std::shared_ptr<SecondOrderSystem> sys = std::make_shared<SecondOrderSystem>();

	SecondOrderSystem::scalar_array_t cntTimeStamp {0, 10};
	SecondOrderSystem::control_vector_array_t uff(2, Eigen::Matrix<double,1,1>::Ones());
	SecondOrderSystem::control_feedback_array_t k(2, Eigen::Matrix<double,1,2>::Zero());

	SecondOrderSystem::controller_t controller;
	controller.time_ = cntTimeStamp;
	controller.uff_ = uff;
	controller.k_ = k;

	sys->setController(controller);

	std::shared_ptr<ControlledSystemBase<2, 1> > sysClone1 = sys->clone();
	std::cout << "The cloned pointer is unique: " << std::boolalpha << sysClone1.unique() << std::noboolalpha << std::endl << std::endl;

	const size_t order = 5;
	IntegratorAdamsBashforth<2, order> integrator(sys);

	std::vector<double> timeTrajectory;
	std::vector<Eigen::Matrix<double,2,1>, Eigen::aligned_allocator<Eigen::Matrix<double,2,1>> > stateTrajectory;
	std::vector<Eigen::Matrix<double,1,1>, Eigen::aligned_allocator<Eigen::Matrix<double,1,1>> > timeEigenTrajectory;

	Eigen::Matrix<double,2,1> x0;
	x0.setZero();

	integrator.integrate(x0, 0.0, 10.0, 0.02, stateTrajectory, timeTrajectory);


//	timeEigenTrajectory.resize(timeTrajectory.size());
//	for (size_t i=0; i<timeTrajectory.size(); i++)
//	{
//		std::cout << "At time " <<  timeTrajectory[i] << "\t state is: " << stateTrajectory[i].transpose() << std::endl;
//		timeEigenTrajectory[i](0) = timeTrajectory[i];
//	}

	bool resultsGood = true;
	if(fabs(timeTrajectory.back() - 10.0) > 1e-6)
		resultsGood = false;

	if(fabs(stateTrajectory.back()[1] - 1.0) > 1e-3)
		resultsGood = false;


	ASSERT_TRUE(resultsGood);
}



//TEST(IntegrationTest, SecondOrderSystem_AdamsBashfortMoulton)
//{
//	std::cout << "INTEGRATION TEST"  << std::endl;
//	std::cout << "========================================" << std::endl;
//	std::cout << "========================================" << std::endl;
//
//	std::shared_ptr<SecondOrderSystem> sys = std::make_shared<SecondOrderSystem>();
//
//	SecondOrderSystem::scalar_array_t cntTimeStamp {0, 10};
//	SecondOrderSystem::control_vector_array_t uff(2, Eigen::Matrix<double,1,1>::Ones());
//	SecondOrderSystem::control_feedback_array_t k(2, Eigen::Matrix<double,1,2>::Zero());
//
//	SecondOrderSystem::controller_t controller;
//	controller.time_ = cntTimeStamp;
//	controller.uff_ = uff;
//	controller.k_ = k;
//
//	sys->setController(controller);
//
//	std::shared_ptr<ControlledSystemBase<2, 1> > sysClone1 = sys->clone();
//	std::cout << "The cloned pointer is unique: " << std::boolalpha << sysClone1.unique() << std::noboolalpha << std::endl << std::endl;
//
//	const size_t order = 5;
//	IntegratorAdamsBashforthMoulton<2, order> integrator(sys);
//
//	std::vector<double> timeTrajectory;
//	std::vector<Eigen::Matrix<double,2,1>, Eigen::aligned_allocator<Eigen::Matrix<double,2,1>> > stateTrajectory;
//	std::vector<Eigen::Matrix<double,1,1>, Eigen::aligned_allocator<Eigen::Matrix<double,1,1>> > timeEigenTrajectory;
//
//	Eigen::Matrix<double,2,1> x0;
//	x0.setZero();
//
//	integrator.integrate(x0, 0.0, 10.0, 0.02, stateTrajectory, timeTrajectory);
//
//
//	timeEigenTrajectory.resize(timeTrajectory.size());
//	for (size_t i=0; i<timeTrajectory.size(); i++)
//	{
//		std::cout << "At time " <<  timeTrajectory[i] << "\t state is: " << stateTrajectory[i].transpose() << std::endl;
//		timeEigenTrajectory[i](0) = timeTrajectory[i];
//	}
//
//	bool resultsGood = true;
//	if(fabs(timeTrajectory.back() - 10.0) > 1e-6)
//		resultsGood = false;
//
//	if(fabs(stateTrajectory.back()[1] - 1.0) > 1e-3)
//		resultsGood = false;
//
//
//	ASSERT_TRUE(resultsGood);
//}

int main(int argc, char** argv)
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

