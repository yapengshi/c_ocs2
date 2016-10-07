/*
 * Integrator.h
 *
 *  Created on: 18.06.2015
 *      Author: neunertm
 */

#ifndef OCS2INTEGRATOR_H_
#define OCS2INTEGRATOR_H_

#include <type_traits>
#include <functional>
#include <cmath>

#include <boost/numeric/odeint.hpp>
#include "eigenIntegration.h"

#include "IntegratorBase.h"


namespace ocs2{

template <size_t STATE_DIM>
using dense_runge_kutta5_t = boost::numeric::odeint::dense_output_runge_kutta <
		boost::numeric::odeint::controlled_runge_kutta <
		boost::numeric::odeint::runge_kutta_dopri5 <
		Eigen::Matrix<double, STATE_DIM, 1>,
		double,
		Eigen::Matrix<double, STATE_DIM, 1>,
		double,
		boost::numeric::odeint::vector_space_algebra > > >;



template <size_t STATE_DIM, class Stepper>
class Integrator : public IntegratorBase<STATE_DIM>
{
public:
	typedef IntegratorBase<STATE_DIM> Base;

	// Constructor
	Integrator(
			const std::shared_ptr<SystemBase<STATE_DIM> >& system,
			const std::shared_ptr<EventHandler<STATE_DIM> >& eventHandler = nullptr
	)
	: IntegratorBase<STATE_DIM>(system, eventHandler)
	  {
		setupSystem();
	  }


	// Equidistant integration based on number of time steps and step length
	//	bool integrate_n_steps(
	//		const typename Base::State_T& initialState,
	//		const double& startTime,
	//		size_t numSteps,
	//		double dt,
	//		typename Base::StateTrajectory_T& stateTrajectory,
	//		typename Base::TimeTrajectory_T& timeTrajectory
	//	 	) override
	//	 {
	//		 typename Base::State_T initialStateInternal = initialState;
	//		 initialize(initialStateInternal, startTime, dt);
	//		 integrate_n_steps(stepper_, systemFunction_, initialStateInternal, startTime, dt, numSteps, Base::observer_.observeWrap);
	//		 Base::retrieveTrajectoriesFromObserver(stateTrajectory, timeTrajectory);
	//		 return true;
	//	 }

	// Equidistant integration based on initial and final time as well as step length
	bool integrate(
			const typename Base::State_T& initialState,
			const double& startTime,
			const double& finalTime,
			double dt,
			typename Base::StateTrajectory_T& stateTrajectory,
			typename Base::TimeTrajectory_T& timeTrajectory
	) override
			{
		typename Base::State_T initialStateInternal = initialState;
		initialize(initialStateInternal, startTime, dt);
		integrate_const(stepper_, systemFunction_, initialStateInternal, startTime, finalTime, dt, Base::observer_.observeWrap);
		Base::retrieveTrajectoriesFromObserver(stateTrajectory, timeTrajectory);
		return true;
			}

	// Adaptive time integration based on start time and final time
	bool integrate(
			const typename Base::State_T& initialState,
			const double& startTime,
			const double& finalTime,
			typename Base::StateTrajectory_T& stateTrajectory,
			typename Base::TimeTrajectory_T& timeTrajectory,
			double dtInitial = 0.01,
			double AbsTol = 1e-9,
			double RelTol = 1e-6,
			size_t maxNumSteps = std::numeric_limits<size_t>::max()
	) override  {


		integrate_adaptive_helper<Stepper>( initialState, startTime, finalTime, stateTrajectory, timeTrajectory, dtInitial, AbsTol, RelTol, maxNumSteps );


		Base::retrieveTrajectoriesFromObserver(stateTrajectory, timeTrajectory);

		return true;
	}

	// Output integration based on a given time trajectory
	bool integrate(
			const typename Base::State_T& initialState,
			const typename Base::TimeTrajectory_T& timeTrajectory,
			typename Base::StateTrajectory_T& stateTrajectory,
			double dtInitial = 0.01,
			double AbsTol = 1e-9,
			double RelTol = 1e-6) override  {

		typename Base::State_T initialStateInternal = initialState;
		initialize(initialStateInternal, timeTrajectory.front(), dtInitial);
		integrate_times(boost::numeric::odeint::make_controlled<Stepper>(AbsTol, RelTol),
				systemFunction_, initialStateInternal, &timeTrajectory.front(), &timeTrajectory.back()+1, dtInitial, Base::observer_.observeWrap);
		Base::retrieveStateTrajectoryFromObserver(stateTrajectory);
		return true;
	}


private:

	void setupSystem()
	{
		systemFunction_ = [this]( const Eigen::Matrix<double, STATE_DIM, 1>& x, Eigen::Matrix<double, STATE_DIM, 1>& dxdt, double t )
	 									{
			const typename Base::State_T& xState(static_cast<const typename Base::State_T& >(x));
			typename Base::State_T& dxdtState(static_cast<typename Base::State_T& >(dxdt));
			this->system_->computeDerivative(t, xState, dxdtState);
	 									};
	}

	void initialize(const typename Base::State_T& initialState, const double& t, double dt)
	{
		initializeStepper(initialState, t, dt);
		Base::observer_.reset();
	}

	template <class S>
	void integrate_adaptive_helper(
			const typename Base::State_T& initialState,
			const double& startTime,
			const double& finalTime,
			typename Base::StateTrajectory_T& stateTrajectory,
			typename Base::TimeTrajectory_T& timeTrajectory,
			double dtInitial = 0.01,
			double AbsTol = 1e-9,
			double RelTol = 1e-6,
			size_t maxNumSteps = 100000,
			bool bla = true);


	// functionality to reset stepper.
	// in the general case, do not reset
	typedef boost::numeric::odeint::dense_output_runge_kutta <
			boost::numeric::odeint::controlled_runge_kutta <
			boost::numeric::odeint::runge_kutta_dopri5 <
			Eigen::Matrix<double, STATE_DIM, 1>,
			double,
			Eigen::Matrix<double, STATE_DIM, 1>,
			double,
			boost::numeric::odeint::vector_space_algebra > > > rk5;


	template <typename S = Stepper>
	typename std::enable_if<std::is_same<S, rk5>::value, void>::type
	initializeStepper(const typename Base::State_T& initialState, const double& t, double dt)
	{
		stepper_.initialize(initialState, t, dt);
	}


	template <typename S = Stepper>
	typename std::enable_if<!std::is_same<S, rk5>::value, void>::type
	initializeStepper(const typename Base::State_T& initialState, const double& t, double dt) {}


	// Member Variables
	std::function<void (const Eigen::Matrix<double, STATE_DIM, 1>&, Eigen::Matrix<double, STATE_DIM, 1>&, double)> systemFunction_;

	Stepper stepper_;

};


template <size_t STATE_DIM, class Stepper>
template <class S>
void Integrator<STATE_DIM, Stepper>::integrate_adaptive_helper<S>(
		const typename Base::State_T& initialState,
		const double& startTime,
		const double& finalTime,
		typename Base::StateTrajectory_T& stateTrajectory,
		typename Base::TimeTrajectory_T& timeTrajectory,
		double dtInitial,
		double AbsTol,
		double RelTol,
		size_t maxNumSteps,
		typename std::enable_if < std::is_same<Stepper, dense_runge_kutta5_t<STATE_DIM>>::value, bool >::type bla){

	typename Base::State_T initialStateInternal = initialState;
	initialize(initialStateInternal, startTime, dtInitial);
	if (maxNumSteps < std::numeric_limits<size_t>::max()){
		Base::observer_.setMaxNumSteps(maxNumSteps, Base::system_);
	}


	integrate_adaptive(boost::numeric::odeint::make_controlled<Stepper>(AbsTol, RelTol),
			systemFunction_, initialStateInternal, startTime, finalTime, dtInitial, Base::observer_.observeWrap);
}

template<class S>
template <size_t STATE_DIM, class Stepper>
void Integrator<STATE_DIM, Stepper>::integrate_adaptive_helper<S>(
		const typename Base::State_T& initialState,
		const double& startTime,
		const double& finalTime,
		typename Base::StateTrajectory_T& stateTrajectory,
		typename Base::TimeTrajectory_T& timeTrajectory,
		double dtInitial,
		double AbsTol,
		double RelTol,
		size_t maxNumSteps,
		typename std::enable_if<!std::is_same<Stepper, dense_runge_kutta5_t<STATE_DIM>>::value, bool>::type bla) {

	typename Base::State_T initialStateInternal = initialState;
	initialize(initialStateInternal, startTime, dtInitial);
	if (maxNumSteps < std::numeric_limits<size_t>::max()){
		Base::observer_.setMaxNumSteps(maxNumSteps, Base::system_);
	}

	integrate_adaptive(Stepper, systemFunction_, initialStateInternal, startTime, finalTime, dtInitial, Base::observer_.observeWrap);
}


template <size_t STATE_DIM>
using IntegratorEuler = Integrator<
		STATE_DIM,
		boost::numeric::odeint::euler<
		Eigen::Matrix<double, STATE_DIM, 1>,
		double,
		Eigen::Matrix<double, STATE_DIM, 1>,
		double,
		boost::numeric::odeint::vector_space_algebra >
>;


template <size_t STATE_DIM>
using IntegratorRK4 = Integrator<
		STATE_DIM,
		boost::numeric::odeint::runge_kutta4<
		Eigen::Matrix<double, STATE_DIM, 1>,
		double,
		Eigen::Matrix<double, STATE_DIM, 1>,
		double,
		boost::numeric::odeint::vector_space_algebra >
>;


template <size_t STATE_DIM>
using IntegratorRK5Variable = Integrator<
		STATE_DIM,
		dense_runge_kutta5_t<STATE_DIM>
>;


template <size_t STATE_DIM>
using ODE45 = Integrator<
		STATE_DIM,
		boost::numeric::odeint::runge_kutta_dopri5<
		Eigen::Matrix<double, STATE_DIM, 1>,
		double,
		Eigen::Matrix<double, STATE_DIM, 1>,
		double,
		boost::numeric::odeint::vector_space_algebra
		>
>;


template <size_t STATE_DIM, size_t ORDER>
using adams_bashforth_uncontrolled =
		boost::numeric::odeint::adams_bashforth<
		ORDER,
		Eigen::Matrix<double, STATE_DIM, 1>,	// state
		double,	// typename value
		Eigen::Matrix<double, STATE_DIM, 1>,	// derivative
		double, 	// typename time
		boost::numeric::odeint::vector_space_algebra
		> ;


template <size_t STATE_DIM, size_t ORDER>
using IntegratorAdamsBashforth = Integrator <
		STATE_DIM,
		adams_bashforth_uncontrolled<STATE_DIM, ORDER>
>;

} // namespace ocs2

#endif /* OCS2INTEGRATOR_H_ */
