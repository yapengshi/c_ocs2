/*
 *
 * Integrator.h
 *
 *  Created on: 07.10.16
 *      Author: mgiftthaler
 *
 */

#ifndef OCS2_INTEGRATOR_H_
#define OCS2_INTEGRATOR_H_

#include <type_traits>
#include <functional>
#include <cmath>

#include <boost/numeric/odeint.hpp>
#include "eigenIntegration.h"

#include "IntegratorBase.h"


namespace ocs2{


/**
 * Defining the steppers
 */
template <size_t STATE_DIM>
using euler_t = boost::numeric::odeint::euler<
		Eigen::Matrix<double, STATE_DIM, 1>,
		double,
		Eigen::Matrix<double, STATE_DIM, 1>,
		double,
		boost::numeric::odeint::vector_space_algebra >;

template <size_t STATE_DIM>
using runge_kutta_4_t = boost::numeric::odeint::runge_kutta4<
		Eigen::Matrix<double, STATE_DIM, 1>,
		double,
		Eigen::Matrix<double, STATE_DIM, 1>,
		double,
		boost::numeric::odeint::vector_space_algebra >;

template <size_t STATE_DIM>
using runge_kutta_dopri5_t = boost::numeric::odeint::runge_kutta_dopri5 <
		Eigen::Matrix<double, STATE_DIM, 1>,
		double,
		Eigen::Matrix<double, STATE_DIM, 1>,
		double,
		boost::numeric::odeint::vector_space_algebra>;

template <size_t STATE_DIM>
using dense_runge_kutta5_t = boost::numeric::odeint::dense_output_runge_kutta <
		boost::numeric::odeint::controlled_runge_kutta <runge_kutta_dopri5_t<STATE_DIM>> >;

template <size_t STATE_DIM, size_t STEPS>
using adams_bashforth_uncontrolled_t =
		boost::numeric::odeint::adams_bashforth<
		STEPS,
		Eigen::Matrix<double, STATE_DIM, 1>,	// state
		double,									// typename value
		Eigen::Matrix<double, STATE_DIM, 1>,	// derivative
		double, 								// typename time
		boost::numeric::odeint::vector_space_algebra> ;

template <size_t STATE_DIM, size_t STEPS>
using adams_bashforth_moulton_uncontrolled_t =
		boost::numeric::odeint::adams_bashforth_moulton<
		STEPS,
		Eigen::Matrix<double, STATE_DIM, 1>,	// state
		double,									// typename value
		Eigen::Matrix<double, STATE_DIM, 1>,	// derivative
		double, 								// typename time
		boost::numeric::odeint::vector_space_algebra> ;


/**
 * \brief Class Integrator
 */
template <size_t STATE_DIM, class Stepper>
class Integrator : public IntegratorBase<STATE_DIM>
{
public:
	typedef IntegratorBase<STATE_DIM> Base;

	/** Constructor
	 * @param system
	 * @param eventHandler
	 */
	Integrator(
			const std::shared_ptr<SystemBase<STATE_DIM> >& system,
			const std::shared_ptr<EventHandler<STATE_DIM> >& eventHandler = nullptr)
	: IntegratorBase<STATE_DIM>(system, eventHandler){
		setupSystem();
	}

	/** Equidistant integration based on initial and final time as well as step length
	 *
	 * @param initialState
	 * @param startTime
	 * @param finalTime
	 * @param dt
	 * @param stateTrajectory
	 * @param timeTrajectory
	 * @return
	 */
	bool integrate(
			const typename Base::State_T& initialState,
			const double& startTime,
			const double& finalTime,
			double dt,
			typename Base::StateTrajectory_T& stateTrajectory,
			typename Base::TimeTrajectory_T& timeTrajectory) override{

		typename Base::State_T initialStateInternal = initialState;
		double startTime_temp = startTime;

		initialize(initialStateInternal, startTime_temp, dt);
		boost::numeric::odeint::integrate_const(stepper_, systemFunction_, initialStateInternal, startTime, finalTime, dt, Base::observer_.observeWrap);
		Base::retrieveTrajectoriesFromObserver(stateTrajectory, timeTrajectory);
		return true;
	}

	/** Adaptive time integration based on start time and final time
	 *
	 * @param initialState
	 * @param startTime
	 * @param finalTime
	 * @param stateTrajectory
	 * @param timeTrajectory
	 * @param dtInitial
	 * @param AbsTol
	 * @param RelTol
	 * @param maxNumSteps
	 * @return
	 */
	bool integrate(
			const typename Base::State_T& initialState,
			const double& startTime,
			const double& finalTime,
			typename Base::StateTrajectory_T& stateTrajectory,
			typename Base::TimeTrajectory_T& timeTrajectory,
			double dtInitial = 0.01,
			double AbsTol = 1e-9,
			double RelTol = 1e-6,
			size_t maxNumSteps = std::numeric_limits<size_t>::max()) override  {

		integrate_adaptive_specialized<Stepper>( initialState, startTime, finalTime, stateTrajectory, timeTrajectory, dtInitial, AbsTol, RelTol, maxNumSteps );

		Base::retrieveTrajectoriesFromObserver(stateTrajectory, timeTrajectory);

		return true;
	}


	/** Output integration based on a given time trajectory
	 *
	 * @param initialState
	 * @param timeTrajectory
	 * @param stateTrajectory
	 * @param dtInitial
	 * @param AbsTol
	 * @param RelTol
	 * @return
	 */
	bool integrate(
			const typename Base::State_T& initialState,
			const typename Base::TimeTrajectory_T& timeTrajectory,
			typename Base::StateTrajectory_T& stateTrajectory,
			double dtInitial = 0.01,
			double AbsTol = 1e-9,
			double RelTol = 1e-6) override  {

		integrate_times_specialized<Stepper>(initialState, timeTrajectory, stateTrajectory, dtInitial, AbsTol, RelTol);

		return true;
	}


private:

	void setupSystem()
	{
		systemFunction_ = [this]( const Eigen::Matrix<double, STATE_DIM, 1>& x, Eigen::Matrix<double, STATE_DIM, 1>& dxdt, double t ){
			const typename Base::State_T& xState(static_cast<const typename Base::State_T& >(x));
			typename Base::State_T& dxdtState(static_cast<typename Base::State_T& >(dxdt));
			this->system_->computeDerivative(t, xState, dxdtState);
		};
	}

	void initialize(typename Base::State_T& initialState, double& t, double dt)
	{
		initializeStepper(initialState, t, dt);
		Base::observer_.reset();
	}


	template <typename S>
	typename std::enable_if<std::is_same<S, runge_kutta_dopri5_t<STATE_DIM>>::value, void>::type
	integrate_adaptive_specialized(
			const typename Base::State_T& initialState,
			const double& startTime,
			const double& finalTime,
			typename Base::StateTrajectory_T& stateTrajectory,
			typename Base::TimeTrajectory_T& timeTrajectory,
			double dtInitial,
			double AbsTol,
			double RelTol,
			size_t maxNumSteps){

		typename Base::State_T initialStateInternal = initialState;
		double startTime_temp = startTime;

		initialize(initialStateInternal, startTime_temp, dtInitial);

		if (maxNumSteps < std::numeric_limits<size_t>::max()){
			Base::observer_.setMaxNumSteps(maxNumSteps, Base::system_);
		}

		boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_controlled<S>(AbsTol, RelTol), systemFunction_, initialStateInternal, startTime, finalTime, dtInitial, Base::observer_.observeWrap);
	}


	template <typename S>
	typename std::enable_if<!std::is_same<S, runge_kutta_dopri5_t<STATE_DIM>>::value, void>::type
	integrate_adaptive_specialized(
			const typename Base::State_T& initialState,
			const double& startTime,
			const double& finalTime,
			typename Base::StateTrajectory_T& stateTrajectory,
			typename Base::TimeTrajectory_T& timeTrajectory,
			double dtInitial,
			double AbsTol,
			double RelTol,
			size_t maxNumSteps) {

		typename Base::State_T initialStateInternal = initialState;
		double startTime_temp = startTime;

		initialize(initialStateInternal, startTime_temp, dtInitial);

		if (maxNumSteps < std::numeric_limits<size_t>::max()){
			Base::observer_.setMaxNumSteps(maxNumSteps, Base::system_);
		}

		boost::numeric::odeint::integrate_adaptive(stepper_, systemFunction_, initialStateInternal, startTime, finalTime, dtInitial, Base::observer_.observeWrap);
	}


	template <typename S = Stepper>
	typename std::enable_if< std::is_same<S, runge_kutta_dopri5_t<STATE_DIM>>::value, void>::type
	integrate_times_specialized(
			const typename Base::State_T& initialState,
			const typename Base::TimeTrajectory_T& timeTrajectory,
			typename Base::StateTrajectory_T& stateTrajectory,
			double dtInitial,
			double AbsTol,
			double RelTol){

		typename Base::State_T initialStateInternal = initialState;
		double startTime_temp = timeTrajectory.front();

		initialize(initialStateInternal, startTime_temp, dtInitial);

		boost::numeric::odeint::integrate_times(boost::numeric::odeint::make_controlled<S>(AbsTol, RelTol),
				systemFunction_, initialStateInternal, &timeTrajectory.front(), &timeTrajectory.back()+1, dtInitial, Base::observer_.observeWrap);

		Base::retrieveStateTrajectoryFromObserver(stateTrajectory);
	}


	template <typename S = Stepper>
	typename std::enable_if<!std::is_same<S, runge_kutta_dopri5_t<STATE_DIM>>::value, void>::type
	integrate_times_specialized(
			const typename Base::State_T& initialState,
			const typename Base::TimeTrajectory_T& timeTrajectory,
			typename Base::StateTrajectory_T& stateTrajectory,
			double dtInitial,
			double AbsTol,
			double RelTol){

		typename Base::State_T initialStateInternal = initialState;
		double startTime_temp = timeTrajectory.front();

		initialize(initialStateInternal, startTime_temp, dtInitial);

		boost::numeric::odeint::integrate_times(stepper_, systemFunction_, initialStateInternal, &timeTrajectory.front(), &timeTrajectory.back()+1, dtInitial, Base::observer_.observeWrap);

		Base::retrieveStateTrajectoryFromObserver(stateTrajectory);
	}


	/**
	 * Functionality to reset stepper. If we integrate with ODE45, we don't need to reset the stepper, hence specialize empty function
	 */
	template <typename S = Stepper>
	typename std::enable_if<std::is_same<S, runge_kutta_dopri5_t<STATE_DIM>>::value, void>::type
	initializeStepper(typename Base::State_T& initialState, double& t, double dt)
	{
		/**do nothing, runge_kutta_5_t does not have a init method */
	}

	/**
	 * Functionality to reset stepper. If we integrate with some other method, eg. adams_bashforth, we need to reset the stepper, hence specialize with initialize call
	 */
	template <typename S = Stepper>
	typename std::enable_if<!(std::is_same<S, runge_kutta_dopri5_t<STATE_DIM>>::value), void>::type
	initializeStepper(typename Base::State_T& initialState, double& t, double dt) {

		stepper_.initialize(systemFunction_, initialState, t, dt);
	}


	std::function<void (const Eigen::Matrix<double, STATE_DIM, 1>&, Eigen::Matrix<double, STATE_DIM, 1>&, double)> systemFunction_;

	Stepper stepper_;
};


/**
 * Defining the integrators
 */
template <size_t STATE_DIM>
using IntegratorEuler = Integrator<STATE_DIM, euler_t<STATE_DIM>>;

template <size_t STATE_DIM>
using IntegratorRK4 = Integrator<STATE_DIM, runge_kutta_4_t<STATE_DIM>>;

template <size_t STATE_DIM>
using IntegratorRK5Variable = Integrator<STATE_DIM, dense_runge_kutta5_t<STATE_DIM>>;

template <size_t STATE_DIM>
using ODE45 = Integrator<STATE_DIM, runge_kutta_dopri5_t<STATE_DIM>>;

template <size_t STATE_DIM, size_t STEPS>
using IntegratorAdamsBashforth = Integrator < STATE_DIM, adams_bashforth_uncontrolled_t<STATE_DIM, STEPS>>;

template <size_t STATE_DIM, size_t STEPS>
using IntegratorAdamsBashforthMoulton = Integrator < STATE_DIM, adams_bashforth_moulton_uncontrolled_t<STATE_DIM, STEPS>>;

} // namespace ocs2

#endif /* OCS2INTEGRATOR_H_ */
