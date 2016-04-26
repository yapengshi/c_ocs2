/*
 * IntegratorBase.h
 *
 *  Created on: 17.12.2015
 *      Author: farbodf
 */

#ifndef OCS2INTEGRATORBASE_H_
#define OCS2INTEGRATORBASE_H_

#include "dynamics/SystemBase.h"
#include "Observer.h"
#include "EventHandler.h"

template <size_t STATE_DIM>
class IntegratorBase
{
public:
	typedef std::vector<double> TimeTrajectory_T;
	typedef Eigen::Matrix<double,STATE_DIM,1> State_T;
	typedef std::vector<State_T, Eigen::aligned_allocator<State_T> > StateTrajectory_T;

	IntegratorBase(
		const std::shared_ptr<SystemBase<STATE_DIM> >& system,
		const std::shared_ptr<EventHandler<STATE_DIM> >& eventHandler = nullptr) :
			observer_(eventHandler),
			system_(system),
			eventHandler_(eventHandler)
	{}

	virtual ~IntegratorBase() {}

	virtual void reset() {
		observer_.reset();
	}

	// Equidistant integration based on number of time steps and step length
	virtual bool integrate(
		const State_T& initialState,
		const double& startTime,
		size_t numSteps,
		double dt,
		StateTrajectory_T& stateTrajectory,
		TimeTrajectory_T& timeTrajectory
	) = 0;

	// Equidistant integration based on initial and final time as well as step length
	virtual bool integrate(
		const State_T& initialState,
		const double& startTime,
		const double& finalTime,
		double dt,
		StateTrajectory_T& stateTrajectory,
		TimeTrajectory_T& timeTrajectory
	) = 0;

	// Adaptive time integration based on start time and final time
	virtual bool integrate(
		const State_T& initialState,
		const double& startTime,
		const double& finalTime,
		StateTrajectory_T& stateTrajectory,
		TimeTrajectory_T& timeTrajectory,
		double dtInitial = 0.01,
		double AbsTol = 1e-6,
		double RelTol = 1e-3
	) = 0;

	// Output integration based on a given time trajectory
	virtual bool integrate(
		const State_T& initialState,
		const TimeTrajectory_T& timeTrajectory,
		StateTrajectory_T& stateTrajectory,
		double dtInitial = 0.01,
		double AbsTol = 1e-9,
		double RelTol = 1e-6
	) = 0;


protected:
	 void retrieveTrajectoriesFromObserver(StateTrajectory_T& stateTrajectory, TimeTrajectory_T& timeTrajectory)
	 {
		 stateTrajectory.swap(observer_.stateTrajectory_);
		 timeTrajectory.swap(observer_.timeTrajectory_);
	 }
	 void retrieveStateTrajectoryFromObserver(StateTrajectory_T& stateTrajectory)
	 {
		 stateTrajectory.swap(observer_.stateTrajectory_);
	 }

	Observer<STATE_DIM> observer_;

	std::shared_ptr<SystemBase<STATE_DIM> > system_;
	std::shared_ptr<EventHandler<STATE_DIM> > eventHandler_;
};



#endif /* OCS2INTEGRATORBASE_H_ */
