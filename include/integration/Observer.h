/*
 * Observer.h
 *
 *  Created on: 18.06.2015
 *      Author: neunertm
 */

#ifndef OCS2OBSERVER_H_
#define OCS2OBSERVER_H_

#include "EventHandler.h"


template <size_t STATE_DIM>
class IntegratorBase;

template <size_t STATE_DIM>
class Observer
{
public:
	typedef std::vector<double> TimeTrajectory_T;
	typedef Eigen::Matrix<double,STATE_DIM,1> State_T;
	typedef std::vector<State_T, Eigen::aligned_allocator<State_T> > StateTrajectory_T;

	friend class IntegratorBase<STATE_DIM>;

	Observer(const std::shared_ptr<EventHandler<STATE_DIM> >& eventHandler = nullptr) :
		observeWrap([this](const State_T& x, const double& t){ this->observe(x,t); }),
		eventHandler_(eventHandler)
	{}

	void reset() {
		stateTrajectory_.clear();
		timeTrajectory_.clear();
	}

	void observe(const State_T& x, const double& t)
	{
		stateTrajectory_.push_back(x);
		timeTrajectory_.push_back(t);

		if (eventHandler_ && eventHandler_->checkEvent(x, t))
		{
			eventHandler_->handleEvent(x, t);
		}
	}

	// Lambda to pass to odeint (odeint takes copies of the observer so we can't pass the class
	std::function<void (const State_T& x, const double& t)> observeWrap;

private:
	std::shared_ptr<EventHandler<STATE_DIM> > eventHandler_;

	StateTrajectory_T stateTrajectory_;
	TimeTrajectory_T timeTrajectory_;

};



#endif /* OCS2OBSERVER_H_ */
