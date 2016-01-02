/*
 * EventHandler.h
 *
 *  Created on: 18.06.2015
 *      Author: neunertm
 */

#ifndef EVENTHANDLER_H_
#define EVENTHANDLER_H_


template <size_t STATE_DIM>
class EventHandler
{
public:
	typedef Eigen::Matrix<double,STATE_DIM,1> State_T;

	EventHandler() {}
	virtual ~EventHandler() {}

	virtual bool checkEvent(const State_T& state, const double& t) = 0;

	virtual void handleEvent(const State_T& state, const double& t) = 0;

private:
};


#endif /* EVENTHANDLER_H_ */
