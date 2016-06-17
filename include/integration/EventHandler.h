/*
 * EventHandler.h
 *
 *  Created on: 18.06.2015
 *      Author: neunertm
 */

#ifndef OCS2_EVENTHANDLER_H_
#define OCS2_EVENTHANDLER_H_

namespace ocs2{

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

} // namespace ocs2

#endif /* OCS2EVENTHANDLER_H_ */
