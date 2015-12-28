/*
 * EventHandler.h
 *
 *  Created on: 18.06.2015
 *      Author: neunertm
 */

#ifndef EVENTHANDLER_H_
#define EVENTHANDLER_H_

#include <core/state/StateBase.h>

namespace ct {
namespace core {
namespace integration {

template <size_t STATE_DIM>
class EventHandler
{
public:
	EventHandler() {}
	virtual ~EventHandler() {}

	virtual bool checkEvent(const state::StateBase<STATE_DIM>& state, const Time& t) = 0;

	virtual void handleEvent(const state::StateBase<STATE_DIM>& state, const Time& t) = 0;

private:
};

}
}
}

#endif /* EVENTHANDLER_H_ */
