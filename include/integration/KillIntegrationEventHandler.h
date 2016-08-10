/*
 * KillIntegrationEventHandler.h
 *
 *  Created on: 09.08.2015
 *      Author: mgiftthaler
 */

#ifndef OCS2_KILL_INTEGRATION_EVENTHANDLER_H_
#define OCS2_KILL_INTEGRATION_EVENTHANDLER_H_

#include <integration/EventHandler.h>

namespace ocs2{

template <size_t STATE_DIM>
class KillIntegrationEventHandler : public EventHandler<STATE_DIM>
{
public:
	typedef Eigen::Matrix<double,STATE_DIM,1> State_T;

	KillIntegrationEventHandler():
		killIntegration_(false)
	{}

	~KillIntegrationEventHandler() {}

	bool checkEvent(const State_T& state, const double& t) override {
		return killIntegration_;
	}

	void handleEvent(const State_T& state, const double& t) override {

		/* throw an exception which stops the integration */
		throw std::runtime_error("Integration terminated due to external event specified by user.");
	}

	void setEvent() {
		killIntegration_ = true;
	}

	void resetEvent() {
		killIntegration_ = false;
	}

private:
	bool killIntegration_;
};

} // namespace ocs2

#endif
