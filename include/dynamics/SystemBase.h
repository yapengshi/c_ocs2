/*
 * SystemBase.h
 *
 *  Created on: 17.06.2015
 *      Author: neunertm
 */

#ifndef SYSTEMBASE_OCS2_H_
#define SYSTEMBASE_OCS2_H_


namespace ocs2{

template <size_t STATE_DIM>
class SystemBase
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	SystemBase()
		: numFunctionCalls_(0) {}

	virtual ~SystemBase() {}

	size_t getNumFunctionCalls() {return numFunctionCalls_;}

	virtual void computeDerivative(
			const double& t,
			const Eigen::Matrix<double,STATE_DIM,1>& state,
			Eigen::Matrix<double,STATE_DIM,1>& derivative) = 0;

protected:
	size_t numFunctionCalls_;

private:
};

} // namespace ocs2

#endif /* SYSTEMBASE_H_ */
