/*
 * SystemBase.h
 *
 *  Created on: 17.06.2015
 *      Author: neunertm
 */

#ifndef SYSTEMBASE_H_
#define SYSTEMBASE_H_


template <size_t STATE_DIM>
class SystemBase
{
public:
	SystemBase() {}
	virtual ~SystemBase() {}

	virtual void computeDerivative(
			const double& t,
			const Eigen::Matrix<double,STATE_DIM,1>& state,
			Eigen::Matrix<double,STATE_DIM,1>& derivative) = 0;

private:
};



#endif /* SYSTEMBASE_H_ */
