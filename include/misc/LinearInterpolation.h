/*
 * LinearInterpolation.h
 *
 *  Created on: Dec 27, 2015
 *      Author: farbod
 */

#ifndef LINEARINTERPOLATION_OCS2_H_
#define LINEARINTERPOLATION_OCS2_H_

#include <Eigen/Dense>
#include <Eigen/StdVector>

#include <memory>
#include <vector>


namespace ocs2{

template <typename Data_T, class Alloc=std::allocator<Data_T> >
class LinearInterpolation
{
public:

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	LinearInterpolation()
		: index_(0),
		  timeStampPtr_(NULL),
		  dataPtr_(NULL)
	{}

	LinearInterpolation(const std::vector<double>* timeStampPtr, const std::vector<Data_T,Alloc>* dataPtr)
		: index_(0),
		  timeStampPtr_(timeStampPtr),
		  dataPtr_(dataPtr)
	{}


	void setData(const std::vector<Data_T,Alloc>* dataPtr)	{

		dataPtr_ = dataPtr;
		reset();
	}

	void setTimeStamp(const std::vector<double>* timeStampPtr)	{

		timeStampPtr_ = timeStampPtr;
		reset();
	}

	void reset(){index_=0;}

	void interpolate(const double& enquiryTime, Data_T& enquiryData, int greatestLessTimeStampIndex = -1) {

		if (timeStampPtr_==NULL)  throw std::runtime_error("timeStampPtr is not initialized.");
		if (dataPtr_==NULL)       throw std::runtime_error("dataPtr is not initialized.");

		if (timeStampPtr_->size()==0)  				  throw std::runtime_error("LinearInterpolation is not initialized.");
		if (timeStampPtr_->size()!=dataPtr_->size())  throw std::runtime_error("The size of timeStamp vector is not equal to the size of data vector.");

		if (timeStampPtr_->size()==1)  {
			enquiryData = dataPtr_->front();
			return;
		}

//		std::cout << "INDEX: " << index_;
		size_t ind;
		if (greatestLessTimeStampIndex == -1)
			ind = find(enquiryTime);
		else {
			ind = greatestLessTimeStampIndex;
			index_ = greatestLessTimeStampIndex;
		}
//		std::cout << " --> " << ind << std::endl;

		if (enquiryTime<timeStampPtr_->front()) {
			enquiryData = dataPtr_->front();
			return;
		}

		if (ind==timeStampPtr_->size()-1) {
			enquiryData = dataPtr_->back();
			return;
		}

//		if (dataPtr_->size() <= ind+1)
//			std::cout << " --> " << ind << std::endl;
		double alpha = (enquiryTime-timeStampPtr_->at(ind+1)) / (timeStampPtr_->at(ind)-timeStampPtr_->at(ind+1));
		enquiryData = alpha*dataPtr_->at(ind) + (1-alpha)*dataPtr_->at(ind+1);
	}

	size_t getGreatestLessTimeStampIndex() { return index_; }

protected:
	size_t find(const double& enquiryTime) {

		size_t index;

		if (timeStampPtr_->at(index_) > enquiryTime) {
			for (int i=index_; i>=0; i--)  {
				index = i;
				if (timeStampPtr_->at(i) <= enquiryTime)
					break;
			}
		} else {
			for (int i=index_; i<timeStampPtr_->size(); i++) {
				index = i;
				if (timeStampPtr_->at(i) > enquiryTime) {
					index = i-1;
					break;
				}
			}
		}

		index_ = index;
		return index;
	}

private:
	size_t index_;

	const std::vector<double>* timeStampPtr_;
	const std::vector<Data_T,Alloc>* dataPtr_;

};


} // namespace ocs2

#endif /* LINEARINTERPOLATION_H_ */
