/*
 * LinearInterpolation.h
 *
 *  Created on: Dec 27, 2015
 *      Author: farbod
 */

#ifndef LINEARINTERPOLATION_H_
#define LINEARINTERPOLATION_H_

#include <Eigen/Dense>
#include <vector>

template <typename Data_T>
class LinearInterpolation
{
public:
	LinearInterpolation(std::vector<double>* const timeStampPtr,
			std::vector<Data_T>* const dataPtr)
		: index_(0),
		  timeStampPtr_(timeStampPtr),
		  dataPtr_(dataPtr)
	{}

	void setData(std::vector<Data_T>* const dataPtr)	{dataPtr_ = dataPtr;}

	void setTimeStamp(std::vector<double>* const timeStampPtr)	{timeStampPtr_ = timeStampPtr;}

	void reset()	{index_=0;}


	size_t find(const double& enquiryTime) {

		if (timeStampPtr_->at(index_) > enquiryTime) {
			for (size_t i=index_-1; i>=0; i--)
				if (timeStampPtr_->at(i) <= enquiryTime) {
					index_ = i;
					break;
				}
		} else {
			for (size_t i=index_; i<timeStampPtr_->size(); i++) {
				if (timeStampPtr_->at(i) > enquiryTime) {
					index_ = i-1;
					break;
				} else if (timeStampPtr_->at(i) == enquiryTime) {
					index_ = i;
					break;
				}
			}
		}

		return index_;
	}

	void interpolate(const double& enquiryTime, Data_T& enquiryData) {

//		std::cout << "INDEX: " << index_;
		size_t ind = find(enquiryTime);
//		std::cout << " --> " << ind << std::endl;

		if (ind==timeStampPtr_->size()-1) {
			enquiryData = dataPtr_->back();
			return;
		}

		double alpha = (enquiryTime-timeStampPtr_->at(ind+1)) / (timeStampPtr_->at(ind)-timeStampPtr_->at(ind+1));
		enquiryData = alpha*dataPtr_->at(ind) + (1-alpha)*dataPtr_->at(ind+1);
	}


private:
	size_t index_;

	std::vector<Data_T>* dataPtr_;
	std::vector<double>* timeStampPtr_;
};



#endif /* LINEARINTERPOLATION_H_ */
