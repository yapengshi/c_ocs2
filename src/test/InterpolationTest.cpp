/*
 * InterpolationTest.cpp
 *
 *  Created on: Dec 27, 2015
 *      Author: farbod
 */

#include <iostream>
#include <cstdlib>

#include "misc/LinearInterpolation.h"


int main(int argc, char* argv[])
{

	std::vector<double> timeStamp(5);
	timeStamp[0] = 0.0;
	timeStamp[1] = 0.5;
	timeStamp[2] = 1.0;
	timeStamp[3] = 1.5;
	timeStamp[4] = 2.0;

	std::vector<Eigen::Vector2d> data(5);
	data[0] << 0.0, 0.0;
	data[1] << 1.0, 0.0;
	data[2] << 2.0, 0.0;
	data[3] << 3.0, 0.0;
	data[4] << 4.0, 0.0;

	for(size_t i=0; i<5; i++) {
		std::cout << "At time " <<timeStamp[i] << "\t data is " <<  data[i].transpose() << std::endl;
	}
	std::cout << std::endl;

	LinearInterpolation<Eigen::Vector2d> linInterpolation(&timeStamp, &data);

//	linInterpolation.setData(&data);
//	linInterpolation.setTimeStamp(&timeStamp);

	for(size_t i=1; i<argc; i++) {
		double enquiryTime;
		enquiryTime = std::atof(argv[i]);

		Eigen::Vector2d enquiryData;
		linInterpolation.interpolate(enquiryTime, enquiryData);
		std::cout << "At time " << enquiryTime << "\t data is " <<  enquiryData.transpose() << std::endl << std::endl;
	}

}

