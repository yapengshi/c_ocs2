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

	LinearInterpolation<Eigen::Vector2d> linInterpolation(&timeStamp, &data);

//	linInterpolation.setData(&data);
//	linInterpolation.setTimeStamp(&timeStamp);

	double enquiryTime;
	enquiryTime = std::atof(argv[1]);

	Eigen::Vector2d enquiryData;
	linInterpolation.interpolate(enquiryTime, enquiryData);
	std::cout << "At time " << enquiryTime << " data is " <<  enquiryData.transpose() << std::endl;

}

