/*
 * InterpolationTest.cpp
 *
 *  Created on: Dec 27, 2015
 *  Updated on: Jun 7, 2016
 *      Author: farbod, mgiftthaler
 */

#include <iostream>
#include <cstdlib>

#include "misc/LinearInterpolation.h"
#include <gtest/gtest.h>

using namespace ocs2;

TEST(InterplationTest, Linear)
{
	std::cout << "INTERPOLATION TEST"  << std::endl;
	std::cout << "========================================" << std::endl;
	std::cout << "========================================" << std::endl;

	std::vector<double> timeStamp(5);
	timeStamp[0] = 0.0;
	timeStamp[1] = 0.5;
	timeStamp[2] = 1.0;
	timeStamp[3] = 1.5;
	timeStamp[4] = 2.0;

	std::vector<Eigen::Vector2d,Eigen::aligned_allocator<Eigen::Vector2d> > data(5);
	data[0] << 0.0, 0.0;
	data[1] << 1.0, 0.0;
	data[2] << 2.0, 0.0;
	data[3] << 3.0, 0.0;
	data[4] << 4.0, 0.0;

	LinearInterpolation<Eigen::Vector2d,Eigen::aligned_allocator<Eigen::Vector2d> > linInterpolation(&timeStamp, &data);


	for(int i=1; i< 5; i++)
	{
		double enquiryTime = 0.5*i;

		Eigen::Vector2d enquiryData;
		linInterpolation.interpolate(enquiryTime, enquiryData);

		std::cout << "At time " << enquiryTime << "\t data is " <<  enquiryData.transpose() << std::endl;

		Eigen::Vector2d nominal_result;
		nominal_result << double(i), 0.0;

		ASSERT_LT((nominal_result - enquiryData).array().abs().maxCoeff(), 1e-6);
	}

}


int main(int argc, char** argv)
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
