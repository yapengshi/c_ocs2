/*
 * PathTweaker.h
 *	A hack that allows users to save data to paths in the catkin root directory
 *	Locates where the executable lies and takes 3 steps back to catkin root directory - returns this absolute path
 *
 *  Created on: June 08, 2016
 *      Author: mgiftthaler
 */

#include <boost/filesystem.hpp>
#include <string>


class PathTweaker
{
public:
	PathTweaker(char *argv[])
	{
		executable_path_= boost::filesystem::system_complete(argv[0]);
		executable_path_.normalize();

		catkin_root_path_ = executable_path_;

		// reduce path to base directory of catkin workspace (3 steps back)
		for (size_t i = 0; i<=3; i++)
		{
			catkin_root_path_ = catkin_root_path_.parent_path();
		}
	}

	std::string getDirectory () {return catkin_root_path_.string();}

private:
	boost::filesystem::path executable_path_;
	boost::filesystem::path catkin_root_path_;
};
