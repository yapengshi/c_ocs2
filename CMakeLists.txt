cmake_minimum_required(VERSION 2.8.3)
project(c_ocs2) 

## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS cmake_modules cereal c_gradient_descent) 
find_package(Eigen3 REQUIRED)
find_package(Ceres REQUIRED)
find_package(Boost REQUIRED system filesystem)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -fopenmp -std=c++11 -Wfatal-errors -Wl,--no-as-needed")
# add flags "-Wl,--no-as-needed" if using gcc 4.8.2 or lower


###################################
## ipopot specific configuration ##
###################################

#SET(IPOPT_SOURCE_DIR "/home/farbod/Programs/catkin_ws/src/ipopt")
SET(IPOPT_SOURCE_DIR "/home/gimarkus/ipopt/CoinIpopt")


link_directories(${IPOPT_SOURCE_DIR}/build/lib)
set(IPOPT_LIBS
    ipopt
    dl
    coinmumps
    coinhsl
    lapack
    blas
    gfortran
    m
    quadmath
    coinmetis
)

###################################
## catkin specific configuration ##
###################################
## The catkin_package macro generates cmake config files for your package
## Declare things to be passed to dependent projects
## INCLUDE_DIRS: uncomment this if you package contains header files
## LIBRARIES: libraries you create in this project that dependent projects also need
## CATKIN_DEPENDS: catkin_packages dependent projects also need
## DEPENDS: system dependencies of this project that dependent projects also need
catkin_package(
  INCLUDE_DIRS include util
      ${EIGEN_INCLUDE_DIRS}
      "${IPOPT_SOURCE_DIR}/build/include/coin"
      ${CERES_INCLUDE_DIRS}
      ${Boost_INCLUDE_DIRS}
  #LIBRARIES c_ocs2 ${Boost_LIBRARIES}
  CATKIN_DEPENDS cereal c_gradient_descent
  DEPENDS ${IPOPT_LIBS} boost eigen
)


########### 
## Build ##
###########

## Specify additional locations of header files
## Your package locations should be listed before other locations
# include_directories(include)
include_directories(
  include
  util
  ${catkin_INCLUDE_DIRS}
  ${EIGEN_INCLUDE_DIRS}
  ${CERES_INCLUDE_DIRS}
  "${IPOPT_SOURCE_DIR}/build/include/coin"  
  ${Boost_INCLUDE_DIRS}
)

## Declare a cpp library
# add_library(c_ocs2
#   src/${PROJECT_NAME}/c_ocs2.cpp
# ) 

## Declare a cpp executable
add_executable(exp1_slqp_mp_test src/test/EXP1_SLQP_MP_Test.cpp)

## Declare a cpp executable
add_executable(exp3_slqp_mp_test src/test/EXP3_SLQP_MP_Test.cpp)

## Declare a cpp executable
add_executable(exp1_gslqp_test src/test/EXP1GSLQPTest.cpp)

## Declare a cpp executable
add_executable(exp3_gslqp_test src/test/EXP3GSLQPTest.cpp)

## Declare a cpp executable
add_executable(exp4_gslqp_test src/test/EXP4GSLQPTest.cpp)

## Declare a cpp executable
add_executable(exp1_test src/test/EXP1Test.cpp)

## Declare a cpp executable
add_executable(exp2_gslqp_test src/test/EXP2GSLQPTest.cpp)

## Declare a cpp executable
add_executable(exp2_test src/test/EXP2Test.cpp)

## Declare a cpp executable 
add_executable(exp1_cost_derivative_test src/test/EXP1DerivativesCheck.cpp)

## Declare a cpp executable 
add_executable(exp3_numerical_derivative_test src/test/EXP3SystemDynamicsDerivativeTest.cpp)


## Add cmake target dependencies of the executable/library
## as an example, message headers may need to be generated before nodes
# add_dependencies(c_ocs2_node c_ocs2_generate_messages_cpp)


## Specify libraries to link a library or executable target against
target_link_libraries(exp1_slqp_mp_test
  ${catkin_LIBRARIES}
  ${Boost_LIBRARIES}
  pthread
)

## Specify libraries to link a library or executable target against
target_link_libraries(exp3_slqp_mp_test
  ${catkin_LIBRARIES}
  ${Boost_LIBRARIES}
  pthread
)

## Specify libraries to link a library or executable target against
target_link_libraries(exp1_gslqp_test
  ${catkin_LIBRARIES}
  ${Boost_LIBRARIES}
)

## Specify libraries to link a library or executable target against
target_link_libraries(exp3_gslqp_test
  ${catkin_LIBRARIES}
  ${Boost_LIBRARIES}
)

## Specify libraries to link a library or executable target against
target_link_libraries(exp4_gslqp_test
  ${catkin_LIBRARIES}  
  ${Boost_LIBRARIES}
)

## Specify libraries to link a library or executable target against
target_link_libraries(exp1_test
  ${catkin_LIBRARIES}
  ${IPOPT_LIBS}
  ${Boost_LIBRARIES}
)

## Specify libraries to link a library or executable target against
target_link_libraries(exp2_gslqp_test
  ${catkin_LIBRARIES}
  ${Boost_LIBRARIES}
)

## Specify libraries to link a library or executable target against
target_link_libraries(exp2_test
  ${catkin_LIBRARIES}
  ${IPOPT_LIBS}
  ${Boost_LIBRARIES}
)

## Specify libraries to link a library or executable target against
target_link_libraries(exp1_cost_derivative_test
  ${catkin_LIBRARIES}
  ${IPOPT_LIBS}
  ${Boost_LIBRARIES}
)

## Specify libraries to link a library or executable target against
target_link_libraries(exp3_numerical_derivative_test
  ${catkin_LIBRARIES}  
  ${Boost_LIBRARIES}
)

#############
## Testing ##
#############

## Info ==============================
## to run tests, cd package/src and run 
## $ catkin run_tests --no-deps --this

## unit tests
catkin_add_gtest(interpolation_test test/InterpolationTest.cpp)

catkin_add_gtest(integration_test test/IntegrationTest.cpp)
target_link_libraries(integration_test ${Boost_LIBRARIES})

catkin_add_gtest(system_dynamics_derivative_test test/SystemDynamicsDerivativeTest.cpp)

catkin_add_gtest(ceres_test test/CeresTest.cpp)

target_link_libraries(ceres_test
  ${catkin_LIBRARIES}
  ${CERES_LIBRARIES}
)
