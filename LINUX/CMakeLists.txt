cmake_minimum_required(VERSION 3.4...3.18)

project(SolverDyns)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

include_directories(../include)
include_directories(../include/fparser)
include_directories(../include/eigen)
include_directories(../include/igl)

add_executable (SolverDyns ../main.cpp ../DynamicSystem.cpp ../DynamicSystem.h ../include/fparser/fparser.cc ../include/fparser/fpoptimizer.cc)