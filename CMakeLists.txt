cmake_minimum_required(VERSION 1.0)

project(overlap)

add_executable(overlap main.cpp overlap.cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g")