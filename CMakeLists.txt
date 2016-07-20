cmake_minimum_required(VERSION 3.5)
project(dpdb)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp DPDB.cpp DPDB.h)
add_executable(dpdb ${SOURCE_FILES})