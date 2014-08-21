cmake_minimum_required(VERSION 2.6)
project(mwSynth) 

set(CMAKE_CXX_FLAGS "-std=c++11 -mavx2")

file(GLOB SRC_FILES "*.cpp")
add_executable(mwSynthDemo ${SRC_FILES})
