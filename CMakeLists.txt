cmake_minimum_required(VERSION 3.10)
project(treeutils)

set(CMAKE_CXX_STANDARD 20)



SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1z -O3 -Wall -march=native")

set(CMAKE_BUILD_TYPE Release)



add_executable(treeutils main.cpp define.h newicklex.h node.h util.h newicklex.cpp BipartiteMWIS.h maxflow.h)

