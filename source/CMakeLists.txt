cmake_minimum_required(VERSION 2.8)

set( CMAKE_VERBOSE_MAKEFILE on )

project(corv_boson)

add_executable(corv_boson SGFBoson.cpp)

include_directories(${CMAKE_SOURCE_DIR}/Library)
include_directories(${CMAKE_SOURCE_DIR}/Corvette)
include_directories(${CMAKE_SOURCE_DIR}/)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/")

find_package(MKL REQUIRED)
include_directories(${MKL_INCLUDE_DIR})
link_directories(${MKL_LIBRARIES})
target_link_libraries(corv_boson
mkl_intel_lp64
mkl_sequential
mkl_core
pthread
m
)


