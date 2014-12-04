# Created by the script cgal_create_cmake_script_with_options
# This is the CMake script for compiling a set of CGAL applications.

project(MA)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "RELWITHDEBINFO")
ENDIF()

ADD_DEFINITIONS(-std=c++11 -g)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(CImg REQUIRED)
include_directories(${CIMG_INCLUDE_DIRS})
link_libraries(${CIMG_PROCESS_LIBS})

find_package(SuiteSparse)
include_directories(${SUITESPARSE_INCLUDE_DIRS})
link_libraries(${SUITESPARSE_LIBRARIES})

if(SUITESPARSEQR_FOUND)
  add_definitions(-DMA_USE_SUITESPARSE_QR)
endif()

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost COMPONENTS timer chrono REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_libraries(${Boost_LIBRARIES})

# CGAL and its components
find_package( CGAL REQUIRED)


include_directories(${MA_SOURCE_DIR}/include)
subdirs(tests)

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# include for local directory

# include for local package


# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )



