# Created by the script cgal_create_cmake_script_with_options
# This is the CMake script for compiling a set of CGAL applications.

project(MA)


cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
 
if ( COMMAND cmake_policy )

  cmake_policy( SET CMP0003 NEW )  

endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(CImg REQUIRED)
include_directories(${CIMG_INCLUDE_DIRS})
link_libraries(${CIMG_PROCESS_LIBS})

find_package(SuiteSparse REQUIRED)
include_directories(${CSPARSE_INCLUDE_DIR})
link_libraries(${SPQR_LIBRARY} ${CHOLMOD_LIBRARIES} lapack blas colamd)
link_libraries(${CSPARSE_LIBRARY})
link_libraries(${SUITESPARSECONFIG_LIBRARY})

find_package(TBB)
if(TBB_FOUND)
  include_directories(${TBB_INCLUDE_DIRS})
  link_libraries(${TBB_LIBRARIES})
endif()


set(Boost_USE_STATIC_LIBS ON)
find_package(Boost COMPONENTS timer chrono REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_libraries(${Boost_LIBRARIES})

# CGAL and its components
find_package( CGAL QUIET COMPONENTS)


add_subdirectory(thirdparty/lbfgspp)
link_libraries(lbfgspp)
include_directories(${MA_SOURCE_DIR}/thirdparty/lbfgspp/include)

include_directories(${MA_SOURCE_DIR}/include)
subdirs(tests)
ADD_DEFINITIONS(-std=c++11 -pg)
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



