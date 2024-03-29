cmake_minimum_required (VERSION 3.5)
project (broomstyx)
list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/Modules")

if(POLICY CMP0074)
	cmake_policy(SET CMP0074 NEW)
endif()

file (GLOB_RECURSE sourcefiles "${PROJECT_SOURCE_DIR}/src/*.cpp")
file (GLOB_RECURSE testfiles "${PROJECT_SOURCE_DIR}/tests/*.cpp")

if (DEFINED ENV{USER_SOURCE_DIR})
    message(STATUS "USER_SOURCE_DIR=$ENV{USER_SOURCE_DIR}")
    file (GLOB_RECURSE userfiles "$ENV{USER_SOURCE_DIR}/*.cpp")
    add_executable (broomstyx ${sourcefiles} ${testfiles} ${userfiles})
else()
    message(STATUS "Environment variable USER_SOURCE_DIR is empty")
    add_executable (broomstyx ${sourcefiles} ${testfiles})
endif()

include_directories (${PROJECT_SOURCE_DIR})
include_directories (${PROJECT_SOURCE_DIR}/src)

# Check for Intel Math Kernel Library
find_package (MKL)

# Check for ViennaCL
find_package (ViennaCL)

# Check for CUDA
include(CheckLanguage)
check_language (CUDA)
if (CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    if (ViennaCL_FOUND)
        target_sources (broomstyx PUBLIC ${PROJECT_SOURCE_DIR}/src/LinearSolvers/ViennaCL_cuda.cu)
    endif()
else()
    message(STATUS "No CUDA support")
endif()

option( ENABLE_OPENMP "Compile with OpenMP" OFF )

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_SOURCE_DIR}/config.h"
)

configure_file (
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
)

# Verbose Makefile
set (CMAKE_VERBOSE_MAKEFILE OFF)

# Set compiler and linker flags
if ( ENABLE_OPENMP )
    set ( OPENMP_CXX_FLAGS "-fopenmp" )
    message(STATUS "Compilation with OpenMP is enabled")
else()
    set ( OPENMP_CXX_FLAGS "" )
    message(STATUS "Compilation with OpenMP is disabled")
endif()

set (CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3 -m64 -std=c++11 -march=native ${OPENMP_CXX_FLAGS} -ffast-math -funroll-loops -Wall -fPIC")
set (CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed -lgomp -lpthread -lm -ldl")
set (CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -std=c++11 -arch=sm_30 -Xcompiler \"${OPENMP_CXX_FLAGS} -ffast-math -funroll-loops -fPIC\"")
