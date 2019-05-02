cmake_minimum_required (VERSION 3.1)
project (broomstyx CXX)
list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/Modules")

# Check for Intel Math Kernel Library
find_package (MKL)

if(NOT (dune-common_DIR OR dune-common_ROOT OR
      "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common)

if( dune-common_FOUND )
  #############################################
  ##  DUNE build system
  #############################################
  list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH})

  #include the dune macros
  include(DuneMacros)

  find_package( BLAS )
  find_package(LAPACK )
  find_package( GSL )

  if( GSL_FOUND )
    dune_register_package_flags(INCLUDE_DIRS ${GSL_INCLUDE_DIRS} LIBRARIES ${GSL_LIBRARIES})
  endif()

  # start a dune project with information from dune.module
  dune_project()

  dune_enable_all_packages()

  #add_subdirectory(doc)
  add_subdirectory(cmake/modules)

  file (GLOB_RECURSE sourcefiles "${PROJECT_SOURCE_DIR}/src/*.cpp")

  add_executable("broomstyx" ${sourcefiles})
  target_link_dune_default_libraries("broomstyx")

  # finalize the dune project, e.g. generating config.h etc.
  finalize_dune_project(GENERATE_CONFIG_H_CMAKE)

else()

  list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules"
    ${dune-common_MODULE_PATH})

  file (GLOB_RECURSE sourcefiles "${PROJECT_SOURCE_DIR}/src/*.cpp")
  add_executable (broomstyx ${sourcefiles})

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
  set (CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3 -m64 -std=c++11 -march=native -fopenmp -ffast-math -funroll-loops -Wall -fPIC")
  set (CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed -lgomp -lpthread -lm -ldl")
  set (CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -std=c++11 -arch=sm_30 -Xcompiler \"-fopenmp -ffast-math -funroll-loops -fPIC\"")
endif()

