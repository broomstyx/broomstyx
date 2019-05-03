# Find the Math Kernel Library from Intel
#
#  MKL_FOUND - System has MKL
#  MKL_INCLUDE_DIRS - MKL include files directories
#  MKL_LIBRARIES - The MKL libraries
#  MKL_INTERFACE_LIBRARY - MKL interface library
#  MKL_THREAD_LIBRARY - MKL thread library
#  MKL_CORE_LIBRARY - MKL core library
#
#  The environment variable MKL_ROOT is used to find the library.
#  If MKL is found "-DUSE_MKL" is added to CMAKE_C_FLAGS and CMAKE_CXX_FLAGS.
#
#  Example usage:
#
#  find_package(MKL)
#  if(MKL_FOUND)
#    target_link_libraries(TARGET ${MKL_LIBRARIES})
#  endif()

option (HAVE_MKL "Using BLAS/LAPACK and PARDISO from Intel MKL" OFF)

# If already in cache, be silent
if (MKL_INCLUDE_DIRS AND MKL_LIBRARIES)
    set (MKL_FIND_QUIETLY TRUE)
endif()

set (INT_LIB "mkl_intel_lp64")
set (THR_LIB "mkl_gnu_thread")
set (COR_LIB "mkl_core")

find_path (MKL_INCLUDE_DIR 
    NAMES mkl.h 
    HINTS $ENV{MKL_ROOT}/include
)

find_library (MKL_INTERFACE_LIBRARY
             NAMES ${INT_LIB}
             PATHS $ENV{MKL_ROOT}/lib
                   $ENV{MKL_ROOT}/lib/intel64
                   $ENV{INTEL}/mkl/lib/intel64
             NO_DEFAULT_PATH
)

find_library (MKL_THREAD_LIBRARY
             NAMES ${THR_LIB}
             PATHS $ENV{MKL_ROOT}/lib
                   $ENV{MKL_ROOT}/lib/intel64
                   $ENV{INTEL}/mkl/lib/intel64
             NO_DEFAULT_PATH
)

find_library (MKL_CORE_LIBRARY
             NAMES ${COR_LIB}
             PATHS $ENV{MKL_ROOT}/lib
                   $ENV{MKL_ROOT}/lib/intel64
                   $ENV{INTEL}/mkl/lib/intel64
             NO_DEFAULT_PATH
)

set (MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR})
set (MKL_LIBRARIES ${MKL_INTERFACE_LIBRARY} ${MKL_THREAD_LIBRARY} ${MKL_CORE_LIBRARY})

if (MKL_INCLUDE_DIR AND
    MKL_INTERFACE_LIBRARY AND
    MKL_THREAD_LIBRARY AND
    MKL_CORE_LIBRARY)

#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_MKL_BLAS -DUSE_MKL_PARDISO")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_MKL_BLAS -DUSE_MKL_PARDISO")
else()
    set (MKL_INCLUDE_DIRS "")
    set (MKL_LIBRARIES "")
    set (MKL_INTERFACE_LIBRARY "")
    set (MKL_THREAD_LIBRARY "")
    set (MKL_CORE_LIBRARY "")
endif()

# Handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE if
# all listed variables are TRUE.
include (FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MKL DEFAULT_MSG MKL_LIBRARIES MKL_INCLUDE_DIRS MKL_INTERFACE_LIBRARY MKL_THREAD_LIBRARY MKL_CORE_LIBRARY)

mark_as_advanced (MKL_INCLUDE_DIRS MKL_LIBRARIES MKL_INTERFACE_LIBRARY MKL_THREAD_LIBRARY MKL_CORE_LIBRARY)

if (MKL_FOUND)
#    include_directories ("${MKL_INCLUDE_DIRS}")
#    target_link_libraries(broomstyx ${MKL_LIBRARIES})
    set (HAVE_MKL ON)
else()
    set (HAVE_MKL OFF)
endif()
