option (HAVE_VIENNACL "Using Iterative solvers from ViennaCL Library" OFF)

find_path (
    VIENNACL_INCLUDE_DIR
    NAMES viennacl/vector.hpp
    HINTS $ENV{VIENNACL_ROOT}
)

set(VIENNACL_INCLUDE_DIRS ${VIENNACL_INCLUDE_DIR})

include (FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(ViennaCL DEFAULT_MSG VIENNACL_INCLUDE_DIRS)

mark_as_advanced (VIENNACL_INCLUDE_DIRS)

if (ViennaCL_FOUND)
    include_directories (${VIENNACL_INCLUDE_DIRS})
    set (HAVE_VIENNACL ON)
else()
    set (HAVE_VIENNACL OFF)
endif()
