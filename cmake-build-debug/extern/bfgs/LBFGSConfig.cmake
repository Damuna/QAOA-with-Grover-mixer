## LBFGS CMake configuration file

# library version information
set(LBFGS_VERSION_STRING "1.10.0")
set(LBFGS_VERSION_MAJOR  1)
set(LBFGS_VERSION_MINOR  10)
set(LBFGS_VERSION_PATCH  0)

# import exported targets
if (NOT TARGET LBFGS::lib)
  include("${CMAKE_CURRENT_LIST_DIR}/LBFGSTargets.cmake")
endif ()

# project libraries
set(LBFGS_LIBRARIES LBFGS::lib)
