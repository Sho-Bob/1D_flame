add_library(solver 
  solver.cpp
  solver_burgers.cpp
  solver_navier_stokes.cpp
  solver_linear_advection.cpp
)

target_include_directories(solver PUBLIC 
  ${PROJECT_INCLUDE_DIRS}
)

target_link_libraries(solver numerics IO grid common physics)
