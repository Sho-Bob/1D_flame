add_library(solver 
  solver.cpp solver.h
  solver_burgers.cpp solver_burgers.h
  solver_navier_stokes.cpp solver_navier_stokes.h
  solver_linear_advection.cpp solver_linear_advection.h
)

target_include_directories(solver PUBLIC 
  ../IO
  ../numerics
  ../Grid
  ../Common/
)

target_link_libraries(solver numerics IO grid)
