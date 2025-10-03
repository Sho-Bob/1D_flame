add_library(numerics
  reconstruction.cpp
  riemann_solver.cpp
)

target_include_directories(numerics PUBLIC
  ${PROJECT_INCLUDE_DIRS}
)


