cmake_minimum_required(VERSION 3.10)
project(Burgers1D)

set(CMAKE_CXX_STANDARD 11)
include_directories(src)

add_executable(Burgers1D 
src/main.cpp 
src/solver.cpp
# src/reconstruction.cpp 
src/vtk_writer.cpp
)
