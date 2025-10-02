# main CMake configuration file for the Burgers1D project
cmake_minimum_required(VERSION 3.10)
project(Burgers1D)

find_package(HDF5 REQUIRED COMPONENTS C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER clang++)

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -w -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -w -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")

set(PROJECT_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/)

add_subdirectory(Math)
add_subdirectory(numerics)
add_subdirectory(Physics)
add_subdirectory(IO)
add_subdirectory(Grid)
add_subdirectory(solver)

add_executable(1D_flame
  main.cpp
)

target_include_directories(1D_flame PUBLIC 
  ${PROJECT_INCLUDE_DIRS}
)

target_link_libraries(1D_flame solver grid IO)
