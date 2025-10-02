add_library(IO
  vtk_writer.cpp
  hdf5_writer.cpp
  input.cpp
)

find_package(tomlplusplus REQUIRED)
message(STATUS "tomlplusplus include dirs: ${tomlplusplus_INCLUDE_DIRS}")

target_link_libraries(IO PRIVATE 
  ${HDF5_LIBRARIES}
  tomlplusplus::tomlplusplus
)
target_include_directories(IO PUBLIC 
  ${PROJECT_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  /usr/local/include
)

