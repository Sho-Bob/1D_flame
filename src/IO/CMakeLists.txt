add_library(IO
  vtk_writer.cpp vtk_writer.h
  hdf5_writer.cpp hdf5_writer.h
)

target_link_libraries(IO PRIVATE ${HDF5_LIBRARIES})
target_include_directories(IO PUBLIC ${HDF5_INCLUDE_DIRS})
