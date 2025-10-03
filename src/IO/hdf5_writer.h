#ifndef HDF5_WRITER_H
#define HDF5_WRITER_H

#include <H5Cpp.h>

using namespace H5;

namespace IO {

  void create_hdf5_file(const std::string& filename);


  void write_structured_mesh_timestep(const std::string& filename,
      int timestep,
      const std::vector<double>& pressure,
      hsize_t nx, hsize_t ny, hsize_t nz,
      const std::string var_name = "u"
      );

  void write_structured_mesh_timestep(const std::string& filename,
      const double* data,
      hsize_t nx, hsize_t ny, hsize_t nz,
      const std::string var_name = "u"
      );
} // namespace IO

#endif
