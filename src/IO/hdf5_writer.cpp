#include "vtk_writer.h"
#include "hdf5_writer.h"
#include <fstream>
#include <sstream>
#include <iomanip>

namespace IO {

  void create_hdf5_file(const std::string& filename) {
    // Create a new file using the default property lists.
    H5File file(filename, H5F_ACC_TRUNC);
    file.close();
  }

void write_structured_mesh_timestep(const std::string& filename,
                                    int timestep,
                                    const std::vector<double>& pressure,
                                    hsize_t nx, hsize_t ny, hsize_t nz)
{
    // Open or create file
    H5File file(filename, H5F_ACC_RDWR);
    // Create a group for this timestep
    std::ostringstream grpname;
    grpname << "/Timestep" << std::setw(4) << std::setfill('0') << timestep;
    Group group = file.createGroup(grpname.str());

    // Data dimensions
    hsize_t dims[3] = {nx, ny, nz};
    DataSpace space(3, dims);

    // Create and write pressure dataset
    DataSet dset = group.createDataSet("pressure",
                                       PredType::NATIVE_DOUBLE,
                                       space);
    dset.write(pressure.data(), PredType::NATIVE_DOUBLE);
} // end of write_structured_mesh_timestep

} // namespace IO

