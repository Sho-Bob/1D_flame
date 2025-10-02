#ifndef VTK_WRITER_H
#define VTK_WRITER_H

#include <vector>
#include <string>

void write_vtk(const std::vector<double>& u, double dx, const std::string& filename);
void write_dat(const std::string& filename, const std::vector<std::vector<double>>& data, const std::vector<std::string>& var_names);

#endif
