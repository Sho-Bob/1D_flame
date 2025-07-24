#ifndef VTK_WRITER_H
#define VTK_WRITER_H

#include <vector>
#include <string>

void write_vtk(const std::vector<double>& u, double dx, const std::string& filename);

#endif