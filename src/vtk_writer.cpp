#include "vtk_writer.h"
#include <fstream>
#include <iomanip>

void write_vtk(const std::vector<double>& u, double dx, const std::string& filename) {
    int N = u.size();
    std::ofstream file(filename);
    file << "# vtk DataFile Version 3.0\n";
    file << "Burgers1D Output\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_GRID\n";
    file << "DIMENSIONS " << N << " 1 1\n";
    file << "POINTS " << N << " float\n";
    for (int i = 0; i < N; ++i) {
        file << std::fixed << std::setprecision(6) << i * dx << " 0 0\n";
    }

    file << "POINT_DATA " << N << "\n";
    file << "SCALARS u float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < N; ++i) {
        file << std::fixed << std::setprecision(6) << u[i] << "\n";
    }

    file.close();
}