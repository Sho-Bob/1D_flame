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

void write_dat(const std::string& filename, const std::vector<std::vector<double>>& data, const std::vector<std::string>& var_names) {
    std::ofstream file(filename);
    int N = data[0].size();
    int num_vars = data.size();

    // Write header
    for (int v = 0; v < num_vars; ++v) {
        file << var_names[v];
        if (v < num_vars - 1) file << "\t";
    }
    file << "\n";

    // Write data
    for (int i = 0; i < N; ++i) {
        for (int v = 0; v < num_vars; ++v) {
            file << std::fixed << std::setprecision(6) << data[v][i];
            if (v < num_vars - 1) file << "\t";
        }
        file << "\n";
    }

    file.close();
}
