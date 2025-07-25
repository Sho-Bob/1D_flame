#ifndef SOLVER_H
#define SOLVER_H

#include <vector>

void burgers_initialize(std::vector<double>& u, int N, double x0, double x1, double dx, int ibd);
void simulate_burgers1d(std::vector<double>& u, double dx, double CFL, double T_end, int ibd);
void Apply_BC(std::vector<double>& u, int N, int ibd);
void reconstruct_MUSCL_minmod(const std::vector<double>& u, std::vector<double>& uL, std::vector<double>& uR, int N, int ibd);
void reconstruct_WENO5(const std::vector<double>& u, std::vector<double>& uL, std::vector<double>& uR, int N, int ibd);
#endif