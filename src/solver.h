#ifndef SOLVER_H
#define SOLVER_H

#include <vector>

void simulate_burgers1d(std::vector<double>& u, double dx, double CFL, double T_end);

#endif