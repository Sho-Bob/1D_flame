#include <iostream>
#include <vector>
#include "grid.h"
#include <solver.h>
#include <solver_burgers.h>
#include <solver_linear_advection.h>
#include <omp.h> 
#include "vtk_writer.h"

int main(){
    const int N = 101;
    const int ibd = 4;
    const int N_total = N + 2*ibd;
    const double x0 = 0.0, x1 = 1.0;
    const double dx = (x1 - x0) / N;
    const double CFL = 0.2;
    const double T_end = 1.0;

    // run solver

    LinearAdvectionSolver solver;
    solver.run_solver();
    return 0;
}
