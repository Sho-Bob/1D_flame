#include <iostream>
#include <vector>
#include <omp.h> 
#include <string>

#include "Solver/solver.h"
#include "Solver/solver_burgers.h"
#include "Solver/solver_linear_advection.h"
#include "Solver/solver_navier_stokes.h"
#include "IO/input.h"

int main(){
    const int N = 101;
    const int ibd = 4;
    const int N_total = N + 2*ibd;
    const double x0 = 0.0, x1 = 1.0;
    const double dx = (x1 - x0) / N;
    const double CFL = 0.2;
    const double T_end = 1.0;

    // run solver
    const std::string filename = "input.toml";
    Input input(filename);

    std::string solver_name = input.getStringParam("solver");

    Solver* solver = nullptr;

    if (solver_name == "burgers") {
        solver = new BurgerSolver;
    } else if (solver_name == "linear-advection") {
        solver = new LinearAdvectionSolver;
    } else if (solver_name == "navier-stokes") {
        solver = new NavierStokesSolver;
    } else {
        std::cerr << "Unknown solver: " << solver_name << std::endl;
        return -1;
    }

    solver->run_solver();
    return 0;
}
