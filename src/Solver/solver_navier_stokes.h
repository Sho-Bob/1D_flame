#ifndef SOLVER_NAVIER_STOKES_H
#define SOLVER_NAVIER_STOKES_H

#include "solver.h"
#include <vector>

class NavierStokesSolver : public Solver { // derived class for BurgerSolver

  public:

    void initialize() override;
    void apply_bc() override;
    void pre_rhs() override;

};

#endif
