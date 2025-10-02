#ifndef SOLVER_BURGER_H
#define SOLVER_BURGER_H

#include "solver.h"
#include <vector>

class BurgerSolver : public Solver { // derived class for BurgerSolver

  public:

    void initialize() override;
    void apply_bc() override;
    void pre_rhs() override;
    void rhs() override;

    void initialize_profile() override;
    void reconstruct_WENO5(const double* u, double* uL, double* uR);

    void output() override;
  private:
    double* u;
    double* uL;
    double* uR;
};

#endif
