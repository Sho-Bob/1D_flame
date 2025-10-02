#ifndef SOLVER_LINEAR_ADVECTION_H
#define SOLVER_LINEAR_ADVECTION_H

#include "solver.h"
#include <vector>

class LinearAdvectionSolver : public Solver { // derived class for LinearAdvectionSolver

  // ===============================================================
  // solving the 1D linear advection equation: 
  // \partial_t phi + U \partial_x phi = 0
  // ===============================================================
  public:

    void initialize() override;
    void apply_bc() override;
    void pre_rhs() override;
    void rhs() override;
    void output() override;

  private:

    double* phi;
    double* phiL;
    double* phiR;
    double* U;
    double* UL;
    double* UR;
};

#endif
