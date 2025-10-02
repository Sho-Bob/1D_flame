#ifndef SOLVER_NAVIER_STOKES_H
#define SOLVER_NAVIER_STOKES_H

#include "solver.h"
#include <vector>

class NavierStokesSolver : public Solver { // derived class for BurgerSolver

  public:

    void initialize() override;
    void apply_bc() override;
    void pre_rhs() override;
    void rhs() override;
    void output() override;


  private:
    // conservatives
    double* rho; // density
    double* rhoE; // total energy
    double** rhoU; // momentum
    // primitive variables
    double** U; // velocity
    double* p; // pressure
    double* T; // temperature
    // additional variables
    double* sos; // speed of sound
};

#endif
