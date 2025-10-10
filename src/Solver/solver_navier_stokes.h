#ifndef SOLVER_NAVIER_STOKES_H
#define SOLVER_NAVIER_STOKES_H

#include "solver.h"
#include <vector>

class NavierStokesSolver : public Solver { // derived class for BurgerSolver

  public:

    void initialize() override;
    void initialize_profile() override;
    void apply_bc() override;
    void pre_rhs() override;
    void rhs() override;
    void output() override;

    // update primitives from conservatives
    void update_primitives() override;

    // update mark_double_flux
    void update_double_flux();

  private:
    // conservatives
    double* rho; // density
    double* rhoE; // total energy
    double** rhoU; // momentum
    double** rhoY; // species density
    // primitive variables
    double** U; // velocity
    double* p; // pressure
    double* T; // temperature
    double** Y; // mass fractions
    // additional variables
    double* sos; // speed of sound

    // physics
    std::vector<std::string> species_names;
    int num_species;

    // numerics for reconstruction
    std::string recon_vars;
    double *rhoL, *rhoR;
    double *rhoEL, *rhoER;
    double *sosL, *sosR;
    double *TL, *TR;
    double *pL, *pR;
    double **UL, **UR;
    double **YL, **YR;

    // double flux method
    int* mark_double_flux;
    double* gammaS;
    double* e0S;
};

#endif
