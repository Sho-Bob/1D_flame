#include "solver_burgers.h"
// #include "limiter.h"
#include "vtk_writer.h"
#include "reconstruction.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <omp.h>
#include <iostream>
#include <cassert>

// ==============================================================
/*
 * Initializion of the 1D Burgers equation with a sine wave
 *
 */
// ==============================================================
void BurgerSolver::initialize(){
  std::cout << "BurgerSolver::initialize()" << std::endl;
  this->num_vars = 1;
  this->conservatives = new double*[this->num_vars];
  this->rhs_conservatives = new double*[this->num_vars];
  for (int v = 0; v < this->num_vars; ++v) {
    this->conservatives[v] = new double[this->nx];
    this->rhs_conservatives[v] = new double[this->nx];
  }

  this->u = this->conservatives[0]; // alias for convenience
  this->uL = new double[this->nx + 1];
  this->uR = new double[this->nx + 1];
  this->initialize_grid();
  this->initialize_profile();
  this->dt = 0.1 * this->dx; // initial time step size
}

void BurgerSolver::initialize_profile() {
#pragma omp parallel for
  for (int idx = 0; idx < this->nx; ++idx) {
    this->x[idx] = this->x0 + idx * this->dx;
    this->u[idx] = std::sin(2.0 * M_PI * this->x[idx]);
  }
//   Apply_BC(u,N,ibd);
} // end initialize_profile

void BurgerSolver::apply_bc(){
  /// Neumann boundary conditions
// #pragma omp parallel for
//   for (int i=0; i<ibd; i++){
//     u[i] = u[ibd];
//     u[N+ibd+i] = u[N+ibd-1];
//   }
}

void BurgerSolver::reconstruct_WENO5(const double* u, double* uL, double* uR) {
  // first order
  const int nfa = this->nx + 1; // number of faces
  for (int ifa = 0 ; ifa < nfa; ifa++) {
    int i = ifa;
    if (i == 0) {
      uL[i] = u[this->nx-1];
    } else {
      uL[i] = u[i - 1];
    }

    if (i == this->nx) {
      uR[i] = u[0];
    } else {
      uR[i] = u[i];
    }
  } // end of loop over icv
}

void BurgerSolver::pre_rhs(){
  this->apply_bc();
  reconstruct_WENO5(this->u, this->uL, this->uR);
}

void BurgerSolver::rhs(){
  //Compute fluxes
#pragma omp parallel for
  for (int i = 0; i < this->nx + 1; i++){
    double qm = 0.5 * (this->uL[i] + std::abs(this->uL[i]));
    double qp = 0.5 * (this->uR[i] - std::abs(this->uR[i]));
    double flux = std::max(0.5 * qm*qm, 0.5 * qp*qp);
    this->rhs_conservatives[0][i] -= flux;
    this->rhs_conservatives[0][i+1] += flux;
  } // end for
} // end rhs

void BurgerSolver::output() {
  std::vector<double> u_vec(this->nx);
  for (int i = 0; i < this->nx; ++i) u_vec[i] = this->u[i];
  write_vtk(u_vec, this->dx,"new" + std::to_string(this->t / this->dt) + ".vtk");
} // end output
