#include "solver_linear_advection.h"
// #include "limiter.h"
#include "vtk_writer.h"
#include "reconstruction.h"
#include <vector>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <cassert>


// ==============================================================
/*
 * Initializion of the 1D Burgers equation with a sine wave
 *
 */
// ==============================================================
void LinearAdvectionSolver::initialize(){
  std::cout << "LinearAdvectionSolver::initialize()" << std::endl;
  this->num_vars = 1; // phi
  this->conservatives = new double*[this->num_vars];
  this->conservatives_old = new double*[this->num_vars];
  this->rhs_conservatives = new double*[this->num_vars];

  this->initialize_grid();

  for (int v = 0; v < this->num_vars; ++v) {
    this->conservatives[v] = new double[this->n_tot[0]];
    this->conservatives_old[v] = new double[this->n_tot[0]];
    this->rhs_conservatives[v] = new double[this->n_tot[0]];
  }

  this->phi = this->conservatives[0]; // alias for convenience
  this->phiL = new double[this->n_tot[0] + 1];
  this->phiR = new double[this->n_tot[0] + 1];
  this->U = new double[this->n_tot[0]];
  this->UL = new double[this->n_tot[0] + 1];
  this->UR = new double[this->n_tot[0] + 1];

  FOR_ICV_G(0) {
    double x = this->coordinates[0][icv];
    if (x < 0.25 || x > 0.75)
      this->phi[icv] = 0.0;
    else {
      this->phi[icv] = 1.0;
    }

    // velocity field
    this->U[icv] = 1.0;

  }
  this->dt = 0.10 * this->Delta[0] / 1.0;
  double cfl = this->U[0] * this->dt / this->Delta[0];
  printf("dt = %g\n", this->dt);
  printf("CFL = %f\n", cfl);
}

void LinearAdvectionSolver::apply_bc(){
  // periodic boundary conditions
  for (int i=0; i< this->num_boundary_points; i++){
    this->phi[i] = this->phi[this->n[0] + i];
    this->phi[this->num_boundary_points + this->n[0] + i] = this->phi[this->num_boundary_points + i];
  }
}

void LinearAdvectionSolver::pre_rhs(){
  // base class pre_rhs
  Solver::pre_rhs();

  this->apply_bc();
  this->reconstruction->first_order(this->phi, this->phiL, this->phiR);
  this->reconstruction->first_order(this->U, this->UL, this->UR);
}

void LinearAdvectionSolver::rhs(){
  //Compute fluxes
#pragma omp parallel for
  FOR_IFA(0) {
    int icv0 = this->mesh->get_icv0(ifa);
    int icv1 = this->mesh->get_icv1(ifa);

    // sanity check for U and phi
    if (this->phiL[ifa] > 1.0 || this->phiL[ifa] < 0.0 ||
        this->phiR[ifa] > 1.0 || this->phiR[ifa] < 0.0) {
      std::cerr << "phi out of bounds at ifa = " << ifa << std::endl;
      assert(false);
    }
    // local lax-Friedrichs flux
    double fL = this->UL[ifa] * this->phiL[ifa];
    double fR = this->UR[ifa] * this->phiR[ifa];
    double alpha = std::max(std::abs(this->UL[ifa]), std::abs(this->UR[ifa]));

    double flux = 0.5 * (fL + fR) - 0.5 * alpha * (this->phiR[ifa] - this->phiL[ifa]);

    // update rhs
    this->rhs_conservatives[0][icv0] -= flux;
    this->rhs_conservatives[0][icv1] += flux;
  } // end for
} // end rhs

void LinearAdvectionSolver::output() {
  if (this->step % 10 != 0) return;
  if (this->step == 0) {
    IO::create_hdf5_file("LinearAdvectionSolver.h5");
  }

  std::vector<double> phi_vec(this->n[0]);
  FOR_ICV(0) phi_vec[icv] = this->phi[this->num_boundary_points + icv];

  IO::write_structured_mesh_timestep("LinearAdvectionSolver.h5",
      this->step,
      phi_vec,
      this->n[0], this->n[1], this->n[2]);
} // end output
