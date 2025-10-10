#include "solver_linear_advection.h"
// #include "limiter.h"
#include "IO/vtk_writer.h"
#include "Numerics/reconstruction.h"
#include "Numerics/riemann_solver.h"
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
  Solver::initialize();

  this->num_vars = 1; // phi
  this->initialize_arrays(); 


  // set pointers
  this->phi = this->conservatives[0]; // alias for convenience
  this->phiL = new double[this->n_tot[0] + 1];
  this->phiR = new double[this->n_tot[0] + 1];
  this->U = new double[this->n_tot[0]];
  this->UL = new double[this->n_tot[0] + 1];
  this->UR = new double[this->n_tot[0] + 1];

  this->initialize_profile();

  this->cfl = this->dt * 1.0 / this->dx; // for linear advection, dt = cfl * dx / U_max, here U_max = 1.0
  std::cout << "CFL = " << this->cfl << std::endl;

}

void LinearAdvectionSolver::initialize_profile(){
  std::string init_case = this->input->getStringParam("init_case");
  if (init_case == "sharp") {
    FOR_ICV_G(0) {
      double x = this->coordinates[0][icv];
      if (x < 0.25 || x > 0.75)
        this->phi[icv] = 1.0;
      else {
        this->phi[icv] = 2.0;
      }

      // velocity field
      this->U[icv] = 1.0;

    }
  } else if (init_case == "smooth") {
    FOR_ICV_G(0) {
      double x = this->coordinates[0][icv];
      this->phi[icv] = 0.5 * (1.0 + sin(2.0 * M_PI * (x - 0.5)));
      // velocity field
      this->U[icv] = 1.0;

    }
  } else if (init_case == "smooth-WENO") { // inspired by Shu's paper on WENO (Efficient implementation of WENO schemes, 1996)
    FOR_ICV_G(0) {
      double x = this->coordinates[0][icv];

      double x_iPlusHalf = x + 0.5 * this->dx;
      double x_iMinusHalf = x - 0.5 * this->dx;

      // phi(x,t=0)=sin(pi*x)
      // cellaveraged of phi
      this->phi[icv] = - 1. / this->dx / M_PI * (cos(M_PI * x_iPlusHalf) - cos(M_PI * x_iMinusHalf));

      // velocity field
      this->U[icv] = 1.0;

    }
  } else if (init_case == "smooth-MUSCL") { // inspired by Pitfall MUSCL paper
    FOR_ICV_G(0) {
      double x = this->coordinates[0][icv];

      double x_iPlusHalf = x + 0.5 * this->dx;
      double x_iMinusHalf = x - 0.5 * this->dx;

      // phi(x,t=0)=sin(pi*x)
      // cellaveraged of phi
      this->phi[icv] = - 1. / this->dx / M_PI * (cos(M_PI * x_iPlusHalf) - cos(M_PI * x_iMinusHalf));

      // velocity field
      this->U[icv] = 1.0;

    }
  } else {
    throw std::runtime_error("Unknown init_case: " + init_case);
  }
}

void LinearAdvectionSolver::apply_bc() {
  // periodic boundary conditions
  for (int i=0; i< this->num_boundary_points; i++){
    this->phi[i] = this->phi[this->n[0] + i];
    this->phi[this->num_boundary_points + this->n[0] + i] = this->phi[this->num_boundary_points + i];
    this->U[i] = this->U[this->n[0] + i];
    this->U[this->num_boundary_points + this->n[0] + i] = this->U[this->num_boundary_points + i];
  }
}

void LinearAdvectionSolver::pre_rhs(){
  // base class pre_rhs
  Solver::pre_rhs();

  this->apply_bc();
  this->reconstruction->reconstruct(this->phi, this->phiL, this->phiR, this->reconstruction_method);
  this->reconstruction->reconstruct(this->U, this->UL, this->UR, this->reconstruction_method);
}

void LinearAdvectionSolver::rhs(){
  //Compute fluxes
#pragma omp parallel for
  FOR_IFA(0) {
    int icv0 = this->mesh->get_icv0(ifa);
    int icv1 = this->mesh->get_icv1(ifa);
    // if (icv0 == this->num_boundary_points || icv1 == this->num_boundary_points) {
    //   printf("Hello world\n");
    // }


    // local lax-Friedrichs flux
    double fL = this->UL[ifa] * this->phiL[ifa];
    double fR = this->UR[ifa] * this->phiR[ifa];
    double alpha = std::max(std::abs(this->UL[ifa]), std::abs(this->UR[ifa]));

    // double flux = 0.5 * (fL + fR) - 0.5 * alpha * (this->phiR[ifa] - this->phiL[ifa]);
    // upwind flux
    double flux = (this->UL[ifa] > 0.0) ? fL : fR;


    // update rhs
    this->rhs_conservatives[0][icv0] -= flux;
    this->rhs_conservatives[0][icv1] += flux;
  } // end for
} // end rhs

void LinearAdvectionSolver::output() {
  if (this->step == 0) {
    // create file
    IO::create_hdf5_file(this->output_file);

    // write mesh
    IO::write_structured_mesh_timestep(this->output_file,
        this->coordinates[0] + this->num_boundary_points,
        this->n[0], 1, 1, "x");
  }

  this->write_cv_data("phi", this->phi);
} // end output
