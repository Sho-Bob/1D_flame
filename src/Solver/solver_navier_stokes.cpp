#include "solver_navier_stokes.h"
// #include "limiter.h"
#include "IO/vtk_writer.h"
#include "IO/hdf5_writer.h"
#include "IO/input.h"
#include "Numerics/reconstruction.h"
#include "Numerics/riemann_solver.h"
#include "Common/common.h"

#include <vector>
#include <string>
#include <omp.h>
#include <iostream>

// ==============================================================
// \brief: initialize multispecies solver
// ==============================================================
void NavierStokesSolver::initialize(){
  std::cout << "NavierStokesSolver::initialize()" << std::endl;
  Solver::initialize();

  // Physics
  this->Mixture = new Physics();
  this->Mixture->set_input_file(this->input);
  this->Mixture->initialize();

  this->species_names = this->input->getStringArrayParam("species");
  int num_species = this->species_names.size();

  this->num_vars = 2 + num_species; // rho, rhoE, rhoU, rhoY(k=1,...,nspecies-1)

  this->conservatives = new double*[this->num_vars];
  this->conservatives_old = new double*[this->num_vars];
  this->rhs_conservatives = new double*[this->num_vars];

  this->num_dim = 1; // 1D solver for now
  this->initialize_grid();

  FOR_VAR {
    this->conservatives[v] = new double[this->n_tot[0]];
    this->conservatives_old[v] = new double[this->n_tot[0]];
    this->rhs_conservatives[v] = new double[this->n_tot[0]];
  }

  this->rho = this->conservatives[0]; // density
  this->rhoE = this->conservatives[1]; // total energy
                                       //
  this->T = new double[this->n_tot[0]]; // temperature
  this->p = new double[this->n_tot[0]]; // pressure
                                        //
  this->rhoU = new double*[this->num_dim]; // momentum
  this->U = new double*[this->num_dim]; // velocity
                                        //
  FOR_IDIM {
    this->rhoU[iDim] = new double[this->n_tot[0]]; // for now
    // this->rhoU[iDim] = this->conservatives[2 + iDim]; // later
    this->U[iDim] = new double[this->n_tot[0]];
  }

  this->initialize_profile();


}

// ==============================================================
// \brief: initialize flow profile
// ==============================================================
void NavierStokesSolver::initialize_profile(){
  std::string init_case = this->input->getStringParam("init_case");
  if (init_case == "contact-discontinuity") {
    std::cout << "Initializing contact discontinuity" << std::endl;
    std::vector<double> rho_vec(this->n[0]);
    FOR_ICV_G(0) {
      double x = this->coordinates[0][icv];

      this->T[icv] = (x < 0.25 || x > 0.75) ? 300.0 : 600.0; // temperature
      this->p[icv] = 1e5; // pressure

      const double Y [2] = {0.0, 1.0}; // mass fractions
      this->Mixture->SetMixture_TPY(this->T[icv], this->p[icv], Y);
      this->rho[icv] = this->Mixture->GetRho();

      if (icv >= this->num_boundary_points && icv < this->num_boundary_points + this->n[0]) {
        rho_vec[icv - this->num_boundary_points] = this->rho[icv];
      }
    }
    IO::create_hdf5_file("navier-stokes.h5");
    IO::write_structured_mesh_timestep("navier-stokes.h5",
        this->step,
        rho_vec,
        this->n[0], 1, 1);
  } else if (init_case == "smooth") {
    std::cout << "Initializing smooth" << std::endl;
    // parameters
    double rho_max = this->input->getDoubleParam("rho_max");
    double rho_min = this->input->getDoubleParam("rho_min");

    std::vector<double> rho_vec(this->n[0]);

    FOR_ICV_G(0) {
      double x = this->coordinates[0][icv];

      this->rho[icv] = 0.5 * (rho_max + rho_min) + 0.5 * (rho_max - rho_min) * sin(2.0 * M_PI * x);
      this->p[icv] = 50e5; // pressure

      const double Y [2] = {0.0, 1.0}; // mass fractions
      this->Mixture->SetMixture_PRY(this->p[icv], this->rho[icv], Y);
      this->T[icv] = this->Mixture->GetT();

      if (icv >= this->num_boundary_points && icv < this->num_boundary_points + this->n[0]) {
        rho_vec[icv - this->num_boundary_points] = this->rho[icv];
      }
    }
    IO::create_hdf5_file("navier-stokes.h5");
    IO::write_structured_mesh_timestep("navier-stokes.h5",
        this->step,
        rho_vec,
        this->n[0], 1, 1);
  } else {
    throw std::runtime_error("Unknown init_case: " + init_case);
  }
} // end initialize_profile

// ==============================================================
// \brief: enforce boundary conditions weakly by filling ghost cells
// ==============================================================
void NavierStokesSolver::apply_bc(){
}

// ==============================================================
// \brief: gradient computation and face reconstruction
// ==============================================================
void NavierStokesSolver::pre_rhs(){
}

// ==============================================================
// \brief: loop over all internal faces to compute the rhs
// ==============================================================
void NavierStokesSolver::rhs(){

} // end rhs

void NavierStokesSolver::output(){
} // end output

