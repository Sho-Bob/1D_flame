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
  this->num_species = this->species_names.size();

  this->num_vars = 2 + this->num_dim +  + this->num_species; // rho, rhoE, rhoU, rhoY(k=1,...,nspecies-1)

  this->initialize_arrays();

  // set pointers
  this->rho = this->conservatives[0]; // density
  this->rhoL = new double[this->n_tot[0] + 1];
  this->rhoR = new double[this->n_tot[0] + 1];

  this->rhoE = this->conservatives[1]; // total energy
  this->rhoEL = new double[this->n_tot[0] + 1];
  this->rhoER = new double[this->n_tot[0] + 1];
                                       //
  this->T = new double[this->n_tot[0]]; // temperature
  this->TL = new double[this->n_tot[0] + 1];
  this->TR = new double[this->n_tot[0] + 1];

  this->p = new double[this->n_tot[0]]; // pressure
  this->pL = new double[this->n_tot[0] + 1];
  this->pR = new double[this->n_tot[0] + 1];

  this->sos = new double[this->n_tot[0]]; // speed of sound
  this->sosL = new double[this->n_tot[0] + 1];
  this->sosR = new double[this->n_tot[0] + 1];

  this->rhoU = new double*[this->num_dim]; // momentum
  this->U = new double*[this->num_dim]; // velocity
  this->UL = new double*[this->num_dim]; // velocity
  this->UR = new double*[this->num_dim]; // velocity

  FOR_IDIM {
    this->rhoU[iDim] = this->conservatives[2 + iDim]; // momentum
    this->U[iDim] = new double[this->n_tot[0]];
    this->UL[iDim] = new double[this->n_tot[0] + 1];
    this->UR[iDim] = new double[this->n_tot[0] + 1];
  }

  this->Y = new double*[num_species]; // mass fractions
  this->YL = new double*[num_species]; // mass fractions
  this->YR = new double*[num_species]; // mass fractions

  this->rhoY = new double*[num_species]; // species density
  LOOP_k_N(num_species) {
    this->rhoY[k] = this->conservatives[2 + this->num_dim + k]; // species density
    this->Y[k] = new double[this->n_tot[0]];
    this->YL[k] = new double[this->n_tot[0] + 1];
    this->YR[k] = new double[this->n_tot[0] + 1];
  }


  // reconstruction variables
  this->recon_vars = this->input->getStringParam("recon_vars", "PRY");
  if (this->recon_vars != "PRY" && this->recon_vars != "TPY") {
    throw std::runtime_error("Unknown recon_vars: " + this->recon_vars);
  }

  // double flux method
  this->mark_double_flux = new int[this->n_tot[0]]; // cell values
  this->gammaS = new double[this->n_tot[0]]; // cell values
  this->e0S = new double[this->n_tot[0]]; // cell values

  this->initialize_profile();

} // end initialize

// ==============================================================
// \brief: initialize flow profile
// ==============================================================
void NavierStokesSolver::initialize_profile() {
  std::string init_case = this->input->getStringParam("init_case");
  if (init_case == "contact-discontinuity") {
    std::cout << "Initializing contact discontinuity" << std::endl;
    std::vector<double> rho_vec(this->n[0]);
    FOR_ICV_G(0) {
      double x = this->coordinates[0][icv];

      this->T[icv] = (x < 0.25 || x > 0.75) ? 300.0 : 600.0; // temperature
      this->p[icv] = 1e5; // pressure
      this->Y[0][icv] = 0.0;
      this->Y[1][icv] = 1.0;

      double Y_ [2] = {this->Y[0][icv], this->Y[1][icv]};

      this->Mixture->SetMixture_TPY(this->T[icv], this->p[icv], Y_);
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
    double rho_max = this->input->getDoubleParam("rho_max", 750.0);
    double rho_min = this->input->getDoubleParam("rho_min", 50.0);
    double u_init = this->input->getDoubleParam("u_init", 100.0);
    double p_init = this->input->getDoubleParam("p_init", 50e5);

    std::vector<double> rho_vec(this->n[0]);

    int iDim = 0; // 1D solver for now

    FOR_ICV_G(iDim) {
      double x = this->coordinates[0][icv];
      double x_right = x + 0.5 * this->Delta[0];
      double x_left = x - 0.5 * this->Delta[0];
      // cellaveraged of rho
      double f_right = (rho_max + rho_min) * x_right - 1.0 / M_PI * (rho_max - rho_min) * cos(M_PI * x_right);
      double f_left = (rho_max + rho_min) * x_left - 1.0 / M_PI * (rho_max - rho_min) * cos(M_PI * x_left);
      this->rho[icv] = 0.5 / this->Delta[0] * (f_right - f_left); // density
                                                                  //
      this->p[icv] = p_init;
      this->U[iDim][icv] = u_init;

      const double Y [2] = {0.0, 1.0}; // mass fractions
      this->Mixture->SetMixture_PRY(this->p[icv], this->rho[icv], Y);
      this->T[icv] = this->Mixture->GetT();
      this->sos[icv] = this->Mixture->GetSos();

      // conservatives
      this->rhoU[iDim][icv] = this->rho[icv] * this->U[iDim][icv];

      double e = this->Mixture->GetE(); // internal energy
      double ekin = 0.5 * this->U[iDim][icv] * this->U[iDim][icv]; // kinetic energy
      this->rhoE[icv] = this->rho[icv] * (e + ekin);

      if (icv >= this->num_boundary_points && icv < this->num_boundary_points + this->n[0]) {
        rho_vec[icv - this->num_boundary_points] = this->rho[icv];
      }
    }
  } else {
    throw std::runtime_error("Unknown init_case: " + init_case);
  }
  // update double flux stufff
  this->update_double_flux();

  std::cout << "Initialization done." << std::endl;
} // end initialize_profile

// ==============================================================
// \brief: enforce boundary conditions weakly by filling ghost cells
// ==============================================================
void NavierStokesSolver::apply_bc(){

  // set the primitives in the ghost cells
  // for now assume periodic

  for (int i=0; i< this->num_boundary_points; i++){
    this->rho[i] = this->rho[this->n[0] + i];
    this->rho[this->num_boundary_points + this->n[0] + i] = this->rho[this->num_boundary_points + i];

    this->T[i] = this->T[this->n[0] + i];
    this->T[this->num_boundary_points + this->n[0] + i] = this->T[this->num_boundary_points + i];

    this->p[i] = this->p[this->n[0] + i];
    this->p[this->num_boundary_points + this->n[0] + i] = this->p[this->num_boundary_points + i];

    this->U[0][i] = this->U[0][this->n[0] + i];
    this->U[0][this->num_boundary_points + this->n[0] + i] = this->U[0][this->num_boundary_points + i];

    LOOP_k_N (this->num_species) {
      this->Y[k][i] = this->Y[k][this->n[0] + i];
      this->Y[k][this->num_boundary_points + this->n[0] + i] = this->Y[k][this->num_boundary_points + i];
    }
  }

}

// ==============================================================
// \brief: gradient computation and face reconstruction
// ==============================================================
void NavierStokesSolver::pre_rhs() {
  // base class pre_rhs
  Solver::pre_rhs(); // zero out the rhs

  this->apply_bc();
  // velocity
  FOR_IDIM this->reconstruction->reconstruct(this->U[iDim], this->UL[iDim], this->UR[iDim], this->reconstruction_method);

  // mass fractions
  LOOP_k_N(this->num_species) 
    this->reconstruction->reconstruct(this->Y[k], this->YL[k], this->YR[k], this->reconstruction_method);

  if (this->recon_vars == "PRY") {
    this->reconstruction->reconstruct(this->rho, this->rhoL, this->rhoR, this->reconstruction_method);
    this->reconstruction->reconstruct(this->p, this->pL, this->pR, this->reconstruction_method);
    // reconstruct speed of sound
  } else if (this->recon_vars == "TPY") {
    this->reconstruction->reconstruct(this->T, this->TL, this->TR, this->reconstruction_method);
    this->reconstruction->reconstruct(this->p, this->pL, this->pR, this->reconstruction_method);
  }

}

// ==============================================================
// \brief: loop over all internal faces to compute the rhs
// ==============================================================
void NavierStokesSolver::rhs(){
  double Frho, FrhoE;
  double Frhou[3]; // for now 3D
  double Frho_scalars[this->num_species];

  double pL_fa, pR_fa;
  double rhoL_fa, rhoR_fa;
  double rhoE_L_fa, rhoE_R_fa;
  double sosL_fa, sosR_fa;

  double fluxes [5 + this->num_species];

  double normal_vec[3] = {1.0, 0.0, 0.0}; // for now 1D
  //Compute fluxes
  FOR_IFA(0) {
    int icv0 = this->mesh->get_icv0(ifa);
    int icv1 = this->mesh->get_icv1(ifa);

    // these are always the same
    pL_fa = this->pL[ifa];
    pR_fa = this->pR[ifa];
    double Y_L[2] = {this->YL[0][ifa], this->YL[1][ifa]};
    double Y_R[2] = {this->YR[0][ifa], this->YR[1][ifa]};

    // determine left and right states at the face, that depends whether
    // double flux is used or not
    if (this->mark_double_flux[icv0] == 0 && this->mark_double_flux[icv1] == 0) { // Fully-conservative
      // compute face values based on reconstructed values
      if (this->recon_vars == "PRY") {
        rhoL_fa = this->rhoL[ifa];
        rhoR_fa = this->rhoR[ifa];
        // set mixture in L
        this->Mixture->SetMixture_PRY(pL_fa, rhoL_fa, Y_L);
        double eL = this->Mixture->GetE();
        double ekinL = 0.5 * this->UL[0][ifa] * this->UL[0][ifa];
        rhoE_L_fa = rhoL_fa * (eL + ekinL);
        sosL_fa = this->Mixture->GetSos();
        // set mixture in R
        this->Mixture->SetMixture_PRY(pR_fa, rhoR_fa, Y_R);
        double eR = this->Mixture->GetE();
        double ekinR = 0.5 * this->UR[0][ifa] * this->UR[0][ifa];
        rhoE_R_fa = rhoR_fa * (eR + ekinR);
        sosR_fa = this->Mixture->GetSos();
      } else if (this->recon_vars == "TPY") {
        this->Mixture->SetMixture_TPY(this->TL[ifa], pL_fa, Y_L);
        rhoL_fa = this->Mixture->GetRho();
        double eL = this->Mixture->GetE();
        double ekinL = 0.5 * this->UL[0][ifa] * this->UL[0][ifa];
        rhoE_L_fa = rhoL_fa * (eL + ekinL);
        sosL_fa = this->Mixture->GetSos();
        this->Mixture->SetMixture_TPY(this->TR[ifa], pR_fa, Y_R);
        rhoR_fa = this->Mixture->GetRho();
        double eR = this->Mixture->GetE();
        double ekinR = 0.5 * this->UR[0][ifa] * this->UR[0][ifa];
        rhoE_R_fa = rhoR_fa * (eR + ekinR);
        sosR_fa = this->Mixture->GetSos();
      }

      // riemann solver
      riemann_solver::HLLC(
          Frho, Frhou, FrhoE, Frho_scalars,
          pL_fa, pR_fa,
          rhoL_fa, rhoR_fa,
          rhoE_L_fa, rhoE_R_fa,
          sosL_fa, sosR_fa,
          &this->YL[0][ifa], &this->YR[0][ifa],
          &this->UL[0][ifa], &this->UR[0][ifa],
          1.0, // lambda
          normal_vec,
          this->num_species
          );

      fluxes[0] = Frho;
      fluxes[1] = FrhoE;
      FOR_IDIM fluxes[2 + iDim] = Frhou[iDim];
      LOOP_k_N(this->num_species) fluxes[2 + this->num_dim + k] = Frho_scalars[k];

      // update rhs
      FOR_VAR {
        this->rhs_conservatives[v][icv0] -= fluxes[v];
        this->rhs_conservatives[v][icv1] += fluxes[v];
      }
    } else { // double flux
      double m_gammaS, m_e0S;
      double Frho0, Frhou0[3], FrhoE0, Frho_scalars0[this->num_species];
      double Frho1, Frhou1[3], FrhoE1, Frho_scalars1[this->num_species];

      // get rho, depending on what was reconstructed
      if (this->recon_vars == "PRY") {
        rhoL_fa = this->rhoL[ifa];
        rhoR_fa = this->rhoR[ifa];
      } else if (this->recon_vars == "TPY") {
        this->Mixture->SetMixture_TPY(this->TL[ifa], pL_fa, Y_L);
        rhoL_fa = this->Mixture->GetRho();
        this->Mixture->SetMixture_TPY(this->TR[ifa], pR_fa, Y_R);
        rhoR_fa = this->Mixture->GetRho();
      }
      // side 0
      m_gammaS = this->gammaS[icv0]; m_e0S = this->e0S[icv0];
      // get sos and rhoE based on reconstructed p, rho, U
      sosL_fa = std::sqrt(m_gammaS * pL_fa / rhoL_fa);
      rhoE_L_fa = pL_fa / (m_gammaS - 1.0) + 0.5 * rhoL_fa * (this->UL[0][ifa] * this->UL[0][ifa] + m_e0S);
      sosR_fa = std::sqrt(m_gammaS * pR_fa / rhoR_fa);
      rhoE_R_fa = pR_fa / (m_gammaS - 1.0) + 0.5 * rhoR_fa * (this->UR[0][ifa] * this->UR[0][ifa] + m_e0S);
      // riemann solver
      riemann_solver::HLLC(
          Frho0, Frhou0, FrhoE0, Frho_scalars0,
          pL_fa, pR_fa,
          rhoL_fa, rhoR_fa,
          rhoE_L_fa, rhoE_R_fa,
          sosL_fa, sosR_fa,
          &this->YL[0][ifa], &this->YR[0][ifa],
          &this->UL[0][ifa], &this->UR[0][ifa],
          1.0, // lambda
          normal_vec,
          this->num_species
          );
      double fluxes0[2 + this->num_dim + this->num_species];
      fluxes0[0] = Frho0;
      fluxes0[1] = FrhoE0;
      FOR_IDIM fluxes0[2 + iDim] = Frhou0[iDim];
      LOOP_k_N(this->num_species) fluxes0[2 + this->num_dim + k] = Frho_scalars0[k];

      // side 1
      m_gammaS = this->gammaS[icv1]; m_e0S = this->e0S[icv1];
      // get sos and rhoE based on reconstructed p, rho, U
      sosL_fa = std::sqrt(m_gammaS * pL_fa / rhoL_fa);
      rhoE_L_fa = pL_fa / (m_gammaS - 1.0) + 0.5 * rhoL_fa * (this->UL[0][ifa] * this->UL[0][ifa] + m_e0S);
      sosR_fa = std::sqrt(m_gammaS * pR_fa / rhoR_fa);
      rhoE_R_fa = pR_fa / (m_gammaS - 1.0) + 0.5 * rhoR_fa * (this->UR[0][ifa] * this->UR[0][ifa] + m_e0S);
      // riemann solver
      riemann_solver::HLLC(
          Frho1, Frhou1, FrhoE1, Frho_scalars1,
          pL_fa, pR_fa,
          rhoL_fa, rhoR_fa,
          rhoE_L_fa, rhoE_R_fa,
          sosL_fa, sosR_fa,
          &this->YL[0][ifa], &this->YR[0][ifa],
          &this->UL[0][ifa], &this->UR[0][ifa],
          1.0, // lambda
          normal_vec,
          this->num_species
          );
      double fluxes1[2 + this->num_dim + this->num_species];

      fluxes1[0] = Frho1;
      fluxes1[1] = FrhoE1;
      FOR_IDIM fluxes1[2 + iDim] = Frhou1[iDim];
      LOOP_k_N(this->num_species) fluxes1[2 + this->num_dim + k] = Frho_scalars1[k];

      // add to rhs, 0 and 1 sides
      FOR_VAR {
        this->rhs_conservatives[v][icv0] -= fluxes0[v];
        this->rhs_conservatives[v][icv1] += fluxes1[v];
      }
      
    } // end of double flux if statement
  } // end for
} // end rhs

// ==============================================================
// \brief: update primitive variables from conservatives
// ==============================================================
void NavierStokesSolver::update_primitives() {
  FOR_ICV(0) { // loop over all icv's
    // velocity
    FOR_IDIM this->U[iDim][icv] = this->rhoU[iDim][icv] / this->rho[icv];
    double ukuk = 0.0;
    FOR_IDIM ukuk += this->U[iDim][icv] * this->U[iDim][icv];

    // mass fractions
    LOOP_k_N(this->num_species) this->Y[k][icv] = this->conservatives[2 + this->num_dim + k][icv] / this->rho[icv];
    // pressure; either fully-conservative or double-flux
    if (this->mark_double_flux[icv] == 0) { // fully-conservative
      double e = this->rhoE[icv] / this->rho[icv] - 0.5 * ukuk;
      this->Mixture->SetMixture_ERY(e, this->rho[icv], &this->Y[0][icv]);
      this->p[icv] = this->Mixture->GetP();

    } else { // double-flux
      this->p[icv] = (this->gammaS[icv] - 1) *
        (this->rhoE[icv] - this->rho[icv] * this->e0S[icv] -
         0.5 * this->rho[icv] * ukuk);
      this->Mixture->SetMixture_PRY(this->p[icv], this->rho[icv], &this->Y[0][icv]);
      double e = this->Mixture->GetE();
      // need to overwrite rhoE if end of time step
      if (this->rk_step == this->order_time - 1) {
        this->rhoE[icv] = this->rho[icv] * (e + 0.5 * ukuk);
      }

    }

    // temperature and etc
    this->T[icv] = this->Mixture->GetT();
    this->sos[icv] = this->Mixture->GetSos();

    // update mark_double_flux, gammaS and e0S if end of time step
    if (this->rk_step == this->order_time - 1) this->update_double_flux();
  } // end for icv
} // end of update_primitives
// ==============================================================

// ==============================================================
void NavierStokesSolver::update_double_flux(){
#pragma omp parallel for
  FOR_ICV_G(0) {
    this->mark_double_flux[icv] = 1; // for now all cells use double flux
                                     //
    double ukuk = 0.0;
    FOR_IDIM ukuk += this->U[iDim][icv] * this->U[iDim][icv];
    double e = (this->rhoE[icv] - 0.5 * this->rho[icv] *
                                      ukuk) / this->rho[icv];
    this->gammaS[icv] =
        this->sos[icv] * this->sos[icv] * this->rho[icv] / this->p[icv];
    this->e0S[icv] = e - this->p[icv] / this->rho[icv] / (this->gammaS[icv] - 1);
  }
} // end of update_double_flux
// =============================================================

void NavierStokesSolver::output(){
  if (this->step == 0) {
    // create file
    IO::create_hdf5_file(this->output_file);

    // write mesh
    IO::write_structured_mesh_timestep(this->output_file,
        this->coordinates[0] + this->num_boundary_points,
        this->n[0], 1, 1, "x");
  }

  // write primitives
  this->write_cv_data("rho", this->rho);
  this->write_cv_data("sos", this->sos);
  this->write_cv_data("p", this->p);
  this->write_cv_data("T", this->T);
  FOR_IDIM this->write_cv_data("U" + std::to_string(iDim), this->U[iDim]);
  LOOP_k_N(this->num_species) this->write_cv_data("Y" + this->species_names[k], this->Y[k]);

} // end output

