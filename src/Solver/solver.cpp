#include "solver.h"
// #include "limiter.h"

#include "IO/vtk_writer.h"
#include "Numerics/reconstruction.h"

#include <vector>
#include <cmath>
#include <algorithm>
#include <omp.h>
#include <iostream>
#include <cassert>

void Solver::run_solver(){
  // initialize solver; set up grid, initial condition
  this->initialize();

  // time-stepping loop
  for (int step = 0; step < 1000; ++step) {
    this->step = step;

    for (int v = 0; v < this->num_vars; ++v) {
      FOR_ICV_G(0) this->conservatives_old[v][icv] = this->conservatives[v][icv];
    }

    for (int idrk = 0; idrk < 3; ++idrk) {
      // \brief: reconstruction, apply BC
      this->pre_rhs();
      // \brief: calculate the rhs
      this->rhs();
      // \brief: update conservatives based on rk step
      this->update_conservatives(idrk);

    }
    // \brief: update conservatives
    this->output();

    // advance time
    this->t += this->dt;
    // printf("Step %d, Time = %f\n", step, this->t);
  }

}

void Solver::initialize(){
    // to be implemented in derived classes
  std::cout << "Base class initialize()" << std::endl;
}

void Solver::apply_bc(){
    // to be implemented in derived classes
  std::cout << "Base class apply_bc()" << std::endl;
}

void Solver::initialize_grid(){
  // TODO: implement TOML reading of these things
  // for now, just hardcode some values
  this->num_dim = 1; 
  this->num_boundary_points = 4; // ghost cells
  this->nx = 101;
  this->x0 = 0.0; this->x1 = 1.0;
  this->mesh = new Grid(this->nx, this->num_boundary_points, this->x0, this->x1);
  // Initialize and set pointers
  this->n = this->mesh->get_n_ptr();
  this->n_tot = this->mesh->get_n_tot_ptr();
  this->coordinates = this->mesh->get_coordinates();
  this->Delta = this->mesh->get_Delta_ptr();
  this->dx = this->Delta[0];

  // pass grid to reconstruction
  this->reconstruction = new Reconstruction(this->mesh);
} // end initialize_grid

void Solver::initialize_profile(){
} // end initialize_profile

// ==============================================================
// \brief: gradient computation and face reconstruction
// ==============================================================
void Solver::pre_rhs() {
  for (int v = 0; v < this->num_vars; ++v) {
    FOR_ICV_G(v) this->rhs_conservatives[v][icv] = 0.0;
  }
} // end pre_rhs


// ==============================================================
// \brief: loop over all cells to compute the rhs
// ==============================================================
void Solver::rhs() {
} // end rhs

// ==============================================================
// \brief: update the conservatives after calculating the rhs
// ==============================================================
void Solver::update_conservatives(const int idrk) {
  for (int v = 0; v < this->num_vars; ++v) {
    FOR_ICV(0) {
      int i = icv;
      if (idrk == 0) { // stage 1
        this->conservatives[v][i] += this->dt / this->dx * this->rhs_conservatives[v][i];
      } else if (idrk == 1) { // stage 2
        this->conservatives[v][i] = 0.75 * this->conservatives_old[v][i] + 0.25 * (this->conservatives[v][i] + this->dt / this->dx * this->rhs_conservatives[v][i]);
      } else if (idrk == 2) {
        this->conservatives[v][i] = (this->conservatives_old[v][i] + 2.0 * (this->conservatives[v][i] + this->dt / this->dx * this->rhs_conservatives[v][i])) / 3.0;
      }
      if (this->conservatives[v][i] != this->conservatives[v][i]) {
        std::cerr << "NaN detected in conservatives at cell " << i << " variable " << v << std::endl;
        std::cout << this->rhs_conservatives[v][i] << std::endl;
        assert( false);
      }
    } // end of loop over icv
  }
} // end update_conservatives

void burgers_initialize(std::vector<double>& u, int N, double x0, double x1, double dx, int ibd){
    #pragma omp parallel for
    for(int i = ibd; i < N+ibd; i++){
        double x = x0 + (i-ibd) * dx;
        u[i] = std::sin(2.0 * M_PI * x);
    }
    Apply_BC(u,N,ibd);
}

void Apply_BC(std::vector<double>& u, int N, int ibd){
    /// Neumann boundary conditions
    #pragma omp parallel for
    for (int i=0; i<ibd; i++){
        u[i] = u[ibd];
        u[N+ibd+i] = u[N+ibd-1];
    }
}

inline double minmod(double a, double b) {
    double kappa = 1.0/3.0;
    double kappa_f = (3.0-kappa)/(1.0+kappa);
    if (a*b*kappa_f > 0.0)
        return std::min(std::abs(a),std::abs(b))*std::copysign(1.0,a);
    else
        return 0.0;
}

void reconstruct_MUSCL_minmod(const std::vector<double>& u, std::vector<double>& uL, std::vector<double>& uR, int N, int ibd) {
    
        #pragma omp parallel for
        for (int i = 0; i <= N; ++i) {
            const double kappa = 1.0/3.0;
            int im2 = i + ibd - 2;
            int im1 = i + ibd - 1;
            int i0  = i + ibd;
            int ip1 = i + ibd + 1;
    
            double duL0 = minmod(u[i0] - u[im1], u[ip1] - u[i0]);
            double duL1 = minmod(u[ip1] - u[i0], u[ip1+1] - u[ip1]);
            uL[i] = u[i0] + 0.25 * ((1.0 - kappa) * duL0 + (1.0 + kappa) * duL1);
    
            double duR0 = minmod(u[ip1] - u[i0], u[ip1+1] - u[ip1]);
            double duR1 = minmod(u[ip1+1] - u[ip1], u[ip1] - u[i0]);
            uR[i] = u[ip1] - 0.25 * ((1.0 + kappa) * duR0 + (1.0 - kappa) * duR1);
        }
    
    
}

void reconstruct_WENO5(
    const std::vector<double>& u,
    std::vector<double>& uL,  // uL[i] = u_{i-1/2}^L (left state at interface i)
    std::vector<double>& uR,  // uR[i] = u_{i-1/2}^R (right state at interface i)
    int N,
    int ibd)
{
    assert(ibd >= 3 && "WENO5 requires at least 3 ghost cells");
    const double eps = 1e-6;
    const double d0 = 0.1, d1 = 0.6, d2 = 0.3;

    // Reconstruct left and right states at each interface i (i = 0 to N)
    #pragma omp parallel for
    for (int i = 0; i <= N; ++i) {
        // For interface i, we need stencil around cells i-1 and i
        
        // Left state uL[i] - extrapolated from cell i-1 using stencil {i-3, i-2, i-1, i, i+1}
        int im3 = (i-1) + ibd - 2;  // u[i-3]
        int im2 = (i-1) + ibd - 1;  // u[i-2] 
        int im1 = (i-1) + ibd;      // u[i-1]
        int i0  = (i-1) + ibd + 1;  // u[i]
        int ip1 = (i-1) + ibd + 2;  // u[i+1]

        // Smoothness indicators for left state
        double beta0 = (13.0/12.0)*std::pow(u[im3] - 2*u[im2] + u[im1], 2)
                     + (1.0/4.0)*std::pow(u[im3] - 4*u[im2] + 3*u[im1], 2);

        double beta1 = (13.0/12.0)*std::pow(u[im2] - 2*u[im1] + u[i0], 2)
                     + (1.0/4.0)*std::pow(u[im2] - u[i0], 2);

        double beta2 = (13.0/12.0)*std::pow(u[im1] - 2*u[i0] + u[ip1], 2)
                     + (1.0/4.0)*std::pow(3*u[im1] - 4*u[i0] + u[ip1], 2);

        // Nonlinear weights for left state
        double alpha0 = d0 / ((eps + beta0)*(eps + beta0));
        double alpha1 = d1 / ((eps + beta1)*(eps + beta1));
        double alpha2 = d2 / ((eps + beta2)*(eps + beta2));

        double sum_alpha = alpha0 + alpha1 + alpha2;

        double w0 = alpha0 / sum_alpha;
        double w1 = alpha1 / sum_alpha;
        double w2 = alpha2 / sum_alpha;

        // Candidate reconstructions for left state
        double q0 = (1.0/3.0)*u[im3] - (7.0/6.0)*u[im2] + (11.0/6.0)*u[im1];
        double q1 = (-1.0/6.0)*u[im2] + (5.0/6.0)*u[im1] + (1.0/3.0)*u[i0];
        double q2 = (1.0/3.0)*u[im1] + (5.0/6.0)*u[i0] - (1.0/6.0)*u[ip1];

        uL[i] = w0*q0 + w1*q1 + w2*q2;

        // Right state uR[i] - extrapolated from cell i using reversed stencil {i+2, i+1, i, i-1, i-2}
        int jm3 = i + ibd + 2;  // u[i+2]
        int jm2 = i + ibd + 1;  // u[i+1]
        int jm1 = i + ibd;      // u[i]
        int j0  = i + ibd - 1;  // u[i-1]
        int jp1 = i + ibd - 2;  // u[i-2]

        // Smoothness indicators for right state (reversed stencil)
        beta0 = (13.0/12.0)*std::pow(u[jm3] - 2*u[jm2] + u[jm1], 2)
              + (1.0/4.0)*std::pow(u[jm3] - 4*u[jm2] + 3*u[jm1], 2);

        beta1 = (13.0/12.0)*std::pow(u[jm2] - 2*u[jm1] + u[j0], 2)
              + (1.0/4.0)*std::pow(u[jm2] - u[j0], 2);

        beta2 = (13.0/12.0)*std::pow(u[jm1] - 2*u[j0] + u[jp1], 2)
              + (1.0/4.0)*std::pow(3*u[jm1] - 4*u[j0] + u[jp1], 2);

        // Nonlinear weights for right state
        alpha0 = d0 / ((eps + beta0)*(eps + beta0));
        alpha1 = d1 / ((eps + beta1)*(eps + beta1));
        alpha2 = d2 / ((eps + beta2)*(eps + beta2));

        sum_alpha = alpha0 + alpha1 + alpha2;

        w0 = alpha0 / sum_alpha;
        w1 = alpha1 / sum_alpha;
        w2 = alpha2 / sum_alpha;

        // Candidate reconstructions for right state
        q0 = (1.0/3.0)*u[jm3] - (7.0/6.0)*u[jm2] + (11.0/6.0)*u[jm1];
        q1 = (-1.0/6.0)*u[jm2] + (5.0/6.0)*u[jm1] + (1.0/3.0)*u[j0];
        q2 = (1.0/3.0)*u[jm1] + (5.0/6.0)*u[j0] - (1.0/6.0)*u[jp1];

        uR[i] = w0*q0 + w1*q1 + w2*q2;
    }
}


void simulate_burgers1d(std::vector<double>& u, double dx, double CFL, double T_end, int ibd){
    const int N = u.size()-2*ibd;
    int n_steps = 0;
    std::vector<double> u_new(N,0.0);
    std::vector<double> u_flux(N+1,0.0);
    std::vector<double> uL(N+1,0.0), uR(N+1,0.0);

    double t = 0.0;

    while(t < T_end){

        n_steps++;
        //Compute dt from CFL condition
        double max_speed = 0.0;
        #pragma omp parallel for reduction(max:max_speed)
        for (int i = ibd;i < N+ibd;i++){
            max_speed = std::max(max_speed,std::abs(u[i]));
        }
        double dt = CFL * dx / max_speed;

        // Reconstruct left and right states (1st order)
        // #pragma omp parallel for
        // for (int i=0;i < N+1;i++){
        // // MUSCL for future work
        //     uL[i] = u[i+ibd-1];
        //     uR[i] = u[i+ibd];
        // }

        // reconstruct_MUSCL_minmod(u, uL, uR, N, ibd);
        reconstruct_WENO5(u, uL, uR, N, ibd);
    
        //Compute fluxes
        #pragma omp parallel for
        for (int i=0; i<N+1; i++){
            double qm = 0.5 * (uL[i] + std::abs(uL[i]));
            double qp = 0.5 * (uR[i] - std::abs(uR[i]));
            u_flux[i] = std::max(0.5 * qm*qm, 0.5 * qp*qp);
        }

        // Update solution
        #pragma omp parallel for
        for (int i=0; i<N; i++){
            u_new[i] = u[i+ibd] - dt/dx * (u_flux[i+1]-u_flux[i]);
        }

        // Update solution
        #pragma omp parallel for
        for (int i=0; i<N; i++){
            u[i+ibd] = u_new[i];
        }

        // Boundary conditions (periodic)
        Apply_BC(u,N,ibd);

        //Write solution to VTK file
        if(n_steps % 10 == 0 || n_steps == 1){
            write_vtk(u,dx,"burgers1d_step_" + std::to_string(n_steps) + ".vtk");
        }

        t += dt;
    }
    
}
