#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include "Grid/grid.h"
#include "Numerics/reconstruction.h"
#include "Common/common.h"
#include "IO/hdf5_writer.h"

class Solver { // base clas for Solver

  public:

    // \brief: run the solver; this is the function called from main
    virtual void run_solver();
    // \brief: initialize the solver; set up grid, initial condition
    virtual void initialize();
    // \brief: enforce boundary conditions weakly
    virtual void apply_bc();
    // \brief: gradient computation; face reconstruction
    virtual void pre_rhs();
    // \brief: rhs computation
    virtual void rhs();
    // \brief: update conservatives
    virtual void update_conservatives(const int idrk);

    // initialize grid and profile
    virtual void initialize_grid();
    virtual void initialize_profile();

    // output solution to file
    virtual void output() {};
  protected:

    // =================================================
    // grid parameters
    // =================================================
    int num_dim; // number of dimensions
    int nx,ny,nz;
    int num_cells;   // number of cells
    int num_boundary_points; // number of ghost cells
    double x0, x1; // domain [x0,x1]
    double y0, y1; // domain [y0,y1]
    double z0, z1; // domain [z0,z1]
    double dx, dy, dz; // grid spacing

    Grid* mesh;
    int* n; // number of grid points in each dimension
    int* n_tot; // number of grid points in each dimension including ghost cells
    double** coordinates; // coordinate locations
    double* Delta; // grid spacing in each dimension

    // coordinates of the grid points
    std::vector<double> x, y, z;

    // =================================================
    // Numerics
    // =================================================
    Reconstruction* reconstruction; // reconstruction object
    // =================================================
    // solution variables
    // =================================================
    int num_vars; // number of variables (e.g., for Navier-Stokes: rho, rho*u, rho*v, rho*w, E)
    double **conservatives;
    double **conservatives_old; // for SSP-RK time-stepping
    double **rhs_conservatives;

    // =================================================
    // time-stepping parameters
    // =================================================
    double dt; // time step size
    double t_end; // end time
    int num_time_steps; // number of time steps
    int step; // current time step
    double t; // current time
};

void burgers_initialize(std::vector<double>& u, int N, double x0, double x1, double dx, int ibd);
void simulate_burgers1d(std::vector<double>& u, double dx, double CFL, double T_end, int ibd);
void Apply_BC(std::vector<double>& u, int N, int ibd);
void reconstruct_MUSCL_minmod(const std::vector<double>& u, std::vector<double>& uL, std::vector<double>& uR, int N, int ibd);
void reconstruct_WENO5(const std::vector<double>& u, std::vector<double>& uL, std::vector<double>& uR, int N, int ibd);
#endif
