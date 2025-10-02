#include "solver_navier_stokes.h"
// #include "limiter.h"
#include "vtk_writer.h"
#include "reconstruction.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <omp.h>
#include <iostream>
#include <cassert>

void NavierStokesSolver::initialize(){
// #pragma omp parallel for
//   for(int i = ibd; i < N+ibd; i++){
//     double x = x0 + (i-ibd) * dx;
//     u[i] = std::sin(2.0 * M_PI * x);
//   }
//   Apply_BC(u,N,ibd);
}

void NavierStokesSolver::apply_bc(){
//   /// Neumann boundary conditions
// #pragma omp parallel for
//   for (int i=0; i<ibd; i++){
//     u[i] = u[ibd];
//     u[N+ibd+i] = u[N+ibd-1];
//   }
}

void NavierStokesSolver::pre_rhs(){
}
