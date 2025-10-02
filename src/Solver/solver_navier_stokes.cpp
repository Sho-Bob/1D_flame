#include "solver_navier_stokes.h"
// #include "limiter.h"
#include "IO/vtk_writer.h"
#include "Numerics/reconstruction.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <omp.h>
#include <iostream>
#include <cassert>

// ==============================================================
// \brief: initialize multispecies solver
// ==============================================================
void NavierStokesSolver::initialize(){
}

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

