#include "grid.h"
#include <iostream>

Grid::Grid(int nx, int num_boundary_points, double x0, double x1) : num_boundary_points(num_boundary_points), x0(x0), x1(x1) {
  this->xmin = x0; this->xmax = x1;
  this->num_dim = 1;  // hard coded
  this->dx = (this->xmax - this->xmin) / (static_cast<double>(nx) - 1.0);

  this->coordinates = new double*[num_dim];
  this->n = new int[num_dim];
  this->n_tot = new int[num_dim];
  this->Delta = new double[num_dim];
  for (int iDim = 0; iDim < num_dim; ++iDim) {
    this->n[iDim] = nx; // hard coded
    this->n_tot[iDim] = nx + 2 * this->num_boundary_points;
    this->coordinates[iDim] = new double[this->n_tot[iDim]];

    // internal cells
    for (int i = this->num_boundary_points; i < this->num_boundary_points + this->n[iDim]; ++i) {
      coordinates[iDim][i] = this->xmin + (i - this->num_boundary_points) * dx;
    }
    // fill in ghost cells
    // left ghost cells
    for (int i = 0; i < this->num_boundary_points; ++i) {
      coordinates[iDim][i] = this->xmin - dx * (this->num_boundary_points - i);
    }
    // right ghost cells
    for (int i = this->num_boundary_points + this->n[iDim]; i < this->n_tot[iDim]; ++i) {
      coordinates[iDim][i] = this->xmax + dx * (i - (this->num_boundary_points + this->n[iDim] - 1));
    }

    this->Delta[iDim] = this->coordinates[iDim][1] - this->coordinates[iDim][0];

  } // loop over dimensions
} // end Grid constructor


Grid::~Grid() {
  for (int iDim = 0; iDim < num_dim; ++iDim) {
    delete[] coordinates[iDim];
  }
  delete[] coordinates;
  delete[] n;
  delete[] n_tot;
} // end Grid destructor
