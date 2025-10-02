#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H

#include <vector>
#include "grid.h"

class Reconstruction {
  public:
    Reconstruction(Grid* mesh_) : mesh(mesh_) {
      this->n = mesh->get_n_ptr();
      this->nx = n[0];
      this->num_boundary_points = this->mesh->get_num_boundary_points();
    }

    void first_order(const double* u, double* uL, double* uR);
  protected:
    int nx;
    int* n;
    int num_boundary_points;
    Grid* mesh;

};
void reconstruct_FirstOrder(const double* u, double* uL, double* uR);
void reconstruct_WENOJS(const std::vector<double>& u, std::vector<double>& uL,  // uL[i] = u_{i-1/2}^L (left state at interface i)
    std::vector<double>& uR,  // uR[i] = u_{i-1/2}^R (right state at interface i)
    int N,int ibd);


#endif // RECONSTRUCTION_H
