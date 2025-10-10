#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H

#include <vector>
#include "Grid/grid.h"

class Reconstruction {
  public:
    Reconstruction(Grid* mesh_) : mesh(mesh_) {
      this->n = mesh->get_n_ptr();
      this->nx = n[0];
      this->num_boundary_points = this->mesh->get_num_boundary_points();
    }

    void reconstruct(const double* u, double* uL, double* uR, const std::string& method = "FirstOrder") {
      if (method == "FirstOrder") {
        first_order(u, uL, uR);
      } else if (method == "WENOJS") {
        WENOJS(u, uL, uR);
      } else if (method == "WENOZ") {
        WENOZ(u, uL, uR);
      } else if (method == "MUSCL") {
        MUSCL(u, uL, uR);
      } else {
        throw std::runtime_error("Reconstruction method not implemented: " + method);
      }
    }

    void first_order(const double* u, double* uL, double* uR);
    void WENOJS(const double* u, double* uL, double* uR);
    void WENOZ(const double* u, double* uL, double* uR);
    void MUSCL(const double* u, double* uL, double* uR);
    void set_kappa_muscl(double kappa) { this->kappa_muscl = kappa; }


  protected:
    int nx;
    int* n;
    int num_boundary_points;
    Grid* mesh;

    // =================================================
    // MUSCL stuff
    // =================================================
    double kappa_muscl = 1.0 / 3.0;

};
void reconstruct_FirstOrder(const double* u, double* uL, double* uR);
void reconstruct_WENOJS(const std::vector<double>& u, std::vector<double>& uL,  // uL[i] = u_{i-1/2}^L (left state at interface i)
    std::vector<double>& uR,  // uR[i] = u_{i-1/2}^R (right state at interface i)
    int N,int ibd);


#endif // RECONSTRUCTION_H
