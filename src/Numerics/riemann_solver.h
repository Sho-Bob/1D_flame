#ifndef RIEMANN_SOLVER_H
#define RIEMANN_SOLVER_H

#include <vector>

namespace riemann_solver {
  // HLLC Riemann solver for Euler equations
  void HLLC(
      double& Frho, double* Frhou, double& FrhoE, double* Frho_scalars,
      const double P_L, const double P_R,
      const double rho_L, const double rho_R, const double rhoE_L, const double rhoE_R,
      const double sos_L, const double sos_R,
      const double* Y_L, const double* Y_R,
      const double* u_L, const double* u_R,
      const double lambda,
      const double* nVec, const int n_scal
      );
} // namespace riemann_solver

#endif // RIEMANN_SOLVER_H

