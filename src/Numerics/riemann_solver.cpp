#include "riemann_solver.h"

#include "Common/common.h"

#include <vector>
#include <cmath>
#include <iostream>


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
      ) {
    // check positvity of the values first
    if (rho_L <= 0.0 || rho_R <= 0.0 || P_L < 0.0 || P_R < 0.0 || sos_L < 0.0 || sos_R < 0.0) {
      std::cerr << "riemann_solver::HLLC: negative density, sos or pressure" << std::endl;
      std::cout << "rho_L = " << rho_L << ", P_L = " << P_L << std::endl;
      std::cout << "rho_R = " << rho_R << ", P_R = " << P_R << std::endl;
      std::cout << "sos_L = " << sos_L << ", sos_R = " << sos_R << std::endl;
      throw std::runtime_error("riemann_solver::HLLC: negative density, sos or pressure");
    }

    // set mixtures to get necessary quantities
    // Left-state
    double uN_L = DOT_PRODUCT(u_L, nVec);
    double u_norm_L = std::sqrt(DOT_PRODUCT(u_L, u_L));

    // Right-state
    double uN_R = DOT_PRODUCT(u_R, nVec);
    double u_norm_R = std::sqrt(DOT_PRODUCT(u_R, u_R));

    double consv_L[5 + n_scal];
    double Flux_L [5 + n_scal];
    consv_L[0] = rho_L;
    Flux_L[0] = rho_L * uN_L;
    LOOP_I3 {
      consv_L[i + 1] = rho_L * u_L[i];
      // consv_L[i + 1] = rho_L * (u_L[i] + (uN_L - u_L[0]) * nVec[i]); // fix the normal velocity term
      Flux_L[i + 1] = rho_L * uN_L * u_L[i] + P_L * nVec[i];
    }
    consv_L[4] = rhoE_L;
    Flux_L[4] = (rhoE_L + P_L) * uN_L;
    LOOP_l_N(n_scal) {
      consv_L[5 + l] = rho_L * Y_L[l];
      Flux_L[5 + l] = rho_L * uN_L * Y_L[l];
    }

    double consv_R[5 + n_scal];
    double Flux_R [5 + n_scal];
    consv_R[0] = rho_R;
    Flux_R[0] = rho_R * uN_R;
    LOOP_I3 {
      consv_R[i + 1] = rho_R * u_R[i];
      // consv_R[i + 1] = rho_R * (u_R[i] + (uN_R - u_R[0]) * nVec[i]); // fix the normal velocity term
      Flux_R[i + 1] = rho_R * uN_R * u_R[i] + P_R * nVec[i];
    }
    consv_R[4] = rhoE_R;
    Flux_R[4] = (rhoE_R + P_R) * uN_R;
    LOOP_l_N(n_scal) {
      consv_R[5 + l] = rho_R * Y_R[l];
      Flux_R[5 + l] = rho_R * uN_R * Y_R[l];
    }

    // wave-speed estimate
    double SL = std::min(uN_L - sos_L, uN_R - sos_R);
    double SR = std::max(uN_L + sos_L, uN_R + sos_R);
    double S_star = (P_R - P_L + rho_L * uN_L * (SL - uN_L) - rho_R * uN_R * (SR - uN_R))
      / (rho_L * (SL - uN_L) - rho_R * (SR - uN_R));
    if (std::isnan(S_star) || std::isinf(S_star)) {
        std::cout << "riemann_solver::HLLC: S_star is NaN or Inf" << std::endl;
        std::cout << "P_L = " << P_L << ", P_R = " << P_R << std::endl;
        std::cout << "rho_L = " << rho_L << ", rho_R = " << rho_R << std::endl;
        std::cout << "sos_L = " << sos_L << ", sos_R = " << sos_R << std::endl;
        exit(-1);
    }
    double pStar = P_L + rho_L * (SL - uN_L) * (S_star - uN_L);

    double consv_L_star[5 + n_scal]; double consv_R_star[5 + n_scal]; 
    double factor_L = rho_L * (SL - uN_L) / (SL - S_star); double factor_R = rho_R * (SR - uN_R) / (SR - S_star);
    consv_L_star[0] = factor_L * 1.0; consv_R_star[0] = factor_R * 1.0;
    LOOP_I3 {
      // consv_L_star[1+i] = factor_L * (u_L[i] + (S_star - uN_L) * nVec[i]); 
      // consv_R_star[1+i] = factor_R * (u_R[i] + (S_star - uN_R) * nVec[i]);
      consv_L_star[1+i] = factor_L * (u_L[i] + (S_star - uN_L) * nVec[i]); 
      consv_R_star[1+i] = factor_R * (u_R[i] + (S_star - uN_R) * nVec[i]);
    }
    consv_L_star[4] = factor_L * (rhoE_L / rho_L + (S_star - uN_L) * (S_star + P_L/(rho_L*(SL - uN_L))));
    consv_R_star[4] = factor_R * (rhoE_R / rho_R + (S_star - uN_R) * (S_star + P_R/(rho_R*(SR - uN_R))));
    LOOP_l_N(n_scal) {
      consv_L_star[5 + l] = factor_L * Y_L[l];
      consv_R_star[5 + l] = factor_R * Y_R[l];
    }

    // compute fluxes
    double Flux[5 + n_scal];
    if (SL >= 0) {
      // left state is supersonic
      LOOP_l_N (5 + n_scal) Flux[l] = Flux_L[l];
    } else if (SR < 0) {
      // right state is supersonic
      LOOP_l_N (5 + n_scal) Flux[l] = Flux_R[l];
    } else if (SL < 0 && 0 <= S_star) {
      LOOP_l_N (5 + n_scal) Flux[l] = Flux_L[l] + SL * (consv_L_star[l] - consv_L[l]);
    } else if (S_star < 0 && 0 <= SR) {
      LOOP_l_N (5 + n_scal) Flux[l] = Flux_R[l] + SR * (consv_R_star[l] - consv_R[l]);
    } else {
      throw std::runtime_error("riemann_solver::HLLC: "
          "Invalid wave speed configuration. "
          "SL = " + std::to_string(SL) + ", S_star = " + std::to_string(S_star) +
          ", SR = " + std::to_string(SR));
    }

    /*
       printf("SL = %.4e, S_star = %.4e, SR = %.4e\n", SL, S_star, SR);
       printf("Flux_L[0] = %.4e, consv_L_star[0] = %.4e, consv_L[0] = %.4e\n",
       Flux_L[0], consv_L_star[0], consv_L[0]);
       printf("Flux_R[0] = %.4e, consv_R_star[0] = %.4e, consv_R[0] = %.4e\n",
       Flux_R[0], consv_R_star[0], consv_R[0]);
       */

    // set fluxes
    Frho = Flux[0]; Frhou[0] = Flux[1]; Frhou[1] = Flux[2]; Frhou[2] = Flux[3];
    FrhoE = Flux[4];
    LOOP_l_N (n_scal) {
      Frho_scalars[l] = Flux[5 + l];
    }
  } // end of HLLC

} // namespace riemann_solver
