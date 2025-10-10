#include "reconstruction.h"

#include "Common/common.h"
// #include "limiter.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <omp.h>
#include <cassert>

void Reconstruction::first_order(const double* u, double* uL, double* uR) {
  FOR_IFA(0) {
    int icv0 = this->mesh->get_icv0(ifa);
    int icv1 = this->mesh->get_icv1(ifa);
    uL[ifa] = u[icv0];
    uR[ifa] = u[icv1];
  }
}

void Reconstruction::WENOJS(const double* u, double* uL, double* uR) {
  FOR_IFA(0) {
    // =========================================================
    // construct for ifa = i+1/2
    //           i+1/2
    // --o--|--o--|--o--|--o--|--o--|--o--
    //  i-2   i-1    i    i+1   i+2   i+3
    // =========================================================
    int icv_i      = this->mesh->get_icv0(ifa);     // left cell index
    int icv_iPlus1 = this->mesh->get_icv1(ifa);     // right cell index
    int icv_iMinus1 = icv_i - 1;
    int icv_iMinus2 = icv_i - 2;
    int icv_iPlus2  = icv_iPlus1 + 1;
    int icv_iPlus3  = icv_iPlus1 + 2;

    // ---------------------------------------------------------
    // Left state reconstruction at i+1/2 (biased to the left)
    // ---------------------------------------------------------
    double v0 = u[icv_iMinus2];
    double v1 = u[icv_iMinus1];
    double v2 = u[icv_i];
    double v3 = u[icv_iPlus1];
    double v4 = u[icv_iPlus2];

    // candidate polynomials
    double p0 = ( 2.0*v0 - 7.0*v1 + 11.0*v2 ) / 6.0;
    double p1 = ( -1.0*v1 + 5.0*v2 + 2.0*v3 ) / 6.0;
    double p2 = ( 2.0*v2 + 5.0*v3 - 1.0*v4 ) / 6.0;

    // smoothness indicators β_k
    double beta0 = (13.0/12.0)*pow(v0 - 2*v1 + v2, 2) + (1.0/4.0)*pow(v0 - 4*v1 + 3*v2, 2);
    double beta1 = (13.0/12.0)*pow(v1 - 2*v2 + v3, 2) + (1.0/4.0)*pow(v1 - v3, 2);
    double beta2 = (13.0/12.0)*pow(v2 - 2*v3 + v4, 2) + (1.0/4.0)*pow(3*v2 - 4*v3 + v4, 2);

    // linear weights γ_k
    double g0 = 0.1, g1 = 0.6, g2 = 0.3;
    double eps = 1e-6;

    // nonlinear weights ω_k
    double a0 = g0 / pow(eps + beta0, 2);
    double a1 = g1 / pow(eps + beta1, 2);
    double a2 = g2 / pow(eps + beta2, 2);
    double sum = a0 + a1 + a2;

    double w0 = a0 / sum;
    double w1 = a1 / sum;
    double w2 = a2 / sum;

    uL[ifa] = w0*p0 + w1*p1 + w2*p2;

    // ---------------------------------------------------------
    // Right state reconstruction at i+1/2 (biased to the right)
    // mirror stencil
    // ---------------------------------------------------------
    double w0r, w1r, w2r;
    {
      double q0 = ( -1.0*v1 + 5.0*v2 + 2.0*v3 ) / 6.0;
      double q1 = ( 2.0*v2 + 5.0*v3 - 1.0*v4 ) / 6.0;
      double q2 = ( 11.0*v3 - 7.0*v4 + 2.0*u[icv_iPlus3] ) / 6.0;

      double beta0r = (13.0/12.0)*pow(v1 - 2*v2 + v3, 2) + (1.0/4.0)*pow(v1 - 4*v2 + 3*v3, 2);
      double beta1r = (13.0/12.0)*pow(v2 - 2*v3 + v4, 2) + (1.0/4.0)*pow(v2 - v4, 2);
      double beta2r = (13.0/12.0)*pow(v3 - 2*v4 + u[icv_iPlus3], 2) + (1.0/4.0)*pow(3*v3 - 4*v4 + u[icv_iPlus3], 2);

      double a0r = g0 / pow(eps + beta0r, 2);
      double a1r = g1 / pow(eps + beta1r, 2);
      double a2r = g2 / pow(eps + beta2r, 2);
      double sumr = a0r + a1r + a2r;

      w0r = a0r / sumr;
      w1r = a1r / sumr;
      w2r = a2r / sumr;

      uR[ifa] = w0r*q0 + w1r*q1 + w2r*q2;
    }
  } // loop over ifa
} // end WENOJS

void Reconstruction::WENOZ(const double* u, double* uL, double* uR) {
  FOR_IFA(0) {
    // =========================================================
    // construct for ifa = i+1/2
    //           i+1/2
    // --o--|--o--|--o--|--o--|--o--|--o--
    //  i-2   i-1    i    i+1   i+2   i+3
    // =========================================================
    int icv_i      = this->mesh->get_icv0(ifa);     // left cell index
    int icv_iPlus1 = this->mesh->get_icv1(ifa);     // right cell index
    int icv_iMinus1 = icv_i - 1;
    int icv_iMinus2 = icv_i - 2;
    int icv_iPlus2  = icv_iPlus1 + 1;
    int icv_iPlus3  = icv_iPlus1 + 2;

    // ---------------------------------------------------------
    // Left state reconstruction at i+1/2 (biased to the left)
    // ---------------------------------------------------------
    double v0 = u[icv_iMinus2];
    double v1 = u[icv_iMinus1];
    double v2 = u[icv_i];
    double v3 = u[icv_iPlus1];
    double v4 = u[icv_iPlus2];

    // candidate polynomials
    double p0 = ( 2.0*v0 - 7.0*v1 + 11.0*v2 ) / 6.0;
    double p1 = ( -1.0*v1 + 5.0*v2 + 2.0*v3 ) / 6.0;
    double p2 = ( 2.0*v2 + 5.0*v3 - 1.0*v4 ) / 6.0;

    // smoothness indicators β_k
    double beta0 = (13.0/12.0)*pow(v0 - 2*v1 + v2, 2) + (1.0/4.0)*pow(v0 - 4*v1 + 3*v2, 2);
    double beta1 = (13.0/12.0)*pow(v1 - 2*v2 + v3, 2) + (1.0/4.0)*pow(v1 - v3, 2);
    double beta2 = (13.0/12.0)*pow(v2 - 2*v3 + v4, 2) + (1.0/4.0)*pow(3*v2 - 4*v3 + v4, 2);

    // linear weights γ_k
    double g0 = 0.1, g1 = 0.6, g2 = 0.3;
    double eps = 1e-6;

    double tau5 = fabs(beta0 - beta2);

    double a0 = g0 * (1.0 + pow(tau5 / (eps + beta0), 2));
    double a1 = g1 * (1.0 + pow(tau5 / (eps + beta1), 2));
    double a2 = g2 * (1.0 + pow(tau5 / (eps + beta2), 2));
    double sum = a0 + a1 + a2;

    double w0 = a0 / sum;
    double w1 = a1 / sum;
    double w2 = a2 / sum;

    uL[ifa] = w0*p0 + w1*p1 + w2*p2;
    // ---------------------------------------------------------
    // Right state reconstruction at i+1/2 (biased to the right)
    // mirror stencil
    // ---------------------------------------------------------
    double w0r, w1r, w2r;
    {
      double q0 = ( -1.0*v1 + 5.0*v2 + 2.0*v3 ) / 6.0;
      double q1 = ( 2.0*v2 + 5.0*v3 - 1.0*v4 ) / 6.0;
      double q2 = ( 11.0*v3 - 7.0*v4 + 2.0*u[icv_iPlus3] ) / 6.0;

      double beta0r = (13.0/12.0)*pow(v1 - 2*v2 + v3, 2) + (1.0/4.0)*pow(v1 - 4*v2 + 3*v3, 2);
      double beta1r = (13.0/12.0)*pow(v2 - 2*v3 + v4, 2) + (1.0/4.0)*pow(v2 - v4, 2);
      double beta2r = (13.0/12.0)*pow(v3 - 2*v4 + u[icv_iPlus3], 2) + (1.0/4.0)*pow(3*v3 - 4*v4 + u[icv_iPlus3], 2);

      double tau5r = fabs(beta0r - beta2r);

      double a0r = g0 * (1.0 + pow(tau5r / (eps + beta0r), 2));
      double a1r = g1 * (1.0 + pow(tau5r / (eps + beta1r), 2));
      double a2r = g2 * (1.0 + pow(tau5r / (eps + beta2r), 2));
      double sumr = a0r + a1r + a2r;

      w0r = a0r / sumr;
      w1r = a1r / sumr;
      w2r = a2r / sumr;

      uR[ifa] = w0r*q0 + w1r*q1 + w2r*q2;
    }
  } // loop over ifa
} // end WENOZ

// =================================================================
// Monotonic Upstream-centered Scheme for Conservation Laws (MUSCL)
// =================================================================
// Reference:
// 1. B. van Leer, H. Nishikawa, Towards the ultimate understandign of MUSCL: Pitfalls in achieving third-order accuracy, Journal of Computational Physics, 446 (2021)
// look at equations (20-23) for the reconstruction formulas
// =================================================================
// =================================================================
void Reconstruction::MUSCL(const double* u, double* uL, double* uR) {
  FOR_IFA(0) {
    // =========================================================
    // construct for ifa = i+1/2
    //           i+1/2
    // --o--|--o--|--o--|--o--|--o--|--o--
    //  i-2   i-1    i    i+1   i+2   i+3
    // =========================================================
    int icv_i      = this->mesh->get_icv0(ifa);     // left cell index
    int icv_iPlus1 = this->mesh->get_icv1(ifa);     // right cell index
    int icv_iMinus1 = icv_i - 1;
    int icv_iMinus2 = icv_i - 2;
    int icv_iPlus2  = icv_iPlus1 + 1;

    double delta_uj = u[icv_i] - u[icv_iMinus1];         // u_i - u_{i-1}
    double delta_ujPlus1 = u[icv_iPlus1] - u[icv_i];   // u_{i+1} - u_i
    double delta_ujPlus2 = u[icv_iPlus2] - u[icv_iPlus1]; // u_{i+2} - u_{i+1}

    double kappa = this->kappa_muscl; // 1/3 for third-order accuracy

    uL[ifa] = u[icv_i] + 0.25 * ((1. - kappa) * delta_uj + (1. + kappa) * delta_ujPlus1);
    uR[ifa] = u[icv_iPlus1] - 0.25 * ((1. + kappa) * delta_ujPlus1 + (1. - kappa) * delta_ujPlus2);
  } // loop over ifa
} // end MUSCL



void reconstruct_WENOJS(
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
