#ifndef COMMON_MACROS_H
#define COMMON_MACROS_H
// Common macros for array indexing and loops
#define FOR_IDIM for (int iDim = 0; iDim < this->num_dim; ++iDim)
#define FOR_ICV(iDim) for (int icv = this->num_boundary_points; icv < this->n[iDim] + this->num_boundary_points; ++icv)
#define FOR_ICV_G(iDim) for (int icv = 0; icv < this->n[iDim] + 2 * this->num_boundary_points; ++icv)
#define FOR_IFA(iDim) for (int ifa = 0; ifa < this->n[iDim] + 1; ++ifa)
#define LOOP_l_N(N) for (int l = 0; l < N; ++l)
#define LOOP_k_N(N) for (int k = 0; k < N; ++k)
#define LOOP_I3 for (int i = 0; i < 3; ++i)
#define FOR_VAR for (int v = 0; v < this->num_vars; ++v)

#define DOT_PRODUCT(a, b) (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])

#endif // COMMON_MACROS_H
