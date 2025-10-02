#ifndef COMMON_MACROS_H
#define COMMON_MACROS_H
// Common macros for array indexing and loops
#define FOR_ICV(iDim) for (int icv = this->num_boundary_points; icv < this->n[iDim] + this->num_boundary_points; ++icv)
#define FOR_ICV_G(iDim) for (int icv = 0; icv < this->n[iDim] + 2 * this->num_boundary_points; ++icv)
#define FOR_IFA(iDim) for (int ifa = 0; ifa < this->n[iDim] + 1; ++ifa)
#endif // COMMON_MACROS_H
