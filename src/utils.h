#ifndef UTILS_H
#define UTILS_H

#include <scorematchingad_forward.h> //includes RcppEigenForward

// Jupp rotation matrix mapping unit vector y1 to unit vector y2.
// Formula: (y1+y2)(y1+y2)^T / (1 + y1.y2) - I
// This is the rotation of minimum geodesic angle carrying y1 to y2.
// Negating gives parallel_transport_mat() in R/rotatedresiduals.R.
// Undefined when y1 = -y2 (antipodal); callers must avoid this case.
mata1 jupp_Rmat(const veca1 & y1, const veca1 & y2);

#endif
