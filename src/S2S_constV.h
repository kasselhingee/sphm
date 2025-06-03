#ifndef S2S_constV_H
#define S2S_constV_H

#include <scorematchingad_forward.h> //includes RcppEigenForward
#include <utils/pADFun.h>
#include "Omega.h"

// function that accepts all the major arguments in usual (non-vector) format
// for a single observation y with covariates x and return the ull for that observation, parallel transporting axes H*(p1)*Kstar to the mean using Jupp's method from base being the first column of P.
// will include projection of om
veca1 ull_S2S_constV(mata1 y, mata1 xs, mata1 xe, mnlink_Omega_cpp<a1type> om, a1type k, a1type a1, veca1 aremaining, mata1 G0);

// for checking ull_S2S_constV via unit testing
// [[Rcpp::export]]
veca1 ull_S2S_constV_forR(mata1 y, mata1 xs, mata1 xe, veca1 omvec, a1type k, a1type a1, veca1 aremaining, mata1 G0);

// axes parameterised by the Cayley transform of K* in H*(P[,1])K*
// the Cayley transform will work okay so long as K* is reasonably close to the basis given by H*. If P[,1] is reasonably close to north pole (due to standardisation of the data) and the axes reasonably close to the ambient basis axes, then K* will be reasonably close to the identity. This is good for using Cayley transforms.

// make the independent values omvec, k, aremaining and vecCayaxes,
// [[Rcpp::export]]
pADFun tape_ull_S2S_constV_nota1(veca1 omvec, a1type k, a1type a1, veca1 aremaining, mata1 G0star, vecd & p_in, vecd & qe_in, matd & yx, matd referencecoords);
#endif

