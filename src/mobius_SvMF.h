#ifndef MOBIUSSVMF
#define MOBIUSSVMF

#include <scorematchingad_forward.h> //includes RcppEigenForward
#include <utils/pADFun.h>
#include "Omega.h"

// function that accepts all the major arguments in usual (non-vector) format and returns log-density for each row of `y`
// for a single observation y with covariates x and return the ull for that observation, parallel transporting axes H*(p1)*Kstar to the mean using Jupp's method from base being the first column of P.
// will include projection of om
veca1 uld_Mobius_SvMF_partran(mata1 y, mata1 xs, mata1 xe, mnlink_Omega_cpp<a1type> om, a1type k, a1type a1, veca1 aremaining, mata1 G0);

// for checking uld_Mobius_SvMF_partran via unit testing
// [[Rcpp::export]]
veca1 uld_Mobius_SvMF_partran_forR(mata1 y, mata1 xs, mata1 xe, veca1 omvec, a1type k, a1type a1, veca1 aremaining, mata1 G0);

// axes parameterised by the Cayley transform of K* in H*(P[,1])K*
// the Cayley transform will work okay so long as K* is reasonably close to the basis given by H*. If P[,1] is reasonably close to north pole (due to standardisation of the data) and the axes reasonably close to the ambient basis axes, then K* will be reasonably close to the identity. This is good for using Cayley transforms.

// make the independent values omvec, k, aremaining and vecCayaxes,
// [[Rcpp::export]]
pADFun tape_uld_Mobius_SvMF_partran_nota1(veca1 omvec, a1type k, a1type a1, veca1 aremaining, mata1 G0star, vecd & p_in, vecd & qe_in, matd & yx, matd referencecoords, std::string G01behaviour);
#endif

