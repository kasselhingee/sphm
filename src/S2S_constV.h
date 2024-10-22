#ifndef S2S_constV_H
#define S2S_constV_H

#include <scorematchingad_forward.h> //includes RcppEigenForward
#include <utils/pADFun.h>
#include "OmegaS2S.h"

// function that accepts all the major arguments in usual (non-vector) format
// for a single observation y with covariates x and return the ull for that observation, parallel transporting axes H*(p1)*Kstar to the mean using Jupp's method from base being the first column of P.
// will include projection of om
veca1 ull_S2S_constV(mata1 y, mata1 x, OmegaS2Scpp<a1type> om, a1type k, a1type a1, veca1 aremaining, mata1 Kstar);

// for checking ull_S2S_constV via unit testing
// [[Rcpp::export]]
veca1 ull_S2S_constV_forR(mata1 y, mata1 x, veca1 omvec, a1type k, a1type a1, veca1 aremaining, mata1 Kstar);

//tape for mean wrt mean with the base point P[,1] fixed for calculating the axes. Returned tape should have dynamic parameters of: k, a1, aremaining, as.vector(Gstar).
pADFun tape_ull_S2S_constV_mean(veca1 & omvec, a1type k, a1type a1, veca1 aremaining, mata1 Gstar, vecd & p_in, matd & yx);

// tape for axes
// axes parameterised by the Cayley transform of K* in H*(P[,1])K*
// the Cayley transform will work okay so long as K* is reasonably close to the basis given by H*. If P[,1] is reasonably close to north pole (due to standardisation of the data) and the axes reasonably close to the ambient basis axes, then K* will be reasonably close to the identity. This is good for using Cayley transforms.
//returned tape should have as dynamic parameters: P1, as.vector(om), k, a1, aremaining
pADFun tape_ull_S2S_constV_axes(veca1 & vecCayaxes, OmegaS2Scpp<a1type> om, a1type k, a1type a1, veca1 aremaining, vecd & p_in, matd & yx);

// tape for aremaining. Like _mean but everything else constant
pADFun tape_ull_S2S_constV_aremaining(veca1 & aremaining, OmegaS2Scpp<a1type> om, a1type k, a1type a1, mata1 Gstar, vecd & p_in, matd & yx);


// tape for k
pADFun tape_ull_S2S_constV_k(veca1 k, OmegaS2Scpp<a1type> om, a1type a1, veca1 & aremaining, mata1 Gstar, vecd & p_in, matd & yx);

// experimental: tape for everything all together except a1
// make the independent values omvec, k, aremaining and vecCayaxes,
// [[Rcpp::export]]
pADFun tape_ull_S2S_constV_nota1(veca1 & omvec, a1type k, a1type a1, veca1 aremaining, mata1 Kstar, vecd & p_in, matd & yx);
#endif

