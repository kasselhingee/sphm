#ifndef MOBIUSSVMF
#define MOBIUSSVMF

#include <scorematchingad_forward.h> //includes RcppEigenForward
#include <utils/pADFun.h>
#include "Omega.h"

//' @noRd
//' @title Log-density at response `y` for SvMF regression with the scales Mobius link function
//' @description For given response row-vectors y, covariate row-vectors xs, xe and regression parameters, returns the log-density for each row of `y` (and corresponding row of xs and xe). When the length of the response (`p`) is not 3, then the normalising constant from the vMF density is approximated.
//' @param y A matrix of row vectors (each row a unit vector)
//' @param xs A matrix of row vectors (each row a unit vector)
//' @param xe A matrix of row vectors
//' @param om A parameter object specifying the parameters of the mean link
//' @param k Concentration of the SvMF error distribution
//' @param a1 The first scale of the SvMF error distribution, which is typically fixed at `1`.
//' @param aremaining A (p-1) vector of the remaining scales of the SvMF error distribution.
//' @param G0 An orthonormal matrix with first column that are the 'base location' for parallel transport of axes on the sphere.
//' @details
//' For a given mean direction `mu`, the axes of the SvMF error distribution are obtained by parallel transport of `G0[,-1]` along the geodesic from `G0[,1]` to `mu`.
//' 
//' The function projects the parameter object `om` such that the matrix Omega is orthogonal to p1, qs1 and qe1.
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

