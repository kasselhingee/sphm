#ifndef SVMF_LL_ALIGNEDP_H
#define SVMF_LL_ALIGNEDP_H

#include <scorematchingad_forward.h> //includes RcppEigenForward

//' The unnormalised log-likelihood of a SvMF Sphere-Sphere Regression with Mobius Mean Link and Variance Axes Aligned with P.
//' @description Three functions in a format compatible with `tapefun`. A difficulty is that the P matrix needs and SVD to get out of the Omega parameterisation so two functions for alternating between optimising the mean (with `kappa` and `a` fixed, and `G` aligned to a fixed `P`) and optimising a, with the mean fixed and concentration fixed, and updating `P` inbetween.
//' @param vec For `_alignedG_mean`: A parameter vector specifying the mean via the Omega vectorisation.
//' @param dyn For `_alignedG_mean`: A p+1+p*p length vector of kappa, then a1, a2, ..., then P as a vector of stacked columns.
//' @param p_in The dimension p
//' @param yx The observations and covariates cbind together as row vectors
// [[Rcpp::export]]
veca1 ull_S2S_alignedG_mean(veca1 & vec, veca1 & dyn, vecd & p_in, matd & yx);

//' @param vec For `_alignedG_a`: A p-2 vector of log(a3), log(a4), log(a5), ... log(a2) will be calculated as the negative sum of the others to satisfy the prod=1 constraint on a2,...
//' @param dyn For `_alignedG_a`: A vector of kappa then a1
//' @param pOmegavec For `_alignedG_a`: A vector of p then the Omega vectorisation, then `as.vector(P)`. Due to an SVD to extract P from Omega vec, taping the dependence on Omega would be unreliable. Furthermore R's SVD routine seems more reliable than Eigen's.
// [[Rcpp::export]]
veca1 ull_S2S_alignedG_a(veca1 & vec, veca1 & dyn, vecd & pOmegavecP, matd & yx);


//' @param k For `_alignedG_k`: A parameter vector specifying the concentration k
//' @param dyn For `_alignedG_k`: A p*q + p + p*p length vector of Omegavec, then a1, a2, ..., then P as a vector of stacked columns.
//' The P and Omega are provided seperately because of the non-smooth nature of SVD.
// [[Rcpp::export]]
veca1 ull_S2S_alignedG_k(veca1 & k, veca1 & dyn, vecd & p_in, matd & yx);

//' @describeIn ull_S2S_alignedG_mean Aligns the columns of the Mobius-link rotation matrix `P` for the mean to the columns of G. Note that the first column of the returned G is the given mean. Returns the matrix G.
// [[Rcpp::export]]
mata1 alignedGcpp(veca1 m, mata1 P);

#endif

