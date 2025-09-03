#ifndef MOBIUSVMF
#define MOBIUSVMF


#include <scorematchingad_forward.h> //includes <RcppEigenForward.h>
#include <Rcpp.h>

//' Preliminary Objective in the style of the `generalfunction` class:
//' @param yx is the response the spherical covariates, and the Euclidean covariates *cbind* together. Each row an observation.
//' @param dyn is a zero length vector that is ignored
//' @param dims_in A vector of `c(p, qe)`. `p` and `qe` are required to separate yx and omvec into their constituents.
//' @details
//' The return is vector of mu.y values where mu is the mean predicted by omvec and x, and y are the observations.
//' This mu.y is (log(vMF density) - log(vMF norm const(k)))/k.
//' Therefore the gradient of below function, for each return value, is: grad(mu.y) = grad(log(vMF density))/k
//' The Fisher information matrix ignoring concentration is var(grad(log(vMF density))) = k^2 var(grad(mu.y)).
//' So standard errors for the preliminary mean require knowledge of the concentration parameter k.
// [[Rcpp::export]]
veca1 prelimobj_cpp(veca1 & omvec, veca1 & dyn, vecd & dims_in, matd & yx);

#endif
