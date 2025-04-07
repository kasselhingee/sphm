#ifndef PRELIMS2S
#define PRELIMS2S


#include <scorematchingad_forward.h> //includes <RcppEigenForward.h>
#include <Rcpp.h>

//' Preliminary Objective in the style of the `generalfunction` class:
//' @param yx is the concatenation of the response, the spherical covariates, then Euclidean covariates *cbind* together
//' @param dims_in A vector of `c(p, qe)`. `p` and `qe` are required to separate yx and omvec into their constituents.
// [[Rcpp::export]]
veca1 prelimobj_cpp(veca1 & omvec, veca1 & yx, vecd & dims_in, matd & ignore);

#endif
