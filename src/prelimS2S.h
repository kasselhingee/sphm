#ifndef PRELIMS2S
#define PRELIMS2S


#include <scorematchingad_forward.h> //includes <RcppEigenForward.h>
#include <Rcpp.h>

//' Preliminary Objective in the style of the `generalfunction` class:
//' @param yx is the response and covariates *cbind* together. Each row an observation.
//' @param dyn is a zero length vector
//' @param p is required to separate yx and omvec. It is passed as a double for compatiblility witn generalfunction, so will have to be rounded to a integer within the function
//' @param dyn ignored
// [[Rcpp::export]]
veca1 pobjS2Scpp(veca1 & omvec, veca1 & dyn, vecd & p_in, matd & yx);

#endif
