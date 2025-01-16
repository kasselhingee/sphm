#ifndef MEANLINKS2S
#define MEANLINKS2S

// using the forward declarations because only the translational unit for RcppExports.cpp need to have access to the wrappers
#include <scorematchingad_forward.h> //includes <RcppEigenForward.h>

// [[Rcpp::depends(RcppEigen)]]

// @param x is made of *row* vectors of covariates
// @param vec is the vectorised form of an OmegaS2S parameterisation
// @param p is the dimension of the response (in ambient space), which is needed to separate `vec` into p1, q1 and Omega.
// @value is a matrix of *row* vectors of means
// [[Rcpp::export]]
mata1 meanlinkS2Scpp(const mata1 &xs, const mata1 &xe, const veca1 &vec, const int p);

#endif

