#ifndef LDSVMF
#define LDSVMF

// using the forward declarations because only the translational unit for RcppExports.cpp need to have access to the wrappers
#include <scorematchingad_forward.h> //includes RcppEigen_forward.h

// [[Rcpp::depends(RcppEigen)]]

// #' Unnormalised log density of the SvMF using the cannonical parameterisation of k, a and G. For p=3, the normalised form is computed, but otherwise the normalising constant is treated as `1` because no analytic formula exist.
// [[Rcpp::export]]
veca1 uldSvMF_cann(mata1 y, a1type k, veca1 a, mata1 G);

// #' Log density of the SvMF using the muV parameterisation
// [[Rcpp::export]]
veca1 uldSvMF_muV(mata1 y, a1type k, veca1 m, a1type a1, mata1 V);

#endif

