# ifndef LDSVMF
# define LDSVMF

// using the forward declarations because only the translational unit for RcppExports.cpp need to have access to the wrappers
#include <RcppEigenForward.h>
#include <scorematchingad_forward.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppEigen)]]

// #' Log density of the SvMF using the cannonical parameterisation of k, a and G
// [[Rcpp::export]]
veca1 ldSvMF_cann(mata1 y, a1type k, veca1 a, mata1 G);

// #' Log density of the SvMF using the muV parameterisation
// [[Rcpp::export]]
veca1 ldSvMF_muV(mata1 y, a1type k, veca1 m, a1type a1, mata1 V);

# endif

