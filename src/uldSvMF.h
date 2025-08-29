#ifndef LDSVMF
#define LDSVMF

// using the forward declarations because only the translational unit for RcppExports.cpp need to have access to the wrappers
#include <scorematchingad_forward.h> //includes RcppEigen_forward.h

// [[Rcpp::depends(RcppEigen)]]

// @title Log-density of the scaled von Mises-Fisher distribution
// @param y a matrix with each row an observation
// @param k concentration
// @param a the scales of the scaled von Mises-Fisher distribution
// @param G An orthogonal matrix with first column the mean of the distribution and remaining columns the orientation axes in the same order as `a`
// @returns A vector of log-density evaluated at row of `y`
// @details For p=3, the normalised form is computed, but otherwise the normalising constant is approximated using `besselImixed()`.
// [[Rcpp::export]]
veca1 uldSvMF_cann(mata1 y, a1type k, veca1 a, mata1 G);

// #' Log density of the SvMF using the muV parameterisation
// [[Rcpp::export]]
veca1 uldSvMF_muV(mata1 y, a1type k, veca1 m, a1type a1, mata1 V);

// a helper
mata1 getHstar(veca1 m);

//' This function approximates the BesselI function by
//' Using BesselItrunc for small values of x
//' Using BesselIasym for large values of x
//' @param threshold is the location at which the calculation switches
// [[Rcpp::export]]
a1type besselImixed(const a1type & x, const double & nu, double threshold, int order, bool log_result = true);

#endif

