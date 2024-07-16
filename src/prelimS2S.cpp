#include <RcppEigenForward.h>
#include <scorematchingad_forward.h>
#include "OmegaS2S.h"

// Preliminary Objective::

// For a parameter set return quadratic distance to constraints matching
// [[Rcpp::export]]
veca1 OmegaS2S_constraints_quad(veca1 & vec, int p) {
  // Convert vector to a OmegaS2Scpp object
  OmegaS2Scpp<a1type> ompar = OmegaS2Scpp_unvec(vec, p);

  // design so that function returns zero vector when constraints satisfied
  int p1_size = ompar.p1.size();
  int q1_size = ompar.q1.size();
  veca1 out(1 + 1 + p1_size + q1_size);
  out(0) = ompar.p1.squaredNorm() - 1.;
  out(1) = ompar.q1.squaredNorm() - 1.;
  out.segment(2, q1_size) = (ompar.p1.transpose() * ompar.Omega).array().square();
  out.segment(2 + q1_size, p1_size) = (ompar.Omega * ompar.q1).array().square();
  return(out);
}

