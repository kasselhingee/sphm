#include <RcppEigenForward.h>
#include <scorematchingad_forward.h>
#include "OmegaS2S.h"

// For a parameter set return quadratic distance to constraints matching
// [[Rcpp::export]]
veca1 OmegaS2S_constraints_quad(veca1 & vec, int p) {
  // Convert vector to a OmegaS2Scpp object
  OmegaS2Scpp<a1type> ompar = OmegaS2Scpp_unvec(vec, p);

  // design so that function returns zero vector when constraints satisfied
  veca1 out(1 + 1 + ompar.p1.size() + ompar.q1.size());
  out << ompar.p1.squaredNorm() - 1., ompar.q1.squaredNorm() - 1., (ompar.p1.transpose() * ompar.Omega).array().square(), (ompar.Omega * ompar.q1).array().square();
  return(out);
}

