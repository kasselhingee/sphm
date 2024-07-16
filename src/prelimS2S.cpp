#include <RcppEigenForward.h>
#include <scorematchingad_forward.h>
#include "OmegaS2S.h"

// For a parameter set return quadratic distance to constraints matching
veca1 OmegaS2S_constraints_quad(OmegaS2Scpp<a1type>& ompar) {
  veca1 out(1 + 1 + ompar.p1.size() + ompar.q1.size());
  out << ompar.p1.squaredNorm() - 1., ompar.q1.squaredNorm() - 1., (ompar.p1.transpose() * ompar.Omega).array().square(), (ompar.Omega * ompar.q1).array().square();
  return(out);
}

