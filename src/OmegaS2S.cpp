#include "OmegaS2S.h"

veca1 OmegaS2S_constraints(veca1 & vec, int p) {
  // Convert vector to a mnlink_Omega_cpp object
  mnlink_Omega_cpp<a1type> ompar = mnlink_Omega_cpp_unvec(vec, p, 0);

  // design so that function returns zero vector when constraints satisfied
  veca1 out(1 + 1);
  out(0) = ompar.p1.squaredNorm() - 1.;
  out(1) = ompar.qs1.squaredNorm() - 1.;
  return(out);
}

//a wrap around OmegaS2S_constraints for use with tapegeneral
veca1 wrap_OmegaS2S_constraints(veca1 & vec, veca1 & ignore1, vecd & p_in, matd & ignore2) {
  veca1 out;
  int p = int(p_in(0) + 0.1);
  out = OmegaS2S_constraints(vec,p);
  return(out);
}

//Constraints on the singular values of Omega - not exact unfortunately, just on total sum
veca1 OmegaS2S_ineqconstaints(veca1 & vec, veca1 & ignore1, vecd & p_in, matd & ignore2){
  int p = int(p_in(0) + 0.1);
  // Convert vector to a mnlink_Omega_cpp object
  mnlink_Omega_cpp<a1type> ompar = mnlink_Omega_cpp_unvec(vec, p, 0);
  a1type ssq_sv = (ompar.Omega.transpose() * ompar.Omega).diagonal().sum();
  veca1 out(1);
  out(0) = ssq_sv - (p-1.);
  return(out);
}



