#include "OmegaS2S.h"

veca1 Omega_constraints(veca1 & vec, int p, int qe=0) {
  // Convert vector to a mnlink_Omega_cpp object
  mnlink_Omega_cpp<a1type> ompar = mnlink_Omega_cpp_unvec(vec, p, qe);

  // design so that function returns zero vector when constraints satisfied
  veca1 out(1 + (ompar.qs>0) + (ompar.qe>0));
  out(0) = ompar.p1.squaredNorm() - 1.;
  if (ompar.qs>0){
    out(1) = ompar.qs1.squaredNorm() - 1.;
  }
  if (ompar.qe>0){
    out(1 + (ompar.qs>0)) = ompar.qe1.squaredNorm() - 1.
  }
  return(out);
}

//a wrap around Omega_constraints for use with tapegeneral
veca1 Omega_constraints_wrap(veca1 & vec, veca1 & ignore1, vecd & dims_in, matd & ignore2) {
  veca1 out;
  int p = int(dims_in(0) + 0.1);
  int qe = int(dims_in(1) + 0.1);
  out = Omega_constraints(vec,p,qe);
  return(out);
}

//Constraints on the singular values of Omega - not exact unfortunately, just on total sum
veca1 Omega_ineqconstraints(veca1 & vec, veca1 & ignore1, vecd & p_in, matd & ignore2){
  int p = int(p_in(0) + 0.1);
  // Convert vector to a mnlink_Omega_cpp object
  mnlink_Omega_cpp<a1type> ompar = mnlink_Omega_cpp_unvec(vec, p, 0);
  a1type ssq_sv = (ompar.Omega.transpose() * ompar.Omega).diagonal().sum();
  veca1 out(1);
  out(0) = ssq_sv - (p-1.);
  return(out);
}



