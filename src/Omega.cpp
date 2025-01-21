#include "Omega.h"

veca1 Omega_constraints(veca1 & vec, int p, int qe) {
  // Convert vector to a mnlink_Omega_cpp object
  mnlink_Omega_cpp<a1type> ompar = mnlink_Omega_cpp_unvec(vec, p, qe);

  // design so that function returns zero vector when constraints satisfied
  veca1 out(1 + (ompar.qs>0) + (ompar.qe>0));
  out(0) = ompar.p1.squaredNorm() - 1.;
  if (ompar.qs>0){
    out(1) = ompar.qs1.squaredNorm() - 1.;
  }
  if (ompar.qe>0){
    out(1 + (ompar.qs>0)) = ompar.qe1.squaredNorm() - 1.;
  }
  return(out);
}

//a wrap around Omega_constraints for use with tapegeneral
veca1 Omega_constraints_wrap(veca1 & vec, veca1 & ignore1, vecd & dims_in, matd & ignore2) {
  veca1 out;
  if (dims_in.size() != 2){Rcpp::stop("dims_in must have two entries");}
  int p = int(dims_in(0) + 0.1);
  int qe = int(dims_in(1) + 0.1);
  out = Omega_constraints(vec,p,qe);
  return(out);
}

//Constraints on the singular values of Omega - not exact unfortunately, just on total sum
//mirrors singularvalssumsquared check in `parameterisations.R`
veca1 Omega_ineqconstraints(veca1 & vec, veca1 & ignore1, vecd & dims_in, matd & ignore2){
  if (dims_in.size() != 2){Rcpp::stop("dims_in must have two entries");}
  int p = int(dims_in(0) + 0.1);
  int qe = int(dims_in(1) + 0.1);
  // Convert vector to a mnlink_Omega_cpp object
  mnlink_Omega_cpp<a1type> ompar = mnlink_Omega_cpp_unvec(vec, p, qe);
  a1type ssq_sv = (ompar.Omega.transpose() * ompar.Omega).diagonal().sum();
  veca1 out(1);
  out(0) = ssq_sv - (ompar.qs>0 + ompar.qe>0) * (p-1.);
  return(out);
}



