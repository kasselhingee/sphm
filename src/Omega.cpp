#include "Omega.h"

veca1 Omega_constraints(veca1 & vec, int p, int qe) {
  // Convert vector to a mnlink_Omega_cpp object
  mnlink_Omega_cpp<a1type> ompar = mnlink_Omega_cpp_unvec(vec, p, qe);

  // design so that function returns zero vector when constraints satisfied
  
  // Compute Omega * Omega^T for commutivity constraint
  mata1 OmOm = ompar.Omega * ompar.Omega.transpose();
  veca1 sphcheck(0);
  veca1 Euccheck(0);
  if (ompar.qs>0){
    sphcheck.resize(1 + ompar.p * (ompar.p - 1)/2);
    sphcheck(0) = ompar.qs1.squaredNorm() - 1.;
    // commutivity constraint
    mata1 Is_tilde = mata1::Zero(ompar.qs + ompar.qe, ompar.qs);
    Is_tilde.topRows(ompar.qs) = mata1::Identity(ompar.qs, ompar.qs);
    mata1 OmpartOmpart = ompar.Omega * Is_tilde * Is_tilde.transpose() * ompar.Omega.transpose();
    mata1 commutediff = OmOm * OmpartOmpart - OmpartOmpart * OmOm; //OmOm etc are always symmetric, so commutediff is always antisymmetric
    // place lower triangular elements into sphcheck
    int idx = 1;
    for (int i = 0; i < p; ++i) {
        for (int j = 0; j < i; ++j) {
            sphcheck(idx) = commutediff(i, j);
            idx++;
        }
    }
  }
  if (ompar.qe>0){
    Euccheck.resize(1 + ompar.p * (ompar.p - 1)/2);//ompar.p * (ompar.p - 1)/2);
    Euccheck(0) = ompar.qe1.squaredNorm() - 1.;
    // commutivity constraint
    mata1 Ie_tilde = mata1::Zero(ompar.qs + ompar.qe, ompar.qe);
    Ie_tilde.bottomRows(ompar.qe) = mata1::Identity(ompar.qe, ompar.qe);
    mata1 OmpartOmpart = ompar.Omega * Ie_tilde * Ie_tilde.transpose() * ompar.Omega.transpose();
    mata1 commutediff = OmOm * OmpartOmpart - OmpartOmpart * OmOm; //OmOm etc are always symmetric, so commutediff is always antisymmetric
    // place lower triangular elements into sphcheck
    int idx = 1;
    for (int i = 0; i < ompar.p; ++i) {
        for (int j = 0; j < i; ++j) {
            Euccheck(idx) = commutediff(i, j);
            idx++;
        }
    }
  }
  veca1 out(1 + sphcheck.size() + Euccheck.size());
  out << ompar.p1.squaredNorm() - 1., sphcheck, Euccheck;
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



