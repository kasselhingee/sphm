#include "Omega.h"
#include "utils.h"

veca1 Omega_constraints(veca1 & vec, int p, int qe) {
  // Convert vector to a mnlink_Omega_cpp object
  mnlink_Omega_cpp<a1type> ompar = mnlink_Omega_cpp_unvec(vec, p, qe);

  // design so that function returns zero vector when constraints satisfied
  
  veca1 sphcheck(0);
  veca1 Euccheck(0);
  if (ompar.qs>0){
    sphcheck.resize(1);
    sphcheck(0) = ompar.qs1.squaredNorm() - 1.;
  }
  if (ompar.qe>0){
    Euccheck.resize(1);
    Euccheck(0) = ompar.qe1.squaredNorm() - 1.;
  }

  // commutivity constraint check // because OmOm = OmpartOmpart(Euc) + OmpartOmpart(Sph), checking both is redundant
  // only needed when BOTH qe > 0 and qs > 0
  // require that commutivity of the *projected* Omega holds
  // Since projected, there are only (p-1) vectors to be orthogonal to each other -->
  // I suspect that means only (p-1)*(p-2)/2 unique constraints, which ones though!?? a sum would avoid this
  veca1 commutecheck(0);
  if ((ompar.qs > 0) && (ompar.qe > 0)){
    mnlink_Omega_cpp<a1type> ompar_proj = Omega_proj_cpp(ompar);
    commutecheck.resize(1);
    // Compute Omega * Omega^T for commutivity constraint
    mata1 OmOm = ompar_proj.Omega * ompar_proj.Omega.transpose();
    mata1 Is_tilde;
    Is_tilde = mata1::Zero(ompar_proj.qs + ompar_proj.qe, ompar_proj.qs);
    Is_tilde.topRows(ompar_proj.qs) = mata1::Identity(ompar_proj.qs, ompar_proj.qs);
    mata1 OmpartOmpart = ompar_proj.Omega * Is_tilde * Is_tilde.transpose() * ompar_proj.Omega.transpose();
    mata1 commutediff = OmOm * OmpartOmpart - OmpartOmpart * OmOm; //OmOm etc are always symmetric, so commutediff is always antisymmetric
    //rotate commutediff so that p1 is the northpole
    veca1 nthpole = vecd::Unit(ompar_proj.p, 0);
    mata1 rotmat = JuppRmat(ompar_proj.p1, nthpole);
    Rcpp::Rcout << "commutediff:" << std::endl << commutediff << std::endl;
    commutecheck(0) = commutediff.norm();
  }
    
  veca1 out(1 + sphcheck.size() + Euccheck.size() + commutecheck.size());
  out << ompar.p1.squaredNorm() - 1., sphcheck, Euccheck, commutecheck;
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
  out(0) = ssq_sv - ((ompar.qs>0) + (ompar.qe>0)) * (p-1.);
  return(out);
}



