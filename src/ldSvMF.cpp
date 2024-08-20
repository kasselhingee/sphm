#include <RcppEigenForward.h>
#include <scorematchingad_forward.h>
#include <Rcpp.h>

// Helper function vMFnormconst
a1type vMFnormconst(a1type kappa, int p) {
  if (p == 3) {
    return 2 * M_PI * (CppAD::exp(kappa) - CppAD::exp(-kappa)) / kappa;
  } else {
    Rcpp::stop("vMF normalising constant for p != 3 not implemented yet");
  }
}

mata1 getHstar(veca1 m) {
  a1type m1 = m(0);  // First element of m
  veca1 mL = m.tail(m.size() - 1);  // Vector m without the first element
  
  // Compute the matrix (1/(1+m1)) * mL %*% t(mL) - diag(1, length(mL))
  mata1 mL_outer = mL * mL.transpose();  // mL %*% t(mL)
  mata1 identity = mata1::Identity(mL.size(), mL.size());
  
  mata1 Hstar = mata1(mL.size() + 1, mL.size());
  Hstar.row(0) = mL.transpose();  // First row is mL
  Hstar.block(1, 0, mL.size(), mL.size()) = (1 / (1 + m1)) * mL_outer - identity;  // Remaining rows
  
  return Hstar;
}

// [[Rcpp::export]]
veca1 ldSvMF(mata1 y, Rcpp::List obj) {
  a1type k = Rcpp::as<a1type>(obj["k"]);
  veca1 m = Rcpp::as<veca1>(obj["m"]);
  a1type a1 = Rcpp::as<a1type>(obj["a1"]);
  mata1 V = Rcpp::as<mata1>(obj["V"]);
  
  int p = m.size();
  a1type lconst = - CppAD::log(vMFnormconst(k, p)) - CppAD::log(a1);
  
  // Assuming getHstar is implemented elsewhere
  mata1 Hstar = getHstar(m); // You'll need to define this
  
  mata1 ystarstarL = y * Hstar;
  veca1 denom = ((y * m / a1).array().square()
                             + ((ystarstarL * V.inverse()).array() * ystarstarL.array()).rowwise().sum().square()).sqrt();
  
  veca1 ll = lconst - (p - 1) * denom.array().log() 
    + (k * (y * m).array()) / (a1 * denom).array();
  
  return ll;
}
