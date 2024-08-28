#include "ldSvMF.h"

// Helper function vMFnormconst
a1type vMFnormconst(a1type kappa, int p) {
  if (p == 3) {
    return 2 * M_PI * (CppAD::exp(kappa) - CppAD::exp(-kappa)) / kappa;
  } else {
    Rcpp::stop("vMF normalising constant for p != 3 not implemented yet");
  }
}

// Helper getHstar
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

veca1 ldSvMF_cann(mata1 y, a1type k, veca1 a, mata1 G) {
  int p = a.size();
  a1type lconst = - CppAD::log(vMFnormconst(k, p)) - CppAD::log(a.coeff(0));
 
  // Scale columns of G by the corresponding elements of a
  mata1 Gscal = G.array().rowwise() / a.transpose().array();
 
  // Compute the denominator sqrt(rowSums((y %*% Gscal)^2))
  veca1 denom = (y * Gscal).rowwise().squaredNorm().cwiseSqrt();
  
  veca1 ll = lconst - (p - 1) * denom.array().log() 
    + (k * (y * Gscal.col(0)).array()) / denom.array();
  
  return ll;
}

veca1 ldSvMF_muV(mata1 y, a1type k, veca1 m, a1type a1, mata1 V) {
  int p = m.size();
  a1type lconst = - CppAD::log(vMFnormconst(k, p)) - CppAD::log(a1);
  
  mata1 Hstar = getHstar(m);
  
  mata1 ystarstarL = y * Hstar;
  veca1 denomA = (y * m / a1).array().square();
  veca1 denomB = ((ystarstarL * V.inverse()).array() * ystarstarL.array()).rowwise().sum();
  veca1 denom = (denomA + denomB).array().sqrt();
  
  veca1 ll = lconst - (p - 1) * denom.array().log() 
    + (k * (y * m).array()) / (a1 * denom).array();
  
  return ll;
}

