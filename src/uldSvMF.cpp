#include "uldSvMF.h"
#include <Rcpp.h>

//' Helper function Bessel I approximation from BesselI::besselIasym()
//' which should be from Asymptotic expansion of Bessel I_nu(x) function   x -> oo
//'       by Abramowitz & Stegun (9.7.1), p.377 
//' I_a(z) = exp(z) / sqrt(2*pi*z) * f(z,..)  where
//'   f(z,..) = 1 - (mu-1)/ (8*z) + (mu-1)(mu-9)/(2! (8z)^2) - ...
//'           = 1- (mu-1)/(8z)*(1- (mu-9)/(2(8z))*(1-(mu-25)/(3(8z))*..))
//' where  mu = 4*a^2  *and*  |arg(z)| < pi/2
// [[Rcpp::export]]
a1type besselIasym(const a1type& x, const a1type& nu, int k_max, bool log_result = true) {
  // Constants
  a1type pi = CppAD::atan(1.0) * 4.0;
  a1type pi2 = 2.0 * pi;
  
  if (x <= pi/2.0) {
    Rcpp::stop("x must be larger than pi/2")
    a1type::abort_recording();
  }
  
  
  // Precompute 8*x for efficiency
  a1type x8 = 8.0 * x;
  
  // Compute the asymptotic series for f(x, nu) up to order k_max
  a1type d = 0.0; // Initialize d
  if (k_max >= 1) {
    for (int k = k_max; k >= 1; --k) {
      // mu = (2*(nu-k)+1)*(2*(nu+k)-1)
      a1type term1 = 2.0 * (nu - k) + 1.0;
      a1type term2 = 2.0 * (nu + k) - 1.0;
      a1type mu = term1 * term2;
      d = (1.0 - d) * mu / (k * x8);
    }
  }
  
  // Compute the result
  if (log_result) {
    // log(f) = log1p(-d)
    a1type log_f = CppAD::log1p(-d);
    return x + log_f - 0.5 * CppAD::log(pi2 * x);
  } else {
    a1type f = 1.0 - d;
    a1type scaling_factor = CppAD::exp(x);
    return scaling_factor * f / CppAD::sqrt(pi2 * x);
  }
}

// Helper function lvMFnormconst
a1type lvMFnormconst(a1type kappa, int p) {
  if (p == 3) {
    return CppAD::log(2 * M_PI) + CppAD::log(CppAD::exp(kappa) - CppAD::exp(-kappa)) - CppAD::log(kappa);
  } else {
    return 1.;
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

veca1 uldSvMF_cann(mata1 y, a1type k, veca1 a, mata1 G) {
  int p = a.size();
  a1type lconst = - lvMFnormconst(k, p) - CppAD::log(a.coeff(0));
 
  // Scale columns of G by the corresponding elements of a
  mata1 Gscal = G.array().rowwise() / a.transpose().array();
 
  // Compute the denominator sqrt(rowSums((y %*% Gscal)^2))
  veca1 denom = (y * Gscal).rowwise().squaredNorm().cwiseSqrt();
  
  veca1 ll = lconst - (p - 1) * denom.array().log() 
    + (k * (y * Gscal.col(0)).array()) / denom.array();
  
  return ll;
}

veca1 uldSvMF_muV(mata1 y, a1type k, veca1 m, a1type a1, mata1 V) {
  int p = m.size();
  a1type lconst = - lvMFnormconst(k, p) - CppAD::log(a1);
  
  mata1 Hstar = getHstar(m);
  
  mata1 ystarstarL = y * Hstar;
  veca1 denomA = (y * m / a1).array().square();
  veca1 denomB = ((ystarstarL * V.inverse()).array() * ystarstarL.array()).rowwise().sum();
  veca1 denom = (denomA + denomB).array().sqrt();
  
  veca1 ll = lconst - (p - 1) * denom.array().log() 
    + (k * (y * m).array()) / (a1 * denom).array();
  
  return ll;
}

