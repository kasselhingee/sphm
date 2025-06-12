#include "uldSvMF.h"
#include <Rcpp.h>

//' Helper function Bessel I approximation from BesselI::besselIasym()
//' which should be from Asymptotic expansion of Bessel I_nu(x) function   x -> oo
//'       by Abramowitz & Stegun (9.7.1), p.377 
//' I_a(z) = exp(z) / sqrt(2*pi*z) * f(z,..)  where
//'   f(z,..) = 1 - (mu-1)/ (8*z) + (mu-1)(mu-9)/(2! (8z)^2) - ...
//'           = 1- (mu-1)/(8z)*(1- (mu-9)/(2(8z))*(1-(mu-25)/(3(8z))*..))
//' where  mu = 4*a^2  *and*  |arg(z)| < pi/2
//' This is useful for large x
//' @param x is the vMF concentration parameter
//' @param nu is such that nu + 1 = d/2, where d is the ambient dimension of the sphere.
// [[Rcpp::export]]
a1type besselIasym(const a1type& x, const double & nu, int order, bool log_result = true) {
  // Constants
  a1type pi = CppAD::atan(1.0) * 4.0;
  a1type pi2 = 2.0 * pi;
  
  // Precompute 8*x for efficiency
  a1type x8 = 8.0 * x;
  
  // Compute the asymptotic series for f(x, nu) up to order
  a1type d = 0.0; // Initialize d
  if (order >= 1) {
    for (int k = order; k >= 1; --k) {
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

//' For small x (i.e. concentration) Hornik and Grun use simple relation by Schou 1979 (and others)
//' for approximating the derivative of log(const(k)) but I want a coarse idea of the value here.
//' I'm going to use the series 10.25.2 from Nist: `https://dlmf.nist.gov/10.25#E2`
//' This looks actually to be just a solution to equation defining the modified Bessel function.
//' (x/2)^nu sum_i{1/i! 1/gamma(nu + i + 1) (x/2)^(2i)}.
//' nu and order are NOT differentiable
//' @param order Maximum order of series to compute
// [[Rcpp::export]]
a1type besselItrunc(const a1type& x, const double & nu, int order, bool log_result = true) {
    a1type sum = 0.0;
    a1type x2 = x / 2.0;
    double gamma_val, ifact;
    for (int i = 0; i <= order; ++i) {
        gamma_val = std::tgamma(i + nu + 1.0);
        ifact = std::tgamma(i + 1.0);  // == i!

        a1type term = CppAD::pow(x2, 2 * i + nu) / (ifact * gamma_val);
        sum = sum + term;
    }

    if (log_result)
        return CppAD::log(sum);
    else
        return sum;
}

//' This function approximates the BesselI function by
//' Using BesselItrunc for small values of x
//' Using BesselIasym for large values of x
//' @param threshold is the location at which the calculation switches
a1type besselImixed(const a1type & x, const double & nu, double threshold, int order, bool log_result) {
  // CppAD::CondExpLe returns one of two T‐typed branches
  // depending on x <= threshold
  a1type thresh = threshold;
  return CppAD::CondExpLe(
    x,                          // condition: x <= threshold?
    thresh,               // compare against this
    besselItrunc(x, nu, order, log_result),   // if true
    besselIasym(x, nu, order, log_result)     // if false
  );
}


//' Helper function lvMFnormconst_approx
//' For p == 3 using an exact formula
//' Otherwise uses *approximations* of the modified Bessel function of the first order.
//' The normalising constant is \eqn{(2 * \pi)^{p/2} besselI(k, p/2 - 1)/k^{p/2 -1}}
//' where \eqn{p} is the dimension of the ambient space of the sphere (i.e. vectors have \eqn{p} entries)
//' @details
//' Returns the log of the normalising constant.
//' The approximation uses a threshold of 10 to choose between the small concentration and large concentration regime,
//' and in each regime the series order used is 15.
// [[Rcpp::export]]
a1type lvMFnormconst_approx(a1type kappa, int p) {
  if (p == 3) {
    return CppAD::log(2 * M_PI) + CppAD::log(CppAD::exp(kappa) - CppAD::exp(-kappa)) - CppAD::log(kappa);
  } else {
    double nu = p/2 - 1.0;
    a1type log_bIval = besselImixed(kappa, nu, 10, 15, true);
    a1type out = (p/2) * CppAD::log(2 * M_PI) + log_bIval - nu * CppAD::log(kappa);
    return(out);
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
  a1type lconst = - lvMFnormconst_approx(k, p) - CppAD::log(a.coeff(0));
 
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
  a1type lconst = - lvMFnormconst_approx(k, p) - CppAD::log(a1);
  
  mata1 Hstar = getHstar(m);
  
  mata1 ystarstarL = y * Hstar;
  veca1 denomA = (y * m / a1).array().square();
  veca1 denomB = ((ystarstarL * V.inverse()).array() * ystarstarL.array()).rowwise().sum();
  veca1 denom = (denomA + denomB).array().sqrt();
  
  veca1 ll = lconst - (p - 1) * denom.array().log() 
    + (k * (y * m).array()) / (a1 * denom).array();
  
  return ll;
}

