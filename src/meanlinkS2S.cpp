
#include "meanlinkS2S.h"
#include "OmegaS2S.h"

mata1 meanlinkS2Scpp(const mata1 &x, const veca1 &vec, const int p) {
  // Convert vector to a mnlink_Omega_cpp object
  mnlink_Omega_cpp<a1type> ompar = mnlink_Omega_cpp_unvec(vec, p);

  mata1 x_t = x.transpose();

  // Extract p1, q1, and Omega
  veca1 p1 = ompar.p1;
  veca1 q1 = ompar.q1;
  mata1 Omega = ompar.Omega;

  // Compute the amount of p1 and Omegax corresponding to each column of x
  veca1 BQx2 = (Omega.transpose() * Omega * x_t).cwiseProduct(x_t).colwise().sum();
  veca1 q1x = q1.transpose() * x_t;
  veca1 p1_mult = ((q1x.array() + 1).square() - BQx2.array()) / ((q1x.array() + 1).square() + BQx2.array());
  veca1 Omegax_mult = (2 * (1 + q1x.array())) / ((q1x.array() + 1).square() + BQx2.array());

  // For each column of x combine p1 and Omega * x according to the above amounts. 
  mata1 p1_part = (p1_mult * p1.transpose()).transpose(); //each column should be p1 scaled by the entry of p1_mult
  mata1 Omegax_part = Omega * x_t;
  for (int col = 0; col < Omegax_part.cols(); ++col){
    Omegax_part.col(col) = Omegax_part.col(col) * Omegax_mult(col);
  }
  mata1 out = p1_part + Omegax_part;

  return out.transpose();
}
