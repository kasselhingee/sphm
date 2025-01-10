
#include "meanlinkS2S.h"
#include "OmegaS2S.h"

mata1 meanlinkS2Scpp(const mata1 &xs, const mata1 &xe, const veca1 &vec, const int p, const int qe = 0) {
  // Convert vector to a mnlink_Omega_cpp object
  mnlink_Omega_cpp<a1type> ompar = mnlink_Omega_cpp_unvec(vec, p, qe);

  mata1 xs_t = xs.transpose();
  mata1 xe_t = xe.transpose();

  // Extract p1, q1, and Omega
  veca1 p1 = ompar.p1;
  veca1 qs1 = ompar.qs1;
  veca1 qe1 = ompar.qe1;
  veca1 ce = ompar.ce;
  mata1 Omega_s = ompar.Omega.leftCols(ompar.qs);
  mata1 Omega_e = ompar.Omega.rightCols(ompar.qe);

  //following proposition 2, but with extension for denominator below Omega_e with offset given by first element of ce.
  veca1 ytilde(xe_t.rows());
  ytilde.setZero();

  //ytilde contribution from spherical covar
  if (ompar.qs > 0){
    veca1 numerator = (Omega_s * xs_t);
    veca1 denominator = (qs1.transpose() * xs_t).array() + 1.0;
    ytilde = ytilde + (numerator.array()/denominator.array()).matrix().transpose();
  }
  if (ompar.qe > 0){
    veca1 numerator = (Omega_e * xe_t).colwise() + cstar; //this is something called broadcasting in Eigen
    veca1 denominator = (qe1.transpose() * xs_t).array() + c1[0];
    ytilde = ytilde + (numerator.array()/denominator.array()).matrix().transpose();
  }
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
