
#include "meanlinkS2S.h"
#include "OmegaS2S.h"

mata1 meanlinkS2Scpp(const mata1 &xs, const mata1 &xe, const veca1 &vec, const int p) {
  int qe = xe.cols();
  // Convert vector to a mnlink_Omega_cpp object
  mnlink_Omega_cpp<a1type> ompar = mnlink_Omega_cpp_unvec(vec, p, qe);

  mata1 xs_t = xs.transpose(); //since xs, xe are matrices of row vectors, xs_t etc are matrics of column vectors
  mata1 xe_t = xe.transpose();

  // Extract p1, q1, and Omega
  veca1 p1 = ompar.p1;
  veca1 qs1 = ompar.qs1;
  veca1 qe1 = ompar.qe1;
  veca1 ce1 = ompar.ce1;
  veca1 PBce = ompar.PBce;
  mata1 Omega_s = ompar.Omega.leftCols(ompar.qs);
  mata1 Omega_e = ompar.Omega.rightCols(ompar.qe);

  //following proposition 2, but with extension for denominator below Omega_e with offset given by first element of ce.
  mata1 ytilde(p, xs_t.cols());
  ytilde.setZero();

  //ytilde contribution from spherical covar
  if (ompar.qs > 0){
    mata1 numerator = (Omega_s * xs_t); //p x xs_t.cols()
    veca1 denominator = (qs1.transpose() * xs_t).array() + 1.0;
    mata1 sph_res = numerator.array().rowwise()/denominator.transpose().array();//broadcast the denominator along each row
    ytilde = ytilde + sph_res;
  }
  if (ompar.qe > 0){
    mata1 numerator = (Omega_e * xe_t).colwise() + PBce; //this is something called broadcasting in Eigen
    veca1 denominator = (qe1.transpose() * xs_t).array() + ce1[0];
    mata1 Euc_res = numerator.array().rowwise()/denominator.transpose().array(); //broadcast the denominator along each row
    ytilde = ytilde + Euc_res; 
  }
  veca1 ytildesizesq = ytilde.colwise().squaredNorm();  //will this produce a row vector?
  veca1 totdenom = ytildesizesq.array() + 1.0;
  mata1 numerator = p1 * (1.0 - ytildesizesq.array()).matrix().transpose() + 2.0 * ytilde;
  mata1 out = numerator.array().rowwise()/totdenom.transpose().array();
  return out.transpose();
}
