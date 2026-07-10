
#include "mobius_link_cpp.h"
#include "Omega.h"

mata1 mobius_link_cpp(const mata1 &xs, const mata1 &xe, const veca1 &vec, const int p) {
  int qe = xe.cols();
  // Convert vector to a mobius_link_Omega_cpp object
  mobius_link_Omega_cpp<a1type> ompar = mobius_link_Omega_cpp_unvec(vec, p, qe);

  mata1 xs_t = xs.transpose(); //since xs, xe are matrices of row vectors, xs_t etc are matrics of column vectors
  mata1 xe_t = xe.transpose();

  // Extract p1, q1, and Omega
  veca1 p1 = ompar.p1;
  veca1 qs1 = ompar.qs1;
  veca1 qe1 = ompar.qe1;
  veca1 ce = ompar.ce;
  mata1 Omega_s = ompar.Omega.leftCols(ompar.qs);
  mata1 Omega_e = ompar.Omega.rightCols(ompar.qe);

  // Implements Proposition 2 of the paper: accumulate ytilde = Omega*xs_stereo + Omega_e*xe_proj,
  // where xs_stereo = Sp(Qs^T xs) is the stereographic projection of the spherical covariate
  // and xe_proj is the stereographic-like Euclidean term. Omega maps the tangent displacement
  // to the response tangent space.
  mata1 ytilde(p, xs_t.cols());
  ytilde.setZero();

  if (ompar.qs > 0){
    // Stereographic projection step: Omega_s * xs / (qs1^T xs + 1).
    // qs1 is the "north-pole" direction in spherical covariate space.
    mata1 numerator = (Omega_s * xs_t); //p x xs_t.cols()
    veca1 denominator = (qs1.transpose() * xs_t).array() + 1.0;
    mata1 sph_res = numerator.array().rowwise()/denominator.transpose().array();//broadcast the denominator along each row
    ytilde = ytilde + sph_res;
  }
  if (ompar.qe > 0){
    // Euclidean stereographic step: Omega_e * xe_perp / (qe1^T xe + ce).
    // ce is a scalar offset ensuring the denominator stays positive for all training points.
    mata1 numerator = (Omega_e * xe_t);
    veca1 denominator = (qe1.transpose() * xe_t).array() + ce[0];
    mata1 Euc_res = numerator.array().rowwise()/denominator.transpose().array(); //broadcast the denominator along each row
    ytilde = ytilde + Euc_res;
  }
  // Inverse stereographic projection back to S^{p-1}: (1-||ytilde||^2)*p1 + 2*ytilde) / (1+||ytilde||^2).
  // This maps the accumulated tangent displacement to a unit vector in the response space.
  veca1 ytildesizesq = ytilde.colwise().squaredNorm();
  veca1 totdenom = ytildesizesq.array() + 1.0;
  mata1 numerator = p1 * (1.0 - ytildesizesq.array()).matrix().transpose() + 2.0 * ytilde;
  mata1 out = numerator.array().rowwise()/totdenom.transpose().array();
  return out.transpose();
}
