
#include "meanlinkS2S.h"
#include "OmegaS2S.h"

mata1 meanlinkS2Scpp(const mata1 &xs, const mata1 &xe, const veca1 &vec, const int p) {
  int qe = xe.cols();
  // Convert vector to a mnlink_Omega_cpp object
  mnlink_Omega_cpp<a1type> ompar = mnlink_Omega_cpp_unvec(vec, p, qe);
  Rcpp::Rcout << "ompar unveced" << std::endl;

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
  Rcpp::Rcout << p1 << std::endl;
  Rcpp::Rcout << qs1 << std::endl;
  Rcpp::Rcout << qe1 << std::endl;
  Rcpp::Rcout << ce1 << std::endl;
  Rcpp::Rcout << PBce << std::endl;
  Rcpp::Rcout << "Omega_s:" << std::endl;
  Rcpp::Rcout << Omega_s << std::endl;
  Rcpp::Rcout << "Omega_e:" << std::endl;
  Rcpp::Rcout << Omega_e << std::endl;

  //following proposition 2, but with extension for denominator below Omega_e with offset given by first element of ce.
  mata1 ytilde(p, xs_t.cols());
  ytilde.setZero();

  //ytilde contribution from spherical covar
  if (ompar.qs > 0){
    mata1 numerator = (Omega_s * xs_t); //p x xs_t.cols()
//    Rcpp::Rcout << "Sph numerator finished" << std::endl;
//    Rcpp::Rcout << numerator << std::endl;
    veca1 denominator = (qs1.transpose() * xs_t).array() + 1.0;
//    Rcpp::Rcout << "Sph denominator finished" << std::endl;
//    Rcpp::Rcout << (qs1.transpose() * xs_t) << std::endl;
//    Rcpp::Rcout << denominator << std::endl;
    mata1 sph_res = numerator.array().rowwise()/denominator.transpose().array();//broadcast the denominator along each row
//    Rcpp::Rcout << "Sph part calculated" << std::endl;
//    Rcpp::Rcout << sph_res << std::endl;
    ytilde = ytilde + sph_res;
//    Rcpp::Rcout << ytilde << std::endl;
  }
  if (ompar.qe > 0){
    mata1 numerator = (Omega_e * xe_t).colwise() + PBce; //this is something called broadcasting in Eigen
    Rcpp::Rcout << "Euc numerator finished" << std::endl;
    Rcpp::Rcout << xe_t.rows() << std::endl;
    Rcpp::Rcout << xe_t.cols() << std::endl;
    Rcpp::Rcout << qe1.rows() << std::endl;
    Rcpp::Rcout << qe1.transpose() * xe_t << std::endl;
    veca1 denominator = (qe1.transpose() * xe_t).array() + ce1[0];
    Rcpp::Rcout << "Euc denominator finished" << std::endl;
    mata1 Euc_res = numerator.array().rowwise()/denominator.transpose().array(); //broadcast the denominator along each row
    Rcpp::Rcout << "Euc part calculated" << std::endl;
    ytilde = ytilde + Euc_res; 
  }
  Rcpp::Rcout << ytilde << std::endl;
  veca1 ytildesizesq = ytilde.colwise().squaredNorm();  //will this produce a row vector?
  Rcpp::Rcout << "ytildesizesq finished" << std::endl;
  Rcpp::Rcout << ytildesizesq << std::endl;
  veca1 totdenom = ytildesizesq.array() + 1.0;
  mata1 numerator = p1 * (1.0 - ytildesizesq.array()).matrix().transpose() + 2.0 * ytilde;
  Rcpp::Rcout << "Part A" << std::endl;
  Rcpp::Rcout << p1 * (1.0 - ytildesizesq.array()).matrix().transpose() << std::endl;
  Rcpp::Rcout << "Part B" << std::endl;
  Rcpp::Rcout << 2.0 * ytilde << std::endl;
  Rcpp::Rcout << "total numerator finished" << std::endl;
  Rcpp::Rcout << numerator << std::endl;
  mata1 out = numerator.array().rowwise()/totdenom.transpose().array();
  Rcpp::Rcout << "out finished" << std::endl;
  Rcpp::Rcout << out << std::endl;
  return out.transpose();
}
