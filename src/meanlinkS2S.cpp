#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::MatrixXd meanlinkS2Scpp(const Eigen::VectorXd &x, const NumericVector &vec, int p) {
  int n = vec.size();
  int q = (n - p) / (1 + p);

  if (q != (int)q) {
    stop("q is not an integer");
  }

  if (q <= p - 1) {
    stop("q must be greater than p - 1");
  }

  // Extract p1, q1, and Omega
  Eigen::VectorXd p1 = Eigen::Map<Eigen::VectorXd>(vec.begin(), p);
  Eigen::VectorXd q1 = Eigen::Map<Eigen::VectorXd>(vec.begin() + p, q);
  Eigen::MatrixXd Omega = Eigen::Map<Eigen::MatrixXd>(vec.begin() + p + q, p, q);

  // Perform matrix calculations
  Eigen::VectorXd BQx2 = (Omega.transpose() * Omega * x).cwiseProduct(x).colwise().sum();
  double q1x = q1.dot(x);
  Eigen::VectorXd p1_mult = ((1 + q1x) * (1 + q1x) - BQx2).array() / ((1 + q1x) * (1 + q1x) + BQx2).array();
  Eigen::VectorXd Omegax_mult = (2 * (1 + q1x)) / ((1 + q1x) * (1 + q1x) + BQx2).array();

  // Compute the outer product
  Eigen::MatrixXd out = p1_mult * p1.transpose() + Omegax_mult.asDiagonal() * Omega * x;

  return out;
}
