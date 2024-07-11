#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
MatrixXd meanlinkS2S(const Eigen::VectorXd &x, const NumericVector &vec, int p) {
  int n = vec.size();
  int q = (n - p) / (1 + p);

  if (q != (int)q) {
    stop("q is not an integer");
  }

  if (q <= p - 1) {
    stop("q must be greater than p - 1");
  }

  // Extract p1, q1, and Omega
  VectorXd p1 = Map<VectorXd>(vec.begin(), p);
  VectorXd q1 = Map<VectorXd>(vec.begin() + p, q);
  MatrixXd Omega = Map<MatrixXd>(vec.begin() + p + q, p, q);

  // Perform matrix calculations
  VectorXd BQx2 = (Omega.transpose() * Omega * x).cwiseProduct(x).colwise().sum();
  double q1x = q1.dot(x);
  VectorXd p1_mult = ((1 + q1x) * (1 + q1x) - BQx2).array() / ((1 + q1x) * (1 + q1x) + BQx2).array();
  VectorXd Omegax_mult = (2 * (1 + q1x)) / ((1 + q1x) * (1 + q1x) + BQx2).array();

  // Compute the outer product
  MatrixXd out = p1_mult * p1.transpose() + Omegax_mult.asDiagonal() * Omega * x;

  return out;
}
