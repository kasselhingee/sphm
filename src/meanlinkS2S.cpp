#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

// @param x is made of column vectors of covariates
// @param vec is the vectorised form of an OmegaS2S parameterisation
// @param p is the dimension of the response (in ambient space), which is needed to separate `vec` into p1, q1 and Omega.
// @value is a matrix of column vectors of means
// [[Rcpp::export]]
Eigen::MatrixXd meanlinkS2Scpp(const Eigen::MatrixXd &x, const Eigen::VectorXd &vec, int p) {
  int n = vec.size();
  int q = (n - p) / (1 + p);

  if (q != static_cast<int>(q)) {
    Rcpp::stop("q is not an integer");
  }

  if (q <= p - 1) {
    Rcpp::stop("q must be greater than p - 1");
  }

  Rcout << "x dimensions: " << x.rows() << " x " << x.cols() << std::endl;

  // Extract p1, q1, and Omega
  Eigen::VectorXd p1 = vec.segment(0, p);
  Eigen::VectorXd q1 = vec.segment(p, q);
  Rcout << "p1 dimensions: " << p1.rows() << " x " << p1.cols() << std::endl;
  Rcout << "q1 dimensions: " << q1.rows() << " x " << q1.cols() << std::endl;
  Eigen::MatrixXd Omega = Eigen::Map<const Eigen::MatrixXd>(vec.segment(p + q, p * q).data(), p, q);
  Rcout << "Omega dimensions: " << Omega.rows() << " x " << Omega.cols() << std::endl;

  // Compute the amount of p1 and Omegax corresponding to each column of x
  Eigen::VectorXd BQx2 = (Omega.transpose() * Omega * x).cwiseProduct(x).colwise().sum();
  Rcout << "BQx2 dimensions: " << BQx2.rows() << " x " << BQx2.cols() << std::endl;
  Eigen::VectorXd q1x = q1.transpose() * x;
  Rcout << "q1x dimensions: " << q1x.rows() << " x " << q1x.cols() << std::endl;
  Eigen::VectorXd p1_mult = ((q1x.array() + 1).square() - BQx2.array()) / ((q1x.array() + 1).square() + BQx2.array());
  Rcout << "p1_mult dimensions: " << p1_mult.rows() << " x " << p1_mult.cols() << std::endl;
  Eigen::VectorXd Omegax_mult = (2 * (1 + q1x.array())) / ((q1x.array() + 1).square() + BQx2.array());
  Rcout << "Omegax_mult dimensions: " << Omegax_mult.rows() << " x " << Omegax_mult.cols() << std::endl;

  // For each column of x combine p1 and Omega * x according to the above amounts. 
  Eigen::MatrixXd p1_part = (p1_mult * p1.transpose()).transpose(); //each column should be p1 scaled by the entry of p1_mult
  Rcout << "p1_part dimensions: " << p1_part.rows() << " x " << p1_part.cols() << std::endl;
  Eigen::MatrixXd Omegax_part = Omega * x;
  Rcout << "Omegax_part dimensions: " << Omegax_part.rows() << " x " << Omegax_part.cols() << std::endl;
  for (int col = 0; col < Omegax_part.cols(); ++col){
    Omegax_part.col(col) = Omegax_part.col(col) * Omegax_mult(col);
  }
  Eigen::MatrixXd out = p1_part + Omegax_part;

  return out;
}
