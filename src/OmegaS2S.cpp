# include "OmegaS2S.h"

OmegaS2S<double> listToOmegaS2S_double(List obj) {
    Eigen::Matrix<T, Eigen::Dynamic, 1> p1 = Rcpp::as<Eigen::Matrix<T, Dynamic, 1>>(obj["p1"]);
    Eigen::Matrix<T, Eigen::Dynamic, 1> q1 = Rcpp::as<Eigen::Matrix<T, Eigen::Dynamic, 1>>(obj["q1"]);
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Omega = Rcpp::as<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>(obj["Omega"]);
    return OmegaS2S<T>(p1, q1, Omega);
}


