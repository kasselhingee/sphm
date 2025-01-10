# include "OmegaS2S.h"

mnlink_Omega_cpp<double> R2OmegaS2S_double(Rcpp::List obj) {
    Eigen::Matrix<double, Eigen::Dynamic, 1> p1 = Rcpp::as<Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>>>(obj["p1"]);
    Eigen::Matrix<double, Eigen::Dynamic, 1> q1 = Rcpp::as<Eigen::Matrix<double, Eigen::Dynamic, 1>>(obj["q1"]);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Omega = Rcpp::as<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(obj["Omega"]);
    return mnlink_Omega_cpp<double>(p1, q1, Omega);
}


mnlink_Omega_cpp<a1type> R2OmegaS2S_a1type(Rcpp::List obj) {
    Eigen::Matrix<a1type, Eigen::Dynamic, 1> p1 = Rcpp::as<Eigen::Matrix<a1type, Eigen::Dynamic, 1>>(obj["p1"]);
    Eigen::Matrix<a1type, Eigen::Dynamic, 1> q1 = Rcpp::as<Eigen::Matrix<a1type, Eigen::Dynamic, 1>>(obj["q1"]);
    Eigen::Matrix<a1type, Eigen::Dynamic, Eigen::Dynamic> Omega = Rcpp::as<Eigen::Matrix<a1type, Eigen::Dynamic, Eigen::Dynamic>>(obj["Omega"]);
    return mnlink_Omega_cpp<a1type>(p1, q1, Omega);
}


