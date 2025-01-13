# include "OmegaS2S.h"

mnlink_Omega_cpp<double> R2mnlink_Omega_cpp_double(Rcpp::List obj) {
    Eigen::Matrix<double, Eigen::Dynamic, 1> p1 = Rcpp::as<Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>>>(obj["p1"]);
    Eigen::Matrix<double, Eigen::Dynamic, 1> qs1 = Rcpp::as<Eigen::Matrix<double, Eigen::Dynamic, 1>>(obj["qs1"]);
    Eigen::Matrix<double, Eigen::Dynamic, 1> qe1 = Rcpp::as<Eigen::Matrix<double, Eigen::Dynamic, 1>>(obj["qe1"]);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Omega = Rcpp::as<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(obj["Omega"]);
    Eigen::Matrix<double, Eigen::Dynamic, 1> ce1 = Rcpp::as<Eigen::Matrix<double, Eigen::Dynamic, 1>>(obj["ce1"]);
    Eigen::Matrix<double, Eigen::Dynamic, 1> PBce = Rcpp::as<Eigen::Matrix<double, Eigen::Dynamic, 1>>(obj["PBce"]);
    return mnlink_Omega_cpp<double>(p1, qs1, qe1, Omega, ce1, PBce);
}


mnlink_Omega_cpp<a1type> R2mnlink_Omega_cpp_a1type(Rcpp::List obj) {
    Rcpp::CharacterVector names = obj.names();
    veca1 p1 = Rcpp::as<veca1>(obj["p1"]);
    veca1 qs1 = Rcpp::as<veca1>(obj["qs1"]);
    veca1 qe1 = Rcpp::as<veca1>(obj["qe1"]);
    mata1 Omega = Rcpp::as<mata1>(obj["Omega"]);
    veca1 ce1 = Rcpp::as<veca1>(obj["ce1"]);
    veca1 PBce = Rcpp::as<veca1>(obj["PBce"]);
    return mnlink_Omega_cpp<a1type>(p1, qs1, qe1, Omega, ce1, PBce);
}


