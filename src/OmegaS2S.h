#ifndef OMEGAS2S_H
#define OMEGAS2S_H

# include <sphm_forward.h>
# include <scorematchingad_forward.h> // includes <RcppEigenForward.h>
# include <RcppEigen.h> //for access to Rcpp::as and Rcpp::wrap for Eigen objects
# include <Rcpp.h>

// [[Rcpp::depends(RcppEigen)]]


// Define the templated C++ struct for the OmegaS2S parametetrisation
template <typename T>
struct OmegaS2Scpp {
    Eigen::Matrix<T, Eigen::Dynamic, 1> p1;
    Eigen::Matrix<T, Eigen::Dynamic, 1> q1;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Omega;

    OmegaS2Scpp(Eigen::Matrix<T, Eigen::Dynamic, 1> p1_, Eigen::Matrix<T, Eigen::Dynamic, 1> q1_, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Omega_) 
        : p1(p1_), q1(q1_), Omega(Omega_) {}
};


// Function to vectorize an OmegaS2Scpp object
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> OmegaS2Scpp_vec(const OmegaS2Scpp<T>& obj) {
    int p1_size = obj.p1.size();
    int q1_size = obj.q1.size();
    int Omega_size = obj.Omega.size();

    // explicitly copy the parameters - need to do this because Eigen will be doing passing of pointers if it can

    Eigen::Matrix<T, Eigen::Dynamic, 1> p1 = obj.p1;
    Eigen::Matrix<T, Eigen::Dynamic, 1> q1 = obj.q1;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Omega = obj.Omega;

    //vectorise
    Eigen::Matrix<T, Eigen::Dynamic, 1> out(p1_size + q1_size + Omega_size);
    out << p1, q1, Eigen::Map< Eigen::Matrix<T, Eigen::Dynamic, 1> >(Omega.data(), Omega_size);

    return out;
}

// Function to unvectorize into an OmegaS2Scpp object
template <typename T>
OmegaS2Scpp<T> OmegaS2Scpp_unvec(const Eigen::Matrix<T, Eigen::Dynamic, 1>& vec, const int p) {
    int total_length = vec.size();
    int q = (total_length - p) / (1 + p);
    
    if (q <= p - 1) {
      Rcpp::stop("q must be greater than p - 1");
    }
    
    Eigen::Matrix<T, Eigen::Dynamic, 1> p1 = vec.segment(0, p);
    Eigen::Matrix<T, Eigen::Dynamic, 1> q1 = vec.segment(p, q);
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Omega = Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >(vec.segment(p + q, p * q).data(), p, q);

    return OmegaS2Scpp<T>(p1, q1, Omega);
}

// Function to project the Omega in an OmegaS2S object to be perpendicular to p1 and q1
template <typename T>
OmegaS2Scpp<T> OmegaS2Sproj(const OmegaS2Scpp<T>& obj) {
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Omegaperpp1;
    Omegaperpp1 = obj.Omega - obj.p1 * (obj.p1.transpose()) * obj.Omega; //remove p1 component
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Omegaperpq1;
    Omegaperpq1 = Omegaperpp1 - Omegaperpp1 * obj.q1 * (obj.q1.transpose()); //remove q1 component
    return OmegaS2Scpp<T>(obj.p1, obj.q1, Omegaperpq1);
}

// Convert an R list to an OmegaS2Scpp struct
OmegaS2Scpp<double> R2OmegaS2Scpp_double(Rcpp::List obj);

OmegaS2Scpp<a1type> R2OmegaS2Scpp_a1type(Rcpp::List obj);


#endif
