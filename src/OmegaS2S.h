# ifndef OMEGAS2S_H
# define OMEGAS2S_H

# include <sphm_forward.h>
# include <Rcpp.h>

// Define the templated C++ struct for the OmegaS2S parametetrisation
template <typename T>
struct OmegaS2Scpp {
    Eigen::Matrix<T, Eigen::Dynamic, 1> p1;
    Eigen::Matrix<T, Eigen::Dynamic, 1> q1;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Omega;

    OmegaS2S(Eigen::Matrix<T, Eigen::Dynamic, 1> p1_, Eigen::Matrix<T, Eigen::Dynamic, 1> q1_, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Omega_) 
        : p1(p1_), q1(q1_), Omega(Omega_) {}
};


// Function to vectorize an OmegaS2Scpp object
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> OmegaS2Scpp_vec(OmegaS2Scpp<T>& obj) {
    int p1_size = obj.p1.size();
    int q1_size = obj.q1.size();
    int Omega_size = obj.Omega.size();

    Eigen::Matrix<T, Eigen::Dynamic, 1> out(p1_size + q1_size + Omega_size);
    out << obj.p1, obj.q1, Map< Eigen::Matrix<T, Eigen::Dynamic, 1> >(obj.Omega.data(), Omega_size);

    return out;
}

// Function to unvectorize into an OmegaS2Scpp object
template <typename T>
OmegaS2Scpp<T> OmegaS2Scpp_unvec(Eigen::Matrix<T, Eigen::Dynamic, 1>, int p) {
    int total_length = vec.size();
    int q = (total_length - p) / (1 + p);
    
    if (q <= p - 1) {
      Rcpp::stop("q must be greater than p - 1");
    }
    
    Eigen::Matrix<T, Dynamic, 1> p1 = vec.segment(0, p);
    Eigen::Matrix<T, Dynamic, 1> q1 = vec.segment(p, q);
    Eigen::Matrix<T, Dynamic, Dynamic> Omega(p, q);
    Omega = Eigen::Map< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >(vec.segment(p + q, p * q).data(), p, q);

    return OmegaS2Scpp<T>(p1, q1, Omega);
}

// Convert an R list to an OmegaS2Scpp struct
// [[Rcpp::export]]
OmegaS2Scpp<double> R2OmegaS2Scpp_double(List obj);

// [[Rcpp::export]]
OmegaS2Scpp<a1type> R2OmegaS2Scpp_a1type(List obj);


# endif
