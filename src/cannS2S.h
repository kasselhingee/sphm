# ifndef CANNS2S_H
# define CANNS2S_H

# include <sphm_forward.h>
# include <scorematchingad_forward.h>
# include <Rcpp.h>

// Define the templated C++ struct for the cannS2S parametetrisation
template <typename T>
struct cannS2Scpp {
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> P;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Q;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> B;

    cannS2Scpp(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> P_, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Q_, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> B_)
        : P(P_), Q(Q_), B(B_) {}
};

// convert from OmegaS2S object to cannS2S object
// Convert OmegaS2Scpp to cannS2Scpp
template <typename T>
cannS2Scpp<T> Omega2cann(const OmegaS2Scpp<T>& obj) {
    //try using isfinite overloading from cppad
    CppAD::AD<double> test_var = 1.0;
    // Explicitly call the functions to force instantiation
    bool test_isfinite = std::isfinite(test_var);


    // Perform SVD on the Omega matrix
    Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> svd(obj.Omega, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Rcpp::Rcout << svd.matrixU() << std::endl;
    Rcpp::Rcout << svd.matrixV() << std::endl;
    
    // Create Q by combining q1 with the V matrix from SVD
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Q(obj.q1.size(), 1 + svd.matrixV().cols());
    Q << obj.q1, svd.matrixV();

    // Create P by combining p1 with the U matrix from SVD
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> P(obj.p1.size(), 1 + svd.matrixU().cols());
    P << obj.p1, svd.matrixU();

    // Create B as a diagonal matrix with the singular values, excluding the last one
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> B = svd.singularValues().head(svd.singularValues().size() - 1).asDiagonal();

    // Return a cannS2Scpp object with P, Q, and B
    return cannS2Scpp<T>(P, Q, B);
}

# endif
