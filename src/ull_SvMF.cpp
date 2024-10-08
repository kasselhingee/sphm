// unnormalised log-likelihood for iid SvMF data
// functions to match tapegeneral()
// results should match (11) and (12) of Scealy and Wood 2019

#include <scorematchingad_forward.h> //includes RcppEigen_forward.h
#include "uldSvMF.h"

//' @param Vvec Vectorised form of matrix V ala vech
//' @param yk vector of y values, then the k. Will be dynamic parameters because the 'y' may be residuals and be updated frequently
//' @param a1m The tuning parameter a1 and the mean, which will often be the northpole
// [[Rcpp::export]]
veca1 ull_SvMF_V(veca1 & Vvec, veca1 & yk, vecd & a1m){
    // Calculate the value of p
    int p = 1 + (-1 + sqrt(8 * Vvec.size() + 1)) / 2;

    // Reshape Vvec into a symmetric matrix V
    mata1 V = mata1::Zero(p - 1, p - 1);
    // Manually fill the lower triangular part - could get Eigen::Map to work
    int idx = 0;
    for (int i = 0; i < p - 1; ++i) {
        for (int j = 0; j <= i; ++j) {
            V(i, j) = Vvec(idx++);
        }
    }
    // Reflect the lower triangular part into the upper triangular part
    V = V.selfadjointView<Eigen::Lower>();


    // Extract y and k from yk
    mata1 y = yk.segment(0, p).transpose();
    a1type k = yk(p);

    //convert a1 to CppAD type
    a1type a1_a1type = a1m(0);
    veca1 m = a1m.segment(1, p);

    // Call uldSvMF_muV with the required arguments
    return uldSvMF_muV(y, k, m, a1_a1type, V);

}

