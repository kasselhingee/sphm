#ifndef OMEGA_H
#define OMEGA_H

# include <sphm_forward.h>
# include <scorematchingad_forward.h> // includes <RcppEigenForward.h>
# include <RcppEigen.h> //for access to Rcpp::as and Rcpp::wrap for Eigen objects
# include <Rcpp.h>

// [[Rcpp::depends(RcppEigen)]]


// Möbius mean-link parameters in the Omega parameterisation (C++ counterpart of
// mobius_link_Omega in R). Fields:
//   p1    — base point on S^{p-1}: the mean direction when covariates are zero.
//   qs1   — "north-pole" direction in spherical covariate space: maps to the origin
//            of the projected space; the stereographic projection runs from -qs1.
//            Must satisfy qs1^T Omega_s = 0 and ||qs1|| = 1. Empty if no spherical covariate.
//   qe1   — "north-pole" reference direction in Euclidean covariate space: plays the
//            same role as qs1 in the Euclidean denominator term; singular at -qe1.
//            Must satisfy qe1^T Omega_e = 0 and ||qe1|| = 1. Empty if no Euclidean covariate.
//   Omega — p x (qs+qe) coefficient matrix, orthogonal to p1: p1^T Omega = 0.
//   ce    — scalar denominator offset for the Euclidean term (length 1 when qe>0, else empty).
// The flat vector layout used by mobius_link_Omega_cpp_vec is:
//   [p1 (p) | qs1 (qs) | qe1 (qe) | vec(Omega) (p*(qs+qe)) | ce (qe>0)]
template <typename T>
struct mobius_link_Omega_cpp {
    Eigen::Matrix<T, Eigen::Dynamic, 1> p1;
    Eigen::Matrix<T, Eigen::Dynamic, 1> qs1; //uninitialised these vectors have 0 length
    Eigen::Matrix<T, Eigen::Dynamic, 1> qe1; //uninitialised these vectors have 0 length
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Omega;
    Eigen::Matrix<T, Eigen::Dynamic, 1> ce;  //uninitialised these vectors have 0 length
    int p = 0;
    int qs = 0;
    int qe = 0;

    mobius_link_Omega_cpp(Eigen::Matrix<T, Eigen::Dynamic, 1> p1_,
                     Eigen::Matrix<T, Eigen::Dynamic, 1> qs1_,
                     Eigen::Matrix<T, Eigen::Dynamic, 1> qe1_,
                     Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Omega_, 
                     Eigen::Matrix<T, Eigen::Dynamic, 1> ce_) 
        : p1(p1_), 
          qs1(qs1_), 
          qe1(qe1_), 
          Omega(Omega_), 
          ce(ce_), 
          p(p1_.size()),
          qs(qs1_.size()), 
          qe(qe1_.size()) {
        if (qe == 0 && (ce.size() > 0)) {
            Rcpp::stop("ce must be empty when qe is 0");
        }
    }
};


// Function to vectorize a mobius_link_Omega_cpp object - to match mobius_link_Omega_vec
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> mobius_link_Omega_cpp_vec(const mobius_link_Omega_cpp<T>& obj) {
    int p1_size = obj.p1.size();
    int Omega_size = obj.Omega.size();
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Omega = obj.Omega;

    //vectorise
    Eigen::Matrix<T, Eigen::Dynamic, 1> out(obj.p + obj.qs + obj.qe + obj.p * (obj.qs + obj.qe) + obj.ce.size());
    out << obj.p1, obj.qs1, obj.qe1, Eigen::Map< Eigen::Matrix<T, Eigen::Dynamic, 1> >(Omega.data(), Omega.size()), obj.ce;

    return out;
}

// Function to unvectorize into a mobius_link_Omega_cpp object.
// The flat vector layout is [p1 (p) | qs1 (qs) | qe1 (qe) | vec(Omega) (p*(qs+qe)) | ce (qe>0)].
// Given total length = p + qs + qe + p*(qs+qe) + (qe>0), solve for qs:
//   qs = (vec.size() - p - (qe>0) - qe - p*qe) / (1 + p)
template <typename T>
mobius_link_Omega_cpp<T> mobius_link_Omega_cpp_unvec(const Eigen::Matrix<T, Eigen::Dynamic, 1>& vec, const int p, const int qe = 0) {
    int qs = (vec.size() - p - (qe > 0) - qe - p * qe) / (1 + p);

    return mobius_link_Omega_cpp<T>(vec.segment(0, p), //p1
                        vec.segment(p, qs), //qs1
                        vec.segment(p + qs, qe), //qe1
                        Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >(vec.segment(p + qs + qe, p * (qs + qe)).data(), p, qs + qe), //Omega
                        vec.segment(p + qs + qe + p * (qs + qe), (qe>0)) //ce
                        );
}

// Function to project the Omega in an Omega parameterisation to be perpendicular to p1 and qs1 and qe1
template <typename T>
mobius_link_Omega_cpp<T> Omega_proj_cpp(const mobius_link_Omega_cpp<T>& inobj) {
    mobius_link_Omega_cpp<T> obj = inobj;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> newOmega;

    // First project orthogonal to p1 (needs p1 as a unit vector)
    obj.p1 = obj.p1 / obj.p1.norm();
    newOmega = obj.Omega - (obj.p1 * obj.p1.transpose()) * obj.Omega;

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Omega_s, Omega_e;

    // Project Omega_s perpendicular to qs1 if qs1 is not empty
    if (obj.qs > 0) {
        obj.qs1 = obj.qs1 / obj.qs1.norm();
        Omega_s = newOmega.leftCols(obj.qs);
        Omega_s = Omega_s - Omega_s * obj.qs1 * obj.qs1.transpose();
    }

    // Project Omega_e perpendicular to qe1 if qe1 is not empty
    if (obj.qe > 0) {
        obj.qe1 = obj.qe1 / obj.qe1.norm();
        Omega_e = newOmega.rightCols(obj.qe);
        Omega_e = Omega_e - Omega_e * obj.qe1 * obj.qe1.transpose();
    }

    // Combine Omega_s and Omega_e into obj.Omega
    if (obj.qs > 0 && obj.qe > 0) {
        newOmega << Omega_s, Omega_e;
    } else if (obj.qs > 0) {
        newOmega = Omega_s;
    } else if (obj.qe > 0) {
        newOmega = Omega_e;
    }
    obj.Omega = newOmega;

    return obj;
}

// For a parameter set return quadratic distance to constraints matching
// [[Rcpp::export]]
veca1 Omega_constraints(veca1 & vec, int p, int qe=0);

//a wrap around Omega_constraints for use with tapegeneral
veca1 Omega_constraints_wrap(veca1 & vec, veca1 & ignore1, vecd & dims_in, matd & ignore2);

//Constraints on the singular values of Omega - not exact unfortunately, just on total sum
veca1 Omega_ineqconstraints(veca1 & vec, veca1 & ignore1, vecd & dims_in, matd & ignore2);

int Omega_veclength(int p, int qs, int qe);

#endif
