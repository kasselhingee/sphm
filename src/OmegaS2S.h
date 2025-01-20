#ifndef OMEGAS2S_H
#define OMEGAS2S_H

# include <sphm_forward.h>
# include <scorematchingad_forward.h> // includes <RcppEigenForward.h>
# include <RcppEigen.h> //for access to Rcpp::as and Rcpp::wrap for Eigen objects
# include <Rcpp.h>

// [[Rcpp::depends(RcppEigen)]]


// Define the templated C++ struct for the mnlink_Omega parametetrisation
template <typename T>
struct mnlink_Omega_cpp {
    Eigen::Matrix<T, Eigen::Dynamic, 1> p1;
    Eigen::Matrix<T, Eigen::Dynamic, 1> qs1; //uninitialised these vectors have 0 length
    Eigen::Matrix<T, Eigen::Dynamic, 1> qe1; //uninitialised these vectors have 0 length
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Omega;
    Eigen::Matrix<T, Eigen::Dynamic, 1> ce1;  //uninitialised these vectors have 0 length
    Eigen::Matrix<T, Eigen::Dynamic, 1> PBce;  //uninitialised these vectors have 0 length
    int p = 0;
    int qs = 0;
    int qe = 0;

    mnlink_Omega_cpp(Eigen::Matrix<T, Eigen::Dynamic, 1> p1_,
                     Eigen::Matrix<T, Eigen::Dynamic, 1> qs1_,
                     Eigen::Matrix<T, Eigen::Dynamic, 1> qe1_,
                     Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Omega_, 
                     Eigen::Matrix<T, Eigen::Dynamic, 1> ce1_,
                     Eigen::Matrix<T, Eigen::Dynamic, 1> PBce_) 
        : p1(p1_), 
          qs1(qs1_), 
          qe1(qe1_), 
          Omega(Omega_), 
          ce1(ce1_), 
          PBce(PBce_), 
          p(p1_.size()),
          qs(qs1_.size()), 
          qe(qe1_.size()) {
        if (qe == 0 && (ce1.size() + PBce.size()) > 0) {
            Rcpp::stop("ce1 and PBce must be empty when qe is 0");
        }
    }
};


// Function to vectorize an mnlink_Omega_cpp object - to match mnlink_Omega_vec
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> mnlink_Omega_cpp_vec(const mnlink_Omega_cpp<T>& obj) {
    int p1_size = obj.p1.size();
    int Omega_size = obj.Omega.size();
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Omega = obj.Omega;

    //vectorise
    Eigen::Matrix<T, Eigen::Dynamic, 1> out(obj.p + obj.qs + obj.qe + obj.p * (obj.qs + obj.qe) + obj.ce1.size() + obj.PBce.size());
    out << obj.p1, obj.qs1, obj.qe1, Eigen::Map< Eigen::Matrix<T, Eigen::Dynamic, 1> >(Omega.data(), Omega.size()), obj.ce1, obj.PBce;

    return out;
}

// Function to unvectorize into an mnlink_Omega_cpp object
template <typename T>
mnlink_Omega_cpp<T> mnlink_Omega_cpp_unvec(const Eigen::Matrix<T, Eigen::Dynamic, 1>& vec, const int p, const int qe = 0) {
    int qs = (vec.size() - p - (p + 1) * (qe > 0) - qe - p * qe) / (1 + p);
   
    return mnlink_Omega_cpp<T>(vec.segment(0, p),
                        vec.segment(p, qs),
                        vec.segment(p + qs, qe),
                        Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >(vec.segment(p + qs + qe, p * (qs + qe)).data(), p, qs + qe),
                        vec.segment(p + qs + qe + p * (qs + qe), (qe>0)),
                        vec.segment(p + qs + qe + p * (qs + qe) + (qe>0), p * (qe>0))
                        );
}

// Function to project the Omega in an OmegaS2S object to be perpendicular to p1 and q1
template <typename T>
mnlink_Omega_cpp<T> Omega_proj_cpp(const mnlink_Omega_cpp<T>& inobj) {
    mnlink_Omega_cpp<T> obj = inobj;
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
        obj.PBce = obj.PBce - (obj.p1 * obj.p1.transpose()) * obj.PBce;
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


#endif
