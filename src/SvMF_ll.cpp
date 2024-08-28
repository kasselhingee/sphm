#include <RcppEigenForward.h>
#include <scorematchingad_forward.h>
#include <Rcpp.h>
#include "OmegaS2S.h"
#include "meanlinkS2S.h"
#include "ldSvMF.h"


//' The log-likelihood of a SvMF Sphere-Sphere Regression with Mobius Mean Link and Variance Axes Aligned with P.
//' @description Three functions in a format compatible with `tapefun`. A difficulty is that the P matrix needs and SVD to get out of the Omega parameterisation so two functions for alternating between optimising the mean (with `kappa` and `a` fixed, and `G` aligned to a fixed `P`) and optimising a, with the mean fixed and concentration fixed, and updating `P` inbetween.
//' @param vec For `_aligned_mean`: A parameter vector specifying the mean via the Omega vectorisation.
//' @param dyn For `_aligned_mean`: A p+1+p*p length vector of kappa, then a1, a2, ..., then P as a vector of stacked columns.
//' @param p_in The dimension p
//' @param yx The observations and covariates cbind together as row vectors
// [[Rcpp::export]]
veca1 ll_SvMF_S2S_aligned_mean(veca1 & vec, veca1 & dyn, vecd & p_in, matd & yx){
  int p = int(p_in(0) + 0.1); //0.1 to make sure p_in is above the integer it represents

  // separate the response the covariates
  mata1 y = yx.leftCols(p);
  mata1 x = yx.rightCols(yx.cols() - p);

  // extract the Omega vector for the mean link, and project it to satisfy p1, q1 orthogonality constraints
  veca1 omvec = vec.block(0,0, p + x.cols() + p*x.cols(), 1);
  OmegaS2Scpp<a1type> om = OmegaS2Scpp_unvec(omvec, p);
  OmegaS2Scpp<a1type> om_projected = OmegaS2Sproj(om);
  veca1 omvec_projected;
  omvec_projected = OmegaS2Scpp_vec(om_projected);

  //get mean
  mata1 ypred;
  ypred = meanlinkS2Scpp(x, omvec_projected, p);

  //evaluate SvMF density using ldSvMF_cann for each row.
  //first get G without its first column and other constant parameters
  a1type k = dyn(0);
  Rcpp::Rcout << k << std::endl;
  veca1 a = dyn.segment(1, p);
  Rcpp::Rcout << a << std::endl;
  mata1 P = Eigen::Map< mata1 > (dyn.segment(1+p, p * p).data(), p, p);
  Rcpp::Rcout << P << std::endl;
  veca1 ld(y.rows());
  mata1 G(p, p);
  for (int i = 0; i < y.rows(); ++i){
    G = alignedGcpp(ypred.row(i).transpose(), P);
    ld(i) = ldSvMF_cann(y.row(i), k, a, G)(0);
  }
  
  return ld;
}

