#include "SvMF_ll_alignedG.h"
#include <Rcpp.h>
#include "OmegaS2S.h"
#include "meanlinkS2S.h"
#include "ldSvMF.h"
#include "cannS2S.h"


veca1 ll_SvMF_S2S_alignedG_mean(veca1 & vec, veca1 & dyn, vecd & p_in, matd & yx){
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
  veca1 a = dyn.segment(1, p);
  mata1 P = Eigen::Map< mata1 > (dyn.segment(1+p, p * p).data(), p, p);
  veca1 ld(y.rows());
  mata1 G(p, p);
  for (int i = 0; i < y.rows(); ++i){
    G = alignedGcpp(ypred.row(i).transpose(), P);
    ld(i) = ldSvMF_cann(y.row(i), k, a, G)(0);
  }
  
  return ld;
}

veca1 ll_SvMF_S2S_alignedG_a(veca1 & vec, veca1 & dyn, vecd & pOmegavecP, matd & yx){
  int p = int(pOmegavecP(0) + 0.1); //0.1 to make sure p_in is above the integer it represents
  //
  veca1 aremaining(p-1);
  aremaining << CppAD::exp(-vec.sum()), vec.array().exp();

  //convert to parameterisation of _alignedG_mean()
  veca1 newvec = pOmegavecP.segment(1, p + (yx.cols() - p) + p * (yx.cols() - p));
  // extract P matrix
  veca1 Pvec = pOmegavecP.tail(p*p);
  veca1 newdyn(p+1+p*p);
  newdyn << dyn(0), dyn(1), aremaining, Pvec;
  vecd p_in = pOmegavecP.head(1);
  veca1 ld = ll_SvMF_S2S_alignedG_mean(newvec, newdyn, p_in, yx);
  return ld;
}


veca1 ll_SvMF_S2S_alignedG_k(veca1 & k, veca1 & dyn, vecd & p_in, matd & yx){
  int p = int(p_in(0) + 0.1); //0.1 to make sure p_in is above the integer it represents

  //convert to parameterisation of _alignedG_mean()
  veca1 newvec = dyn.head(dyn.size() - p - p*p); //should just be the Omegavec
  veca1 newdyn(p + 1 + p*p);
  newdyn << k, dyn.tail(p + p*p);

  veca1 ld = ll_SvMF_S2S_alignedG_mean(newvec, newdyn, p_in, yx);
  return ld;
}


mata1 alignedGcpp(veca1 m, mata1 P) {
    int p = m.size();
    
    // Compute m %*% t(m) which is the outer product of m with itself
    mata1 mproj = m * m.transpose();

    // first remove the mean direction from all directions in P. P 'no m'
    mata1 Pnom = (Eigen::MatrixXd::Identity(p, p) - mproj) * P;

    // Remove the mean direction from all directions in P
    for (int j = 1; j < p - 1; ++j) {
        // Normalize Pnom[, j]
        Pnom.col(j) /= Pnom.col(j).norm();

        // Update Pnom by removing the jth direction from subsequent columns
        Pnom.block(0, j + 1, p, p - j - 1) -= Pnom.col(j) * Pnom.col(j).transpose() * Pnom.block(0, j + 1, p, p - j - 1);
    }

    // Normalize the last column of Pnom since it misses out in the above loop
    Pnom.col(p - 1) /= Pnom.col(p - 1).norm();

    // Replace the first column of Pnom with m
    Pnom.col(0) = m;

    return Pnom;
}

