#include "SvMF_ll_alignedG_.h"
#include <Rcpp.h>
#include "OmegaS2S.h"
#include "meanlinkS2S.h"
#include "ldSvMF.h"
#include "SvMF_links.h"
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
    G = alignedG_cpp(ypred.row(i).transpose(), P);
    ld(i) = ldSvMF_cann(y.row(i), k, a, G)(0);
  }
  
  return ld;
}

veca1 ll_SvMF_S2S_alignedG_a(veca1 & vec, veca1 & dyn, vecd & pOmegavec, matd & yx){
  int p = int(pOmegavec(0) + 0.1); //0.1 to make sure p_in is above the integer it represents

  //convert to parameterisation of _alignedG_mean()
  veca1 newvec = pOmegavec.tail(pOmegavec.size() - 1);
  // extract P matrix
  mata1 P = Omega2cann(OmegaS2Scpp_unvec(newvec, p)).P;
  veca1 Pvec = Eigen::Map<veca1>(P.data(), P.size());
  veca1 newdyn(p+1+p*p);
  newdyn << dyn(0), dyn(1), vec, Pvec;
  vecd p_in = pOmegavec.segment(0,1);
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

