#include "S2S_constV.h"
#include "meanlinkS2S.h"
#include "uldSvMF.h"

mata1 JuppRmat(veca1 const & y1, veca1 const & y2){
  veca1 sum = y1 + y2;
  a1type denom = 1.0 + y1.dot(y2);//(y1.transpose() * y2).coeff(0,0)
  mata1 ident = mata1::Identity(y1.size(), y1.size());
  mata1 out = ((sum * sum.transpose()).array() / denom).matrix() - ident;
  return out;
}


veca1 ull_S2S_constV(mata1 y, mata1 x, OmegaS2Scpp<a1type> om, a1type k, a1type a1, veca1 aremaining, mata1 Gstar){
  int p = om.p1.size();
  //check that ncol(y) == p
  if (y.cols() != p){Rcpp::stop("width of y does not equal length of p1");}

  // project Omega matrix to satisfy orthogonality to p1 and q1
  OmegaS2Scpp<a1type> om_projected = OmegaS2Sproj(om);
  veca1 omvec_projected = OmegaS2Scpp_vec(om_projected);

  //get mean
  mata1 ypred;
  ypred = meanlinkS2Scpp(x, omvec_projected, p);

  //evaluate SvMF density of each observation
  veca1 ld(y.rows());
  mata1 G(p, p);
  G.col(0) = om_projected.p1;
  veca1 a(p);
  a(0) = a1;
  a.segment(1, p-1) = aremaining;
  for (int i = 0; i < y.rows(); ++i){
    G.block(0, 1, p, p-1) = JuppRmat(om_projected.p1, y.row(i)) * Gstar;
    ld(i) = uldSvMF_cann(y.row(i), k, a, G)(0);
  }
  return ld;
}



