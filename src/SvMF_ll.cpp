#include <RcppEigenForward.h>
#include <scorematchingad_forward.h>
#include <Rcpp.h>
#include "OmegaS2S.h"
#include "meanlinkS2S.h"

//' The log-likelihood of a SvMF Sphere-Sphere Regression with Mobius Mean Link and Variance Axes Aligned with P.
//' @description In a format compatible with `tapefun`. A huge difficulty is that the P matrix needs and SVD to get out of the Omega parameterisation.
//' @param vec A parameter vector. The first p + q + p*q parameters specify the mean via the Omega vectorisation. The remaining p-1 parameters are the concentration k then the vector a2, a3...
//' @param dyn The a1 parameter 
veca1 ll_SvMF_S2S_aligned(veca1 & vec, veca1 & a1, vecd & p_in, matd & yx){
  int p = int(p_in(0) + 0.1); //0.1 to make sure p_in is above the integer it represents

  // separate the response the covariates
  mata1 y = yx.leftCols(p);
  mata1 x = yx.block(0, p, yx.rows(), yx.cols() - p);

  // extract the Omega vector for the mean link, and project it to satisfy p1, q1 orthogonality constraints
  veca1 omvec = vec.block(0,0, p + q.rows() + p*q, 1);
  OmegaS2Scpp<a1type> om = OmegaS2Scpp_unvec(omvec, p);
  OmegaS2Scpp<a1type> om_projected = OmegaS2Sproj(om);
  veca1 omvec_projected;
  omvec_projected = OmegaS2Scpp_vec(om_projected);

  //get mean
  mata1 ypred;
  ypred = meanlinkS2Scpp(x, omvec_projected, p);

  //get the P matrix from the omvec_projected THIS ISNT DIFFERENTIABLE STABLY
  //the order of the a matters here too now 

  //for each row of covariates, get G assuming aligned association with P of mean link
  //and apply ldSvMF_cann().
  for (int i = 0; i < x.rows() - 1; ++i){
    alignedGcpp(ypred.row(i), 

  }

}

