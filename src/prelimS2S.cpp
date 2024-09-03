#include <scorematchingad_forward.h> //includes <RcppEigenForward.h>
#include "OmegaS2S.h"
#include "meanlinkS2S.h"
#include "tapegeneral.h"

//' Preliminary Objective in the style of the `generalfunction` class:
//' @param yx is the response and covariates *cbind* together. Each row an observation.
//' @param dyn is a zero length vector
//' @param p is required to separate yx and omvec. It is passed as a double for compatiblility witn generalfunction, so will have to be rounded to a integer within the function
//' @param dyn ignored
// [[Rcpp::export]]
veca1 pobjS2Scpp(veca1 & omvec, veca1 & dyn, vecd & p_in, matd & yx){
  int p = int(p_in(0) + 0.1); //0.1 to make sure p_in is above the integer it represents
  mata1 y = yx.leftCols(p);
  mata1 x = yx.block(0, p, yx.rows(), yx.cols() - p);
 
  OmegaS2Scpp<a1type> om = OmegaS2Scpp_unvec(omvec, p);
  OmegaS2Scpp<a1type> om_projected = OmegaS2Sproj(om);
  veca1 omvec_projected;
  omvec_projected = OmegaS2Scpp_vec(om_projected);  

  mata1 ypred;
  ypred = meanlinkS2Scpp(x, omvec_projected, p);
  veca1 obj(1);
  obj(0) = -1 * (ypred.array() * y.array()).sum()/y.rows();
  return(obj);
}

//' Tape the preliminary objective
// [[Rcpp::export]]
Rcpp::XPtr< CppAD::ADFun<double> > pobjS2Stape(veca1 & omvec, vecd & p_in, matd & yx) {
  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer
  veca1 dyn(0); //empty vector for passing to general taping function
  *out = tapefun(*pobjS2Scpp, omvec, dyn, p_in, yx);
  Rcpp::XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
}


// For a parameter set return quadratic distance to constraints matching
// [[Rcpp::export]]
veca1 OmegaS2S_constraints(veca1 & vec, int p) {
  // Convert vector to a OmegaS2Scpp object
  OmegaS2Scpp<a1type> ompar = OmegaS2Scpp_unvec(vec, p);

  // design so that function returns zero vector when constraints satisfied
  int p1_size = ompar.p1.size();
  int q1_size = ompar.q1.size();
  veca1 out(1 + 1);
  out(0) = ompar.p1.squaredNorm() - 1.;
  out(1) = ompar.q1.squaredNorm() - 1.;
  return(out);
}

//a wrap around OmegaS2S_constraints for use with tapegeneral
veca1 wrap_OmegaS2S_constraints(veca1 & vec, veca1 & ignore1, vecd & p_in, matd & ignore2) {
  veca1 out;
  int p = int(p_in(0) + 0.1);
  out = OmegaS2S_constraints(vec,p);
  return(out);
}


//' Tape the constraint
// [[Rcpp::export]]
Rcpp::XPtr< CppAD::ADFun<double> > OmegaS2S_constraintstape(veca1 & omvec, vecd & p_in) {
  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer
  veca1 dyn(0); //empty vector for passing to general taping function
  matd constmat; //empty matrix for passing to general taping function
  *out = tapefun(*wrap_OmegaS2S_constraints, omvec, dyn, p_in, constmat);
  Rcpp::XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
}

