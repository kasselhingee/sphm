#include "prelimS2S.h"
#include "OmegaS2S.h"
#include "meanlinkS2S.h"
#include "tapegeneral.h"

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

veca1 OmegaS2S_constraints(veca1 & vec, int p) {
  // Convert vector to a OmegaS2Scpp object
  OmegaS2Scpp<a1type> ompar = OmegaS2Scpp_unvec(vec, p);

  // design so that function returns zero vector when constraints satisfied
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

//Constraints on the singular values of Omega - not exact unfortunately, just on total sum
veca1 OmegaS2S_ineqconstaints(veca1 & vec, veca1 & ignore1, vecd & p_in, matd & ignore2){
  int p = int(p_in(0) + 0.1);
  // Convert vector to a OmegaS2Scpp object
  OmegaS2Scpp<a1type> ompar = OmegaS2Scpp_unvec(vec, p);
  a1type ssq_sv = (ompar.Omega.transpose() * ompar.Omega).diagonal().sum();
  veca1 out(1);
  out(0) = ssq_sv - (p-1.);
  return(out);
}


