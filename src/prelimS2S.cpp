#include "prelimS2S.h"
#include "OmegaS2S.h"
#include "mnlink_cpp.h"
#include "tapegeneral.h"

veca1 pobjS2Scpp(veca1 & omvec, veca1 & dyn, vecd & p_in, matd & yx){
  int p = int(p_in(0) + 0.1); //0.1 to make sure p_in is above the integer it represents
  mata1 y = yx.leftCols(p);
  mata1 x = yx.block(0, p, yx.rows(), yx.cols() - p);
 
  mnlink_Omega_cpp<a1type> om = mnlink_Omega_cpp_unvec(omvec, p, 0);
  mnlink_Omega_cpp<a1type> om_projected = Omega_proj_cpp(om);
  veca1 omvec_projected;
  omvec_projected = mnlink_Omega_cpp_vec(om_projected);  

  mata1 ypred;
  ypred = mnlink_cpp(x, mata1(x.rows(), 0), omvec_projected, p);
  veca1 obj(1);
  obj(0) = -1 * (ypred.array() * y.array()).sum()/y.rows();
  return(obj);
}


