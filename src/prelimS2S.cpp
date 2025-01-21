#include "prelimS2S.h"
#include "Omega.h"
#include "mnlink_cpp.h"
#include "tapegeneral.h"

veca1 prelimobj_cpp(veca1 & omvec, veca1 & dyn, vecd & dims_in, matd & yx){
  int p = int(dims_in(0) + 0.1); //0.1 to make sure p_in is above the integer it represents
  int qe = int(dims_in(1) + 0.1); //0.1 to make sure p_in is above the integer it represents
  mnlink_Omega_cpp<a1type> om = mnlink_Omega_cpp_unvec(omvec, p, 0);

  mata1 y = yx.leftCols(p);
  mata1 xs = yx.rightCols(om.qs + om.qe).leftCols(om.qs);
  mata1 xe = yx.rightCols(om.qe);
 
  mnlink_Omega_cpp<a1type> om_projected = Omega_proj_cpp(om);
  veca1 omvec_projected;
  omvec_projected = mnlink_Omega_cpp_vec(om_projected);  

  mata1 ypred;
  ypred = mnlink_cpp(xs, xe, omvec_projected, p, qe);
  veca1 obj(1);
  obj(0) = -1 * (ypred.array() * y.array()).sum()/y.rows();
  return(obj);
}


