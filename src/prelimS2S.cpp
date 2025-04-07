#include "prelimS2S.h"
#include "Omega.h"
#include "mnlink_cpp.h"
#include "tapegeneral.h"

veca1 prelimobj_cpp(veca1 & omvec, veca1 & yx, vecd & dims_in, matd & ignore){
  int p = int(dims_in(0) + 0.1); //0.1 to make sure p_in is above the integer it represents
  int qe = int(dims_in(1) + 0.1); //0.1 to make sure p_in is above the integer it represents
  mnlink_Omega_cpp<a1type> om = mnlink_Omega_cpp_unvec(omvec, p, qe);
  
  veca1 y = yx.head(p);
  veca1 xs = yx.tail(om.qs + om.qe).head(om.qs);
  veca1 xe = yx.tail(om.qe);
 
  mnlink_Omega_cpp<a1type> om_projected = Omega_proj_cpp(om);
  veca1 omvec_projected;
  omvec_projected = mnlink_Omega_cpp_vec(om_projected);  

  mata1 ypred;
  ypred = mnlink_cpp(xs.transpose(), xe.transpose(), omvec_projected, p);
  veca1 obj(1);
  obj(0) = -1 * (ypred.array() * y.transpose().array()).sum();
  return(obj);
}


