#include "mobius_vMF.h"
#include "Omega.h"
#include "mobius_link_cpp.h"
#include "tapegeneral.h"

veca1 prelimobj_cpp(veca1 & omvec, veca1 & dyn, vecd & dims_in, matd & yx){
  // +0.1 guards against floating-point truncation when R passes an exact integer as a double
  int p = int(dims_in(0) + 0.1);
  int qe = int(dims_in(1) + 0.1);
  mobius_link_Omega_cpp<a1type> om = mobius_link_Omega_cpp_unvec(omvec, p, qe);

  mata1 y = yx.leftCols(p);
  mata1 xs = yx.rightCols(om.qs + om.qe).leftCols(om.qs);
  mata1 xe = yx.rightCols(om.qe);

  // Project Omega onto the Stiefel manifold before evaluating the link: the optimiser
  // satisfies Omega constraints only approximately, so projection ensures p1, qs1, qe1
  // are unit vectors and Omega is orthogonal to them.
  mobius_link_Omega_cpp<a1type> om_projected = Omega_proj_cpp(om);
  veca1 omvec_projected;
  omvec_projected = mobius_link_Omega_cpp_vec(om_projected);

  mata1 ypred;
  ypred = mobius_link_cpp(xs, xe, omvec_projected, p);
  veca1 obj;
  // obj = mu . y (dot product of predicted mean with observed y), which equals
  // (log vMF density - log vMF norm const) / k. Its gradient w.r.t. mu equals
  // grad(log vMF density) / k, making this the natural vMF surrogate objective.
  obj = (ypred.array() * y.array()).rowwise().sum();
  return(obj);
}


