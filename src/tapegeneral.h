#ifndef TAPEGENERAL
#define TAPEGENERAL

# include <scorematchingad_forward.h> // includes <RcppEigenForward.h>
# include <utils/pADFun.h>
# include "sphm_forward.h"
# include <Rcpp.h>


//' @noRd
//' @title Create a CppAD automatic differentiation tape of a named C++ function
//' @description
//' Records the computation of a named C++ function using CppAD operator overloading, producing
//' a `pADFun` tape object that supports fast evaluation and exact differentiation.
//' The returned tape supports `$forward(0, x)` (evaluate), `$Jacobian(x)`, and `$Hessian0(x)`.
//'
//' The `func_name` must be a key in the internal `function_map` (see `src/function_map.h`):
//' - `"prelimobj_cpp"` — vMF preliminary objective
//' - `"Omega_constraints_wrap"` — Omega orthonormality equality constraints
//' - `"Omega_ineqconstraints"` — Omega inequality constraints
//'
//' Eventually this function is expected to be replaced by bespoke tape functions (as was done
//' for the SvMF objective via `tape_ld_mobius_SvMF_partransport_nota1`).
//'
//' @param func_name Name of the C++ function to tape.
//' @param ind_t Independent variables at their taping values. Differentiation is with respect
//'   to these.
//' @param dyn_t Dynamic variables at their taping values. These can be updated after taping
//'   without re-taping (use `$new_dynamic()`).
//' @param constvec Constant vector baked into the tape. Re-taping is needed to change these.
//' @param constmat Constant matrix baked into the tape. Re-taping is needed to change these.
//' @param check_for_nan If `TRUE`, the tape detects NaN values during evaluation (useful for
//'   debugging but slower).
// [[Rcpp::export]]
pADFun tape_namedfun(std::string func_name, veca1 & ind_t, veca1 & dyn_t, vecd & constvec, matd & constmat, bool check_for_nan);

// tape purely within C++ using generalfunction class (which isn't exported yet)
CppAD::ADFun<double> tapefun(generalfunction fun, veca1 & ind_t, veca1 & dyn_t, vecd & constvec, matd & constmat, bool check_for_nan);


#endif
