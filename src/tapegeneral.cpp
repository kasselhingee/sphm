# include <RcppEigenForward.h>
# include <scorematchingad_forward.h>
# include "sphm_forward.h"


//' Function for taping a general function. The function must have signature
//' `veca1 fun(const veca1 & independent, const veca2 & dynamic, const vecd & constvec, const matd & constmat)`.
//' Differentiation of `fun` will occur with respect to the independent arguments. The taping will keep of dependence on the dynamic arguments so that the value of the dynamic arguments can be changed in the tape. The constants (constvec and constmat) will be baked into the tape (to change these constants `tapefun` will have to be called again.
//' 
//' If some of the arguments are unwanted, I'm hoping I can write the function so that vectors of length zero are okay.
//'
//' @param fun A function with the correct signature.
//' @param ind_t The value of the independent argument to use for taping.
//' @param dyn_t The value of the dynamic argument to use for taping.
//' @param constants The value of the constants argument.

CppAD::ADFun<double> tapefun(generalfunction fun, veca1 & ind_t, veca1 & dyn_t, const vecd & constvec, matd & constmat) {
  CppAD::Independent(ind_t, dyn_t);
  veca1 y;
  y = fun(ind_t, dyn_t, constvec, constmat);
  CppAD::ADFun<double> tape;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
  tape.Dependent(ind_t, y);
  return(tape);
}


//' Tape using a pointer to a function created by RcppXPtrUtils::cppXPtr
//' [[Rcpp::Export]]
Rcpp::XPtr< CppAD::ADFun<double> > tapefun(Rcpp::XPtr<generalfunction> funptr, veca1 & ind_t, veca1 & dyn_t, const vecd & constvec, matd & constmat) {
  generalfunction fun = *Rcpp::XPtr<generalfunction>(funptr);

  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer
  *out = tapefun(fun,
                 ind_t,
                 dyn_t,
                 constvec,
                 constmat);

  Rcpp::XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
}
