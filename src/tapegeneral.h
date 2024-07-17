# ifndef TAPEGENERAL
# define TAPEGENERAL

# include <RcppEigenForward.h>
# include <scorematchingad_forward.h>
# include "sphm_forward.h"
# include <Rcpp.h>


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
CppAD::ADFun<double> tapefun(generalfunction fun, veca1 & ind_t, veca1 & dyn_t, const vecd & constvec, matd & constmat);

//' Tape using a pointer to a function created by RcppXPtrUtils::cppXPtr
//' [[Rcpp::Export]]
Rcpp::XPtr< CppAD::ADFun<double> > tapefun(Rcpp::XPtr<generalfunction> funptr, veca1 & ind_t, veca1 & dyn_t, const vecd & constvec, matd & constmat);
# endif
