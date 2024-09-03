#ifndef TAPEGENERAL
#define TAPEGENERAL

# include <scorematchingad_forward.h> // includes <RcppEigenForward.h>
# include "sphm_forward.h"
# include <Rcpp.h>
# include "function_map.h"


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

//' Tape using a function name in function_map 
//' @param func_name Name of function to tape. Name must be in the internal `function_map` object.
// [[Rcpp::export]]
Rcpp::XPtr< CppAD::ADFun<double> > tape_namedfun(std::string func_name, veca1 & ind_t, veca1 & dyn_t, vecd & constvec, matd & constmat);

// tape purely within C++ using generalfunction class (which isn't exported yet)
CppAD::ADFun<double> tapefun(generalfunction fun, veca1 & ind_t, veca1 & dyn_t, vecd & constvec, matd & constmat);

//' @describeIn tape_namedfun Tape using a pointer to a function created by RcppXPtrUtils::cppXPtr
//' @param funptr A pointer to a function created by RcppXPtrUtils::cppXPtr
// [[Rcpp::export]]
Rcpp::XPtr< CppAD::ADFun<double> > tape_funptr(Rcpp::XPtr<generalfunction> funptr, veca1 & ind_t, veca1 & dyn_t, vecd & constvec, matd & constmat);


#endif
