# include "tapegeneral.h"

CppAD::ADFun<double> tapefun(generalfunction fun, veca1 & ind_t, veca1 & dyn_t, vecd & constvec, matd & constmat) {
  CppAD::Independent(ind_t, dyn_t);
  veca1 y;
  y = fun(ind_t, dyn_t, constvec, constmat);
  CppAD::ADFun<double> tape;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
  tape.Dependent(ind_t, y);
  return(tape);
}


Rcpp::XPtr< CppAD::ADFun<double> > tapefun(Rcpp::XPtr<generalfunction> funptr, veca1 & ind_t, veca1 & dyn_t, vecd & constvec, matd & constmat) {
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


//function_map definition, the ONLY definition of this global variable
//see sphm_forward.h for more explanation on use
std::map<std::string, generalfunction> function_map;


