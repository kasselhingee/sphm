# include "tapegeneral.h"
# include "function_map.h"

CppAD::ADFun<double> tapefun(generalfunction fun, veca1 & ind_t, veca1 & dyn_t, vecd & constvec, matd & constmat) {
  CppAD::Independent(ind_t, dyn_t);
  veca1 y;
  y = fun(ind_t, dyn_t, constvec, constmat);
  CppAD::ADFun<double> tape;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
  tape.Dependent(ind_t, y);
  return(tape);
}


Rcpp::XPtr< CppAD::ADFun<double> > tape_funptr(Rcpp::XPtr<generalfunction> funptr, veca1 & ind_t, veca1 & dyn_t, vecd & constvec, matd & constmat) {
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


Rcpp::XPtr< CppAD::ADFun<double> > tape_namedfun(std::string func_name, veca1 & ind_t, veca1 & dyn_t, vecd & constvec, matd & constmat) {
  //look up the function
  if (function_map.find(func_name) == function_map.end()) {
      Rcpp::stop("Function not found: " + func_name);
  }
  generalfunction fun = function_map[func_name];

  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer
  *out = tapefun(fun,
                 ind_t,
                 dyn_t,
                 constvec,
                 constmat);

  Rcpp::XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
}

