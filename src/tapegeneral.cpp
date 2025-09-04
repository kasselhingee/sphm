# include "tapegeneral.h"
# include "function_map.h"
# include "ldSvMF.h" // purely for tape_besselImixed() below

CppAD::ADFun<double> tapefun(generalfunction fun, veca1 & ind_t, veca1 & dyn_t, vecd & constvec, matd & constmat, bool check_for_nan) {
  CppAD::Independent(ind_t, dyn_t);
  veca1 y;
  y = fun(ind_t, dyn_t, constvec, constmat);
  CppAD::ADFun<double> tape;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
  tape.Dependent(ind_t, y);
  tape.check_for_nan(check_for_nan);
  return(tape);
}


pADFun tape_namedfun(std::string func_name, veca1 & ind_t, veca1 & dyn_t, vecd & constvec, matd & constmat, bool check_for_nan) {
  //look up the function
  if (function_map.find(func_name) == function_map.end()) {
      Rcpp::stop("Function not found: " + func_name);
  }
  generalfunction fun = function_map[func_name];

  CppAD::ADFun<double> tape;
  tape = tapefun(fun,
                 ind_t,
                 dyn_t,
                 constvec,
                 constmat,
                 check_for_nan);

  pADFun out(tape, ind_t, dyn_t, func_name);
  return(out);
}

//' Function to create tapes of besselImixed() from ldSvMF purely for testing differentiation
// [[Rcpp::export]]
pADFun tape_besselImixed(veca1 & x, const double & nu, double threshold, int order, bool log_result = true) {
  CppAD::Independent(x);
  veca1 y(1);
  y(0) = besselImixed(x(0), nu, threshold, order, log_result);
  CppAD::ADFun<double> tape;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
  tape.Dependent(x, y);
  tape.check_for_nan(false);
  veca1 dyn_t(0);
  pADFun out(tape, x, dyn_t, "besselImixed");
  return(out);
}
 