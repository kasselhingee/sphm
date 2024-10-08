# include "tapegeneral.h"
# include "function_map.h"

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

