#ifndef FUNCTION_MAP
#define FUNCTION_MAP

#include <sphm_forward.h>
#include "prelimS2S.h"
#include "Omega.h"

// defining a function map that will convert string names to a pointer to the function
std::map<std::string, generalfunction> function_map = {
   {"prelimobj_cpp", &prelimobj_cpp},
   {"Omega_constraints_wrap", &Omega_constraints_wrap},
   {"Omega_ineqconstraints", &Omega_ineqconstraints}
};


#endif
