#ifndef FUNCTION_MAP
#define FUNCTION_MAP

#include <sphm_forward.h>
#include "prelimS2S.h"
#include "Omega.h"
#include "S2S_alignedG.h"

// defining a function map that will convert string names to a pointer to the function
std::map<std::string, generalfunction> function_map = {
   {"prelimobj_cpp", &prelimobj_cpp},
   {"Omega_constraints_wrap", &Omega_constraints_wrap},
   {"Omega_ineqconstraints", &Omega_ineqconstraints},
   {"ull_S2S_alignedG_mean", &ull_S2S_alignedG_mean},
   {"ull_S2S_alignedG_a", &ull_S2S_alignedG_a},
   {"ull_S2S_alignedG_k", &ull_S2S_alignedG_k},
};


#endif
