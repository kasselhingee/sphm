#ifndef FUNCTION_MAP
#define FUNCTION_MAP

#include <sphm_forward.h>
#include "prelimS2S.h"
#include "S2S_alignedG.h"

// defining a function map that will convert string names to a pointer to the function
std::map<std::string, generalfunction> function_map = {
   {"pobjS2Scpp", &pobjS2Scpp},
   {"Omega_constraints_wrap", &Omega_constraints_wrap},
   {"OmegaS2S_ineqconstaints", &OmegaS2S_ineqconstaints},
   {"ull_S2S_alignedG_mean", &ull_S2S_alignedG_mean},
   {"ull_S2S_alignedG_a", &ull_S2S_alignedG_a},
   {"ull_S2S_alignedG_k", &ull_S2S_alignedG_k},
};


#endif
