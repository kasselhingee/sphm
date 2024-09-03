#ifndef FUNCTION_MAP
#define FUNCTION_MAP

#include <sphm_forward.h>
#include "prelimS2S.h"
#include "SvMF_ll_alignedG.h"

// defining a function map that will convert string names to a pointer to the function
std::map<std::string, generalfunction> function_map = {
   {"pobjS2Scpp", &pobjS2Scpp},
   {"wrap_OmegaS2S_constraints", &wrap_OmegaS2S_constraints},
   {"ll_SvMF_S2S_alignedG_mean", &ll_SvMF_S2S_alignedG_mean},
   {"ll_SvMF_S2S_alignedG_a", &ll_SvMF_S2S_alignedG_a},
   {"ll_SvMF_S2S_alignedG_k", &ll_SvMF_S2S_alignedG_k},
};


#endif
