#ifndef FUNCTION_MAP
#define FUNCTION_MAP

#include <sphm_forward.h>
#include "mobius_vMF.h"
#include "Omega.h"

// Maps string function names to C++ function pointers for use by tape_namedfun() (tapegeneral.h).
// R-side code passes a function name as a string at runtime; tape_namedfun looks it up here
// to select the correct objective/constraint function to tape with CppAD.
// All entries must match the generalfunction signature: veca1(veca1&, veca1&, vecd&, matd&).
// Currently only used for the vMF optimisation; the SvMF objective uses a bespoke tape
// (tape_ld_mobius_SvMF_partransport_nota1) that bypasses this map.
// To add a new entry: implement the function in its own .h/.cpp, include the header here,
// and add a {"name", &function} line below.
std::map<std::string, generalfunction> function_map = {
   {"prelimobj_cpp", &prelimobj_cpp},
   {"Omega_constraints_wrap", &Omega_constraints_wrap},
   {"Omega_ineqconstraints", &Omega_ineqconstraints}
};


#endif
