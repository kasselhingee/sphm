# ifndef SPHM_FORWARD
# define SPHM_FORWARD

# include <scorematchingad_forward.h> // includes <RcppEigenForward.h>

// the following is for passing around general functions for taping by tapefun
// Each object of this class *points* to function with the required signature below
typedef veca1 (*generalfunction)(const veca1 &, const veca1 &, const vecd &, const matd &);

// defining a function map that will convert string names to a pointer to the function
// it will declared as extern here so that when the compiler will look elsewhere for the definition
// a definition of function_map is in the central cpp file and defined exactly ONCE: tapegeneral.cpp
// adding to the function_map can be done in any file by (e.g. ???):
// 1. including this header
// 2. declaring a static bool initialisation function that modifies function_map to have more entries
// WARNING: this creates a global variable, which has many difficulties with use
extern std::map<std::string, generalfunction> function_map;

# endif
