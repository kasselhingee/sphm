# ifndef SPHM_FORWARD
# define SPHM_FORWARD

# include <scorematchingad_forward.h> // includes <RcppEigenForward.h>

// the following is for passing around general functions for taping by tapefun
// Each object of this class *points* to function with the required signature below
typedef veca1 (*generalfunction)(const veca1 &, const veca1 &, const vecd &, const matd &);


# endif
