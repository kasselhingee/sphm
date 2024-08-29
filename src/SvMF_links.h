# ifndef SVMF_LINKS
# define SVMF_LINKS

// using the forward declarations because only the translational unit for RcppExports.cpp need to have access to the wrappers
#include <RcppEigenForward.h>
#include <scorematchingad_forward.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppEigen)]]

//' @describeIn SvMF_varlinks Aligns the columns of the Mobius-link rotation matrix `P` for the mean to the columns of G. Note that the first column of the returned G is the given mean. Returns the matrix G.
// [[Rcpp::export]]
mata1 alignedGcpp(veca1 m, mata1 P);

# endif

