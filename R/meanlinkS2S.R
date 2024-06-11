# The link function for sphere to sphere.

# From Shogo's script


vnorm=function(x) sqrt(sum(x^2))

# stereographic projection
Sp=function(x) {
  if (all(x==-e1)){rep(1e+9,d-1)}
  else{2/vnorm(x+e1)^2*x[2:d]}
}

# inverse stereographic projection
iSp=function(y) 1/(1+vnorm(y)^2)*c(1-vnorm(y)^2,2*y[1:(d-1)])

#' The mean link for spherical covariates
#' @param x a vector of covariate values
#' @param P P matrix: a p x p (?orthonormal) matrix
#' @param B B matrix: a (p-1) x (p-1) diagonal matrix with elements between zero and one ordered in decreasing size.
#' @param Q The rotation-like matrix `Q` for rotating the covariate vector `x`.
#' @export
meanlinkS2S <- function(x,P,Q,B) P%*%iSp(B%*%Sp(t(Q)%*%x))

mu <- meanlinkS2S

# want meanlinkS2S with more indenpendent paramaters: p1, q1, Omega
# would be nice to have the diagonal elements separate too?
