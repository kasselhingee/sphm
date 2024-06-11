# The link function for sphere to sphere.

# From Shogo's script

# stereographic projection
Sp=function(x, e1 = c(1,rep(0,length(x)-1))) {
  if (all(x==-e1)){rep(1e+9,length(x)-1)}
  else{2/vnorm(x+e1)^2*x[2:length(x)]}
}

# inverse stereographic projection
iSp=function(y) 1/(1+vnorm(y)^2)*c(1-vnorm(y)^2,2*y)

#' The mean link for spherical covariates
#' @param x a vector of covariate values
#' @param P P matrix: a p x p (?orthonormal) matrix
#' @param B B matrix: a (p-1) x (p-1) diagonal matrix with elements between zero and one ordered in decreasing size.
#' @param Q The rotation-like matrix `Q` for rotating the covariate vector `x`.
#' @examples
#' P <- Q <- diag(1, 3)
#' B <- diag(c(0.9, 0.2))
#' x <- c(0.1, 0.2, 0.3)
#' x <- x / sqrt(sum(x^2))
#' meanlinkS2S(x, P, Q, B)
#' @export
meanlinkS2S <- function(x,P,Q,B){
  stopifnot(abs(sum(x^2) - 1) < sqrt(.Machine$double.eps))
  stopifnot(max(abs(B-diag(diag(B)))) < sqrt(.Machine$double.eps))
  return(P%*%iSp(B%*%Sp(t(Q)%*%x)))
}


# want meanlinkS2S with more indenpendent paramaters: p1, q1, Omega
#' The mean link for spherical covariates using the Omega parameterisation
#' @param x a vector of covariate values
#' @param P P matrix: a p x p (?orthonormal) matrix
#' @param B B matrix: a (p-1) x (p-1) diagonal matrix with elements between zero and one ordered in decreasing size.
#' @param Q The rotation-like matrix `Q` for rotating the covariate vector `x`.
#' @export
meanlinkS2S_Omega <- function(x, p1, q1, Omega){


}


# would be nice to have the diagonal elements separate too?
