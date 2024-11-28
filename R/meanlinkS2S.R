# The link function for sphere to sphere.
# stereographic projection

Sp=function(x, e1 = c(1,rep(0,length(x)-1))) {
  if (all(x==-e1)){rep(1e+9,length(x)-1)}
  else{1/(1+x[1])*x[2:length(x)]}
}

# inverse stereographic projection
iSp=function(y) 1/(1+vnorm(y)^2)*c(1-vnorm(y)^2,2*y)

#' The mean link for spherical covariates
#' @param x a vector of covariate values or tidy-style array of covariate values (each row a vector of covariates)
#' @param P P matrix: a p x p (?orthonormal) matrix
#' @param B B matrix: a (p-1) x (p-1) diagonal matrix with elements between zero and one ordered in decreasing size.
#' @param Q The rotation-like matrix `Q` for rotating the covariate vector `x`.
#' @examples
#' P <- Q <- diag(1, 3)
#' B <- diag(c(0.9, 0.2))
#' x <- c(0.1, 0.2, 0.3)
#' x <- x / sqrt(sum(x^2))
#' meanlinkS2S(x, P, Q, B)
#' @param paramobj A parameterisation object for either the cannonical PQB or Omega parameterisation. See [`OmegaS2S()`] and [`cannS2S()`]. Used if P, Q, and B are NULL.
#' @details If `x` is a single vector then the mean is computed directly.
#' If `x` is an array then the parameters are converted to `OmegaS2S`.
#' @useDynLib sphm, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
meanlinkS2S <- function(x,P = NULL,Q = NULL,B = NULL, paramobj = NULL, check = TRUE){
  if (is.array(x)){
    if (is.null(P) & is.null(Q) & is.null(B)){paramobj <- as_OmegaS2S(paramobj)}
    else {paramobj <- as_OmegaS2S(cannS2S(P, Q, B))}
    return(meanlinkS2S_Omega(x, paramobj, check = check))
  } else {
    if (is.null(P) & is.null(Q) & is.null(B)){paramobj <- as_mnlink_cann(paramobj)}
    else {paramobj <- cannS2S(P, Q, B)}
    if (check){cannS2S_check(paramobj)}
    return(meanlinkS2S_cann(x, paramobj))
  }
}

meanlinkS2S_cann <- function(x, paramobj){
  stopifnot(inherits(paramobj, "mnlink_cann"))
  list2env(paramobj, envir = environment())
  stopifnot(abs(sum(x^2) - 1) < sqrt(.Machine$double.eps))
  return(drop(P%*%iSp(B%*%Sp(t(Q)%*%x))))
}


# want meanlinkS2S with more indenpendent paramaters: p1, q1, Omega
#' The mean link for spherical covariates using the Omega parameterisation
#' @description Computes mean from p1, q1, Omega avoiding the svd in [`param_omega2cann()`]. It uses the expression for the link for the Omega parameterisation.
#' @param x a vector of covariate values or tidy-style array of covariate values (each row a vector of covariates)
#' @noRd
meanlinkS2S_Omega <- function(x, paramobj, check = TRUE){
  stopifnot(inherits(paramobj, "OmegaS2S"))
  if (check){OmegaS2S_check(paramobj)}
  list2env(paramobj, envir = environment())
  if (inherits(x, "array")){x <- t(x)} #so rows of x become column vectors
  
  BQx2 <- rowSums(t(x) * t(t(Omega) %*% Omega %*% x))
  q1x <- drop(t(q1) %*% x)
  p1_mult <- ( (1+q1x)^2 - BQx2 ) / ( (1+q1x)^2 + BQx2 )
  Omegax_mult <- (2 * (1+q1x)) / ( (1+q1x)^2 + BQx2 )
  
  out <- outer(p1_mult, p1) + #row i of this corresponds to a p1mult[i] %*% p1
          Omegax_mult * t(Omega %*% x) #row i of this corresponds to Omegaxmult[i] * x[i]
  if (!inherits(x, "array")){out <- drop(out)}
  return(out)
}


# would be nice to have the diagonal elements separate too?
