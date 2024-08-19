#' SvMF Parameterisations
#' The SvMF canonical parameterisation is a set of positive real numbers (kappa, `a`) and a matrix of orthonormal column vectors `Gamma`, where the product of the `a` is 1. `a_1` and `kappa` are not jointly estimatable by MLE, so I'm expecting them to be closely related.
#' There are three other parametrisations.
#' Alternatively this set of numbers are the eigenvalues and eigenvectors of the symmetric matrix E(yy^T). The determinant of this matrix would be kappa * a_1.
#' A third parametrisation is available: using kappa, mu and a matrix V that is symmetric positive definite with det(V) = 1.
#' A fourth represents the above V matrix as an orthogonal matrix and a set of positive real values.

SvMFcann <- function(kappa, a, Gamma){
  obj <- list(k = kappa, a = a, G = Gamma)
  class(obj) <- c("SvMFcann", class(obj))
  return(obj)
}
SvMFcann_check <- function(obj){
  list2env(obj, envir = environment())
  checks <- c(
  `a length` = length(a) == ncol(G),
  `Gamma square` = ncol(G) == nrow(G),
  `kappa positive` = k > 0,
  `a positive` = a > 0,
  `a product` = (prod(a[-1]) - 1)^2 < sqrt(.Matchine$double.eps),
  `Gamma orthogonal` = (t(G) %*% G - diag(1, ncol(G)))^2 <  sqrt(.Matchine$double.eps)
  )
  if (any(!checks)){
    stop(sprinf("Parameters fail SvMF checks: %s", paste(names(checks)[!checks], collapse = ", ")))
  }
  return(NULL)
}


