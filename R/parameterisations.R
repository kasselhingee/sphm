#' Parameterisation Classes
#' @description Parameterisations of the link functions are stored as lists.
#' @param P P matrix: a p x p (?orthonormal) matrix
#' @param B B matrix: a (p-1) x (p-1) diagonal matrix with elements between zero and one ordered in decreasing size.
#' @param Q The rotation-like matrix `Q` for rotating the covariate vector `x`.
cannS2S <- function(P, Q, B, check = TRUE){
  obj <- list(P = P, Q = Q, B = B)
  class(obj) <- c("cannS2S", class(obj))
  cannS2S_check(obj)
  return(obj)
}
as_cannS2S <- function(obj){
  if (inherits(obj, "cannS2S")){return(obj)}
  if (inherits(obj, "OmegaS2S")){return(Omega2cann(obj, check = FALSE))}
  if (!inherits(obj, "list")){stop("obj isn't a cannS2S, OmegaS2S or a list.")}
  spotifnot(all(c("P", "Q", "B") %in% names(obj)))
  class(obj) <- c("cannS2S", class(obj))
  return(obj)
}

#' @param p1 First column of the P matrix (vector of length `p`)
#' @param q1 First column of the Q matrix (vector of length `q`)
#' @param Omega A `p` by `q` matrix representing `P* B t(Q*)`.
OmegaS2S <- function(p1, q1, Omega, check = TRUE){
  obj <- list(
    p1 = p1,
    q1 = q1,
    Omega = Omega
  )
  class(obj) <- c("OmegaS2S", class(obj))
  if (check) {OmegaS2S_check(obj)}
  return(obj)
}
as_OmegaS2S <- function(obj){
  if (inherits(obj, "cannS2S")){return(cann2Omega(obj, check = FALSE))}
  if (inherits(obj, "OmegaS2S")){return(obj)}
  if (!inherits(obj, "list")){stop("obj must be either a cannS2S, OmegaS2S or list.")}
  stopifnot(all(c("p1", "q1", "Omega") %in% names(obj)))
  class(obj) <- c("OmegaS2S", class(obj))
  return(obj)
}
#' @noRd
#' Vectorsises and inverse of this vectorisation for the OmegaS2S parameterisation.
#' @param obj An OmegaS2S parameter object.
OmegaS2S_vec <- function(obj){
  stopifnot(inherits(obj, "OmegaS2S"))
  list2env(obj, envir = environment())
  out <- c(p1, q1, as.vector(Omega))
  class(out) <- "OmegaS2S_vec"
  return(out)
}

#' @noRd
#' @param vec is a vector like `OmegaS2S_vec()`
#' @param p The dimension of the response (The dimension of covariates will be infered from `p`).
OmegaS2S_unvec <- function(vec, p, check = TRUE){
  # length of vec = p + q + p*q
  # (l - p)/(1+p) = q
  q <- (length(vec) - p)/(1+p)
  stopifnot(q == as.integer(q))
  stopifnot(q > p -1)
  
  OmegaS2S(p1 = vec[1:p],
           q1 = vec[p + 1:q],
           Omega = matrix(vec[p+q+1:(p*q)], nrow = p, ncol = q, byrow = FALSE),
           check = check)
}

#' For converting between parameterisations of the link function
#' @examples
#' P <- diag(1, 3)
#' Q <- rbind(diag(1, 3), 0, 0)
#' B <- diag(c(0.9, 0.2))
#' param_cann2omega(P, Q, B)
#' @details
#' 
#' # Warning
#' Sign of columns of P and Q are lost by this transformation.
#' @export
cann2Omega <- function(obj, check = TRUE){
  if (check){cannS2S_check(obj)}
  list2env(obj, envir = environment())
  p1 <- P[, 1]
  q1 <- Q[, 1]
  Omega <- P[,-1] %*% B %*% t(Q[, -1])
  return(OmegaS2S(p1, q1, Omega, check = FALSE))
}

#' # Warning
#' Apart from p1 and q1, sign of columns of P and Q cannot be recovered from Omega.
Omega2cann <- function(obj, check = TRUE){
  list2env(obj, envir = environment())
  if (check){check_omega(p1, q1, Omega)}
  svdres <- svd(Omega, nu = nrow(Omega) - 1, nv = nrow(Omega) - 1)

  Q <- cbind(q1, svdres$v)
  P <- cbind(p1, svdres$u)
  B <- diag(svdres$d[-length(svdres$d)])
  if (check){check_cann(P, Q, B)}
  return(cannS2S(P, Q, B))
}

cannS2S_check <- function(obj){
  stopifnot(inherits(obj, "cannS2S"))
  list2env(obj, envir = environment())
  stopifnot(max(abs(B-diag(diag(B)))) < sqrt(.Machine$double.eps))
  stopifnot(max(abs(P %*% t(P) - diag(1, ncol(P)))) < sqrt(.Machine$double.eps))
  stopifnot(max(abs(t(Q) %*% Q - diag(1, ncol(Q)))) < sqrt(.Machine$double.eps))
  stopifnot(ncol(P) == ncol(Q))
  stopifnot(ncol(B) == ncol(P) - 1)
}

OmegaS2S_check <- function(obj){
  vals <- OmegaS2S_check_internal(obj)
  good <- (vals < sqrt(.Machine$double.eps))
  if (!all(good)){
    stop(paste("The following checks failed.", 
               paste0(names(vals)[!good], ": ", format(sqrt(vals[!good]), digits = 2), collapse = ", ") #sqrt here converts squared sizes to actual sizes
    ))
  }
}
OmegaS2S_check_internal <- function(obj){ #uses squared values for smoothness
  stopifnot(inherits(obj, "OmegaS2S"))
  list2env(obj, envir = environment())
  return(c(
    p1sizediff = (vnorm(p1) - 1)^2,
    q1sizediff = (vnorm(q1) - 1)^2,
    p1Omega = vnorm2(t(p1) %*% Omega),
    Omegaq1 = vnorm2(Omega %*% q1))
  )
}
