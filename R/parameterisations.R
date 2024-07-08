#' Parameterisation Classes
#' @description Parameterisations of the link functions are stored as lists.
cannS2S <- function(P, Q, B){
  obj <- list(P = P, Q = Q, B = B)
  class(obj) <- c("cannS2S", class(obj))
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

OmegaS2S <- function(p1, q1, Omega){
  obj <- list(
    p1 = p1,
    q1 = q1,
    Omega = Omega
  )
  class(obj) <- c("OmegaS2S", class(obj))
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
  list2env(obj, envir = environment())
  if (check){check_cann(P, Q, B)}
  p1 <- P[, 1]
  q1 <- Q[, 1]
  Omega <- P[,-1] %*% B %*% t(Q[, -1])
  if (check){check_omega(p1, q1, Omega)}
  return(OmegaS2S(p1, q1, Omega))
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
  stopifnot(inherits(obj, "OmegaS2S"))
  list2env(obj, envir = environment())
  stopifnot(max(abs(t(p1) %*% Omega)) < sqrt(.Machine$double.eps))
  stopifnot(max(abs(Omega %*% q1)) < sqrt(.Machine$double.eps))
}
