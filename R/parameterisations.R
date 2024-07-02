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
param_cann2omega <- function(P,Q,B, check = TRUE){
  if (check){check_cann(P, Q, B)}
  p1 <- P[, 1]
  q1 <- Q[, 1]
  Omega <- P[,-1] %*% B %*% t(Q[, -1])
  if (check){check_omega(p1, q1, Omega)}
  return(list(
    p1 = p1,
    q1 = q1,
    Omega = Omega
  ))
}

#' # Warning
#' Apart from p1 and q1, sign of columns of P and Q cannot be recovered from Omega.
param_omega2cann <- function(p1, q1, Omega, check = TRUE){
  if (check){check_omega(p1, q1, Omega)}
  svdres <- svd(Omega, nu = nrow(Omega) - 1, nv = nrow(Omega) - 1)

  Q <- cbind(q1, svdres$v)
  P <- cbind(p1, svdres$u)
  B <- diag(svdres$d[-length(svdres$d)])
  if (check){check_cann(P, Q, B)}
  return(list(
    P = P,
    Q = Q,
    B = B
  ))
}

check_cann <- function(P, Q, B){
  stopifnot(max(abs(B-diag(diag(B)))) < sqrt(.Machine$double.eps))
  stopifnot(max(abs(P %*% t(P) - diag(1, ncol(P)))) < sqrt(.Machine$double.eps))
  stopifnot(max(abs(t(Q) %*% Q - diag(1, ncol(Q)))) < sqrt(.Machine$double.eps))
  stopifnot(ncol(P) == ncol(Q))
  stopifnot(ncol(B) == ncol(P) - 1)
}

check_omega <- function(p1, q1, Omega){
  stopifnot(max(abs(t(p1) %*% Omega)) < sqrt(.Machine$double.eps))
  stopifnot(max(abs(Omega %*% q1)) < sqrt(.Machine$double.eps))
}
