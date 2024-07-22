#' Parameterisation Classes
#' @description Parameterisations of the link functions are stored as lists.
#' @param P P matrix: a p x p (?orthonormal) matrix
#' @param B B matrix: a (p-1) x (p-1) diagonal matrix with elements between zero and one ordered in decreasing size.
#' @param Q The rotation-like matrix `Q` for rotating the covariate vector `x`.
cannS2S <- function(P, Q, B, check = TRUE){
  obj <- list(P = P, Q = Q, B = B)
  class(obj) <- c("cannS2S", class(obj))
  if (check){cannS2S_check(obj)}
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
  names(p1) <- paste0("p1_", 1:length(p1))
  names(q1) <- paste0("q1_", 1:length(q1))
  out <- c(p1, q1, as.vector(Omega))
  names(out)[(1 + length(p1) + length(q1)):length(out)] <- 
    paste0("Omega_", as.vector(outer(1:nrow(Omega), 1:ncol(Omega), function(x,y){paste0(x, ",", y)})))
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
  names(vec) <- NULL
  
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
  if (check){Omega2cann_check(p1, q1, Omega)}
  svdres <- svd(Omega, nu = nrow(Omega) - 1, nv = nrow(Omega) - 1)

  Q <- cbind(q1, svdres$v)
  P <- cbind(p1, svdres$u)
  B <- diag(svdres$d[-length(svdres$d)])
  return(cannS2S(P, Q, B, check = check))
}

cannS2S_check <- function(obj){
  stopifnot(inherits(obj, "cannS2S"))
  list2env(obj, envir = environment())
  stopifnot(max(abs(B-diag(diag(B)))) < sqrt(.Machine$double.eps))
  stopifnot(max(abs(P %*% t(P) - diag(1, ncol(P)))) < sqrt(.Machine$double.eps))
  stopifnot(max(abs(t(Q) %*% Q - diag(1, ncol(Q)))) < sqrt(.Machine$double.eps))
  stopifnot(ncol(P) == ncol(Q))
  stopifnot(ncol(B) == ncol(P) - 1)
  if (any(diag(B) > 1)){warning("Elements of B are larger than 1")}
  if (any(diag(B) < 0)){warning("Elements of B are negative")}
  return(NULL)
}

OmegaS2S_check <- function(obj){
  vals <- OmegaS2S_check_internal(obj)
  good <- (vals < sqrt(.Machine$double.eps))
  if (!all(good)){
    stop(paste("The following checks failed.", 
               paste0(names(vals)[!good], ": ", format(sqrt(vals[!good]), digits = 2), collapse = ", ") #sqrt here converts squared sizes to actual sizes
    ))
  }
  if (sum(diag(t(obj$Omega) %*% obj$Omega)) > nrow(obj$Omega) - 1){warning("The sum of squared singular values of Omega is greater than p - 1, which means that there are singular values of Omega with size greater than 1.")}
  return(NULL)
}
OmegaS2S_check_internal <- function(obj){ #uses squared values for smoothness
  stopifnot(inherits(obj, "OmegaS2S"))
  list2env(obj, envir = environment())
  return(c(
    p1sizediff = (vnorm(p1) - 1)^2,
    q1sizediff = (vnorm(q1) - 1)^2,
    p1Omega = (t(p1) %*% Omega)^2,
    Omegaq1 = (Omega %*% q1)^2
  ))
}
# if the values are close to satisfying the constraints, it might make sense to project and scale p1 and q1 to satisfy the constraints
# will use canonical parameterisation to do this because orthogonality of the columns will make for easier projections
OmegaS2S_proj <- function(obj, method = "Omega"){
  stopifnot(inherits(obj, "OmegaS2S"))
  stopifnot(method %in% c("Omega", "p1q1"))
  if (method == "Omega") {
    list2env(obj, envir = environment())
    # first project orthogonal to p1 (needs p1 unit vector)
    p1 <- p1/vnorm(p1)
    Omegaperpp1 <- Omega -  (p1 %*% t(p1)) %*% Omega
    # similarly to q1
    q1 <- q1/vnorm(q1)
    Omegaperpq1 <- Omegaperpp1 - Omegaperpp1 %*% q1 %*% t(q1)
    return(OmegaS2S(p1, q1, Omegaperpq1, check = FALSE))
  }
  if (method == "p1q1") {
    cann <- Omega2cann(obj, check = FALSE)
    newp1 <- cann$P[,1] - cann$P[, -1] %*% t(cann$P[, -1]) %*% cann$P[,1]
    newp1 <- newp1 / vnorm(newp1)
    cann$P[,1] <- newp1
  
    newq1 <- cann$Q[,1] - cann$Q[, -1] %*% t(cann$Q[, -1]) %*% cann$Q[,1]
    newq1 <- newq1 / vnorm(newq1)
    cann$Q[,1] <- newq1
    return(as_OmegaS2S(cann))
  }
}
