#' @name mnlink_params
#' @title Parameterisation Classes
#' @description Parameterisations of the link functions.
#' These methods check and convert between parameterisations.
#' Actual mean link calculations are performed by other functions.
#' @param P Final rotation matrix on the response sphere: a p x p orthonormal matrix with positive determinant.
#' @param Bs Scaling matrix for spherical covariates: a (p-1) x (p-1) diagonal matrix with elements between zero and one ordered in decreasing size. `NULL` if no Sph covariates.
#' @param Be Scaling matrix for Euclidean covariates: a (p-1) x (p-1) diagonal matrix with elements between zero and one ordered in decreasing size. `NULL` if no Euc covariates.
#' @param Qs The qs x p rotation-like matrix `R_s` for rotating the spherical covariate vector. `NULL` if no Sph covariates.
#' @param Qe The qe x p rotation-like matrix `R_e` for rotating the Euclidean covariate vector. `NULL` if no Euc covariates.
#' @param ce The additive offset \eqn{c_e} for Euclidean covariates only. Vector of length qe. `NULL` if no Euc covariates.
#' @details
#' # Cannonical Parameterisation
#' The `P`, `Bs`, `Be`, `Qs`, `Qe` and `ce` is slightly more flexible than Shogo's link function with both Euclidean covariates and a spherical covariate that matches Remark 2 of Manuscript (Nov 29, 2024).
#' The link (1) from that manuscript can be obtained by including an extra zero-valued Euclidean covariate as the first covariate and forcing \eqn{q_{1e}} to be `(1, 0, ...)` to match the index of the constant covariate and setting `ce[1]` to be zero. I think these changes will not affect the estimation method as both \eqn{q_{1e}} and `ce[1]` separate out of the "Omega" parameterisation.
#' 
#' Andy's link function for Euclidean covariates needs an additional scaling \eqn{b_{im}} parameter for the imaginary component Andy's link function to be parameterised. It will also need `ce = 0` and `Bs` and `Qs` will be ignored since spherical covariates not incorporated yet.
#' 
#' # Omega Parameterisation
#' The link functions are simplified by writing \eqn{\Omega_s = P^* B_s {Q_s^*}^T} and \eqn{\Omega_e = P^* B_e {Q_e^*}^T}, and finally \eqn{\Omega = [\Omega_s \,\, \Omega_e]}.
#' This parameterisation helps optimisation as optimisation in Stiefel manifolds is harder than other spaces, and also reflects the sign ambiguity of columns of P with the matching columns of `Qe` and `Qs`.
#' 
NULL

cannS2S <- function(P, Q, B, check = TRUE){
  mnlink_cann(P = P, Bs = B, Qs = Q)
}
mnlink_cann <- function(P, Bs = NULL, Qs = NULL, Be = NULL, Qe = NULL, ce = NULL, check = TRUE){
  stopifnot(is.matrix(P))
  obj <- list(P = P, Bs = Bs, Qs = Qs, Be = Be, Qe = Qe)
  class(obj) <- c("mnlink_cann", class(obj))
  if (check){mnlink_cann_check(obj)}
  return(obj)
}
as_mnlink_cann <- function(obj){
  if (inherits(obj, "mnlink_cann")){return(obj)}
  if (inherits(obj, "mnlink_Omega")){return(Omega2cann(obj, check = FALSE))}
  if (!inherits(obj, "list")){stop("obj isn't a cannS2S, OmegaS2S or a list.")}
  if ("P" %in% names(obj)){return(mnlink_cann(obj, check = FALSE))}
  if ("p1" %in% names(obj)){return(Omega2cann(mnlink_Omega(obj, check = FALSE), check = FALSE))}
  return(obj)
}

#' @name mnlink_params
#' @param p1 First column of the P matrix (vector of length `p`)
#' @param qe1 First column of the Qe matrix (vector of length `qe`). `NULL` if no Euc covariates.
#' @param qs1 First column of the Qs matrix (vector of length `qs`). `NULL` if no Sph covariates.
#' @param Omega A `p` by `qe + qs` matrix representing 
#' \deqn{\Omega = [\Omega_s \Omega_e] = [P^* B_s {Q_s^*}^T \,   P^* B_e {Q_e^*}^T]}
NULL

OmegaS2S <- function(p1, q1, Omega, check = TRUE){
  mnlink_Omega(p1 = p1, qs1 = q1, Omega = Omega, check = check)
}

mnlink_Omega <- function(p1, qs1 = NULL, Omega, qe1 = NULL, ce = NULL, check = TRUE){
  obj <- list(
    p1 = p1,
    qs1 = qs1,
    qe1 = qe1,
    Omega = Omega,
    ce = ce
  )
  class(obj) <- c("mnlink_Omega", class(obj))
  if (check) {mnlink_Omega_check(obj)}
  return(obj)
}

as_mnlink_Omega <- function(obj){
  if (inherits(obj, "mnlink_cann")){return(cann2Omega(obj, check = FALSE))}
  if (inherits(obj, "mnlink_Omega")){return(obj)}
  if (!inherits(obj, "list")){stop("obj must be either a cannS2S, OmegaS2S or list.")}
  if ("P" %in% names(obj)){return(cann2Omega(mnlink_cann(obj, check = FALSE)))}
  if ("p1" %in% names(obj)){return(mnlink_Omega(obj, check = FALSE))}
  return(obj)
}
#' @noRd
#' Vectorsises and inverse of this vectorisation for the OmegaS2S parameterisation.
#' @param obj An OmegaS2S parameter object.
OmegaS2S_vec <- function(obj){
  stopifnot(inherits(obj, "mnlink_Omega"))
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
  if (check){mnlink_cann_check(obj)}
  p1 <- obj$P[, 1]
  qs1 <- obj$Qs[, 1]
  qe1 <- obj$Qe[, 1]
  Omega_s <- Omega_e <- NULL
  if (!is.null(qs1)){
    Omega_s <- obj$P[,-1] %*% obj$Bs %*% t(obj$Qs[, -1])
  }
  if (!is.null(qe1)){
    Omega_e <- obj$P[,-1] %*% obj$Be %*% t(obj$Qe[, -1])
  }
  Omega <- cbind(Omega_s, Omega_e)
  return(mnlink_Omega(p1, qs1 = qs1, Omega = Omega, qe1 = qe1, ce = obj$ce, check = FALSE))
}

#' # Warning
#' Apart from p1 and q1, sign of columns of P and Q cannot be recovered from Omega.
Omega2cann <- function(obj, check = TRUE){
  if (check){mnlink_Omega_check(obj)}
  list2env(obj, envir = environment())
  svdres <- svd(Omega, nu = nrow(Omega) - 1, nv = nrow(Omega) - 1)

  Q <- cbind(q1, svdres$v)
  P <- cbind(p1, svdres$u)
  B <- diag(svdres$d[-length(svdres$d)])
  return(cannS2S(P, Q, B, check = check))
}

mnlink_cann_check <- function(obj){
  stopifnot(inherits(obj, "mnlink_cann"))
  
  #check P matrix
  p <- nrow(obj$P)
  stopifnot(p == ncol(obj$P))
  stopifnot(max(abs(obj$P %*% t(obj$P) - diag(1, ncol(obj$P)))) < sqrt(.Machine$double.eps))
  
  #check spherical covariate parameters
  if (!is.null(obj$Qs)){
    stopifnot(!is.null(obj$Bs))
    stopifnot(ncol(obj$Bs) == p - 1)
    stopifnot(nrow(obj$Bs) == p - 1)
    stopifnot(ncol(obj$Qs) == p)
    
    stopifnot(max(abs(obj$Bs-diag(diag(obj$Bs)))) < sqrt(.Machine$double.eps))
    stopifnot(max(abs(t(obj$Qs) %*% obj$Qs - diag(1, ncol(obj$Qs)))) < sqrt(.Machine$double.eps))
    if (any(diag(obj$Bs) > 1)){warning("Elements of Bs are larger than 1")}
    if (any(diag(obj$Bs) < 0)){warning("Elements of Bs are negative")}
  } else {
    stopifnot(is.null(obj$Bs))
    stopifnot(is.null(obj$ce))
  }
  
  #check Euc covariate parameters
  if (!is.null(obj$Qe)){
    stopifnot(!is.null(obj$Be))
    stopifnot(ncol(obj$Be) == p - 1)
    stopifnot(nrow(obj$Be) == p - 1)
    stopifnot(ncol(obj$Qe) == p)
    
    stopifnot(max(abs(obj$Be-diag(diag(obj$Be)))) < sqrt(.Machine$double.eps))
    stopifnot(max(abs(t(obj$Qe) %*% obj$Qe - diag(1, ncol(obj$Qe)))) < sqrt(.Machine$double.eps))
    stopifnot(is.vector(obj$ce))
    stopifnot(length(obj$ce) == p-1)
    if (any(diag(obj$Be) > 1)){warning("Elements of Be are larger than 1")}
    if (any(diag(obj$Be) < 0)){warning("Elements of Be are negative")}
  } else {
    stopifnot(is.null(obj$Be))
    stopifnot(is.null(obj$ce))
  }
  return(NULL)
}

mnlink_Omega_check <- function(obj){
  # Check dimensions and nullness of elements
  stopifnot(length(obj$p1) == nrow(obj$Omega))
  stopifnot(length(obj$qs1) + length(obj$qe1) == ncol(obj$Omega))
  stopifnot(is.null(obj$qe1) + is.null(obj$ce) %in% c(0, 2))
  stopifnot(length(obj$ce) == length(obj$qe1))
  
  vals <- mnlink_Omega_check_numerical(obj)
  good <- (vals < sqrt(.Machine$double.eps))
  if (!all(good)){
    stop(paste("The following checks failed.", 
               paste0(names(vals)[!good], ": ", format(sqrt(vals[!good]), digits = 2), collapse = ", ") #sqrt here converts squared sizes to actual sizes
    ))
  }
  # sum of squared singular values is sum(Bs^2 + Be^2).
  # If all of Bs and Be are less than or equal to 1 then sum of squared singular values is less than 2*(p-1) if there are both Spherical and Euclidean covariates 
  singularvalssumsquared <- sum(diag(t(obj$Omega) %*% obj$Omega))
  if (singularvalssumsquared > (!is.null(obj$qe1) + !is.null(obj$qs1)) * (nrow(obj$Omega) - 1)){warning(sprintf("The sum of squared singular values of Omega is %0.2f, which means that there are scales in either Be or Bs that are greater than 1.", singularvalssumsquared))}
  
  return(NULL)
}
mnlink_Omega_check_numerical <- function(obj){ #uses squared values for smoothness
  stopifnot(inherits(obj, "mnlink_Omega"))
  # list2env(obj, envir = environment())
  qs <- switch(is.null(obj$qs1), 0, length(obj$qs1))
  qe <- switch(is.null(obj$qe1), 0, length(obj$qe1))
  checkvals <- c(
    p1sizediff = (vnorm(obj$p1) - 1)^2,
    p1Omega = (t(obj$p1) %*% obj$Omega)^2
  )
  if (qs > 0){
    checkvals <- c(
      checkvals,
      qs1sizediff = (vnorm(obj$qs1) - 1)^2,
      Omegaqs1 = (obj$Omega[, seq.int(1, qs)] %*% obj$qs1)^2
    )
  }
  if (qe > 0){
    checkvals <- c(
      qe1sizediff = (vnorm(obj$qe1) - 1)^2,
      Omegaqe1 = (obj$Omega[, qs + seq.int(1, qe)] %*% obj$qe1)^2
    )
  }
  return(checkvals)
}
# if the values are close to satisfying the constraints, it might make sense to project and scale p1 and q1 to satisfy the constraints
# will use canonical parameterisation to do this because orthogonality of the columns will make for easier projections
OmegaS2S_proj <- function(obj, method = "Omega"){
  stopifnot(inherits(obj, "mnlink_Omega"))
  stopifnot(method %in% c("Omega", "p1q1"))
  if (method == "Omega") {
    list2env(obj, envir = environment())
    # first project orthogonal to p1 (needs p1 unit vector)
    p1 <- p1/vnorm(p1)
    Omegaperpp1 <- Omega -  (p1 %*% t(p1)) %*% Omega
    # now t(p1) %*% Omegaperpp1 = 0
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
    return(as_mnlink_Omega(cann))
  }
}
