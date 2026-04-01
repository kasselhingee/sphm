#' SvMF Parameterisations
#' @name SvMFparams
#' @description The SvMF canonical parameterisation is a set of positive real numbers (kappa, `a`) and a matrix of orthonormal column vectors `Gamma`, where the product of the `a` is 1. `a_1` and `kappa` are not jointly estimatable by MLE, so I'm expecting them to be closely related.
#' There are three other parametrisations.
#' Alternatively this set of numbers are the eigenvalues and eigenvectors of the symmetric matrix E(yy^T). The determinant of this matrix would be kappa * a_1.
#' A third parametrisation is available: using kappa, mu and a matrix V that is symmetric positive definite matrix with det(V) = 1.
#' A fourth represents the above V matrix as an orthogonal matrix and a set of positive real values.
#' 
#' @examples 
#' SvMF_cann(k = 0.5, a = c(1, rep(1, 5-1)), G = diag(1, 5))
NULL

#' @describeIn SvMFparams The canonical parameterisation of a SvMF using a kappa, a vector and Gamma, as used in equation 7 of Scealy and Wood 2017.
#' @family SvMF-distribution
#' @export
SvMF_cann <- function(k, a, G){
  obj <- list(k = k, a = a, G = G)
  class(obj) <- c("SvMF_cann", class(obj))
  return(obj)
}
SvMF_cann_check <- function(obj){
  list2env(obj, envir = environment())
  checks <- c(
  `a length` = length(a) == ncol(G),
  `Gamma square` = ncol(G) == nrow(G),
  `kappa positive` = k > 0,
  `a positive` = a > 0,
  `a product` = (prod(a[-1]) - 1)^2 < sqrt(.Machine$double.eps),
  `Gamma orthogonal` = (t(G) %*% G - diag(1, ncol(G)))^2 <  sqrt(.Machine$double.eps)
  )
  SvMF_unimodalcriterion(obj)
  if (any(!checks)){
    stop(sprintf("Parameters fail SvMF checks: %s", paste(names(checks)[!checks], collapse = ", ")))
  }
  return(NULL)
}

#' @describeIn SvMFparams The parameterisation of a SvMF that uses a kappa, mean, a_1, and a symmetric positive definite matrix `V` with `det(V) = 1`. This parameterisation is defined in equations 8 and 9 of Scealy and Wood 2019.
#' The log-likelihood using this parameterisation is in equation 11 of Scealy and Wood 2019.
SvMF_muV <- function(k, m, a1, V){
  obj <- list(k = k, m = m, a1 = a1, V = V)
  class(obj) <- c("SvMF_muV", class(obj))
  return(obj)
}
SvMF_muV_check <- function(obj, tol = sqrt(.Machine$double.eps)){
  list2env(obj, envir = environment())
  checks <- c(
    `kappa positive` = k > 0,
    `a1 positive` = a1 > 0,
    `m size` = (vnorm(m) - 1)^2 <  tol,
    `V square` = ncol(V) == nrow(V),
    `V size` = ncol(V) == length(m) - 1,
    `V symmetric` = isSymmetric(V, tol = tol),
    `det(V)` = (det(V) - 1)^2 < tol
  )
  if (any(!checks)){
    stop(sprintf("Parameters fail checks: %s", paste(names(checks)[!checks], collapse = ", ")))
  }
  return(NULL)
}
as_SvMF_muV <- function(obj){
  if (inherits(obj, "SvMF_muV")){return(obj)}
  if (inherits(obj, "SvMF_cann")){return(SvMF_cann2muV(obj))}
  if (!inherits(obj, "list")){stop("obj isn't a SvMF_muV, SvMF_cann or a list.")}
  if ("V" %in% names(obj)){return(do.call(SvMF_muV, obj))}
  if ("G" %in% names(obj)){return(SvMF_muV(do.call(SvMF_cann, obj)))}
}

# conversion
SvMF_cann2muV <- function(obj){
  list2env(obj, envir = environment())
  m <- G[, 1]
  Hstar <- getHstar(m)
  Kstar <- t(Hstar) %*% G[, -1] #this solves Hstar %*% Kstar = G[, -1] because by orthogonality of H, t(Hstar) %*% Hstar = I.
  V <- Kstar %*% diag(a[-1]^2) %*% t(Kstar)
  SvMF_muV(k, m, a[1], V)
}
SvMF_muV2cann <- function(obj){
  list2env(obj, envir = environment())
  Hstar <- getHstar(m)
  es <- eigen(V)
  Kstar <- es$vectors
  G <- cbind(m, Hstar %*% Kstar)
  rownames(G) <- NULL
  SvMF_cann(k, c(a1, sqrt(es$values)), G = G)
}

#' @noRd 
#' @description Get the Hstar matrix from Scealy and Wood (2019, Section 3) for a mean vector m
getHstar <- function(m){
  m1 <- m[1]
  mL <- m[-1]
  Hstar <- rbind(mL, (1/(1+m1)) * mL %*% t(mL) - diag(1, length(mL)))
  return(Hstar)
}

getH <- function(m){
  cbind(m, getHstar(m))
}

#Scealy and Wood (2019) Proposition 1 check for unimodality
SvMF_unimodalcriterion <- function(cannparam){
  a1 <- cannparam$a[1]
  a2 <- cannparam$a[2]
  shapecalc <- a1*(length(cannparam$a)-1)*((a2/a1)^2 - 1)
  if (a1 >= 1-sqrt(.Machine$double.eps)){
    if (a2 < a1){warning("a2 is smaller than a1 and SvMF may be multimodal.")}
    if (cannparam$k < shapecalc){warning("Given the axial scales, concentration is small so the SvMF may be multimodal.")}
  }
  invisible(shapecalc)
}
