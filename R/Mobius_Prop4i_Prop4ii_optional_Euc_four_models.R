#' @name mobius_link_params
#' @title Parameterisation Classes for the Scaled Mobius Mean Link
#' @description Parameterisations of the link functions.
#' These methods, create, check and convert between parameterisations of the mean link.
#' Actual mean link calculations are performed by other functions.
#' @param P Final rotation matrix on the response sphere: a p x p orthonormal matrix with positive determinant.
#' @param Bs Scaling matrix for spherical covariates: a (p-1) x (p-1) diagonal matrix with elements between zero and one ordered in decreasing size. `NULL` if no Sph covariates.
#' @param Be Scaling matrix for Euclidean covariates: a (p-1) x (p-1) diagonal matrix with elements between zero and one ordered in decreasing size. `NULL` if no Euc covariates.
#' @param Qs The qs x p rotation-like matrix `R_s` for rotating the spherical covariate vector. `NULL` if no Sph covariates.
#' @param Qe The qe x p rotation-like matrix `R_e` for rotating the Euclidean covariate vector. `NULL` if no Euc covariates.
#' @param ce A single real value. The additive offset \eqn{c_e} in the denominator below Euclidean covariates in the 'H' link. `NULL` if no Euc covariates.
#' @details
#' # Cannonical Parameterisation
#' The parameters here are for a mean link defined as
#' \deqn{\mu_{SpEuc}(x) = P\mathcal{S}^{[-1]}\left(B_s \mathcal{S}(Q_s^\top x_s)  +  \frac{B_e(Q_e[,-1]^\top x_e}{Q_e[,1]^\top x_e + c_e}\right)}
#' where \eqn{\mu(x)} is a p-length unit vector,
#' \eqn{x_s} is a `qs`-length unit vector,
#' \eqn{x_e} is a `qe`-length vector,
#' `P` is a p x p rotation matrix,
#' `Qs` is a `qs x p` orthonormal matrix (`t(Qs) %*% Qs == diag(p)`),
#' `Bs` is a diagonal matrix (p-1) x (p-1) matrix,
#' `Be` is a diagonal matrix (p-1) x (p-1) matrix,
#' `Qs` is a `qe x p` orthonormal matrix,
#' `ce` is a real-valued number.
#' The Euclidean covariate in this link is essentially transformed by \eqn{B_e\mathcal{S}\left(Q_e^\top x/c_1 \right)}, thus we call this link the `SpEuc` link.
#'
#' The `SpEuc` link above is slightly more flexible than the link function defined at the start of "Regression for spherical responses with linear and spherical covariates using a scaled link function" manuscript.
#' The link in Equation (1) from that manuscript is 
#' \deqn{\mu_{LinEuc}(x) = P\mathcal{S}^{[-1]}\left(B_s \mathcal{S}(Q_s^\top x_s)  +  B_e(Q_e[,-1]^\top x_e\right),}
#' which can be written as \eqn{\mu_{SpEuc}} with an extra zero-valued Euclidean covariate as the first covariate and fixing \eqn{Q_e[,1]=q_{e1}} to be `(1, 0, ...)` to match the index of the constant covariate and setting `ce=1`.
#' This is the `LinEuc` (for linear Euclidean) link form.
#' 
#' # Omega Parameterisation
#' The link functions are simplified by writing \eqn{\Omega_s = P^* B_s {Q_s^*}^T} and \eqn{\Omega_e = P^* B_e {Q_e^*}^T}, \eqn{\Omega = [\Omega_s \,\, \Omega_e]}.
#' This parameterisation helps optimisation as optimisation in Stiefel manifolds is harder than other spaces, and also reflects the sign ambiguity of columns of P with the matching columns of `Qe` and `Qs`.
#' 
NULL

#' @describeIn mobius_link_params Represent the parameters of the Mobius mean link [`mobius_link()`] using the canonical form.
#' @family link-function
#' @export
mobius_link_cann <- function(P, Bs = NULL, Qs = NULL, Be = NULL, Qe = NULL, ce = NULL, check = TRUE){
  stopifnot(is.matrix(P))
  obj <- list(P = P, Bs = Bs, Qs = Qs, Be = Be, Qe = Qe, ce = ce)
  obj <- lapply(obj, function(element){if (length(element) == 0){return(NULL)}else{return(element)}})

  # export
  class(obj) <- c("mobius_link_cann", class(obj))
  if (check){mobius_link_cann_check(obj)}
  return(obj)
}

#' @describeIn mobius_link_params Represent the parameters of the Mobius mean link [`mobius_link()`] using the canonical form.
#' @family link-function
#' @export
as_mobius_link_cann <- function(obj){
  if (inherits(obj, "mobius_link_cann")){return(obj)}
  if (inherits(obj, "mobius_link_Omega")){return(Omega2cann(obj, check = FALSE))}
  if (!inherits(obj, "list")){stop("obj isn't a mobius_link_cann or mobius_link_Omega class, or a list.")}
  if ("P" %in% names(obj)){return(do.call(mobius_link_cann, c(obj, list(check = FALSE))))}
  if ("p1" %in% names(obj)){return(Omega2cann(do.call(mobius_link_Omega, c(obj, list(check = FALSE))), check = FALSE))}
  return(obj)
}

# Flatten a canonical mean link object to a named numeric vector. Used for
# passing parameters to/from C++ optimisers and numerical checks.
mobius_link_cann_vec <- function(obj){
  stopifnot(inherits(obj, "mobius_link_cann"))
  Pvec <- as.vector(obj$P)
  names(Pvec) <- as.vector(outer(seq.int(1, nrow(obj$P)), seq.int(1, nrow(obj$P)), 
                                 FUN = paste, sep = ","))
  names(Pvec) <- paste0("P", names(Pvec))
  Qsvec <- as.vector(obj$Qs)
  if (!is.null(Qsvec)){
    names(Qsvec) <- as.vector(outer(seq.int(1, nrow(obj$Qs)), seq.int(1, ncol(obj$Qs)), 
                                    FUN = paste, sep = ","))
    names(Qsvec) <- paste0("Qs", names(Qsvec))
  }
  Qevec <- as.vector(obj$Qe)
  if (!is.null(Qevec)){
    names(Qevec) <- as.vector(outer(seq.int(1, nrow(obj$Qe)), seq.int(1, ncol(obj$Qe)), 
                                    FUN = paste, sep = ","))
    names(Qevec) <- paste0("Qe", names(Qevec))
  }
  Bsvec <- as.vector(diag(obj$Bs))
  if (length(Bsvec) > 0){names(Bsvec) <- paste0("Bs", 1:(nrow(obj$Bs)))}
  Bevec <- as.vector(diag(obj$Be))
  if (length(Bevec) > 0){names(Bevec) <- paste0("Be", 1:(nrow(obj$Be)))}
  cevec <- obj$ce
  if (!is.null(cevec)){names(cevec) <- paste0("ce", c("",1:(length(cevec) - 1)))}
  return(c(Pvec, Bsvec, Qsvec, Bevec, Qevec, cevec))
}




#' @name mobius_link_params
#' @param p1 First column of the P matrix (vector of length `p`)
#' @param qe1 First column of the Qe matrix (vector of length `qe`). `NULL` if no Euc covariates.
#' @param qs1 First column of the Qs matrix (vector of length `qs`). `NULL` if no Sph covariates.
#' @param Omega A `p` by `qe + qs` matrix representing 
#' \deqn{\Omega = [\Omega_s \Omega_e] = [P^* B_s {Q_s^*}^T \,   P^* B_e {Q_e^*}^T]}
#' @param ce The value of \eqn{c_e}.
NULL

#' @describeIn mobius_link_params Represent the parameters of the Mobius mean link [`mobius_link()`] using the Omega form.
#' @family link-function
#' @export
mobius_link_Omega <- function(p1, qs1 = vector("numeric", 0), Omega, qe1 = vector("numeric", 0), ce = vector("numeric", 0), check = TRUE){
  if (is.null(qs1)){qs1 <- vector("numeric", 0)}
  if (is.null(qe1)){qe1 <- vector("numeric", 0)}
  if (is.null(ce)){ce <- vector("numeric", 0)}
  obj <- list(
    p1 = p1,
    qs1 = qs1,
    qe1 = qe1,
    Omega = Omega,
    ce = ce
  )

  #export
  class(obj) <- c("mobius_link_Omega", class(obj))
  if (check) {mobius_link_Omega_check(obj)}
  return(obj)
}

#' @describeIn mobius_link_params Represent the parameters of the Mobius mean link [`mobius_link()`] using the Omega form.
#' @family link-function
#' @export
as_mobius_link_Omega <- function(obj){
  if (inherits(obj, "mobius_link_cann")){return(cann2Omega(obj, check = FALSE))}
  if (inherits(obj, "mobius_link_Omega")){return(obj)}
  if (!inherits(obj, "list")){stop("obj must be either a mobius_link_cann or mobius_link_Omega class, or a list.")}
  if ("P" %in% names(obj)){return(cann2Omega(do.call(mobius_link_cann, c(obj, list(check = FALSE))), check = FALSE))}
  if ("p1" %in% names(obj)){return(do.call(mobius_link_Omega, c(obj, list(check = FALSE))))}
  return(obj)
}

#' @noRd
#' @description Converts Omega parameterisation to a vector.
#' @param obj An `mobius_link_Omega` parameter object.
mobius_link_Omega_vec <- function(obj){
  stopifnot(inherits(obj, "mobius_link_Omega"))
  p1 <- obj$p1
  qs1 <- obj$qs1
  qe1 <- obj$qe1
  Omega <- obj$Omega
  ce <- obj$ce
  Omegavec <- as.vector(Omega)
  names(Omegavec) <- as.vector(
      outer(seq.int(1, length.out = length(p1)), 
        c( paste0("s", seq.int(1, length.out = length(qs1)), recycle0 = TRUE),
           paste0("e", seq.int(1, length.out = length(qe1)), recycle0 = TRUE)),
        FUN = paste, sep = ","))
  names(Omegavec) <- paste0("Omega_", names(Omegavec))
  names(p1) <- paste0("p1_", seq.int(1, length.out = length(p1)), recycle0 = TRUE)
  names(qs1) <- paste0("qs1_", seq.int(1, length.out = length(qs1)), recycle0 = TRUE)
  names(qe1) <- paste0("qe1_", seq.int(1, length.out = length(qe1)), recycle0 = TRUE)
  names(ce) <- paste0("ce", seq.int(1, length.out = length(ce)), recycle0 = TRUE)
  out <- c(p1, qs1, qe1, Omegavec, ce)
  class(out) <- "mobius_link_Omega_vec"
  return(out)
}

#' @noRd
#' @description Inverse of `mobius_link_Omega_vec()`
#' @param vec is a vector created by `mobius_link_Omega_vec()`
#' @param p The dimension of the response (The dimension of covariates will be infered from `p`).
#' @param qe Number of Euclidean covariates
mobius_link_Omega_unvec <- function(vec, p, qe = 0, check = TRUE){
  # length of vec = p + qs + qe + p*(qs + qe) + (qe>0)
  # l - 1*(qe>0) - qe - p*qe = qs + p*qs
  # (l - (qe>0) - qe - p*qe)/(1+p) = qs
  qs <- (length(vec) - p - (qe>0) - qe - p*qe)/(1+p)
  stopifnot(qs == as.integer(qs))
  names(vec) <- NULL
  
  mobius_link_Omega(p1 = vec[1:p],
           qs1 = vec[p + seq.int(1, length.out = qs)],
           qe1 = vec[p + qs + seq.int(1, length.out = qe)],
           Omega = matrix(vec[p + qs + qe + seq.int(1, length.out = p*(qe + qs))], nrow = p, ncol = qs + qe, byrow = FALSE),
           ce = vec[p + qs + qe + p*(qs + qe) + seq.int(1, length.out = (qe > 0))],
           check = check)
}

#' @describeIn mobius_link_params For converting between parameterisations of the link function.
#' Sign of columns of P and Q are lost by the `cann2Omega()` transformation.
#' @family link-function
#' @export
cann2Omega <- function(obj, check = TRUE){
  if (check){mobius_link_cann_check(obj)}
  p1 <- obj$P[, 1]
  qs1 <- obj$Qs[, 1]
  qe1 <- obj$Qe[, 1]
  Omega_s <- Omega_e <- ce <- NULL
  if (!is.null(qs1)){
    # P[,-1] and Qs[,-1] are the P* and Qs* sub-matrices (columns 2:p of P and Qs);
    # column 1 of each is the base-point direction and does not appear in Omega.
    Omega_s <- obj$P[,-1] %*% obj$Bs %*% t(obj$Qs[, -1])
  }
  if (!is.null(qe1)){
    Omega_e <- obj$P[,-1] %*% obj$Be %*% t(obj$Qe[, -1])
    ce <- obj$ce
  }
  Omega <- cbind(Omega_s, Omega_e)
  return(mobius_link_Omega(p1, qs1 = qs1, Omega = Omega, qe1 = qe1, ce = ce, check = FALSE))
}

#' @noRd
#' @description Converts an Omega parameterisation to the cannonical parameterisation.
#' Apart from p1 and q1, sign of columns of P and Q cannot be recovered from Omega.
Omega2cann <- function(obj, check = TRUE){
  if (check){mobius_link_Omega_check(obj)}
  svdres <- svd(obj$Omega, nu = nrow(obj$Omega) - 1, nv = nrow(obj$Omega) - 1)

  P <- cbind(obj$p1, svdres$u)
  
  # much of the rest uses the SVD of Omega as written in the Euclidean link document (not ce)
  Qs <- Qe <- Bs <- Be <- ce <- NULL
  if (length(obj$qs1) > 0){
    # svdres$v rows carry unnormalised Qs* columns; their norms encode the singular values
    # of Omega, so Bs = diag(norms * singular_values) recovers the diagonal scale matrix.
    Qs_unnorm <- svdres$v[seq.int(1, length.out = length(obj$qs1)), , drop = FALSE]
    Qs_norms <- sqrt(colSums(Qs_unnorm^2))
    Qsstar <- t(t(Qs_unnorm)/ Qs_norms)
    Qs <- cbind(obj$qs1, Qsstar)
    Bs <- diag(Qs_norms * svdres$d[-nrow(obj$Omega)])
  }
  if (length(obj$qe1) > 0){
    Qe_unnorm <- svdres$v[length(obj$qs1) + seq.int(1, length.out = length(obj$qe1)), , drop = FALSE]
    Qe_norms <- sqrt(colSums(Qe_unnorm^2))
    Qestar <- t(t(Qe_unnorm)/ Qe_norms)
    Qe <- cbind(obj$qe1, Qestar)
    Be <- diag(Qe_norms * svdres$d[-nrow(obj$Omega)])
    ce <- obj$ce
  }
  
  out <- mobius_link_cann(P, Bs = Bs, Qs = Qs, Be = Be, Qe = Qe, ce = ce, check = check)
  if (det(out$P) < 0){#if negative determinant, flip an axis to get positive determinant of P (i.e. a rotation matrix)
    out <- P_signswitch(out, ncol(out$P))
  }
  return(out)
}

mobius_link_cann_check <- function(obj){
  stopifnot(inherits(obj, "mobius_link_cann"))
  
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
    
    row(obj$Bs)!=col(obj$Bs)
    stopifnot(max(abs(obj$Bs[row(obj$Bs)!=col(obj$Bs)]), 0) < sqrt(.Machine$double.eps))
    stopifnot(max(abs(t(obj$Qs) %*% obj$Qs - diag(1, ncol(obj$Qs)))) < sqrt(.Machine$double.eps))
    if (any(diag(obj$Bs) < 0)){warning("Elements of Bs are negative")}
  } else {
    stopifnot(is.null(obj$Bs))
  }
  
  #check Euc covariate parameters
  if (!is.null(obj$Qe)){
    stopifnot(!is.null(obj$Be))
    stopifnot(ncol(obj$Be) == p - 1)
    stopifnot(nrow(obj$Be) == p - 1)
    stopifnot(ncol(obj$Qe) == p)
    
    stopifnot(max(abs(obj$Be[row(obj$Be)!=col(obj$Be)]), 0) < sqrt(.Machine$double.eps))
    stopifnot(max(abs(t(obj$Qe) %*% obj$Qe - diag(1, ncol(obj$Qe)))) < sqrt(.Machine$double.eps))
    stopifnot(is.vector(obj$ce))
    stopifnot(length(obj$ce) == 1)
    if (any(diag(obj$Be) < 0)){warning("Elements of Be are negative")}
  } else {
    stopifnot(is.null(obj$Be))
    stopifnot(is.null(obj$ce))
  }
  return(NULL)
}

mobius_link_Omega_check <- function(obj){
  # Check dimensions and nullness of elements
  if (is.null(obj$p1)){stop("p1 should be non-null")}
  if (is.null(obj$qs1)){stop("qs1 should be non-null")}
  if (is.null(obj$qe1)){stop("qe1 should be non-null")}
  if (is.null(obj$ce)){stop("ce should be non-null")}
  
  elementnames <- c("p1", "qs1", "qe1", "ce", "Omega")
  nullelements <- vapply(obj[elementnames], is.null, FUN.VALUE = FALSE) 
  names(nullelements) <- elementnames #needed in case an element is completely missing
  if (any(nullelements)){stop("The following elements are null: ", paste(names(which(nullelements)), collapse = ", "))}
  isvecs <- vapply(obj[names(obj)!="Omega"], is.vector, FUN.VALUE = FALSE)
  if (any(!isvecs)){stop("The following elements should be vectors: ", paste(names(which(!isvecs)), collapse = ", "))}

  stopifnot(length(obj$p1) == nrow(obj$Omega))
  stopifnot(length(obj$qs1) + length(obj$qe1) == ncol(obj$Omega))
  stopifnot( ( (length(obj$qe1) > 0) + (length(obj$ce) > 0)) %in% c(0, 2))
  if(length(obj$qe1) > 0){
    stopifnot(length(obj$ce) == 1)
  }
  
  vals <- mobius_link_Omega_check_numerical(obj)
  good <- (vals < sqrt(.Machine$double.eps))
  if (!all(good)){
    stop(paste("The following checks failed.", 
               paste0(names(vals)[!good], ": ", format(sqrt(vals[!good]), digits = 2), collapse = ", ") #sqrt here converts squared sizes to actual sizes
    ))
  }
  return(NULL)
}
mobius_link_Omega_check_numerical <- function(obj){ #uses squared values for smoothness
  stopifnot(inherits(obj, "mobius_link_Omega"))
  # list2env(obj, envir = environment())
  qs <- length(obj$qs1)
  qe <- length(obj$qe1)
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
      checkvals,
      qe1sizediff = (vnorm(obj$qe1) - 1)^2,
      Omegaqe1 = (obj$Omega[, qs + seq.int(1, qe)] %*% obj$qe1)^2
    )
  }
  if ((qs > 0) & (qe > 0)) {
    # Commutativity check: the Möbius link is well-defined only when the spherical and
    # Euclidean sub-matrices of Omega commute under the products OmOm and OmpartOmpart.
    # This ensures the link's stereographic and Euclidean components compose consistently.
    # See ExtraOmegaConstraint.pdf for the derivation; the C++ version in Omega.cpp gives
    # a more refined vector-valued constraint used in the optimiser.
    OmOm <- obj$Omega %*% t(obj$Omega)
    Is_tilde <- diag(1, qs + qe, qs)
    OmpartOmpart <- obj$Omega %*% (Is_tilde %*% t(Is_tilde)) %*% t(obj$Omega)
    Omega_comm = (sum((OmOm %*% OmpartOmpart - OmpartOmpart %*% OmOm)^2)^2)
    # second ^2 is a scaling hack to match the magnitude of other checks; see C++ for refined version
    checkvals <- c(checkvals, Omega_comm = Omega_comm)
  }
  return(checkvals)
}

# Project an Omega parameter object to exactly satisfy the Omega constraints:
# normalise p1, qs1, qe1 to unit length, then project Omega to be orthogonal to p1,
# Omega_s to be orthogonal to qs1, and Omega_e to be orthogonal to qe1.
# Used after small numerical violations during optimisation.
Omega_proj <- function(obj){
  stopifnot(inherits(obj, "mobius_link_Omega"))
  # first project orthogonal to p1 (needs p1 unit vector)
  obj$p1 <- obj$p1/vnorm(obj$p1)
  newOmega <- obj$Omega -  (obj$p1 %*% t(obj$p1)) %*% obj$Omega
  # now t(p1) %*% newOmega = 0
  
  Omega_s <- Omega_e <- NULL
  if (length(obj$qs1) > 0){# project Omega_s perpendicular to qs1
    obj$qs1 <- obj$qs1/vnorm(obj$qs1)
    Omega_s <- newOmega[, seq.int(1, length.out = length(obj$qs1)), drop = FALSE]
    Omega_s <- Omega_s - Omega_s %*% obj$qs1 %*% t(obj$qs1)
  }
  if (length(obj$qe1) > 0){# project Omega_e perpendicular to qe1
    obj$qe1 <- obj$qe1/vnorm(obj$qe1)
    Omega_e <- newOmega[, length(obj$qs1) + seq.int(1, length.out = length(obj$qe1)), drop = FALSE]
    Omega_e <- Omega_e - Omega_e %*% obj$qe1 %*% t(obj$qe1)
  }
  obj$Omega <- cbind(Omega_s, Omega_e)
  return(obj)
}

# Switch the sign of the Euclidean part of the mean link (qe1, ce, and Omega_e columns).
# The Euclidean part of the Omega parameterisation is invariant to this sign flip because
# ce is a free parameter; used to standardise sign conventions when recovering canonical params.
Euc_signswitch <- function(obj){
  if (inherits(obj, "mobius_link_Omega")){
    obj$qe1 <- -1 * obj$qe1
    obj$ce <- -1 * obj$ce
    obj$Omega[, length(obj$qs1) + (1:length(obj$qe1))] <- -1 * obj$Omega[, length(obj$qs1) + (1:length(obj$qe1))]
    return(obj)
  }
  if (inherits(obj, "mobius_link_cann")){
    obj$ce <- -1 * obj$ce
    obj$Qe <- -1 * obj$Qe
    return(obj)
  }
}
# Flip the sign of selected columns of P (and corresponding columns of Qs, Qe)
# in a canonical parameter object. Used to enforce a canonical sign convention
# (e.g. making P a rotation matrix after SVD). Column 1 cannot be switched.
P_signswitch <- function(obj, cols){
  stopifnot(inherits(obj, "mobius_link_cann"))
  if (is.logical(cols)){cols <- which(cols)}
  if (1 %in% cols){stop("Sign of first column of P cannot be easily swapped")}
  obj$P[, cols] <- -1 * obj$P[,cols]
  if (!is.null(obj$Qs)){obj$Qs[,cols] <- -1 * obj$Qs[,cols]}
  if (!is.null(obj$Qe)){
    obj$Qe[,cols] <- -1 * obj$Qe[,cols]
  }
  return(obj)
}

#' @family link-function
#' @title Randomly generate mean link parameters
#' @description Uniformly generates orthogonal matrices `Qs`, `Qe` and `P`.
#' The diagonal elements of `Bs` and `Be` are drawn uniformly from between `0` and `1`.
#' The parameter `ce` is generated uniformly from between `0` and `1`.
#' @param preseed Before the generation of each parameter matrix, the `seed` is incremented by `1`, starting from `preseed`.
#' @param p The length of response vectors
#' @param qs The length of spherical covariate vectors
#' @param qe The length of Euclidean covariate vectors
#' @return A `mobius_link_cann` object.
#' @export
rand_mobius_link_cann <- function(p = 3, qs = 5, qe = 4, preseed = 0){
  stopifnot((qe==0) | (qe >= p-1))
  stopifnot((qs==0) | (qs >= p-1))
  set.seed(preseed + 1)
  P <- mclust::randomOrthogonalMatrix(p, p)
  Qs <- Qe <- Bs <- Be <- ce <- NULL
  if (qs > 0){
    set.seed(preseed + 2)
    Qs <- mclust::randomOrthogonalMatrix(qs, p)
    set.seed(preseed + 3)
    Bs <- diag(x = sort(runif(p-1), decreasing = TRUE), nrow = p-1)
  }
  if (qe > 0){
    set.seed(preseed + 2 + 10)
    Qe <- mclust::randomOrthogonalMatrix(qe, p)
    set.seed(preseed + 3 + 10)
    Be <-diag(x = sort(runif(p-1), decreasing = TRUE), nrow = p-1) 
    set.seed(preseed + 4 + 10)
    ce <- runif(1)
  }
  paramobj <- mobius_link_cann(P, Bs = Bs, Qs = Qs, Be = Be, Qe = Qe, ce = ce, check = FALSE)
  return(paramobj)
}

#' @noRd
#' @description To simplify testing places result of `rand_mobius_link_cann()` into the calling environment.
#' @details
#' Places all the parameters in the environment by running list2env(obj, envir = environment())
rand_mobius_link_cann__place_in_env <- function(p = 3, qs = 5, qe = 4, preseed = 0){
  paramobj <- rand_mobius_link_cann(p = p, qs = qs, qe = qe, preseed = preseed)
  target_env <- if (is.null(environment(-1))) .GlobalEnv else environment(-1)
  list2env(c(paramobj, list(paramobj = paramobj, qs = qs, qe = qe, p = p)), envir = target_env)
  return(NULL)
}

# Check whether a mean link parameter object uses the "LinEuc" form (linear Euclidean link).
# In LinEuc form, qe1 = (1,0,...) and ce = 1, so the Euclidean part acts as a linear predictor.
is_LinEuc <- function(obj, tol = sqrt(.Machine$double.eps)){
  if (inherits(obj, "mobius_link_Omega")){
    if (length(obj$qe1) > 0){
      checks <- c((obj$qe1 - c(1, rep(0, length(obj$qe1) - 1)))^2, (obj$ce - 1)^2)
      if (all(checks < tol)){return(TRUE)}else{return(FALSE)}
    }
  }
  if (inherits(obj, "mobius_link_cann")){
    if (!is.null(obj$Qe)){
      checks <-  c((obj$Qe[,1] - c(1, rep(0, nrow(obj$Qe) - 1)))^2, (obj$ce[1] - 1)^2)
      if (all(checks < tol)){return(TRUE)}else{return(FALSE)}
    }
  }
  warning("obj doesn't have Euclidean component")
  return(NULL)
}

#' @title Obtain dimensions corresponding to a mean link parameter set
#' @param x A mean link parameter object of class `mobius_link_Omega` or `mobius_link_cann`.
#' @return An integer vector of p (length of response unit vectors), qs (length of spherical covariate unit vectors) and qe (length of Euclidean covariate vectors).
#' @name dim_mobius_link_params
NULL

#' @rdname dim_mobius_link_params
#' @export
#' @method dim mobius_link_cann
dim.mobius_link_cann <- function(x){
  c(p = ncol(x$P), 
    qs = switch(1 + is.null(x$Qs), nrow(x$Qs), 0),
    qe = switch(1 + is.null(x$Qe), nrow(x$Qe), 0))
}

#' @rdname dim_mobius_link_params
#' @export
#' @method dim mobius_link_Omega
dim.mobius_link_Omega <- function(x){
  c(p = length(x$p1), 
    qs = length(x$qs1),
    qe = length(x$qe1))
}
#' Calculate the Mean Given Covariates
#' @description
#' Implements mean link:
#' \deqn{\mu(x) = P\mathcal{S}^{-1}\left(B_s \mathcal{S}(Q_s^\top x_s)  +  \frac{B_e\left(Q_e[,-1]^\top x_e\right)}{Qe[,1]^\top x_e + c_e}\right).}
#' @param xs A matrix of row-vectors of the spherical covariate.
#' @param xe A matrix of row-vectors of the Euclidean covariates.
#' @param param Parameters of the mean link. As an object of class "mobius_link_Omega" or "mobius_link_cann". See [`mobius_link_params`]. 
#' @details
#' This general form of the mean link encompases the primary form of the mean link in "Regression for spherical responses with linear and spherical covariates using a scaled link function" and a more general form that uses the a stereographic-like projection of the Euclidean covariates.
#' See [`mobius_link_params`] for further details.
#'
#' If `param` is of class "mobius_link_Omega" then means are computed as
#' \deqn{\mu(x) = \frac{(1-\|\tilde{y}(x)\|^2) P[,1] + 2 \tilde{y}(x)}{1+\|\tilde{y}(x)\|^2}}
#' where
#' \deqn{\tilde{y}(x) = \frac{\Omega_s x_s}{1+Q_s[,1]^\top x_s} + \frac{\Omega_e x_e}{c_e+{Q_e[,1]}^\top x_e}}
#' and \eqn{x_s} and \eqn{x_e} are the spherical and Euclidean covariate.
#' @family link-function
#' @export
mobius_link <- function(xs = NULL, xe = NULL, param = NULL, check = TRUE){
  if (!is.null(xs)){
    if (inherits(xs, "mobius_link_Omega") | inherits(xs, "mobius_link_cann")){
      stop("xs is a parameter object (mobius_link_Omega or mobius_link_cann), but should be a matrix of covariate values.")
    }
    stopifnot(inherits(xs, "matrix"))}
  if (!is.null(xe)){
    if (inherits(xe, "mobius_link_Omega") | inherits(xe, "mobius_link_cann")){
      stop("xe is a parameter object (mobius_link_Omega or mobius_link_cann), but should be a matrix of covariate values.")
    }
    stopifnot(inherits(xe, "matrix"))}
  if (inherits(param, "mobius_link_Omega")){
    # by C++
    if (is.null(xs)){xs <- matrix(ncol = 0, nrow = nrow(xe))}
    if (is.null(xe)){xe <- matrix(ncol = 0, nrow = nrow(xs))}
    # Checks that xs and xe compatible with param
    stopifnot(ncol(xs) == length(param$qs1))
    stopifnot(ncol(xe) == length(param$qe1))
    # Evaluate
    out <- mobius_link_cpp(xs, xe, mobius_link_Omega_vec(param), length(param$p1))
  } else if (inherits(param, "mobius_link_cann")){
    out <- mobius_link_pred_cann(xs, xe, param)
  } else {
    stop("param is not of the correct class")
  }
  return(out)
}

# Evaluate the mean link using the canonical parameterisation. Follows the SpEuc form
# (Remark 1 of the manuscript), which is slightly more general than the primary form.
mobius_link_pred_cann <- function(xs = NULL, xe = NULL, paramobj){
  stopifnot(inherits(paramobj, "mobius_link_cann"))
  if (!is.null(xs) && !is.null(xe)){stopifnot(nrow(xs) == nrow(xe))}
  y <- matrix(0, nrow = max(nrow(xs), nrow(xe)), ncol = ncol(paramobj$P) - 1)
  
  if (!is.null(paramobj$Qs)){
    # xs stored as row vectors, so xs %*% Qs right-multiplies; Sp is applied row-wise.
    y <- y + Sp(xs %*% paramobj$Qs) %*% paramobj$Bs
  }
  if (!is.null(paramobj$Qe)){
    xetilde <- xe %*% paramobj$Qe #first column is used in denominator
    # xetilde[,1] = Qe[,1]^T xe acts as the denominator in the stereographic-like Euclidean
    # projection; Qe[,1] is the "north-pole" reference direction in covariate space (the
    # singular direction is -Qe[,1]), and ce is an offset ensuring the denominator stays
    # positive for all training points.
    numerator <- (xetilde[, -1, drop = FALSE]) %*% paramobj$Be
    denominator <- xetilde[,1] + paramobj$ce
    y <- y + (numerator/denominator)
  }
  # iSp maps the accumulated displacement y back to S^{p-1}; t(P) rotates the result
  # to the target orientation (transpose because the output is stored as a row vector).
  out <- iSp(y) %*% t(paramobj$P)
  return(out)
}
#' @title Objective function for vMF regression
#' @description
#' Internal objective function used by [`mobius_vMF()`].
#' For given parameters of [`mobius_link()`], computes \deqn{\frac{1}{n}\sum_{i=1}^n y_i^T \mu(x_i)}
#' where \eqn{y_i} are observed unit vectors and \eqn{\mu(x_i)} is the corresponding predicted unit vector.
#' Maximising this is equivalent to maximising the vMF log-likelihood over the mean link parameters.
#' @inheritParams mobius_link
#' @keywords internal
prelimobj <- function(y, xs = NULL, xe = NULL, param){
  predictedmeans <- mobius_link(xs = xs, xe = xe, param = param, check = FALSE)
  stopifnot(nrow(y) == nrow(predictedmeans))
  stopifnot(ncol(y) == ncol(predictedmeans))
  return(-mean(rowSums(y * predictedmeans)))
}

#' @title Regression with vMF Error and Scaled Mobius Link
#' @importClassesFrom scorematchingad Rcpp_ADFun
#' @description Regression of spherical response with vMF error and scaled Mobius link [`mobius_link()`].
#' @details Optimisation of link parameters is by maximising 
#' \deqn{\sum_i=1^n y_i^T \mu(x_i)}
#' where \eqn{y_i} are observed unit vectors and \eqn{\mu(x_i)} is the corresponding predicted unit vector.
#' The concentration of the vMF is then estimated from the residuals.
#' 
#' Uses `nloptr`.
#' 
#' Before fitting, rotates y, xs and xe into standard form. If supplied, `start`, is updated accordingly.
#' 
#' If `type == "LinEuc"` a column of zeros called `'dummyzero'` is added to the front of `xe`.
#' 
#' @param y A matrix of unit vectors, each row corresponds to a single unit vector.
#' @param xs A matrix of spherical covariate vectors (also unit vectors), each row corresponds to the same row in `y`.
#' @param xe A matrix of Euclidean covariate vectors, each row corresponds to the same row in `y`.
#' @param start is a starting parameter object. For LinEuc mean link the Qe matrix must have an extra row and column that at the front/top, with 1 in the first entry (and zero elsewhere).
#' @param ... Passed as options to [`nloptr()`]. 
#' @param intercept `TRUE` to include a Euclidean intercept term using a covariate that is always `1`. This is needed for centering of Euclidean covariates, which is part of standardising the covariates. If `intercept = FALSE` then the Euclidean covariates will not be standardised.
#' @param lb Passed to [`nloptr()`]. Usually not used.
#' @param ub Passed to [`nloptr()`]. Usually not used.
#' @inheritParams mobius_SvMF
#' @family regression
#' @export
mobius_vMF <- function(y, xs = NULL, xe = NULL, start = NULL, type = "SpEuc", fix_qs1 = FALSE, fix_qe1 = (type == "LinEuc"), intercept = TRUE, lb = NULL, ub = NULL, ...){
  p <- ncol(y)
  preplist <- list(y = y, xs = xs, xe = xe, start = start)
  # For LinEuc: prepend a dummy-zero column to xe so the first covariate direction is always
  # fixed at the north pole (required by the LinEuc parameterisation); also append an
  # intercept column of ones if intercept = TRUE. Updates start to match the augmented xe.
  preplist <- addEuccovars(preplist, type = type, intercept = intercept)
  # Rotate y so its empirical mean direction aligns with the north pole; center and whiten xe
  # (if intercept = TRUE). Update start to match the new coordinates. This makes the default
  # starting point (near-identity matrices) a sensible initial guess.
  preplist <- standardise_data(preplist, intercept)
  # If start not supplied, choose start close to identities since data standardised
  preplist <- defaultstart(preplist, type)

  # For LinEuc: verify that qe1 equals the north pole (is_LinEuc) and that the dummy column
  # of xe is numerically zero, both required by the LinEuc parameterisation.
  if ((type == "LinEuc") && (!is.null(preplist$xe))){
    stopifnot(is_LinEuc(preplist$start))
    stopifnot(all(preplist$xe[, 1]^2 < sqrt(.Machine$double.eps)))
  }

  ### More detailed preparation ###
  om0 <- as_mobius_link_Omega(preplist$start)
  # check inputs:
  check_meanlink(preplist$y, preplist$xs, preplist$xe, om0)

  # Set up constraints: tape, fixed-parameter mask, and (possibly perturbed) starting point
  conprep <- estprep_meanconstraints(om0, fix_qs1, fix_qe1)
  # conprep$om0vec          - starting parameters as a vector: c(p1, qs1, qe1, Omega, ce)
  # conprep$isfixed         - logical mask: which positions in om0vec are held fixed
  # conprep$x0              - free-parameter starting values (om0vec[!isfixed]), possibly
  #                           perturbed to avoid a singular constraint Jacobian
  # conprep$constraint_tape - AD tape for orthonormality constraints on P, Qs, Qe
  # Rebuild om0vec incorporating any perturbation in conprep$x0. om0vec then serves as the
  # fixed-value carrier in the later t_sfi2u call: t_sfi2u(nlopt$solution, om0vec, isfixed)
  # reconstructs the full parameter vector after optimisation using om0vec's fixed positions.
  om0vec <- scorematchingad:::t_sfi2u(conprep$x0, conprep$om0vec, conprep$isfixed)

  # Prepare objective tape.
  # tape_namedfun records prelimobj_cpp under CppAD, producing a tape that evaluates the
  # objective and its gradient at any parameter value via $forward(0, x) and $Jacobian(x).
  # om0vec is the independent (differentiable) variable; the data (y, xs, xe) and dimensions
  # (p, qe) are baked in as constants.
  objtape <- tape_namedfun("prelimobj_cpp", om0vec, vector(mode = "numeric"), c(p, length(om0$qe1)), cbind(preplist$y,preplist$xs,preplist$xe), check_for_nan = FALSE)
  objtape <- scorematchingad::avgrange(objtape) #objtape initially returns a value for each measurement. Average here to get average over all data.
  
  # Re-tape the objective with fixed parameters baked in as constants, reducing the tape's
  # domain from length(om0vec) to sum(!isfixed) (free parameters only).
  objtape <- scorematchingad::fixindependent(objtape, objtape$xtape, conprep$isfixed)
  # Use the taping values as the starting point for optimisation. These are the (possibly
  # perturbed) starting values baked into the tape by estprep_meanconstraints.
  x0 <- objtape$xtape
  
  # prepare nloptr options
  # NLOPT_LD_SLSQP: Sequential Least Squares Programming; gradient-based constrained
  # optimiser suited here because the AD tapes supply exact gradients at low cost.
  # tol_constraints_eq = 1E-1: constraints (orthonormality) are satisfied only approximately
  # during optimisation; exact orthonormality is restored afterwards by Omega_proj.
  default_opts <- list(algorithm = "NLOPT_LD_SLSQP",
                xtol_rel = 1E-10, #1E-04,
                tol_constraints_eq = rep(1E-1, conprep$constraint_tape$range),
                # check_derivatives = TRUE, check_derivatives_print = 'errors', check_derivatives_tol = 1E-3,
                # print_level = 3,
                maxeval = 1E4)
  ellipsis_args <- list(...)
  combined_opts <- utils::modifyList(default_opts, ellipsis_args)
  
  # activate a progress bar
  pb <- progress::progress_bar$new(total = combined_opts$maxeval + 5, format = ":bar :percent :current :tick_rate elapsed::elapsedfull eta::eta")
  
  # Optimisation
  # current dynamic parameter values of tapes will be used
  # nloptr minimises; the objective sum_i y_i . mu(x_i) is to be maximised, hence the negation.
  nlopt <- nloptr::nloptr(
    x0 = x0,
    eval_f = function(theta){
      if (!pb$finished) pb$tick()
      list(objective = -objtape$forward(0, theta), gradient = -objtape$Jacobian(theta))
      },
    eval_g_eq =  function(theta){
      list(constraints = conprep$constraint_tape$forward(0, theta),
           jacobian = matrix(conprep$constraint_tape$Jacobian(theta), byrow = TRUE, ncol = length(theta)))},
    opts = combined_opts,
    lb = lb,
    ub = ub
  )
  if (!(nlopt$status %in% 1:4)){warning(nlopt$message)}
  
  # Estimate concentration k by a separate 1D search.
  # Given the optimised mean link, the per-observation vMF log-likelihood is
  #   -log_norm_const(k, p) + k * mean(y . ypred).
  # The first term depends only on k and p; the second is k times the already-maximised
  # objective (-nlopt$objective = mean(y . ypred)). So k can be found by a scalar optimise
  # over k, treating mean(y . ypred) as fixed.
  res <- optimise(function(k){
    -vMF_log_norm_const_exact(k, p) + k * (-nlopt$objective) #full vMF log-likelihood (standardised by number of observations)
  }, lower = 1E-8, upper = 1E5, maximum = TRUE)
  k <- res$maximum
  
  #output some diagnostics - vector names would be nice here
  nlopt$solution_grad_f <- -objtape$Jacobian(nlopt$solution)
  nlopt$solution_g_eq <- conprep$constraint_tape$forward(0, nlopt$solution)
  nlopt$solution_jac_g_eq <- matrix(conprep$constraint_tape$Jacobian(nlopt$solution),
                                     byrow = TRUE, ncol = length(nlopt$solution))
  nlopt$solution_Hes_f <- matrix(-objtape$Hessian0(nlopt$solution),
         nrow = objtape$domain,
         byrow = TRUE)

  # remove the tapes from the return to save on memory
  nlopt$eval_f <- nlopt$eval_g_eq <- nlopt$nloptr_environment <- NULL
   
  # insert any fixed values of mean parameters
  fullparam <- scorematchingad:::t_sfi2u(nlopt$solution, om0vec, conprep$isfixed)
 
  # SLSQP satisfies orthonormality constraints only up to tol_constraints_eq = 1E-1;
  # Omega_proj projects P, Qs, Qe back onto the exact Stiefel manifold.
  projectedom <- Omega_proj(mobius_link_Omega_unvec(fullparam, p, length(om0$qe1), check = FALSE))
  try({mobius_link_Omega_check(projectedom)})

  ### Making nicer return objects ###
  # Residuals and distances are computed in the standardised coordinate system (before
  # reverting), since they are invariant to the choice of coordinates.
  pred <- mobius_link(xs = preplist$xs, xe = preplist$xe, param = projectedom)
  dists <- acos(rowSums(pred * preplist$y))
  # rotated_resid parallel-transports each observation to the north pole. The returned matrix
  # has p columns: the first is the (normalised) response direction (always ~1 for small
  # residuals). [, -1] drops it, keeping only the p-1 tangent components.
  rresids_tmp <- rotated_resid(preplist$y, pred, north_pole(ncol(preplist$y)))
  rresids <- rresids_tmp[, -1]
  attr(rresids, "samehemisphere") <-  attr(rresids_tmp, "samehemisphere")
  colnames(rresids) <- paste0("r", 1:ncol(rresids))
  
  # Invert the rotations and centering applied by standardise_data, expressing the estimated
  # parameters in the original user-supplied coordinates.
  est <- undo_recoordinate_Omega(projectedom,
                          yrot = attr(preplist$y, "std_rotation"),
                          xsrot = attr(preplist$xs, "std_rotation"), #if xs/xe is NULL then attr(xs/xe, ..) is NULL too
                          xerot = attr(preplist$xe, "std_rotation"),
                          xecenter = attr(preplist$xe, "std_center"),
                          onescovaridx = preplist$onescovaridx)
  # The Euclidean link has an inherent sign ambiguity: replacing (ce, Qe, Be) with
  # (-ce, -Qe, -Be) leaves the mean link unchanged. Convention is ce >= 0; Euc_signswitch
  # flips all three if ce < 0.
  if (isTRUE(est$ce < 0)){est <- Euc_signswitch(est)}
  
  # DoF
  DoF <- mobius_dof(p, length(est$qs1), length(est$qe1), fix_qs1 = fix_qs1, fix_qe1 = fix_qe1) + 
    1 #concentration
  # AIC using result from concentration search
  AIC <- 2*DoF - 2 * res$objective * nrow(y)
  lLik <- res$objective * nrow(y)
  
  niceout <- list(
    est = est,
    mean = est,
    k = k,
    obj = nlopt$objective,
    solution = projectedom, #non-standardised solution
    nlopt = nlopt,
    y = y,
    xs = xs,
    # destandardise_Euc inverts the centering and whitening from standardise_data; using
    # preplist$xe also recovers any columns added by addEuccovars (dummy-zero and intercept).
    xe = if (!is.null(xe)){if (intercept){destandardise_Euc(preplist$xe, attr(preplist$xe, "std_center"), attr(preplist$xe, "std_rotation"))} else {xe}},
    pred = destandardise_sph(pred, rotation = attr(preplist$y, "std_rotation")),
    rresids = rresids,
    dists = dists,
    DoF = DoF,
    AIC = AIC,
    lLik = lLik,
    start = start,
    linktype = list(type = type, fix_qs1 = fix_qs1, fix_qe1 = fix_qe1, intercept = intercept)
  )
  return(niceout)
}

#' @title Degrees of freedom of the Mobius mean link function
#' @description The parameters of the Mobius mean link function [`mobius_link()`] have a number of constraints.
#' This function incorporates these constraints to obtain the total degrees of freedom of the parameters.
#' @param p The length of response vectors
#' @param qs The length of spherical covariate vectors
#' @param qe The length of Euclidean covariate vectors
#' @param fix_qs1 Whether the first column of `Qs` is fixed (i.e. not estimated and not free).
#' @param fix_qe1 Whether `ce` and the first column of `Qe` is fixed (i.e. not estimated and not free).
#' @return An integer
#' @family regression
#' @export
mobius_dof <- function(p, qs = 0, qe = 0, fix_qs1 = FALSE, fix_qe1 = FALSE){
  if (qs == 0){fix_qs1 <- FALSE} #ignore fix_qs1
  if (qe == 0){fix_qe1 <- FALSE} #ignore fix_qs1
  DoF <- DoF_Stiefel(p, p) + #P
    DoF_Stiefel(qs-fix_qs1, p-fix_qs1) + #Qs
    DoF_Stiefel(qe-fix_qe1, p-fix_qe1) + #Qe
    (qs>0)*(p-1) + #Bs
    (qe>0)*(p-1) + #Be
    1*((qe>0) & (!fix_qe1)) #ce
}

check_meanlink <- function(y, xs, xe, om0){
  try(mobius_link_Omega_check(om0))
  p <- ncol(y)
  stopifnot(p == length(om0$p1))
  if (!is.null(xs)){
    stopifnot(ncol(xs) == length(om0$qs1))
  } else {
    stopifnot(length(om0$qs1) == 0)
  }
  # check Euc info if Euc is part of the link
  if (!is.null(xe)){
    stopifnot(ncol(xe) == length(om0$qe1))
  } else {
    stopifnot(length(om0$qe1) == 0)
  }
}

#' @title Prepare mean-link constraint setup for optimisation
#' @description
#' Prepares the constraint tape, fixed-parameter mask, and starting point for constrained
#' optimisation of the Mobius mean link parameters in [`mobius_vMF()`].
#' Specifically, the function:
#' 1. Vectorises `om0` to `om0vec`, the starting parameters as `c(p1, qs1, qe1, Omega, ce)`.
#' 2. Creates an AD constraint tape evaluating the orthonormality constraints on `P`, `Qs`,
#'    `Qe` and their Jacobians.
#' 3. Determines `isfixed` — a logical mask of which positions in `om0vec` are held constant,
#'    based on `fix_qs1` and `fix_qe1`.
#' 4. Reduces the tape via `fixindependent` and `keeprange`, dropping constraint outputs that
#'    become constant once fixed parameters are removed.
#' 5. Computes `x0 = om0vec[!isfixed]`, the free-parameter starting vector.
#' 6. Checks the constraint Jacobian at `x0` for singularity (near-zero singular values).
#' 7. If singular, perturbs `x0` by adding small noise (from a pregenerated matrix in
#'    `R/sysdata.rda`) to the off-diagonal columns of `Qs` and `Qe`, then recomputes `x0`.
#' @param om0 Starting [`mobius_link_Omega`] object.
#' @param fix_qs1 Logical; whether to hold the first column of `Qs` fixed.
#' @param fix_qe1 Logical; whether to hold `ce` and the first column of `Qe` fixed.
#' @return A list:
#' \describe{
#'   \item{om0vec}{Starting parameters as a flat vector: `c(p1, qs1, qe1, Omega, ce)`.}
#'   \item{isfixed}{Logical vector marking fixed positions in `om0vec`.}
#'   \item{x0}{Starting values for the free (non-fixed) parameters, possibly perturbed to
#'     avoid a singular constraint Jacobian.}
#'   \item{constraint_tape}{AD tape evaluating orthonormality constraints on `P`, `Qs`, `Qe`
#'     and their Jacobians, with fixed-parameter dimensions already reduced.}
#' }
#' @keywords internal
estprep_meanconstraints <- function(om0, fix_qs1, fix_qe1){
  dims_in <- c(length(om0$p1), length(om0$qe1))
  om0vec <- mobius_link_Omega_vec(om0)
 
  # generate tapes
  # tape_namedfun records Omega_constraints_wrap under CppAD: evaluates the orthonormality
  # constraints and their Jacobian at any om0vec via $forward(0, x) and $Jacobian(x).
  constraint_tape <- tape_namedfun("Omega_constraints_wrap", om0vec, vector(mode = "numeric"), dims_in, matrix(nrow = 0, ncol = 0), check_for_nan = FALSE)

  # fix mean link parameters depending on arguments
  # use the starting parameters om0 to detect whether we have xs and xe as their form is more predictable due to the Omega class
  omfixed <- lapply(om0, function(x) x * 0)
  if (fix_qe1 && (length(om0$qe1) > 0)){
    # if qe1 is fixed and qe1 is a non-empty vector, then record the corresponding omfixed elements as 1
    omfixed$qe1 <- omfixed$qe1 + 1
    omfixed$ce <- omfixed$ce + 1
  }
  
  # fix qs1 to the starting values if qs1 exists
  # use the starting parameters to check as their form is more constrained by the Omega class
  if (fix_qs1 && (length(om0$qs1) > 0)){
    omfixed$qs1 <- omfixed$qs1 + 1
  }
  
  # Vectorise the omfixed marker through the Omega class so the logical mask aligns exactly
  # with the positions in om0vec.
  isfixed <- mobius_link_Omega_vec(as_mobius_link_Omega(omfixed)) > 0.5
  # Bake the fixed-parameter values into the constraint tape, reducing its domain to free
  # parameters only (so the constraint Jacobian is over free parameters only).
  constraint_tape <- scorematchingad::fixindependent(constraint_tape, om0vec, isfixed)
  # After fixing parameters, some constraint equations become identically zero (e.g. the
  # orthonormality constraint on a fixed column is pre-satisfied). Drop those to keep the
  # constraint Jacobian full-rank.
  keep <- which(vapply((1:constraint_tape$range)-1, function(i){!constraint_tape$parameter(i)}, FUN.VALUE = FALSE))
  constraint_tape <- scorematchingad::keeprange(constraint_tape, keep)
  
  # check Jacobians of constraints are non-singular for the starting parameters.
  # For pathological params (e.g. the default starting params of no rotations), it can be zero.
  # If it is singular, perturb starting omega very slightly
  x0 <- om0vec[!isfixed]
  Jac_eq <- matrix(constraint_tape$Jacobian(x0), byrow = TRUE, ncol = constraint_tape$domain)
  if (any(abs(svd(Jac_eq)$d) < sqrt(.Machine$double.eps))){
    # modify Qs and Qe so that they dont individually form orthogonal vectors by adding a little bit of noise
    cann <- as_mobius_link_cann(om0)
    scalemodifier <- max(abs(rbind(cann$Qs, cann$Qe)))/1E3
    cann$Qs[,-1] <- cann$Qs[,-1] + scalemodifier*savednoisemat[1:nrow(cann$Qs), 1:(ncol(cann$Qs)-1)] #use pregenerated 100 x 100 noise matrix in R/sysdata.rda, built by data-raw/generate-savednoisemat.R
    cann$Qe[,-1] <- cann$Qe[,-1] + scalemodifier*savednoisemat[1:nrow(cann$Qe), 1:(ncol(cann$Qe)-1)] #use pregenerated 100 x 100 noise matrix in R/sysdata.rda, built by data-raw/generate-savednoisemat.R
    om0new <- as_mobius_link_Omega(cann)
    x0 <- mobius_link_Omega_vec(om0new)[!isfixed]
    }
  Jac_eq <- matrix(constraint_tape$Jacobian(x0), byrow = TRUE, ncol = constraint_tape$domain)
  if (any(abs(svd(Jac_eq)$d) < sqrt(.Machine$double.eps))){
    warning("Initial parameters lead to a singular constraint Jacobian")
  }
  
  return(list(
    om0vec = om0vec,
    x0 = x0, #x0 may be perturbed to avoid singular Jac_eq 
    isfixed = isfixed,
    constraint_tape = constraint_tape
    ))
}

#' @title Try other initial starting parameters for a given regression
#' @description Given a vMF regression, repeat the optimisation from initial parameters randomly generated by [`rand_mobius_link_cann()`].
#' @param mod_vMF Result of [`mobius_vMF()`]
#' @param preseed Passed to [`rand_mobius_link_cann()`]
#' @details `fix_qe1` and `fix_qs1` of `mod_vMF` will be respected.
#' @return An vMF regression using the same data as `mod_vMF`.
#' @family regression
#' @export
mobius_vMF_refit <- function(mod_vMF, preseed = 1){
  inparam <- as_mobius_link_Omega(mod_vMF$est)
  dims <- dim.mobius_link_Omega(inparam)
  start <- rand_mobius_link_cann(p = dims["p"], 
                           qs = dims[["qs"]] - (dims[["qs"]] > 0 & mod_vMF$linktype$fix_qs1), 
                           qe = dims[["qe"]] - (dims[["qe"]] > 0 & mod_vMF$linktype$fix_qe1), 
                           preseed = preseed)
  # Goal: for each fixed first column (qe1 or qs1), construct a random starting matrix whose
  # first column equals the fixed value and whose remaining columns are random in the
  # orthogonal complement. Strategy:
  #   1. rand_mobius_link_cann already generated a random matrix for the free columns only.
  #   2. Prepend a zero placeholder column to form the full p x qe (or p x qs) matrix.
  #   3. Apply the inverse parallel transport from the north pole to qe1 (or qs1), mapping
  #      the placeholder first column to qe1 and rotating the free columns accordingly.
  #   4. Overwrite the first column with the exact fixed value (guards against the edge case
  #      where qe1 = -north_pole and parallel transport is a reflection, not a rotation).
  if (dims[["qe"]] > 0 & mod_vMF$linktype$fix_qe1){
    start$ce <- inparam$ce
    rotmat <- parallel_transport_mat(inparam$qe1, north_pole(dims[["qe"]]))
    rotQe <- rbind(0, start$Qe)
    # Two cases depending on whether rand_mobius_link_cann returned a square or tall matrix
    if (ncol(rotQe) == dims["p"]){
      rotQe <- cbind(0, rotQe[,-1, drop = FALSE])
    } else if (ncol(rotQe) == dims["p"] - 1) { #case when qe=p and fix_qe1=TRUE
      rotQe <- cbind(0, rotQe)
    }
    rotQe[1,1] <- 1
    start$Qe <- t(rotmat) %*% rotQe
    # Overwrite first column: parallel transport may be a reflection when qe1 = -north_pole
    start$Qe[,1] <- inparam$qe1
  }
  if (dims[["qs"]] > 0 & mod_vMF$linktype$fix_qs1){
    rotmat <- parallel_transport_mat(inparam$qs1, north_pole(dims[["qs"]]))
    rotQs <- rbind(0, start$Qs)
    if (ncol(rotQs) == dims["p"]){
      rotQs <- cbind(0, rotQs[,-1, drop = FALSE])
    } else if (ncol(rotQs) == dims["p"] - 1) { #case when qs=p and fix_qs1=TRUE
      rotQs <- cbind(0, rotQs)
    }
    rotQs[1,1] <- 1
    start$Qs <- t(rotmat) %*% rotQs
    # Overwrite first column: parallel transport may be a reflection when qs1 = -north_pole
    start$Qs[,1] <- inparam$qs1
  }
  start <- as_mobius_link_Omega(start)
  stopifnot(all(dim.mobius_link_Omega(start) == dims))
  mobius_vMF(y = mod_vMF$y,
             xs = mod_vMF$xs,
             xe = mod_vMF$xe,
             fix_qs1 = mod_vMF$linktype$fix_qs1,
             fix_qe1 = mod_vMF$linktype$fix_qe1,
             type = mod_vMF$linktype$type,
             intercept = mod_vMF$linktype$intercept,
             start = start)
}
#' SvMF Regression with Scaled Mobius Link and Parallel Transported Axes
#' @description
#' Performs regression for spherical response data with a scaled Mobius mean link and scale von Mises-Fisher error distribution (SvMF).
#' @details
#' The mean is assumed to follow [`mobius_link()`].
#' The concentration and scaling in the SvMF is assumed constant across observations.
#' The symmetry axes of the SvMF at location \eqn{\mu} are assumed to be the result of parallel transport along the geodesic of axes `G0[,-1]` from  `G0[,1]` to \eqn{\mu}. The axis at this base location `G0[,1]` are estimated.
#'
#' Unless requested by `doprelim = FALSE`, a preliminary regression using a vMF error distribution is performed.
#'
#' When p!=3, the optimisation first approximates the normalising constant of the vMF distribution, then the concentration is optimised seperately.
#' @param y Response data on a sphere
#' @param xs Covariate data on a sphere
#' @param xe Covariate data in Euclidean Space
#' @param k Starting concentration. If omitted, then concentration of the preliminary vMF regression is used.
#' @param a The scaling vector `a`. `a[1]` is a fixed tuning parameter and the remainining is used as an initial value for the optimisation. If omitted, then a moment estimate using residuals of from the preliminary vMF regression is used.
#' @param mean Parameters for the mean link, used as a starting guess for the preliminary vMF regression if `doprelim = TRUE`.
#' @param G0 A `p x p` rotation matrix specifying the starting guess of the axes of the SvMF distribution. G0 should have positive determinant because in the estimation routine G0 or parts of G0 are representented using Cayley transforms. If omitted, then a moment estimate using residuals of from the preliminary vMF regression is used.
#' @param G0reference A `p x p` rotation matrix specifying a set of coordinates with which to represent G0 for the estimation. This is because internally `G0` or `G0[,-1]` is represented using a Cayley transform, which performs better the closer the matrix is to the identity. Ideally the columns of `G0reference` will be close to the best `G0`.
#' @param G01behaviour "p1" identifies `G0[,1]` with `p1`. "fixed" fixes `G0[,1]` to its initial value. "free" allows `G0[,1]` to be estimated freely.
#' @param doprelim When `FALSE` the preliminary von Mises-Fisher regression and subsequent moment estimation of `a` and `G0` is omitted. The provided parameters are used as the initial values for an optimisation of all parameters of the SvMF regression all together using [`nloptr::nloptr()`].
#' @param ... Named optional arguments passed as a list to the `opts` argument of [`nloptr::nloptr()`].
#' @param type Specify the link type ("SpEuc" or "LinEuc")
#' @param fix_qs1 Fix qs1 to the starting values.
#' @param fix_qe1 Fix qe1 to the starting values - typically used in when "LinEuc" link used.
#' @param intercept Include an intercept (Euclidean covariate that is always `1`).
#' @return A list:
#' \describe{
#'   \item{mean}{Estimated mean parameters}
#'   \item{k}{Estimated concentration parameter k}
#'   \item{a}{Estimated vector of scales `a` (except `a[1]` which is fixed)}
#'   \item{G0}{Estimated base location and symmetry axes.}
#'   \item{obj}{Final objective value from NLopt}
#'   \item{nlopt}{NLopt optimization result object}
#'   \item{y}{Response data}
#'   \item{xs}{Spherical covariate values}
#'   \item{xe}{Euclidean covariate values (including any automatically added columns like dummy zeros and an intercept)}
#'   \item{pred}{Predicted values on sphere}
#'   \item{rresids_I}{Residuals rotated (by parallel transport) to the north pole}
#'   \item{rresids_G0, rresids_std}{Residuals rotated to `G0[,1]` and expressed according to the axes given by `G0[,-1]`. `rresids_std` additionaly scales the rotated residuals by \eqn{\sqrt{k} a[1]/a[j]}.}
#'   \item{dists}{Residual geodesic distances (geodesic distance between predicted value and observed value)}
#'   \item{DoF}{Degrees of freedom of the optimisation}
#'   \item{AIC}{Akaike Information Criterion}
#'   \item{lLik}{Log-likelihood}
#'   \item{initial}{Initial values}
#'   \item{preest}{Results of the preliminary vMF regression}
#' }
#' @family regression
#' @export
mobius_SvMF <- function(y, xs = NULL, xe = NULL, mean = NULL, k = NULL, a = NULL, G0 = NULL, G0reference = NULL, G01behaviour = "p1", type = "LinEuc", fix_qs1 = FALSE, fix_qe1 = (type == "LinEuc"), intercept = TRUE, doprelim = TRUE, ...){

  # Two-phase estimation:
  # Phase 1 (prelim): vMF regression for the mean link, followed by moment estimation of G0
  # and a. This provides a warm start for the joint optimisation.
  # Phase 2 (finalest): joint maximum likelihood over all parameters (mean, k, a, G0).
  if (doprelim){
    preest <- mobius_SvMF_partransport_prelim(y, xs, xe,
                                              mean = mean,
                                              G0 = G0, G01behaviour = G01behaviour,
                                              type = type, fix_qs1 = fix_qs1, fix_qe1 = fix_qe1,
                                              intercept = intercept, ...)
  } else {
    # doprelim = FALSE: caller must supply all starting values directly.
    stopifnot(!is.null(mean))
    stopifnot(!is.null(k))
    stopifnot(!is.null(a))
    stopifnot(!is.null(G0))
    preest <- list(
      mean = mean,
      k = k,
      a = a,
      G0 = G0)
  }

  # Explicit caller-supplied k, a, G0reference take precedence over the prelim estimates.
  finalest <- mobius_SvMF_joint_fit(y, xs, xe,
                           mean = preest$mean,
                           k = if(!is.null(k)){k}else{preest$k},
                           a = if(!is.null(a)){a}else{preest$a},
                           G0 = preest$G0,
                           G0reference = if(!is.null(G0reference)){G0reference}else{preest$G0},
                           G01behaviour = G01behaviour,
                           type = type, fix_qs1 = fix_qs1, fix_qe1 = fix_qe1,
                           intercept = intercept, ...)
  return(c(finalest, list(preest = preest)))
}

# Joint maximum likelihood optimisation of all SvMF regression parameters for fixed
# scale structure (a and G0 are estimated, but the constraint prod(a[-1]) = 1 is enforced).
# Called by mobius_SvMF() after the preliminary vMF estimate. Uses automatic differentiation
# (CppAD via scorematchingad) and a gradient-based solver (nloptr/SLSQP).
mobius_SvMF_joint_fit <- function(y, xs, xe, mean, k, a, G0 = NULL, G0reference = NULL, G01behaviour = "p1", type = "LinEuc", fix_qs1 = FALSE, fix_qe1 = (type == "LinEuc"), intercept = TRUE, lb = NULL, ub = NULL, ...){
  initial <-  list(
    mean = mean,
    k = k, 
    a = a, 
    G0 = G0
  )
  SvMF_cann_check(SvMF_cann(k = initial$k, a = initial$a, G = initial$G0))

  p <- ncol(y)
  preplist <- list(y = y, xs = xs, xe = xe, start = mean)
  # if needed, add Euclidean covariates and update start accordingly
  preplist <- addEuccovars(preplist, type = type, intercept = intercept)
  # standardise y, xe and xe and update start accordingly. Dont standardise xe if intercept = FALSE
  preplist <- standardise_data(preplist, intercept)
  # G0 lives in the same space as y, so it must be rotated by the same std_rotation used on y.
  if (!is.null(G0)){preplist$G0 <- attr(preplist$y, "std_rotation") %*% G0}
  if (!is.null(G0reference)){preplist$G0reference <- attr(preplist$y, "std_rotation") %*% G0reference}
  # If start not supplied, choose start close to identities since data standardised
  preplist <- defaultstart(preplist, type)
  # Caller-supplied a and k take precedence over any defaults set by defaultstart.
  if (!is.null(a)){preplist$a <- a}
  if (!is.null(k)){preplist$k <- k}

  # Check LinEuc link initiated properly
  if ((type == "LinEuc") && (!is.null(preplist$xe))){
    stopifnot(is_LinEuc(preplist$start))
    stopifnot(all(preplist$xe[, 1]^2 < sqrt(.Machine$double.eps)))
  }

  ### More detailed preparation ###
  om0 <- as_mobius_link_Omega(preplist$start)
  # check inputs:
  check_meanlink(preplist$y, preplist$xs, preplist$xe, om0)
  stopifnot(length(preplist$a) == p)
  a1 <- preplist$a[1]
  aremaining <- preplist$a[-1]
  stopifnot(isTRUE(all.equal(prod(aremaining), 1)))
  if (is.null(preplist$G0reference)){
    G0reference <- preplist$G0
  } else {
    G0reference <- preplist$G0reference
  }
   
  qs <- length(om0$qs1)
  qe <- length(om0$qe1)


  # Set up constraints: tape, fixed-parameter mask, and (possibly perturbed) starting point
  conprep <- estprep_meanconstraints(om0, fix_qs1, fix_qe1)
  # conprep$om0vec          - starting parameters as a vector: c(p1, qs1, qe1, Omega, ce)
  # conprep$isfixed         - logical mask: which positions in om0vec are held fixed
  # conprep$x0              - free-parameter starting values (om0vec[!isfixed]), possibly
  #                           perturbed to avoid a singular constraint Jacobian
  # conprep$constraint_tape - AD tape for orthonormality constraints on P, Qs, Qe
  # Rebuild om0vec incorporating any perturbation in conprep$x0. om0vec then serves as the
  # fixed-value carrier in the later t_sfi2u call: t_sfi2u(meanpars, om0vec, isfixed)
  # reconstructs the full mean-link parameter vector after optimisation.
  om0vec <- scorematchingad:::t_sfi2u(conprep$x0, conprep$om0vec, conprep$isfixed)
  
  # Prepare objective tape.
  # tape_ld_mobius_SvMF_partransport_nota1 records the full SvMF log-density under CppAD.
  # Independent variables (differentiated w.r.t.): free mean-link parameters from om0vec,
  # followed by k, aremaining, and the Cayley-transform coordinates of G0 — all in one flat
  # vector. The Cayley transform maps the rotation manifold to an unconstrained space, so G0
  # requires no extra constraints in the optimiser.
  # Constants baked in: data (y, xs, xe), dimensions, G0reference, G01behaviour.
  objtape_ind <- tape_ld_mobius_SvMF_partransport_nota1(omvec = om0vec, k = preplist$k,
                                       a1 = a1,
                                       aremaining = aremaining,
                                       G0 = preplist$G0,
                                       p = p, qe = qe,
                                       yx = cbind(preplist$y, preplist$xs, preplist$xe),
                                       referencecoords = G0reference,
                                       G01behaviour = G01behaviour)
  # Re-tape with fixed mean-link parameters baked in. The isfixed vector is padded with zeros
  # for the non-mean-link components (k, aremaining, Cayley G0) since those are always free.
  objtape_ind <- scorematchingad::fixindependent(objtape_ind, objtape_ind$xtape, c(conprep$isfixed, rep(0, length(objtape_ind$xtape) - length(conprep$isfixed))))

  # objtape initially returns a value for each measurement. Average here to get average over all data.
  objtape <- scorematchingad::avgrange(objtape_ind)
  # Use the taping values as the starting point for optimisation.
  x0 <- objtape$xtape

  # prepare nloptr options
  # NLOPT_LD_SLSQP: Sequential Least Squares Programming; gradient-based constrained
  # optimiser suited here because the AD tapes supply exact gradients at low cost.
  # tol_constraints_eq = 1E-1: constraints (orthonormality) are satisfied only approximately
  # during optimisation; exact orthonormality is restored afterwards by Omega_proj.
  default_opts <- list(algorithm = "NLOPT_LD_SLSQP",
                       xtol_rel = 1E-10, #1E-04,
                       tol_constraints_eq = rep(1E-1, conprep$constraint_tape$range),
                       # check_derivatives = TRUE, check_derivatives_print = 'errors', check_derivatives_tol = 1E-3,
                       # print_level = 3,
                       maxeval = 1E4)
  ellipsis_args <- list(...)
  combined_opts <- utils::modifyList(default_opts, ellipsis_args)

  # The free mean-link parameters occupy positions 1:length(conprep$x0) in the flat vector;
  # position length(conprep$x0) + 1 is k. Enforce k > 0 as a lower bound.
  if (is.null(lb)){
    lb <- rep(-Inf, objtape$domain)
    lb[length(conprep$x0) + 1] <- 0
  }

  # activate a progress bar
  pb <- progress::progress_bar$new(total = combined_opts$maxeval + 5, format = ":bar :percent :current :tick_rate elapsed::elapsedfull eta::eta")
  
  
  # Optimisation
  # current dynamic parameter values of tapes will be used
  # nloptr minimises; the SvMF log-likelihood is to be maximised, hence the negation.
  nlopt <- nloptr::nloptr(
    x0 = x0,
    eval_f = function(theta){
      if (!pb$finished) pb$tick()
      list(objective = -objtape$forward(0, theta), gradient = -objtape$Jacobian(theta))
      },
    lb = lb,
    ub = ub,
    eval_g_eq =  function(theta){list(
      # Orthonormality constraints apply only to the mean-link parameters (first
      # constraint_tape$domain elements of theta). The zero padding gives zero Jacobian
      # for k, aremaining, and the Cayley G0 components, which are unconstrained.
      constraints = conprep$constraint_tape$forward(0, theta[1:conprep$constraint_tape$domain]),
      jacobian = cbind(matrix(conprep$constraint_tape$Jacobian(theta[1:conprep$constraint_tape$domain]), byrow = TRUE, ncol = conprep$constraint_tape$domain),
             matrix(0, nrow = conprep$constraint_tape$range, ncol = length(theta) - conprep$constraint_tape$domain))
    )},
    opts = combined_opts
  )
  if (!(nlopt$status %in% 1:4)){warning(nlopt$message)}

  #output some diagnostics - vector names would be nice here
  nlopt$solution_grad_f <- -objtape$Jacobian(nlopt$solution)
  nlopt$solution_g_eq <- conprep$constraint_tape$forward(0, nlopt$solution[1:conprep$constraint_tape$domain])
  nlopt$solution_jac_g_eq <- matrix(conprep$constraint_tape$Jacobian(nlopt$solution[1:conprep$constraint_tape$domain]),
                                     byrow = TRUE, ncol = length(nlopt$solution[1:conprep$constraint_tape$domain]))
  nlopt$solution_Hes_f <- matrix(-objtape$Hessian0(nlopt$solution),
         nrow = objtape$domain,
         byrow = TRUE)
  nlopt$ldens <- drop(objtape_ind$forward(0,nlopt$solution))
  lLik <- sum(nlopt$ldens)

  # remove the tapes from the return to save on memory
  nlopt$eval_f <- nlopt$eval_g_eq <- nlopt$eval_g_ineq <- nlopt$nloptr_environment <- NULL

  # The solution vector is partitioned as [free mean-link | k | aremaining | Cayley G0].
  # Reinsert fixed mean-link values (t_sfi2u uses om0vec's fixed positions as the template).
  meanpars <- nlopt$solution[1:length(conprep$x0)]
  meanpars <- scorematchingad:::t_sfi2u(meanpars, om0vec, conprep$isfixed)
  fullparam <- c(meanpars, nlopt$solution[-(1:length(conprep$x0))])

  # Reconstruct the full typed parameter list (om, k, aremaining, G0) from the flat vector,
  # applying the inverse Cayley transform to recover G0 from its unconstrained coordinates.
  estparamlist <- mobius_SvMF_partransport_nota1_fromvecparams_forR(fullparam, p, qs, qe,
                                                  referencecoords = G0reference,
                                                  G01behaviour = G01behaviour,
                                                  G01 = preplist$G0[,1])
  
  #project mean pars to have correct orthogonality
  projectedom <- Omega_proj(mobius_link_Omega_unvec(estparamlist$omvec, p, length(om0$qe1), check = FALSE))
  try({mobius_link_Omega_check(projectedom)})

  # The model is identifiable only up to permutation of (a_j, G0[:,j]) pairs. By convention
  # the scales are returned in decreasing order with the G0 columns reordered to match.
  aord <- order(estparamlist$aremaining, decreasing = TRUE)
  aremaining <- estparamlist$aremaining[aord]
  estparamlist$G0[,-1] <- estparamlist$G0[,-1][, aord]

  pred <- mobius_link(xs = preplist$xs, xe = preplist$xe, param = projectedom)
  if (p!=3){
    # When p != 3 the joint AD optimisation uses a numerical approximation to the vMF
    # normalising constant (required for AD differentiability). After convergence,
    # mobius_SvMF_konly refines k by an exact univariate search using the more accurate
    # Bessel function implementation in base R, correcting for the approximation error in k.
    result <- mobius_SvMF_konly(y = preplist$y, ymean = pred, a = c(a1, aremaining), G0 = estparamlist$G0)
    estparamlist$k <- result$k
    lLik <- result$lLik
  }
  
  ### Making nicer return objects ###
  # Residuals and distances are computed in the standardised coordinate system (before
  # reverting), since they are invariant to the choice of coordinates.
  dists <- acos(rowSums(pred * preplist$y))
  # rresids_std: parallel-transported to G0[:,1] and scaled by sqrt(k)*a[1]/a[j] —
  # under high concentration these are approximately i.i.d. standard normal.
  rresids_std <- resid_SvMF_partransport(preplist$y, pred, estparamlist$k, c(a1, estparamlist$aremaining), estparamlist$G0, scale = TRUE)
  # rresids_G0: same parallel transport to G0[:,1] but unscaled.
  rresids_G0 <- resid_SvMF_partransport(preplist$y, pred, G0 = estparamlist$G0, scale = FALSE)
  # rresids_I: parallel-transported to the north pole; [, -1] drops the response component
  # (always ~1 for small residuals), keeping only the p-1 tangent components.
  rresids_I_tmp <- rotated_resid(preplist$y, pred, north_pole(ncol(preplist$y)))
  rresids_I <- rresids_I_tmp[, -1]
  attr(rresids_I, "samehemisphere") <-  attr(rresids_I_tmp, "samehemisphere")
  colnames(rresids_I) <- paste0("r", 1:ncol(rresids_I))

  # Invert the rotations and centering applied by standardise_data, expressing the estimated
  # mean-link parameters in the original user-supplied coordinates.
  est <- undo_recoordinate_Omega(projectedom,
                          yrot = attr(preplist$y, "std_rotation"),
                          xsrot = attr(preplist$xs, "std_rotation"), #if xs/xe is NULL then attr(xs/xe, ..) is NULL too
                          xerot = attr(preplist$xe, "std_rotation"),
                          xecenter = attr(preplist$xe, "std_center"),
                          onescovaridx = preplist$onescovaridx)
  # The Euclidean link has an inherent sign ambiguity: replacing (ce, Qe, Be) with
  # (-ce, -Qe, -Be) leaves the mean link unchanged. Convention is ce >= 0; Euc_signswitch
  # flips all three if ce < 0.
  if (isTRUE(est$ce < 0)){est <- Euc_signswitch(est)}

  # G0 is in standardised y coordinates; t(std_rotation) reverts it to the original space.
  G0 <- estparamlist$G0
  G0 <- t(attr(preplist$y, "std_rotation")) %*% G0
  # Sign convention for G0 axes:
  # (1) make first element of each non-first column positive, removing the per-axis sign ambiguity
  G0[,-1] <- standardise_col_signs(G0[,-1])
  # (2) ensure G0 is a proper rotation matrix (det = +1) by flipping the final column if needed
  if (det(G0) < 0){G0[,p] <- -G0[,p]}
  # DoF: mean link + concentration + scales + G0 axes
  # (p-1)-1: aremaining has p-1 components but prod(aremaining) = 1 removes one degree of freedom
  # G0: DoF_Stiefel(p,p) when G0[,1] is free; DoF_Stiefel(p-1, p-1) when G0[,1] is
  #     identified with p1 or fixed (one fewer free column)
  DoF <- mobius_dof(p, length(est$qs1), length(est$qe1), fix_qs1 = fix_qs1, fix_qe1 = fix_qe1) +
    1 + #concentration
    (p-1)-1 + #aremaining given that prod(aremaining) = 1
    if (G01behaviour == "free"){ #G0 freedom
      DoF_Stiefel(p,p)
    } else {
      DoF_Stiefel(p-1, p-1)
    }
  # AIC
  AIC <- 2*DoF - 2 * lLik
  
  if (p==3) {
    if (estparamlist$k < 1E-15){
    	warning("Estimated concentration is very small and computation of the vMF normalising constant may be breaking down.")
    }
  }
  
  #Scealy and Wood (2019) Proposition 1 check for unimodality
  SvMF_cann_check(SvMF_cann(k = estparamlist$k, a = c(a1, aremaining), G = G0))
  
  niceout <- list(
    mean = est,
    k = estparamlist$k,
    a = c(a1, aremaining),
    G0 = G0,
    obj = nlopt$objective,
    nlopt = nlopt,
    y = y,
    xs = xs,
    # destandardise_Euc inverts the centering and whitening from standardise_data; using
    # preplist$xe also recovers any columns added by addEuccovars (dummy-zero and intercept).
    xe = if (!is.null(xe)){if (intercept){destandardise_Euc(preplist$xe, attr(preplist$xe, "std_center"), attr(preplist$xe, "std_rotation"))} else {xe}},
    pred = destandardise_sph(pred, rotation = attr(preplist$y, "std_rotation")),
    rresids_I = rresids_I,
    rresids_G0 = rresids_G0,
    rresids_std = rresids_std,
    dists = dists,
    DoF = DoF,
    AIC = AIC,
    lLik = lLik,
    initial = initial
  )
  return(niceout)
}

# Preliminary estimation step for mobius_SvMF(): runs a vMF regression to get a good
# starting mean link, then estimates G0 axes (parallel transport base) and shape
# parameters a using method-of-moments on the rotated residuals (Scealy & Wood 2019, Sec 4.1).
mobius_SvMF_partransport_prelim <- function(y, xs, xe, mean = NULL, G0 = NULL, G01behaviour = "p1", type = "LinEuc", fix_qs1 = FALSE, fix_qe1 = (type == "LinEuc"), intercept = TRUE, ...){
  prelim <- mobius_vMF(y = y, xs = xs, xe = xe, 
             start = mean, 
             type = type, fix_qs1 = fix_qs1, fix_qe1 = fix_qe1, intercept = intercept, ...)
  # update starting values accordingly
  mean <- as_mobius_link_Omega(prelim$est)
  k <- prelim$k
  p <- ncol(y)
  
  # get/choose G01 depending on behaviour
  if (G01behaviour == "fixed" && is.null(G0)){stop("At least the first column of G0 must be supplied when G01behaviour = 'fixed'.")}
  G01 <- switch(G01behaviour,
         p1 = prelim$est$p1,
         free = if(is.null(G0)){prelim$est$p1}else{G0[,1]},
         fixed = G0[,1])
  # get rotated residuals
  rresid <- rotated_resid(y, prelim$pred, base = G01)
  if (!is.null(G0) && all(!is.na(G0))){
    # G0 is supplied: re-express its axes relative to the new G01 (prelim$est$p1) by parallel-
    # transporting the columns of G0[,-1] from G0[,1] to G01. The minus sign is because
    # -jupp_Rmat(a,b) equals parallel_transport_mat(a,b). Then estimate scales from the rotated residuals.
    if (G01behaviour == "p1"){
      G0 <- cbind(G01, -jupp_Rmat(G0[,1], G01) %*% G0[,-1])
    }
    aremaining <- SvMF_prelim_scales(rresid, G0)
  } else {
    # G0 not supplied: estimate principal-axis directions from the sample second-moment
    # matrix of the rotated residuals (Scealy & Wood 2019, §4.1).
    G0 <- SvMF_moment_axes(rresid, G01)
    # Estimate scale parameters a using the high-concentration approximation.
    aremaining <- SvMF_prelim_scales(rresid, G0)
  }
  
  prelim <- list(
    mean = mean,
    k = k,
    a = c(1, aremaining),
    G0 = G0,
    nlopt = prelim$nlopt
  )
  return(prelim)
}

# Estimate only the concentration parameter k by univariate MLE, holding a and G0 fixed.
# Used after the main joint optimisation converges when p != 3 (because the normalising
# constant can be computed more accurately for a scalar search than during joint AD optimisation).
mobius_SvMF_konly <- function(y, ymean, a, G0){
  yrot <- undo_partransport(y = y, ymean = ymean, G01 = G0[,1])
  res <- optimise(function(k){
    sum(SvMF_log_lik_cann(yrot, SvMF_cann(k = k, a = a, G = G0)))
  }, lower = 1E-8, upper = 1E5, maximum = TRUE)
  SvMF_cann_check(SvMF_cann(k = res$maximum, a = a, G = G0))
  if (res$maximum == 1E-8){warning("Concentration at numerical lower limit of 1E-8")}
  if (res$maximum == 1E5){warning("Concentration at numerical upper limit of 1E5")}
  return(list(
    k = res$maximum,
    a = a,
    G0 = G0,
    lLik = res$objective
  ))
}

#' @title Function for simulating data given mean link and SvMF parameters
#' @inheritParams mobius_SvMF_log_lik
#' @return A matrix of `p+1` columns and the same number of rows as `xs` or `xe`. The final column is the log-density of the simulated response.
#' @family regression
#' @export
rmobius_SvMF <- function(xs, xe, mean, k, a, G0){
  ymean <- mobius_link(xs = xs, xe = xe, param = mean)
  
  # simulate noise
  y_ld <- t(apply(ymean, 1, function(mn){
    # Construct local axes at the predicted mean mn by parallel-transporting the columns of
    # G0[,-1] from G0[,1] to mn. (-jupp_Rmat(a,b) equals parallel_transport_mat(a,b).)
    G <- cbind(mn, -jupp_Rmat(G0[,1], mn) %*% G0[,-1])
    obs <- rSvMF(1, SvMF_cann(k, a, G))
    ld <- ldSvMF_cann(obs, k = k, a = a, G = G)
    return(c(obs, ld))
  }))
  return(y_ld)
}


# @param y matrix of observations
# @param ymean matrix of predicted means
# @param G01 first column of the G0 matrix
# Parallel-transports each observation y[i,] back from ymean[i,] to G01, bringing all
# observations into a common reference frame (as if they all had mean direction G01).
# This is required to evaluate the pooled SvMF log-likelihood via SvMF_log_lik_cann /
# ldSvMF_cann, which assume the mean direction is G01.
undo_partransport <- function(y, ymean, G01){
  #rotate all observations to reverse the transport from G0[,1] to ymean
  yrot <- lapply(1:nrow(y), function(idx){
    drop(t(parallel_transport_mat(G01, ymean[idx, ])) %*% y[idx, ])
  })
  yrot <- do.call(rbind, yrot)
  return(yrot)
}

#' @title log-density of data for a given SvMF regression
#' @param mean Parameter object (see [`mobius_link_params`]) specifying mean link.
#' @param k Concentration of the SvMF error distribution
#' @param a Scales of the SvMF error distribution
#' @param G0 The base location of parallel transport along with axes \eqn{\gamma_{0j}}, \eqn{j=2,...,p}.
#' @inheritParams mobius_SvMF
#' @description The log-density of each row of `y` for a given SvMF regression. Two methods are used. The approximate method used in the optimisation of all parameters (labelled `Cpp`) and an exact log-density using highly accurate Bessel function implementations from base `R` (labelled `R`).
#' @return A matrix with two columns and the same number of rows as `y`.
#' @family regression
#' @export
mobius_SvMF_log_lik <- function(y, xs, xe, mean, k, a, G0){
  ymean <- mobius_link(xs = xs, xe = xe, param = mean)
  # diff: the discrepancy between the approximate normalising constant used in the C++ AD tape
  # and the exact value; stored as an attribute for diagnostic use.
  diff <- vMF_log_norm_const(k, ncol(y)) - vMF_log_norm_const_exact(k, ncol(y))

  yrot <- undo_partransport(y, ymean, G01 = G0[,1])
  # ldCpp: uses the same approximation as the optimiser (fast, suitable for AD).
  ldCpp <- ldSvMF_cann(yrot, k = k, a = a, G = G0)
  # ldR: uses exact Bessel functions from base R (accurate, not AD-differentiable).
  ldR <- SvMF_log_lik_cann(yrot, SvMF_cann(k = k, a = a, G = G0))
  ld <- cbind(Cpp = ldCpp, R = ldR)
  attr(ld, "error") <- diff
  return(ld)
}
# =============================================================================
# Proposition 4(ii) extension for the spherical regression code
#
# Mean link:
#   mu(xs) = Rtilde0 %*% xs,  Rtilde0 in O(p)
#
# The estimation of G0 is kept in the same form as in the original code:
#   * preliminary moment estimation from parallel-transported residuals;
#   * G01behaviour = "p1", "free", or "fixed";
#   * final joint likelihood estimation using Cayley coordinates relative to
#     G0reference;
#   * permutation/sign standardisation of (a_j, G0[,j]), j >= 2.
#
# Source this file AFTER the original package R files have been loaded.
# =============================================================================

# Access non-exported sphm functions used by the Proposition 4(ii) extension.
# They are available in namespace:sphm but are not necessarily attached by
# library(sphm). Import them when needed so the code also works in a clean
# R session.
.prop4ii_import_sphm_internal <- function(fname) {
  if (!exists(fname, mode = "function", inherits = TRUE)) {
    obj <- try(getFromNamespace(fname, "sphm"), silent = TRUE)
    if (!inherits(obj, "try-error")) {
      assign(fname, obj, envir = .GlobalEnv)
    }
  }
  invisible(NULL)
}

invisible(lapply(
  c(
    "ldSvMF_cann",
    "SvMF_log_lik_cann",
    "SvMF_cann_check",
    "resid_SvMF_partransport",
    "rotated_resid",
    "north_pole",
    "parallel_transport_mat",
    "standardise_col_signs",
    "jupp_Rmat",
    "SvMF_prelim_scales",
    "SvMF_moment_axes",
    "Sp",
    "iSp"
  ),
  .prop4ii_import_sphm_internal
))

# Save the original functions only once. This makes repeated source() calls safe.
if (!exists(".mobius_link_definition1", inherits = FALSE)) {
  .mobius_link_definition1 <- mobius_link
}
if (!exists(".mobius_vMF_definition1", inherits = FALSE)) {
  .mobius_vMF_definition1 <- mobius_vMF
}
if (!exists(".mobius_SvMF_definition1", inherits = FALSE)) {
  .mobius_SvMF_definition1 <- mobius_SvMF
}

.is_prop4ii_type <- function(type) {
  if (length(type) != 1L || is.na(type)) return(FALSE)
  key <- tolower(gsub("[^[:alnum:]]", "", as.character(type)))
  key %in% c("prop4ii", "proposition4ii")
}

.prop4ii_check_unit_rows <- function(x, name,
                                     tol = 1000 * sqrt(.Machine$double.eps)) {
  if (!is.matrix(x) || !is.numeric(x)) {
    stop(name, " must be a numeric matrix.")
  }
  if (any(!is.finite(x))) stop(name, " contains non-finite values.")
  nr <- sqrt(rowSums(x^2))
  err <- max(abs(nr - 1))
  if (err > tol) {
    stop(name, " must contain unit vectors as rows; maximum norm error is ",
         format(err, digits = 4), ".")
  }
  invisible(NULL)
}

.prop4ii_nearest_orthogonal <- function(M, det_sign = NULL) {
  M <- as.matrix(M)
  if (nrow(M) != ncol(M)) stop("M must be square.")
  p <- nrow(M)
  z <- svd(M, nu = p, nv = p)
  D <- diag(p)
  R0 <- z$u %*% t(z$v)
  if (!is.null(det_sign)) {
    requested <- if (det_sign >= 0) 1 else -1
    current <- if (det(R0) >= 0) 1 else -1
    D[p, p] <- requested / current
  }
  z$u %*% D %*% t(z$v)
}

# -----------------------------------------------------------------------------
# Proposition 4(ii) parameter class and mean calculation
# -----------------------------------------------------------------------------

mobius_link_prop4ii <- function(Rtilde0, check = TRUE) {
  Rtilde0 <- as.matrix(Rtilde0)
  storage.mode(Rtilde0) <- "double"
  obj <- list(Rtilde0 = Rtilde0)
  class(obj) <- c("mobius_link_prop4ii", "list")
  if (check) mobius_link_prop4ii_check(obj)
  obj
}

mobius_link_prop4ii_check <- function(obj,
                                      tol = 100 * sqrt(.Machine$double.eps)) {
  if (!inherits(obj, "mobius_link_prop4ii")) {
    stop("obj must inherit from 'mobius_link_prop4ii'.")
  }
  R <- obj$Rtilde0
  if (!is.matrix(R) || nrow(R) != ncol(R) || nrow(R) < 2L) {
    stop("Rtilde0 must be a square p by p matrix with p >= 2.")
  }
  if (any(!is.finite(R))) stop("Rtilde0 contains non-finite values.")
  err <- max(abs(crossprod(R) - diag(ncol(R))))
  if (err > tol) {
    stop(sprintf("Rtilde0 is not orthogonal: max|t(R)R-I| = %.3e.", err))
  }
  invisible(NULL)
}

as_mobius_link_prop4ii <- function(obj, check = TRUE) {
  if (inherits(obj, "mobius_link_prop4ii")) {
    if (check) mobius_link_prop4ii_check(obj)
    return(obj)
  }
  if (is.matrix(obj)) return(mobius_link_prop4ii(obj, check = check))
  if (is.list(obj) && !is.null(obj$Rtilde0)) {
    return(mobius_link_prop4ii(obj$Rtilde0, check = check))
  }
  if (inherits(obj, "mobius_link_cann")) {
    p <- nrow(obj$P)
    if (!is.null(obj$Qe) || is.null(obj$Qs) ||
        !all(dim(obj$Qs) == c(p, p))) {
      stop("The canonical object is not a Proposition 4(ii) link.")
    }
    if (max(abs(obj$Bs - diag(p - 1L))) >
        100 * sqrt(.Machine$double.eps)) {
      stop("Proposition 4(ii) requires Bs = I_(p-1).")
    }
    return(mobius_link_prop4ii(obj$P %*% t(obj$Qs), check = check))
  }
  stop("obj must be a mobius_link_prop4ii object, a square matrix, or a compatible canonical object.")
}

print.mobius_link_prop4ii <- function(x, ...) {
  cat("Proposition 4(ii) spherical link\n")
  cat("  p = qs =", nrow(x$Rtilde0), ", qe = 0\n")
  cat("  det(Rtilde0) =", format(det(x$Rtilde0), digits = 6), "\n")
  print(x$Rtilde0, ...)
  invisible(x)
}

dim.mobius_link_prop4ii <- function(x) {
  p <- nrow(x$Rtilde0)
  c(p = p, qs = p, qe = 0L)
}

# Wrapper: retain the original Definition 1 link for its original classes.
mobius_link <- function(xs = NULL, xe = NULL, param = NULL, check = TRUE) {
  if (inherits(param, "mobius_link_prop4ii")) {
    if (check) mobius_link_prop4ii_check(param)
    if (is.null(xs)) stop("The Proposition 4(ii) link requires xs.")
    xs <- as.matrix(xs)
    if (!is.null(xe) && ncol(as.matrix(xe)) > 0L) {
      stop("Proposition 4(ii) requires qe = 0; xe must be NULL.")
    }
    if (ncol(xs) != nrow(param$Rtilde0)) {
      stop("ncol(xs) must equal nrow(Rtilde0).")
    }
    # Rows contain transposed column vectors: mu_i^T = xs_i^T Rtilde0^T.
    return(xs %*% t(param$Rtilde0))
  }
  .mobius_link_definition1(xs = xs, xe = xe, param = param, check = check)
}

# -----------------------------------------------------------------------------
# Direct vMF estimation of Rtilde0 by orthogonal Procrustes
# -----------------------------------------------------------------------------

.prop4ii_procrustes <- function(y, xs, det_sign = NULL) {
  p <- ncol(y)
  C <- crossprod(y, xs)             # Y^T X
  z <- svd(C, nu = p, nv = p)
  D <- diag(p)
  R0 <- z$u %*% t(z$v)
  if (!is.null(det_sign)) {
    requested <- if (det_sign >= 0) 1 else -1
    current <- if (det(R0) >= 0) 1 else -1
    D[p, p] <- requested / current
  }
  list(R = z$u %*% D %*% t(z$v), singular_values = z$d)
}

.prop4ii_estimate_k <- function(rbar, p) {
  rbar <- min(max(as.numeric(rbar), -1 + 1e-12), 1 - 1e-12)
  objective <- function(k) {
    -vMF_log_norm_const_exact(k, p) + k * rbar
  }
  fit <- optimise(objective, lower = 1e-8, upper = 1e5, maximum = TRUE)
  list(k = fit$maximum, loglik_per_obs = fit$objective)
}

mobius_vMF_prop4ii <- function(y, xs, xe = NULL,
                               det_constraint = c("orthogonal", "rotation"),
                               start = NULL, check = TRUE) {
  det_constraint <- match.arg(det_constraint)
  y <- as.matrix(y)
  xs <- as.matrix(xs)
  storage.mode(y) <- storage.mode(xs) <- "double"

  if (!is.null(xe) && ncol(as.matrix(xe)) > 0L) {
    stop("Proposition 4(ii) requires qe = 0; xe must be NULL.")
  }
  if (nrow(y) != nrow(xs) || ncol(y) != ncol(xs)) {
    stop("For Proposition 4(ii), y and xs must have the same n by p dimensions.")
  }
  if (ncol(y) < 2L) stop("p must be at least 2.")
  if (check) {
    .prop4ii_check_unit_rows(y, "y")
    .prop4ii_check_unit_rows(xs, "xs")
  }
  if (!is.null(start)) {
    warning("start is ignored because the vMF estimate of Rtilde0 is available in closed form.")
  }

  det_sign <- if (det_constraint == "rotation") 1 else NULL
  pro <- .prop4ii_procrustes(y, xs, det_sign = det_sign)
  Rhat <- pro$R
  est <- mobius_link_prop4ii(Rhat)
  pred <- xs %*% t(Rhat)
  dots <- pmax(-1, pmin(1, rowSums(y * pred)))
  rbar <- mean(dots)
  kfit <- .prop4ii_estimate_k(rbar, ncol(y))

  p <- ncol(y)
  n <- nrow(y)
  DoF <- p * (p - 1L) / 2L + 1L
  lLik <- n * kfit$loglik_per_obs
  dists <- acos(dots)
  rresids <- NULL
  if (exists("rotated_resid", mode = "function")) {
    rr <- try(rotated_resid(y, pred, north_pole(p)), silent = TRUE)
    if (!inherits(rr, "try-error")) {
      rresids <- rr[, -1, drop = FALSE]
      attr(rresids, "samehemisphere") <- attr(rr, "samehemisphere")
      colnames(rresids) <- paste0("r", seq_len(ncol(rresids)))
    }
  }

  out <- list(
    est = est,
    mean = est,
    Rtilde0 = Rhat,
    k = kfit$k,
    obj = -rbar,
    solution = est,
    nlopt = list(
      status = 1L,
      message = "Closed-form orthogonal Procrustes solution.",
      objective = -rbar,
      solution = as.vector(Rhat)
    ),
    y = y,
    xs = xs,
    xe = NULL,
    pred = pred,
    rresids = rresids,
    dists = dists,
    DoF = DoF,
    AIC = 2 * DoF - 2 * lLik,
    lLik = lLik,
    start = start,
    singular_values = pro$singular_values,
    linktype = list(type = "Prop4ii", det_constraint = det_constraint,
                    intercept = FALSE)
  )
  class(out) <- c("mobius_vMF_prop4ii", "mobius_vMF", "list")
  out
}

print.mobius_vMF_prop4ii <- function(x, ...) {
  cat("vMF regression with Proposition 4(ii) link\n")
  cat("  det(Rtilde0) =", format(det(x$Rtilde0), digits = 6), "\n")
  cat("  k =", format(x$k, digits = 6), "\n")
  cat("  log-likelihood =", format(x$lLik, digits = 7), "\n")
  cat("  AIC =", format(x$AIC, digits = 7), "\n")
  invisible(x)
}

predict.mobius_vMF_prop4ii <- function(object, newdata = NULL,
                                        xs = newdata, ...) {
  if (is.null(xs)) xs <- object$xs
  mobius_link(xs = as.matrix(xs), param = object$mean)
}

# Wrapper retaining the original vMF implementation for the original link types.
mobius_vMF <- function(y, xs = NULL, xe = NULL, start = NULL,
                       type = "SpEuc", fix_qs1 = FALSE,
                       fix_qe1 = (type == "LinEuc"), intercept = TRUE,
                       lb = NULL, ub = NULL,
                       det_constraint = c("orthogonal", "rotation"), ...) {
  if (.is_prop4ii_type(type)) {
    return(mobius_vMF_prop4ii(
      y = y, xs = xs, xe = xe, start = start,
      det_constraint = det_constraint
    ))
  }
  .mobius_vMF_definition1(
    y = y, xs = xs, xe = xe, start = start, type = type,
    fix_qs1 = fix_qs1, fix_qe1 = fix_qe1, intercept = intercept,
    lb = lb, ub = ub, ...
  )
}

# -----------------------------------------------------------------------------
# Cayley-coordinate helpers
# -----------------------------------------------------------------------------

.prop4ii_skew_dim <- function(p) as.integer(p * (p - 1L) / 2L)

.prop4ii_vec_to_skew <- function(theta, p) {
  d <- .prop4ii_skew_dim(p)
  if (length(theta) != d) stop("Incorrect skew-vector length.")
  A <- matrix(0, p, p)
  if (d > 0L) {
    ij <- which(upper.tri(A), arr.ind = TRUE)
    A[ij] <- theta
    A[cbind(ij[, 2], ij[, 1])] <- -theta
  }
  A
}

.prop4ii_cayley <- function(theta, p) {
  A <- .prop4ii_vec_to_skew(theta, p)
  I <- diag(p)
  solve(I - A, I + A)
}

.prop4ii_inverse_cayley <- function(Q, tol = 1e-10) {
  Q <- as.matrix(Q)
  p <- nrow(Q)
  I <- diag(p)
  if (ncol(Q) != p || rcond(Q + I) < tol) {
    return(rep(0, .prop4ii_skew_dim(p)))
  }
  A <- (Q - I) %*% solve(Q + I)
  A <- (A - t(A)) / 2
  A[upper.tri(A)]
}

.prop4ii_unit <- function(x, name = "vector") {
  x <- as.numeric(x)
  nrm <- sqrt(sum(x^2))
  if (!is.finite(nrm) || nrm <= sqrt(.Machine$double.eps)) {
    stop(name, " must be a finite nonzero vector.")
  }
  x / nrm
}

.prop4ii_frame_with_base <- function(base, tangent, tol = 1e-10) {
  base <- .prop4ii_unit(base, "base")
  p <- length(base)
  tangent <- as.matrix(tangent)
  if (!all(dim(tangent) == c(p, p - 1L))) {
    stop("tangent must be p by (p-1).")
  }
  base_col <- matrix(base, ncol = 1L)
  tangent <- tangent - base_col %*% crossprod(base, tangent)
  q <- qr(tangent, tol = tol)
  if (q$rank < p - 1L) {
    seed <- diag(p)[, -which.max(abs(base)), drop = FALSE]
    seed <- seed - base_col %*% crossprod(base, seed)
    tangent <- qr.Q(qr(seed), complete = FALSE)[, seq_len(p - 1L), drop = FALSE]
  } else {
    tangent <- qr.Q(q, complete = FALSE)[, seq_len(p - 1L), drop = FALSE]
  }
  G <- cbind(base, tangent)
  if (det(G) < 0) G[, p] <- -G[, p]
  G
}

# Same parallel-transport re-expression used in the original preliminary code.
.prop4ii_rebase_G0 <- function(G0, newbase, tol = 1e-8) {
  G0 <- .prop4ii_nearest_orthogonal(G0, det_sign = 1)
  oldbase <- G0[, 1]
  newbase <- .prop4ii_unit(newbase, "newbase")
  if (sum(oldbase * newbase) <= -1 + tol) return(NULL)

  tangent <- NULL
  if (exists("jupp_Rmat", mode = "function")) {
    tangent <- try(-jupp_Rmat(oldbase, newbase) %*%
                     G0[, -1, drop = FALSE], silent = TRUE)
  }
  if (is.null(tangent) || inherits(tangent, "try-error")) {
    tangent <- try(parallel_transport_mat(oldbase, newbase) %*%
                     G0[, -1, drop = FALSE], silent = TRUE)
  }
  if (inherits(tangent, "try-error") || any(!is.finite(tangent))) return(NULL)
  .prop4ii_frame_with_base(newbase, tangent)
}

.prop4ii_normalise_scales <- function(a, p, a1 = NULL) {
  a <- as.numeric(a)
  if (length(a) != p || any(!is.finite(a)) || any(a <= 0)) {
    stop("a must be a positive vector of length p.")
  }
  if (!is.null(a1)) a[1] <- a1
  if (p > 1L) a[-1] <- a[-1] / exp(mean(log(a[-1])))
  a
}

.prop4ii_scales_to_eta <- function(a) {
  p <- length(a)
  if (p <= 2L) return(numeric(0))
  log(a[2:(p - 1L)])
}

.prop4ii_eta_to_scales <- function(eta, a1, p) {
  if (p == 2L) return(c(a1, 1))
  if (length(eta) != p - 2L) stop("Incorrect number of scale coordinates.")
  c(a1, exp(eta), exp(-sum(eta)))
}

# -----------------------------------------------------------------------------
# Preliminary SvMF step: G0 estimation is the original procedure
# -----------------------------------------------------------------------------

mobius_SvMF_partransport_prelim_prop4ii <- function(
    y, xs, xe = NULL, mean = NULL, G0 = NULL,
    G01behaviour = "p1",
    det_constraint = c("orthogonal", "rotation"), check = TRUE, ...) {

  det_constraint <- match.arg(det_constraint)
  G01behaviour <- match.arg(G01behaviour, c("p1", "free", "fixed"))
  y <- as.matrix(y)
  xs <- as.matrix(xs)

  if (!is.null(xe) && ncol(as.matrix(xe)) > 0L) {
    stop("Proposition 4(ii) requires qe = 0; xe must be NULL.")
  }
  if (nrow(y) != nrow(xs) || ncol(y) != ncol(xs)) {
    stop("For Proposition 4(ii), y and xs must have identical dimensions.")
  }
  if (check) {
    .prop4ii_check_unit_rows(y, "y")
    .prop4ii_check_unit_rows(xs, "xs")
  }

  prelim <- mobius_vMF_prop4ii(
    y = y, xs = xs, xe = NULL, start = mean,
    det_constraint = det_constraint, check = FALSE
  )
  mean <- prelim$est
  k <- prelim$k
  Rtilde0 <- prelim$Rtilde0

  # This is exactly the original switch, with prelim$est$p1 replaced by the
  # identifiable representative p1 = Rtilde0[,1].
  if (G01behaviour == "fixed" && is.null(G0)) {
    stop("At least the first column of G0 must be supplied when G01behaviour = 'fixed'.")
  }
  G01 <- switch(
    G01behaviour,
    p1 = Rtilde0[, 1],
    free = if (is.null(G0)) Rtilde0[, 1] else G0[, 1],
    fixed = G0[, 1]
  )

  rresid <- rotated_resid(y, prelim$pred, base = G01)

  if (!is.null(G0) && all(!is.na(G0))) {
    if (G01behaviour == "p1") {
      # Original re-expression of supplied axes at the new p1 base.
      G0 <- cbind(G01, -jupp_Rmat(G0[, 1], G01) %*%
                         G0[, -1, drop = FALSE])
    }
    aremaining <- SvMF_prelim_scales(rresid, G0)
  } else {
    # Original moment estimator of the axes and scales.
    G0 <- SvMF_moment_axes(rresid, G01)
    aremaining <- SvMF_prelim_scales(rresid, G0)
  }

  list(
    mean = mean,
    Rtilde0 = Rtilde0,
    k = k,
    a = c(1, aremaining),
    G0 = G0,
    nlopt = prelim$nlopt,
    pred = prelim$pred
  )
}

# -----------------------------------------------------------------------------
# Joint SvMF likelihood. G0 uses the original Cayley-coordinate structure.
# -----------------------------------------------------------------------------

.prop4ii_joint_loglik <- function(y, xs, Rtilde0, k, a, G0,
                                  approximate = TRUE) {
  if (!is.finite(k) || k <= 0 || any(!is.finite(a)) || any(a <= 0)) {
    return(-Inf)
  }
  pred <- xs %*% t(Rtilde0)
  yrot <- try(undo_partransport(y, pred, G01 = G0[, 1]), silent = TRUE)
  if (inherits(yrot, "try-error") || any(!is.finite(yrot))) return(-Inf)

  # Try the same fast approximate density used by the original optimiser first.
  # For some otherwise valid parameter values ldSvMF_cann() can return a
  # non-finite value.  In that case fall back to the exact R implementation
  # rather than treating the whole determinant component as having -Inf
  # likelihood.  This does not change the model; it is only a numerical
  # safeguard for Proposition 4(ii).
  ld <- NULL
  if (approximate && exists("ldSvMF_cann", mode = "function")) {
    ld_try <- try(ldSvMF_cann(yrot, k = k, a = a, G = G0), silent = TRUE)
    if (!inherits(ld_try, "try-error") && all(is.finite(ld_try))) {
      ld <- ld_try
    }
  }
  if (is.null(ld)) {
    ld_try <- try(
      SvMF_log_lik_cann(yrot, SvMF_cann(k = k, a = a, G = G0)),
      silent = TRUE
    )
    if (!inherits(ld_try, "try-error") && all(is.finite(ld_try))) {
      ld <- ld_try
    }
  }
  if (is.null(ld)) return(-Inf)
  sum(ld)
}

.prop4ii_build_joint_spec <- function(Rref, prelim, G0reference,
                                      G01behaviour) {
  p <- nrow(Rref)
  G01behaviour <- match.arg(G01behaviour, c("p1", "free", "fixed"))

  dR <- .prop4ii_skew_dim(p)
  dA <- max(p - 2L, 0L)
  dG <- if (G01behaviour == "free") {
    .prop4ii_skew_dim(p)
  } else {
    .prop4ii_skew_dim(p - 1L)
  }

  pos <- 0L
  take <- function(n) {
    if (n == 0L) return(integer(0))
    out <- pos + seq_len(n)
    pos <<- pos + n
    out
  }
  idxR <- take(dR)
  idxK <- take(1L)
  idxA <- take(dA)
  idxG <- take(dG)

  Gstart <- .prop4ii_nearest_orthogonal(prelim$G0, det_sign = 1)
  Gcoord <- if (is.null(G0reference)) {
    Gstart
  } else {
    .prop4ii_nearest_orthogonal(G0reference, det_sign = 1)
  }

  if (G01behaviour == "free") {
    Gref <- Gcoord
    Q0 <- Gstart %*% t(Gref)
    thetaG0 <- .prop4ii_inverse_cayley(Q0)
    Ganchor <- NULL
    fixed_base <- NULL
  } else {
    target_base <- if (G01behaviour == "p1") Rref[, 1] else Gstart[, 1]

    Ganchor <- .prop4ii_rebase_G0(Gcoord, target_base)
    if (is.null(Ganchor)) Ganchor <- .prop4ii_rebase_G0(Gstart, target_base)
    if (is.null(Ganchor)) stop("Could not construct the G0 reference frame.")

    Gstart2 <- .prop4ii_rebase_G0(Gstart, target_base)
    if (is.null(Gstart2)) stop("Could not rebase the starting G0 frame.")

    H0 <- crossprod(Ganchor[, -1, drop = FALSE],
                    Gstart2[, -1, drop = FALSE])
    H0 <- .prop4ii_nearest_orthogonal(H0, det_sign = 1)
    thetaG0 <- .prop4ii_inverse_cayley(H0)
    Gref <- NULL
    fixed_base <- if (G01behaviour == "fixed") target_base else NULL
  }

  par0 <- numeric(pos)
  par0[idxR] <- 0
  par0[idxK] <- log(max(prelim$k, 1e-8))
  if (dA > 0L) par0[idxA] <- .prop4ii_scales_to_eta(prelim$a)
  if (dG > 0L) par0[idxG] <- thetaG0

  # As in the original Cayley parameterisation, rotation coordinates are
  # unconstrained. Bounds are used only for k and the log-scale coordinates.
  lower <- rep(-Inf, pos)
  upper <- rep( Inf, pos)
  lower[idxK] <- log(1e-8)
  upper[idxK] <- log(1e5)
  if (dA > 0L) {
    lower[idxA] <- -10
    upper[idxA] <- 10
  }

  list(
    p = p,
    Rref = Rref,
    a1 = prelim$a[1],
    G01behaviour = G01behaviour,
    Gref = Gref,
    Ganchor = Ganchor,
    fixed_base = fixed_base,
    idxR = idxR,
    idxK = idxK,
    idxA = idxA,
    idxG = idxG,
    par0 = par0,
    lower = lower,
    upper = upper
  )
}

.prop4ii_unpack_joint <- function(theta, spec) {
  p <- spec$p
  Rtilde0 <- .prop4ii_cayley(theta[spec$idxR], p) %*% spec$Rref
  k <- exp(theta[spec$idxK])
  eta <- if (length(spec$idxA)) theta[spec$idxA] else numeric(0)
  a <- .prop4ii_eta_to_scales(eta, spec$a1, p)

  if (spec$G01behaviour == "free") {
    G0 <- .prop4ii_cayley(theta[spec$idxG], p) %*% spec$Gref
  } else {
    Gbase <- spec$Ganchor
    if (spec$G01behaviour == "p1") {
      # Same identification as the original code: G0[,1] is p1. Here the
      # identifiable representative is p1 = Rtilde0[,1].
      Gbase <- .prop4ii_rebase_G0(spec$Ganchor, Rtilde0[, 1])
      if (is.null(Gbase)) return(NULL)
    }
    H <- if (p == 2L) {
      matrix(1, 1, 1)
    } else {
      .prop4ii_cayley(theta[spec$idxG], p - 1L)
    }
    B <- diag(p)
    B[-1, -1] <- H
    G0 <- Gbase %*% B
  }

  list(Rtilde0 = Rtilde0, k = k, a = a, G0 = G0)
}


.prop4ii_fd_gradient <- function(fn, x, lower, upper,
                                 rel_step = .Machine$double.eps^(1/3)) {
  g <- numeric(length(x))
  f0 <- fn(x)
  for (j in seq_along(x)) {
    h <- rel_step * (1 + abs(x[j]))
    xp <- xm <- x
    xp[j] <- min(x[j] + h, upper[j])
    xm[j] <- max(x[j] - h, lower[j])
    if (xp[j] > x[j] && xm[j] < x[j]) {
      g[j] <- (fn(xp) - fn(xm)) / (xp[j] - xm[j])
    } else if (xp[j] > x[j]) {
      g[j] <- (fn(xp) - f0) / (xp[j] - x[j])
    } else if (xm[j] < x[j]) {
      g[j] <- (f0 - fn(xm)) / (x[j] - xm[j])
    } else {
      g[j] <- 0
    }
  }
  g
}

.prop4ii_fit_component <- function(y, xs, Rref, prelim, G0reference,
                                   G01behaviour, algorithm, opts,
                                   approximate = TRUE) {
  spec <- .prop4ii_build_joint_spec(
    Rref = Rref,
    prelim = prelim,
    G0reference = G0reference,
    G01behaviour = G01behaviour
  )

  eval_f <- function(theta) {
    pars <- .prop4ii_unpack_joint(theta, spec)
    if (is.null(pars)) return(1e100)
    ll <- .prop4ii_joint_loglik(
      y, xs, pars$Rtilde0, pars$k, pars$a, pars$G0,
      approximate = approximate
    )
    if (!is.finite(ll)) return(1e100)
    -ll / nrow(y)
  }

  default_opts <- list(
    algorithm = algorithm,
    xtol_rel = 1e-9,
    maxeval = 10000,
    print_level = 0
  )
  combined_opts <- utils::modifyList(default_opts, opts)

  uses_gradient <- grepl("^NLOPT_LD_", combined_opts$algorithm)
  eval_nloptr <- if (uses_gradient) {
    function(theta) {
      value <- eval_f(theta)
      gradient <- .prop4ii_fd_gradient(
        eval_f, theta, spec$lower, spec$upper
      )
      list(objective = value, gradient = gradient)
    }
  } else {
    eval_f
  }

  fit <- nloptr::nloptr(
    x0 = spec$par0,
    eval_f = eval_nloptr,
    lb = spec$lower,
    ub = spec$upper,
    opts = combined_opts
  )

  pars <- .prop4ii_unpack_joint(fit$solution, spec)
  ll <- if (is.null(pars)) -Inf else
    .prop4ii_joint_loglik(
      y, xs, pars$Rtilde0, pars$k, pars$a, pars$G0,
      approximate = approximate
    )

  list(nlopt = fit, pars = pars, lLik_approx = ll, spec = spec)
}

mobius_SvMF_joint_fit_prop4ii <- function(
    y, xs, xe = NULL, mean, k, a, G0,
    G0reference = NULL, G01behaviour = "p1",
    det_constraint = c("orthogonal", "rotation"),
    algorithm = "NLOPT_LD_SLSQP", opts = list(), check = TRUE) {

  det_constraint <- match.arg(det_constraint)
  G01behaviour <- match.arg(G01behaviour, c("p1", "free", "fixed"))
  y <- as.matrix(y)
  xs <- as.matrix(xs)
  storage.mode(y) <- storage.mode(xs) <- "double"

  if (!is.null(xe) && ncol(as.matrix(xe)) > 0L) {
    stop("Proposition 4(ii) requires qe = 0; xe must be NULL.")
  }
  if (nrow(y) != nrow(xs) || ncol(y) != ncol(xs)) {
    stop("For Proposition 4(ii), y and xs must have identical dimensions.")
  }
  p <- ncol(y)
  if (check) {
    .prop4ii_check_unit_rows(y, "y")
    .prop4ii_check_unit_rows(xs, "xs")
  }

  start_mean <- as_mobius_link_prop4ii(mean)
  start_R <- start_mean$Rtilde0
  start_a <- .prop4ii_normalise_scales(a, p, a1 = a[1])
  start_G0 <- .prop4ii_nearest_orthogonal(G0, det_sign = 1)

  initial <- list(mean = start_mean, Rtilde0 = start_R,
                  k = k, a = start_a, G0 = start_G0)

  target_signs <- if (det_constraint == "rotation") 1 else c(1, -1)
  fits <- vector("list", length(target_signs))

  for (j in seq_along(target_signs)) {
    sgn <- target_signs[j]
    if ((if (det(start_R) >= 0) 1 else -1) == sgn) {
      Rref <- start_R
    } else {
      Rref <- .prop4ii_procrustes(y, xs, det_sign = sgn)$R
    }

    pre_j <- initial
    pre_j$Rtilde0 <- Rref
    pre_j$mean <- mobius_link_prop4ii(Rref)

    # Under G01behaviour="p1", the original G0 base must follow p1.
    if (G01behaviour == "p1") {
      rebased <- .prop4ii_rebase_G0(pre_j$G0, Rref[, 1])
      if (is.null(rebased)) stop("Could not rebase G0 to Rtilde0[,1].")
      pre_j$G0 <- rebased
    }

    fits[[j]] <- .prop4ii_fit_component(
      y = y,
      xs = xs,
      Rref = Rref,
      prelim = pre_j,
      G0reference = if (is.null(G0reference)) pre_j$G0 else G0reference,
      G01behaviour = G01behaviour,
      algorithm = algorithm,
      opts = opts,
      approximate = TRUE
    )
    fits[[j]]$det_sign <- sgn
  }

  ll_approx <- vapply(fits, function(z) z$lLik_approx, numeric(1))
  if (all(!is.finite(ll_approx))) stop("All Proposition 4(ii) optimisations failed.")
  best <- fits[[which.max(ll_approx)]]
  pars <- best$pars

  # Original identifiability convention for (a_j, G0[,j]), j >= 2.
  aord <- order(pars$a[-1], decreasing = TRUE)
  pars$a[-1] <- pars$a[-1][aord]
  pars$G0[, -1] <- pars$G0[, -1, drop = FALSE][, aord, drop = FALSE]
  if (exists("standardise_col_signs", mode = "function")) {
    pars$G0[, -1] <- standardise_col_signs(pars$G0[, -1, drop = FALSE])
  }
  if (det(pars$G0) < 0) pars$G0[, p] <- -pars$G0[, p]

  est <- mobius_link_prop4ii(pars$Rtilde0)
  pred <- xs %*% t(pars$Rtilde0)

  # The original code refines k with the exact likelihood when p != 3.
  if (p != 3L) {
    kres <- mobius_SvMF_konly(y = y, ymean = pred,
                              a = pars$a, G0 = pars$G0)
    pars$k <- kres$k
    lLik <- kres$lLik
  } else {
    yrot <- undo_partransport(y, pred, G01 = pars$G0[, 1])
    lLik <- sum(SvMF_log_lik_cann(
      yrot, SvMF_cann(k = pars$k, a = pars$a, G = pars$G0)
    ))
  }

  # Same final validity/unimodality check as in the original code.
  SvMF_cann_check(SvMF_cann(k = pars$k, a = pars$a, G = pars$G0))

  dots <- pmax(-1, pmin(1, rowSums(y * pred)))
  dists <- acos(dots)

  rrI <- rotated_resid(y, pred, north_pole(p))
  rresids_I <- rrI[, -1, drop = FALSE]
  attr(rresids_I, "samehemisphere") <- attr(rrI, "samehemisphere")
  colnames(rresids_I) <- paste0("r", seq_len(ncol(rresids_I)))

  rresids_G0 <- resid_SvMF_partransport(
    y, pred, G0 = pars$G0, scale = FALSE
  )
  rresids_std <- resid_SvMF_partransport(
    y, pred, pars$k, pars$a, pars$G0, scale = TRUE
  )

  dmean <- .prop4ii_skew_dim(p)
  dscale <- max(p - 2L, 0L)
  dG0 <- if (G01behaviour == "free") {
    .prop4ii_skew_dim(p)
  } else {
    .prop4ii_skew_dim(p - 1L)
  }
  DoF <- dmean + 1L + dscale + dG0

  out <- list(
    mean = est,
    est = est,
    Rtilde0 = pars$Rtilde0,
    k = pars$k,
    a = pars$a,
    G0 = pars$G0,
    obj = -lLik / nrow(y),
    nlopt = best$nlopt,
    y = y,
    xs = xs,
    xe = NULL,
    pred = pred,
    rresids_I = rresids_I,
    rresids_G0 = rresids_G0,
    rresids_std = rresids_std,
    dists = dists,
    DoF = DoF,
    AIC = 2 * DoF - 2 * lLik,
    lLik = lLik,
    initial = initial,
    component_fits = lapply(fits, function(z) {
      list(det_sign = z$det_sign,
           status = z$nlopt$status,
           message = z$nlopt$message,
           lLik_approx = z$lLik_approx)
    }),
    linktype = list(
      type = "Prop4ii",
      det_constraint = det_constraint,
      G01behaviour = G01behaviour,
      intercept = FALSE
    )
  )
  class(out) <- c("mobius_SvMF_prop4ii", "mobius_SvMF", "list")
  out
}

mobius_SvMF_prop4ii <- function(
    y, xs, xe = NULL, mean = NULL, k = NULL, a = NULL,
    G0 = NULL, G0reference = NULL, G01behaviour = "p1",
    doprelim = TRUE,
    det_constraint = c("orthogonal", "rotation"),
    algorithm = "NLOPT_LD_SLSQP", opts = list(), check = TRUE, ...) {

  det_constraint <- match.arg(det_constraint)
  G01behaviour <- match.arg(G01behaviour, c("p1", "free", "fixed"))
  dots <- list(...)
  if (length(dots)) opts <- utils::modifyList(opts, dots)

  if (doprelim) {
    preest <- mobius_SvMF_partransport_prelim_prop4ii(
      y = y, xs = xs, xe = xe, mean = mean, G0 = G0,
      G01behaviour = G01behaviour,
      det_constraint = det_constraint, check = check
    )
  } else {
    if (is.null(mean) || is.null(k) || is.null(a) || is.null(G0)) {
      stop("When doprelim = FALSE, mean, k, a, and G0 must all be supplied.")
    }
    mean <- as_mobius_link_prop4ii(mean)
    preest <- list(
      mean = mean,
      Rtilde0 = mean$Rtilde0,
      k = k,
      a = a,
      G0 = G0
    )
  }

  finalest <- mobius_SvMF_joint_fit_prop4ii(
    y = y,
    xs = xs,
    xe = xe,
    mean = preest$mean,
    k = if (is.null(k)) preest$k else k,
    a = if (is.null(a)) preest$a else a,
    G0 = preest$G0,
    G0reference = if (is.null(G0reference)) preest$G0 else G0reference,
    G01behaviour = G01behaviour,
    det_constraint = det_constraint,
    algorithm = algorithm,
    opts = opts,
    check = check
  )
  finalest$preest <- preest
  finalest
}

print.mobius_SvMF_prop4ii <- function(x, ...) {
  cat("SvMF regression with Proposition 4(ii) link\n")
  cat("  G0 estimation: original parallel-transport/Cayley procedure\n")
  cat("  G01behaviour =", x$linktype$G01behaviour, "\n")
  cat("  det(Rtilde0) =", format(det(x$Rtilde0), digits = 6), "\n")
  cat("  k =", format(x$k, digits = 6), "\n")
  cat("  a =", paste(format(x$a, digits = 5), collapse = ", "), "\n")
  cat("  log-likelihood =", format(x$lLik, digits = 7), "\n")
  cat("  AIC =", format(x$AIC, digits = 7), "\n")
  invisible(x)
}

predict.mobius_SvMF_prop4ii <- function(object, newdata = NULL,
                                         xs = newdata, ...) {
  if (is.null(xs)) xs <- object$xs
  mobius_link(xs = as.matrix(xs), param = object$mean)
}

# Wrapper retaining the original SvMF implementation for the original link types.
mobius_SvMF <- function(
    y, xs = NULL, xe = NULL, mean = NULL, k = NULL, a = NULL,
    G0 = NULL, G0reference = NULL, G01behaviour = "p1",
    type = "LinEuc", fix_qs1 = FALSE,
    fix_qe1 = (type == "LinEuc"), intercept = TRUE,
    doprelim = TRUE,
    det_constraint = c("orthogonal", "rotation"),
    prop4ii_algorithm = "NLOPT_LD_SLSQP",
    prop4ii_opts = list(), ...) {

  if (.is_prop4ii_type(type)) {
    return(mobius_SvMF_prop4ii(
      y = y, xs = xs, xe = xe, mean = mean, k = k, a = a,
      G0 = G0, G0reference = G0reference,
      G01behaviour = G01behaviour, doprelim = doprelim,
      det_constraint = det_constraint,
      algorithm = prop4ii_algorithm,
      opts = prop4ii_opts,
      ...
    ))
  }

  .mobius_SvMF_definition1(
    y = y, xs = xs, xe = xe, mean = mean, k = k, a = a,
    G0 = G0, G0reference = G0reference,
    G01behaviour = G01behaviour,
    type = type, fix_qs1 = fix_qs1, fix_qe1 = fix_qe1,
    intercept = intercept, doprelim = doprelim, ...
  )
}
# =============================================================================
# MINIMAL PATCH: Proposition 4(ii) spherical part + optional LinEuc covariates
#
# The spherical part is constrained by p = qs and Bs = I_(p-1).  With no
# Euclidean covariates, the link is exactly
#
#     mu(xs) = Rtilde0 xs.
#
# If Euclidean covariates are supplied, the original LinEuc construction is
# retained while Bs is fixed at the identity:
#
#     mu(xs, xe) = P iSp{ Sp(Qs^T xs) + Be Qe*^T xe_tilde },
#
# where Rtilde0 = P Qs^T, xe_tilde contains standardised Euclidean covariates
# and (when intercept = TRUE) an intercept, Qe* has orthonormal columns, and Be
# is diagonal.  Thus setting the Euclidean contribution to zero gives the
# Proposition 4(ii) link exactly.
#
# This section is intentionally appended at the end of the combined source file
# so that the original code is not rewritten.  Only the Proposition 4(ii)
# wrappers are overridden.
# =============================================================================

if (!exists(".mobius_vMF_prop4ii_spherical_only", inherits = FALSE)) {
  .mobius_vMF_prop4ii_spherical_only <- mobius_vMF_prop4ii
}
if (!exists(".mobius_SvMF_prelim_prop4ii_spherical_only", inherits = FALSE)) {
  .mobius_SvMF_prelim_prop4ii_spherical_only <-
    mobius_SvMF_partransport_prelim_prop4ii
}
if (!exists(".mobius_SvMF_joint_prop4ii_spherical_only", inherits = FALSE)) {
  .mobius_SvMF_joint_prop4ii_spherical_only <- mobius_SvMF_joint_fit_prop4ii
}
if (!exists(".mobius_SvMF_prop4ii_spherical_only", inherits = FALSE)) {
  .mobius_SvMF_prop4ii_spherical_only <- mobius_SvMF_prop4ii
}

.prop4ii_has_euc <- function(xe) {
  !is.null(xe) && ncol(as.matrix(xe)) > 0L
}

.prop4ii_vmf_log_norm_const_exact <- function(k, p) {
  if (exists("vMF_log_norm_const_exact", mode = "function")) {
    return(vMF_log_norm_const_exact(k, p))
  }
  k <- as.numeric(k)
  nu <- p / 2 - 1
  out <- numeric(length(k))
  small <- k < 1e-7
  out[small] <- log(2) + (p / 2) * log(pi) - lgamma(p / 2)
  if (any(!small)) {
    kk <- k[!small]
    logI <- log(besselI(kk, nu, expon.scaled = TRUE)) + kk
    out[!small] <- (p / 2) * log(2 * pi) + logI - nu * log(kk)
  }
  out
}

# Replace the earlier helper so the pure-spherical and optional-Euclidean paths
# both work even when the non-exported package helper is unavailable.
.prop4ii_estimate_k <- function(rbar, p) {
  rbar <- min(max(as.numeric(rbar), -1 + 1e-12), 1 - 1e-12)
  objective <- function(k) {
    -.prop4ii_vmf_log_norm_const_exact(k, p) + k * rbar
  }
  fit <- optimise(objective, lower = 1e-8, upper = 1e5, maximum = TRUE)
  list(k = fit$maximum, loglik_per_obs = fit$objective)
}

.prop4ii_prepare_euc <- function(xe, intercept = TRUE,
                                 center = NULL, scale = NULL,
                                 fitting = FALSE) {
  xe <- as.matrix(xe)
  storage.mode(xe) <- "double"
  if (any(!is.finite(xe))) stop("xe contains non-finite values.")
  q <- ncol(xe)

  if (fitting) {
    if (intercept) {
      center <- if (q) colMeans(xe) else numeric(0)
      scale <- if (q) apply(xe, 2, stats::sd) else numeric(0)
      scale[!is.finite(scale) | scale <= sqrt(.Machine$double.eps)] <- 1
    } else {
      center <- rep(0, q)
      scale <- rep(1, q)
    }
  } else {
    if (is.null(center) || length(center) != q) {
      stop("The number of columns of xe does not match the fitted model.")
    }
    if (is.null(scale) || length(scale) != q) {
      stop("The Euclidean scaling information is invalid.")
    }
  }

  xstd <- if (q) {
    sweep(sweep(xe, 2, center, "-"), 2, scale, "/")
  } else {
    matrix(numeric(0), nrow(xe), 0L)
  }
  xmodel <- if (intercept) cbind(xstd, `(Intercept)` = 1) else xstd

  list(
    raw = xe,
    standardised = xstd,
    model = xmodel,
    center = as.numeric(center),
    scale = as.numeric(scale),
    q = q,
    m = ncol(xmodel),
    intercept = isTRUE(intercept)
  )
}

.prop4ii_stiefel_dim <- function(m, k) {
  as.integer(m * k - k * (k + 1L) / 2L)
}

.prop4ii_complete_orthogonal <- function(V) {
  V <- as.matrix(V)
  m <- nrow(V)
  k <- ncol(V)
  if (k > m) stop("A Stiefel matrix cannot have more columns than rows.")
  V <- qr.Q(qr(V), complete = FALSE)[, seq_len(k), drop = FALSE]
  if (k == m) return(V)
  z <- svd(t(V), nu = 0, nv = m)$v
  W <- z[, (k + 1L):m, drop = FALSE]
  Q <- cbind(V, W)
  if (det(Q) < 0) Q[, m] <- -Q[, m]
  Q
}

.prop4ii_stiefel_cayley <- function(theta, Vref) {
  Vref <- as.matrix(Vref)
  m <- nrow(Vref)
  k <- ncol(Vref)
  dA <- k * (k - 1L) / 2L
  dB <- (m - k) * k
  if (length(theta) != dA + dB) stop("Incorrect Stiefel-coordinate length.")

  A <- matrix(0, k, k)
  if (dA > 0L) {
    ij <- which(upper.tri(A), arr.ind = TRUE)
    A[ij] <- theta[seq_len(dA)]
    A[cbind(ij[, 2], ij[, 1])] <- -theta[seq_len(dA)]
  }
  B <- if (dB > 0L) {
    matrix(theta[dA + seq_len(dB)], nrow = m - k, ncol = k)
  } else {
    matrix(numeric(0), 0L, k)
  }

  K <- matrix(0, m, m)
  K[seq_len(k), seq_len(k)] <- A
  if (m > k) {
    K[k + seq_len(m - k), seq_len(k)] <- B
    K[seq_len(k), k + seq_len(m - k)] <- -t(B)
  }

  Qref <- .prop4ii_complete_orthogonal(Vref)
  Qref %*% .prop4ii_cayley(K[upper.tri(K)], m)[, seq_len(k), drop = FALSE]
}

# Extended parameter class.  The original call mobius_link_prop4ii(Rtilde0)
# remains valid; Euclidean fields are added only when xe is used.
mobius_link_prop4ii <- function(Rtilde0, P = NULL, Qe_star = NULL,
                                Be = NULL, xe_center = NULL,
                                xe_scale = NULL, intercept = TRUE,
                                check = TRUE) {
  Rtilde0 <- as.matrix(Rtilde0)
  storage.mode(Rtilde0) <- "double"
  has_euc <- !is.null(P) || !is.null(Qe_star) || !is.null(Be)

  obj <- list(
    Rtilde0 = Rtilde0,
    P = if (has_euc) as.matrix(P) else NULL,
    Qe_star = if (has_euc) as.matrix(Qe_star) else NULL,
    Be = if (has_euc) as.numeric(Be) else NULL,
    xe_center = if (has_euc) as.numeric(xe_center) else NULL,
    xe_scale = if (has_euc) as.numeric(xe_scale) else NULL,
    intercept = if (has_euc) isTRUE(intercept) else FALSE
  )
  class(obj) <- c("mobius_link_prop4ii", "list")
  if (check) mobius_link_prop4ii_check(obj)
  obj
}

mobius_link_prop4ii_check <- function(obj,
                                      tol = 100 * sqrt(.Machine$double.eps)) {
  if (!inherits(obj, "mobius_link_prop4ii")) {
    stop("obj must inherit from 'mobius_link_prop4ii'.")
  }
  R <- obj$Rtilde0
  if (!is.matrix(R) || nrow(R) != ncol(R) || nrow(R) < 2L) {
    stop("Rtilde0 must be a square p by p matrix with p >= 2.")
  }
  if (any(!is.finite(R))) stop("Rtilde0 contains non-finite values.")
  if (max(abs(crossprod(R) - diag(ncol(R)))) > tol) {
    stop("Rtilde0 is not orthogonal.")
  }

  has_euc <- !is.null(obj$P)
  if (!has_euc) {
    if (!is.null(obj$Qe_star) || !is.null(obj$Be)) {
      stop("P, Qe_star, and Be must either all be supplied or all be NULL.")
    }
    return(invisible(NULL))
  }

  p <- nrow(R)
  k <- p - 1L
  P <- obj$P
  V <- obj$Qe_star
  b <- obj$Be
  if (!all(dim(P) == c(p, p))) stop("P must be p by p.")
  if (max(abs(crossprod(P) - diag(p))) > tol || det(P) < 0) {
    stop("P must be a proper orthogonal p by p matrix.")
  }
  if (!is.matrix(V) || ncol(V) != k || nrow(V) < k) {
    stop("Qe_star must be m by (p-1), with m >= p-1.")
  }
  if (max(abs(crossprod(V) - diag(k))) > tol) {
    stop("Qe_star must have orthonormal columns.")
  }
  if (length(b) != k || any(!is.finite(b)) || any(b <= 0) || any(b >= 1)) {
    stop("Be must contain p-1 values strictly between zero and one.")
  }
  q <- nrow(V) - as.integer(isTRUE(obj$intercept))
  if (q < 0 || length(obj$xe_center) != q || length(obj$xe_scale) != q) {
    stop("Euclidean centering/scaling information has the wrong dimension.")
  }
  invisible(NULL)
}

as_mobius_link_prop4ii <- function(obj, check = TRUE) {
  if (inherits(obj, "mobius_link_prop4ii")) {
    if (check) mobius_link_prop4ii_check(obj)
    return(obj)
  }
  if (is.matrix(obj)) return(mobius_link_prop4ii(obj, check = check))
  if (is.list(obj) && !is.null(obj$Rtilde0)) {
    args <- obj[intersect(
      c("Rtilde0", "P", "Qe_star", "Be", "xe_center", "xe_scale", "intercept"),
      names(obj)
    )]
    args$check <- check
    return(do.call(mobius_link_prop4ii, args))
  }
  if (inherits(obj, "mobius_link_cann")) {
    p <- nrow(obj$P)
    if (is.null(obj$Qs) || !all(dim(obj$Qs) == c(p, p))) {
      stop("The canonical object does not have p = qs.")
    }
    if (max(abs(obj$Bs - diag(p - 1L))) >
        100 * sqrt(.Machine$double.eps)) {
      stop("The spherical part requires Bs = I_(p-1).")
    }
    R <- obj$P %*% t(obj$Qs)
    if (is.null(obj$Qe)) return(mobius_link_prop4ii(R, check = check))

    if (!is_LinEuc(obj)) {
      stop("Only the LinEuc Euclidean form is supported with type='Prop4ii'.")
    }
    V <- obj$Qe[-1, -1, drop = FALSE]
    return(mobius_link_prop4ii(
      Rtilde0 = R,
      P = obj$P,
      Qe_star = V,
      Be = diag(obj$Be),
      xe_center = rep(0, nrow(V) - 1L),
      xe_scale = rep(1, nrow(V) - 1L),
      intercept = TRUE,
      check = check
    ))
  }
  stop("obj is not a compatible Proposition 4(ii) parameter object.")
}

print.mobius_link_prop4ii <- function(x, ...) {
  cat("Proposition 4(ii) spherical link")
  if (!is.null(x$P)) cat(" with optional LinEuc covariates")
  cat("\n")
  cat("  p = qs =", nrow(x$Rtilde0), "\n")
  cat("  qe =", if (is.null(x$P)) 0L else length(x$xe_center), "\n")
  cat("  det(Rtilde0) =", format(det(x$Rtilde0), digits = 6), "\n")
  if (!is.null(x$Be)) {
    cat("  Be =", paste(format(x$Be, digits = 5), collapse = ", "), "\n")
  }
  invisible(x)
}

dim.mobius_link_prop4ii <- function(x) {
  p <- nrow(x$Rtilde0)
  c(p = p, qs = p, qe = if (is.null(x$P)) 0L else length(x$xe_center))
}

# Wrapper: Definition 1 classes continue to use the original implementation.
mobius_link <- function(xs = NULL, xe = NULL, param = NULL, check = TRUE) {
  if (!inherits(param, "mobius_link_prop4ii")) {
    return(.mobius_link_definition1(xs = xs, xe = xe,
                                    param = param, check = check))
  }
  if (check) mobius_link_prop4ii_check(param)
  if (is.null(xs)) stop("The Proposition 4(ii) link requires xs.")
  xs <- as.matrix(xs)
  if (ncol(xs) != nrow(param$Rtilde0)) {
    stop("ncol(xs) must equal nrow(Rtilde0).")
  }

  if (is.null(param$P)) {
    if (.prop4ii_has_euc(xe)) {
      stop("This fitted Proposition 4(ii) link has no Euclidean component.")
    }
    return(xs %*% t(param$Rtilde0))
  }

  if (is.null(xe)) stop("xe must be supplied for this fitted model.")
  ep <- .prop4ii_prepare_euc(
    xe, intercept = param$intercept,
    center = param$xe_center, scale = param$xe_scale,
    fitting = FALSE
  )
  if (nrow(xs) != nrow(ep$model)) stop("xs and xe must have the same number of rows.")

  # Rtilde0 = P Qs^T, hence Qs = Rtilde0^T P.
  Qs <- t(param$Rtilde0) %*% param$P
  spherical_part <- Sp(xs %*% Qs)
  euclidean_part <- sweep(ep$model %*% param$Qe_star, 2, param$Be, "*")
  iSp(spherical_part + euclidean_part) %*% t(param$P)
}

.prop4ii_mean_spec_euc <- function(Rref, Pref, Vref, bstart,
                                   xe_center, xe_scale, intercept) {
  p <- nrow(Rref)
  k <- p - 1L
  m <- nrow(Vref)
  dR <- .prop4ii_skew_dim(p)
  dP <- .prop4ii_skew_dim(p)
  dV <- .prop4ii_stiefel_dim(m, k)
  dB <- k

  pos <- 0L
  take <- function(n) {
    if (n == 0L) return(integer(0))
    z <- pos + seq_len(n)
    pos <<- pos + n
    z
  }
  idxR <- take(dR)
  idxP <- take(dP)
  idxV <- take(dV)
  idxB <- take(dB)

  bstart <- pmin(pmax(as.numeric(bstart), 1e-5), 1 - 1e-5)
  par0 <- numeric(pos)
  par0[idxB] <- qlogis(bstart)
  lower <- rep(-Inf, pos)
  upper <- rep( Inf, pos)
  lower[idxB] <- -12
  upper[idxB] <- 12

  list(
    p = p, k = k, m = m,
    Rref = Rref,
    Pref = .prop4ii_nearest_orthogonal(Pref, det_sign = 1),
    Vref = qr.Q(qr(Vref), complete = FALSE)[, seq_len(k), drop = FALSE],
    xe_center = xe_center,
    xe_scale = xe_scale,
    intercept = intercept,
    idxR = idxR, idxP = idxP, idxV = idxV, idxB = idxB,
    par0 = par0, lower = lower, upper = upper
  )
}

.prop4ii_unpack_mean_euc <- function(theta, spec) {
  R <- .prop4ii_cayley(theta[spec$idxR], spec$p) %*% spec$Rref
  P <- .prop4ii_cayley(theta[spec$idxP], spec$p) %*% spec$Pref
  V <- .prop4ii_stiefel_cayley(theta[spec$idxV], spec$Vref)
  b <- plogis(theta[spec$idxB])
  mean <- mobius_link_prop4ii(
    Rtilde0 = R, P = P, Qe_star = V, Be = b,
    xe_center = spec$xe_center, xe_scale = spec$xe_scale,
    intercept = spec$intercept, check = FALSE
  )
  list(mean = mean, Rtilde0 = R, P = P, Qe_star = V, Be = b)
}

.prop4ii_fit_vmf_euc_component <- function(y, xs, xe, spec,
                                           algorithm, opts) {
  eval_f <- function(theta) {
    pars <- .prop4ii_unpack_mean_euc(theta, spec)
    pred <- try(mobius_link(xs = xs, xe = xe, param = pars$mean,
                            check = FALSE), silent = TRUE)
    if (inherits(pred, "try-error") || any(!is.finite(pred))) return(1e100)
    -mean(rowSums(y * pred))
  }

  default_opts <- list(
    algorithm = algorithm,
    xtol_rel = 1e-9,
    maxeval = 10000,
    print_level = 0
  )
  combined_opts <- utils::modifyList(default_opts, opts)
  uses_gradient <- grepl("^NLOPT_LD_", combined_opts$algorithm)
  eval_nloptr <- if (uses_gradient) {
    function(theta) {
      list(
        objective = eval_f(theta),
        gradient = .prop4ii_fd_gradient(
          eval_f, theta, spec$lower, spec$upper
        )
      )
    }
  } else {
    eval_f
  }

  fit <- nloptr::nloptr(
    x0 = spec$par0,
    eval_f = eval_nloptr,
    lb = spec$lower,
    ub = spec$upper,
    opts = combined_opts
  )
  pars <- .prop4ii_unpack_mean_euc(fit$solution, spec)
  pred <- mobius_link(xs = xs, xe = xe, param = pars$mean, check = FALSE)
  list(nlopt = fit, pars = pars, pred = pred,
       objective = -mean(rowSums(y * pred)))
}

mobius_vMF_prop4ii <- function(
    y, xs, xe = NULL,
    det_constraint = c("orthogonal", "rotation"),
    start = NULL, check = TRUE, intercept = TRUE,
    algorithm = "NLOPT_LD_SLSQP", opts = list(), ...) {

  if (!.prop4ii_has_euc(xe)) {
    return(.mobius_vMF_prop4ii_spherical_only(
      y = y, xs = xs, xe = NULL,
      det_constraint = det_constraint, start = start, check = check
    ))
  }

  det_constraint <- match.arg(det_constraint)
  y <- as.matrix(y)
  xs <- as.matrix(xs)
  xe <- as.matrix(xe)
  storage.mode(y) <- storage.mode(xs) <- storage.mode(xe) <- "double"
  if (nrow(y) != nrow(xs) || nrow(y) != nrow(xe) || ncol(y) != ncol(xs)) {
    stop("For type='Prop4ii', y and xs must be n by p and xe must have n rows.")
  }
  p <- ncol(y)
  k <- p - 1L
  if (check) {
    .prop4ii_check_unit_rows(y, "y")
    .prop4ii_check_unit_rows(xs, "xs")
  }
  ep <- .prop4ii_prepare_euc(xe, intercept = intercept, fitting = TRUE)
  if (ep$m < k) {
    stop("With Euclidean covariates, ncol(xe) + intercept must be at least p-1.")
  }

  start_obj <- if (is.null(start)) NULL else as_mobius_link_prop4ii(start)
  target_signs <- if (det_constraint == "rotation") 1 else c(1, -1)
  fits <- vector("list", length(target_signs))

  for (j in seq_along(target_signs)) {
    sgn <- target_signs[j]
    Rref <- if (!is.null(start_obj) &&
                (if (det(start_obj$Rtilde0) >= 0) 1 else -1) == sgn) {
      start_obj$Rtilde0
    } else {
      .prop4ii_procrustes(y, xs, det_sign = sgn)$R
    }
    Pref <- if (!is.null(start_obj) && !is.null(start_obj$P)) {
      start_obj$P
    } else {
      .prop4ii_nearest_orthogonal(Rref, det_sign = 1)
    }
    Vref <- if (!is.null(start_obj) && !is.null(start_obj$Qe_star) &&
                all(dim(start_obj$Qe_star) == c(ep$m, k))) {
      start_obj$Qe_star
    } else {
      diag(ep$m)[, seq_len(k), drop = FALSE]
    }
    bstart <- if (!is.null(start_obj) && !is.null(start_obj$Be)) {
      start_obj$Be
    } else {
      rep(0.05, k)
    }

    spec <- .prop4ii_mean_spec_euc(
      Rref = Rref, Pref = Pref, Vref = Vref, bstart = bstart,
      xe_center = ep$center, xe_scale = ep$scale,
      intercept = intercept
    )
    fits[[j]] <- .prop4ii_fit_vmf_euc_component(
      y = y, xs = xs, xe = xe, spec = spec,
      algorithm = algorithm, opts = opts
    )
    fits[[j]]$det_sign <- sgn
  }

  values <- vapply(fits, function(z) z$objective, numeric(1))
  best <- fits[[which.min(values)]]
  pars <- best$pars
  pred <- best$pred
  dots <- pmax(-1, pmin(1, rowSums(y * pred)))
  kfit <- .prop4ii_estimate_k(mean(dots), p)
  n <- nrow(y)

  dR <- .prop4ii_skew_dim(p)
  dP <- .prop4ii_skew_dim(p)
  dV <- .prop4ii_stiefel_dim(ep$m, k)
  dB <- k
  DoF <- dR + dP + dV + dB + 1L
  lLik <- n * kfit$loglik_per_obs

  rr <- try(rotated_resid(y, pred, north_pole(p)), silent = TRUE)
  rresids <- if (inherits(rr, "try-error")) NULL else rr[, -1, drop = FALSE]
  if (!is.null(rresids)) {
    attr(rresids, "samehemisphere") <- attr(rr, "samehemisphere")
    colnames(rresids) <- paste0("r", seq_len(ncol(rresids)))
  }

  out <- list(
    est = pars$mean,
    mean = pars$mean,
    Rtilde0 = pars$Rtilde0,
    P = pars$P,
    Qe_star = pars$Qe_star,
    Be = pars$Be,
    k = kfit$k,
    obj = best$objective,
    solution = best$nlopt$solution,
    nlopt = best$nlopt,
    y = y, xs = xs, xe = xe,
    pred = pred,
    rresids = rresids,
    dists = acos(dots),
    DoF = DoF,
    AIC = 2 * DoF - 2 * lLik,
    lLik = lLik,
    start = start,
    component_fits = lapply(fits, function(z) {
      list(det_sign = z$det_sign, status = z$nlopt$status,
           message = z$nlopt$message, objective = z$objective)
    }),
    linktype = list(
      type = "Prop4ii", det_constraint = det_constraint,
      intercept = intercept, euclidean = TRUE
    )
  )
  class(out) <- c("mobius_vMF_prop4ii", "mobius_vMF", "list")
  out
}

# Wrapper retaining the original vMF implementation for its original link types.
mobius_vMF <- function(
    y, xs = NULL, xe = NULL, start = NULL,
    type = "SpEuc", fix_qs1 = FALSE,
    fix_qe1 = (type == "LinEuc"), intercept = TRUE,
    lb = NULL, ub = NULL,
    det_constraint = c("orthogonal", "rotation"),
    prop4ii_algorithm = "NLOPT_LD_SLSQP",
    prop4ii_opts = list(), ...) {

  if (.is_prop4ii_type(type)) {
    dots <- list(...)
    if (length(dots)) prop4ii_opts <- utils::modifyList(prop4ii_opts, dots)
    return(mobius_vMF_prop4ii(
      y = y, xs = xs, xe = xe, start = start,
      det_constraint = det_constraint, intercept = intercept,
      algorithm = prop4ii_algorithm, opts = prop4ii_opts
    ))
  }
  .mobius_vMF_definition1(
    y = y, xs = xs, xe = xe, start = start, type = type,
    fix_qs1 = fix_qs1, fix_qe1 = fix_qe1, intercept = intercept,
    lb = lb, ub = ub, ...
  )
}

mobius_SvMF_partransport_prelim_prop4ii <- function(
    y, xs, xe = NULL, mean = NULL, G0 = NULL,
    G01behaviour = "p1",
    det_constraint = c("orthogonal", "rotation"), check = TRUE,
    intercept = TRUE, algorithm = "NLOPT_LD_SLSQP",
    opts = list(), ...) {

  if (!.prop4ii_has_euc(xe)) {
    return(.mobius_SvMF_prelim_prop4ii_spherical_only(
      y = y, xs = xs, xe = NULL, mean = mean, G0 = G0,
      G01behaviour = G01behaviour,
      det_constraint = det_constraint, check = check
    ))
  }

  det_constraint <- match.arg(det_constraint)
  G01behaviour <- match.arg(G01behaviour, c("p1", "free", "fixed"))
  prelim <- mobius_vMF_prop4ii(
    y = y, xs = xs, xe = xe, start = mean,
    det_constraint = det_constraint, check = check,
    intercept = intercept, algorithm = algorithm, opts = opts
  )

  # With Euclidean covariates, p1 is identified as the first column of P.
  p1 <- prelim$mean$P[, 1]
  if (G01behaviour == "fixed" && is.null(G0)) {
    stop("At least the first column of G0 must be supplied when G01behaviour='fixed'.")
  }
  G01 <- switch(
    G01behaviour,
    p1 = p1,
    free = if (is.null(G0)) p1 else G0[, 1],
    fixed = G0[, 1]
  )

  rresid <- rotated_resid(as.matrix(y), prelim$pred, base = G01)
  if (!is.null(G0) && all(!is.na(G0))) {
    if (G01behaviour == "p1") {
      G0 <- cbind(G01, -jupp_Rmat(G0[, 1], G01) %*%
                         G0[, -1, drop = FALSE])
    }
    aremaining <- SvMF_prelim_scales(rresid, G0)
  } else {
    G0 <- SvMF_moment_axes(rresid, G01)
    aremaining <- SvMF_prelim_scales(rresid, G0)
  }

  list(
    mean = prelim$mean,
    Rtilde0 = prelim$Rtilde0,
    P = prelim$P,
    Qe_star = prelim$Qe_star,
    Be = prelim$Be,
    k = prelim$k,
    a = c(1, aremaining),
    G0 = G0,
    nlopt = prelim$nlopt,
    pred = prelim$pred
  )
}

.prop4ii_joint_spec_euc <- function(mean_ref, prelim, G0reference,
                                    G01behaviour) {
  mean_ref <- as_mobius_link_prop4ii(mean_ref)
  p <- nrow(mean_ref$Rtilde0)
  k <- p - 1L
  mspec <- .prop4ii_mean_spec_euc(
    Rref = mean_ref$Rtilde0,
    Pref = mean_ref$P,
    Vref = mean_ref$Qe_star,
    bstart = mean_ref$Be,
    xe_center = mean_ref$xe_center,
    xe_scale = mean_ref$xe_scale,
    intercept = mean_ref$intercept
  )

  dA <- max(p - 2L, 0L)
  dG <- if (G01behaviour == "free") {
    .prop4ii_skew_dim(p)
  } else {
    .prop4ii_skew_dim(p - 1L)
  }
  pos <- length(mspec$par0)
  take <- function(n) {
    if (n == 0L) return(integer(0))
    z <- pos + seq_len(n)
    pos <<- pos + n
    z
  }
  idxK <- take(1L)
  idxA <- take(dA)
  idxG <- take(dG)

  Gstart <- .prop4ii_nearest_orthogonal(prelim$G0, det_sign = 1)
  Gcoord <- if (is.null(G0reference)) Gstart else
    .prop4ii_nearest_orthogonal(G0reference, det_sign = 1)

  if (G01behaviour == "free") {
    Gref <- Gcoord
    thetaG0 <- .prop4ii_inverse_cayley(Gstart %*% t(Gref))
    Ganchor <- NULL
  } else {
    target <- if (G01behaviour == "p1") mean_ref$P[, 1] else Gstart[, 1]
    Ganchor <- .prop4ii_rebase_G0(Gcoord, target)
    if (is.null(Ganchor)) Ganchor <- .prop4ii_rebase_G0(Gstart, target)
    if (is.null(Ganchor)) stop("Could not construct the G0 reference frame.")
    Gstart2 <- .prop4ii_rebase_G0(Gstart, target)
    H0 <- crossprod(Ganchor[, -1, drop = FALSE],
                    Gstart2[, -1, drop = FALSE])
    H0 <- .prop4ii_nearest_orthogonal(H0, det_sign = 1)
    thetaG0 <- .prop4ii_inverse_cayley(H0)
    Gref <- NULL
  }

  par0 <- c(mspec$par0, rep(0, pos - length(mspec$par0)))
  par0[idxK] <- log(max(prelim$k, 1e-8))
  if (dA > 0L) par0[idxA] <- .prop4ii_scales_to_eta(prelim$a)
  if (dG > 0L) par0[idxG] <- thetaG0
  lower <- c(mspec$lower, rep(-Inf, pos - length(mspec$lower)))
  upper <- c(mspec$upper, rep( Inf, pos - length(mspec$upper)))
  lower[idxK] <- log(1e-8)
  upper[idxK] <- log(1e5)
  if (dA > 0L) {
    lower[idxA] <- -10
    upper[idxA] <- 10
  }

  list(
    mean_spec = mspec,
    mean_length = length(mspec$par0),
    p = p,
    a1 = prelim$a[1],
    G01behaviour = G01behaviour,
    Gref = Gref,
    Ganchor = Ganchor,
    idxK = idxK, idxA = idxA, idxG = idxG,
    par0 = par0, lower = lower, upper = upper
  )
}

.prop4ii_unpack_joint_euc <- function(theta, spec) {
  mp <- .prop4ii_unpack_mean_euc(
    theta[seq_len(spec$mean_length)], spec$mean_spec
  )
  kappa <- exp(theta[spec$idxK])
  eta <- if (length(spec$idxA)) theta[spec$idxA] else numeric(0)
  a <- .prop4ii_eta_to_scales(eta, spec$a1, spec$p)

  if (spec$G01behaviour == "free") {
    G0 <- .prop4ii_cayley(theta[spec$idxG], spec$p) %*% spec$Gref
  } else {
    Gbase <- spec$Ganchor
    if (spec$G01behaviour == "p1") {
      Gbase <- .prop4ii_rebase_G0(spec$Ganchor, mp$P[, 1])
      if (is.null(Gbase)) return(NULL)
    }
    H <- if (spec$p == 2L) matrix(1, 1, 1) else
      .prop4ii_cayley(theta[spec$idxG], spec$p - 1L)
    B <- diag(spec$p)
    B[-1, -1] <- H
    G0 <- Gbase %*% B
  }

  c(mp, list(k = kappa, a = a, G0 = G0))
}

.prop4ii_joint_loglik_euc <- function(y, xs, xe, mean, k, a, G0,
                                      approximate = TRUE) {
  if (!is.finite(k) || k <= 0 || any(!is.finite(a)) || any(a <= 0)) {
    return(-Inf)
  }
  pred <- try(mobius_link(xs = xs, xe = xe, param = mean,
                          check = FALSE), silent = TRUE)
  if (inherits(pred, "try-error") || any(!is.finite(pred))) return(-Inf)
  yrot <- try(undo_partransport(y, pred, G01 = G0[, 1]), silent = TRUE)
  if (inherits(yrot, "try-error") || any(!is.finite(yrot))) return(-Inf)
  # Same numerical safeguard as in the spherical-only Proposition 4(ii)
  # likelihood: use the fast approximation when finite and otherwise fall back
  # to the exact R density.
  ld <- NULL
  if (approximate && exists("ldSvMF_cann", mode = "function")) {
    ld_try <- try(ldSvMF_cann(yrot, k = k, a = a, G = G0), silent = TRUE)
    if (!inherits(ld_try, "try-error") && all(is.finite(ld_try))) ld <- ld_try
  }
  if (is.null(ld)) {
    ld_try <- try(
      SvMF_log_lik_cann(yrot, SvMF_cann(k = k, a = a, G = G0)),
      silent = TRUE
    )
    if (!inherits(ld_try, "try-error") && all(is.finite(ld_try))) ld <- ld_try
  }
  if (is.null(ld)) return(-Inf)
  sum(ld)
}

.prop4ii_fit_joint_euc_component <- function(
    y, xs, xe, mean_ref, prelim, G0reference, G01behaviour,
    algorithm, opts) {

  spec <- .prop4ii_joint_spec_euc(
    mean_ref = mean_ref, prelim = prelim,
    G0reference = G0reference,
    G01behaviour = G01behaviour
  )
  eval_f <- function(theta) {
    pars <- .prop4ii_unpack_joint_euc(theta, spec)
    if (is.null(pars)) return(1e100)
    ll <- .prop4ii_joint_loglik_euc(
      y, xs, xe, pars$mean, pars$k, pars$a, pars$G0,
      approximate = TRUE
    )
    if (!is.finite(ll)) return(1e100)
    -ll / nrow(y)
  }

  default_opts <- list(
    algorithm = algorithm,
    xtol_rel = 1e-9,
    maxeval = 10000,
    print_level = 0
  )
  combined_opts <- utils::modifyList(default_opts, opts)
  uses_gradient <- grepl("^NLOPT_LD_", combined_opts$algorithm)
  eval_nloptr <- if (uses_gradient) {
    function(theta) {
      list(
        objective = eval_f(theta),
        gradient = .prop4ii_fd_gradient(
          eval_f, theta, spec$lower, spec$upper
        )
      )
    }
  } else {
    eval_f
  }

  fit <- nloptr::nloptr(
    x0 = spec$par0,
    eval_f = eval_nloptr,
    lb = spec$lower,
    ub = spec$upper,
    opts = combined_opts
  )
  pars <- .prop4ii_unpack_joint_euc(fit$solution, spec)
  ll <- if (is.null(pars)) -Inf else .prop4ii_joint_loglik_euc(
    y, xs, xe, pars$mean, pars$k, pars$a, pars$G0,
    approximate = TRUE
  )
  list(nlopt = fit, pars = pars, lLik_approx = ll, spec = spec)
}

mobius_SvMF_joint_fit_prop4ii <- function(
    y, xs, xe = NULL, mean, k, a, G0,
    G0reference = NULL, G01behaviour = "p1",
    det_constraint = c("orthogonal", "rotation"),
    algorithm = "NLOPT_LD_SLSQP", opts = list(), check = TRUE,
    intercept = TRUE) {

  if (!.prop4ii_has_euc(xe)) {
    return(.mobius_SvMF_joint_prop4ii_spherical_only(
      y = y, xs = xs, xe = NULL, mean = mean, k = k, a = a, G0 = G0,
      G0reference = G0reference, G01behaviour = G01behaviour,
      det_constraint = det_constraint, algorithm = algorithm,
      opts = opts, check = check
    ))
  }

  det_constraint <- match.arg(det_constraint)
  G01behaviour <- match.arg(G01behaviour, c("p1", "free", "fixed"))
  y <- as.matrix(y)
  xs <- as.matrix(xs)
  xe <- as.matrix(xe)
  p <- ncol(y)
  if (nrow(y) != nrow(xs) || nrow(y) != nrow(xe) || ncol(xs) != p) {
    stop("For type='Prop4ii', y and xs must be n by p and xe must have n rows.")
  }
  if (check) {
    .prop4ii_check_unit_rows(y, "y")
    .prop4ii_check_unit_rows(xs, "xs")
  }

  start_mean <- as_mobius_link_prop4ii(mean)
  if (is.null(start_mean$P)) stop("The starting mean has no Euclidean component.")
  start_a <- .prop4ii_normalise_scales(a, p, a1 = a[1])
  start_G0 <- .prop4ii_nearest_orthogonal(G0, det_sign = 1)
  initial <- list(mean = start_mean, Rtilde0 = start_mean$Rtilde0,
                  P = start_mean$P, Qe_star = start_mean$Qe_star,
                  Be = start_mean$Be, k = k, a = start_a, G0 = start_G0)

  target_signs <- if (det_constraint == "rotation") 1 else c(1, -1)
  fits <- vector("list", length(target_signs))
  for (j in seq_along(target_signs)) {
    sgn <- target_signs[j]
    Rref <- if ((if (det(start_mean$Rtilde0) >= 0) 1 else -1) == sgn) {
      start_mean$Rtilde0
    } else {
      .prop4ii_procrustes(y, xs, det_sign = sgn)$R
    }
    mean_ref <- mobius_link_prop4ii(
      Rtilde0 = Rref,
      P = start_mean$P,
      Qe_star = start_mean$Qe_star,
      Be = start_mean$Be,
      xe_center = start_mean$xe_center,
      xe_scale = start_mean$xe_scale,
      intercept = start_mean$intercept,
      check = FALSE
    )
    pre_j <- initial
    pre_j$mean <- mean_ref
    pre_j$Rtilde0 <- Rref
    if (G01behaviour == "p1") {
      rebased <- .prop4ii_rebase_G0(pre_j$G0, mean_ref$P[, 1])
      if (is.null(rebased)) stop("Could not rebase G0 to P[,1].")
      pre_j$G0 <- rebased
    }
    fits[[j]] <- .prop4ii_fit_joint_euc_component(
      y = y, xs = xs, xe = xe,
      mean_ref = mean_ref, prelim = pre_j,
      G0reference = if (is.null(G0reference)) pre_j$G0 else G0reference,
      G01behaviour = G01behaviour,
      algorithm = algorithm, opts = opts
    )
    fits[[j]]$det_sign <- sgn
  }

  lla <- vapply(fits, function(z) z$lLik_approx, numeric(1))
  if (all(!is.finite(lla))) stop("All Proposition 4(ii)+Euclidean optimisations failed.")
  best <- fits[[which.max(lla)]]
  pars <- best$pars

  # Original residual-axis identifiability convention.
  aord <- order(pars$a[-1], decreasing = TRUE)
  pars$a[-1] <- pars$a[-1][aord]
  pars$G0[, -1] <- pars$G0[, -1, drop = FALSE][, aord, drop = FALSE]
  if (exists("standardise_col_signs", mode = "function")) {
    pars$G0[, -1] <- standardise_col_signs(pars$G0[, -1, drop = FALSE])
  }
  if (det(pars$G0) < 0) pars$G0[, p] <- -pars$G0[, p]

  pred <- mobius_link(xs = xs, xe = xe, param = pars$mean, check = FALSE)
  if (p != 3L) {
    kres <- mobius_SvMF_konly(y = y, ymean = pred,
                              a = pars$a, G0 = pars$G0)
    pars$k <- kres$k
    lLik <- kres$lLik
  } else {
    yrot <- undo_partransport(y, pred, G01 = pars$G0[, 1])
    lLik <- sum(SvMF_log_lik_cann(
      yrot, SvMF_cann(k = pars$k, a = pars$a, G = pars$G0)
    ))
  }
  SvMF_cann_check(SvMF_cann(k = pars$k, a = pars$a, G = pars$G0))

  dots <- pmax(-1, pmin(1, rowSums(y * pred)))
  rrI <- rotated_resid(y, pred, north_pole(p))
  rresids_I <- rrI[, -1, drop = FALSE]
  attr(rresids_I, "samehemisphere") <- attr(rrI, "samehemisphere")
  colnames(rresids_I) <- paste0("r", seq_len(ncol(rresids_I)))
  rresids_G0 <- resid_SvMF_partransport(y, pred, G0 = pars$G0, scale = FALSE)
  rresids_std <- resid_SvMF_partransport(
    y, pred, pars$k, pars$a, pars$G0, scale = TRUE
  )

  ep_m <- nrow(pars$mean$Qe_star)
  kmean <- p - 1L
  dmean <- .prop4ii_skew_dim(p) +       # Rtilde0
    .prop4ii_skew_dim(p) +              # P
    .prop4ii_stiefel_dim(ep_m, kmean) + # Qe_star
    kmean                               # Be
  dscale <- max(p - 2L, 0L)
  dG0 <- if (G01behaviour == "free") {
    .prop4ii_skew_dim(p)
  } else {
    .prop4ii_skew_dim(p - 1L)
  }
  DoF <- dmean + 1L + dscale + dG0

  out <- list(
    mean = pars$mean,
    est = pars$mean,
    Rtilde0 = pars$Rtilde0,
    P = pars$P,
    Qe_star = pars$Qe_star,
    Be = pars$Be,
    k = pars$k,
    a = pars$a,
    G0 = pars$G0,
    obj = -lLik / nrow(y),
    nlopt = best$nlopt,
    y = y, xs = xs, xe = xe,
    pred = pred,
    rresids_I = rresids_I,
    rresids_G0 = rresids_G0,
    rresids_std = rresids_std,
    dists = acos(dots),
    DoF = DoF,
    AIC = 2 * DoF - 2 * lLik,
    lLik = lLik,
    initial = initial,
    component_fits = lapply(fits, function(z) {
      list(det_sign = z$det_sign, status = z$nlopt$status,
           message = z$nlopt$message, lLik_approx = z$lLik_approx)
    }),
    linktype = list(
      type = "Prop4ii", det_constraint = det_constraint,
      G01behaviour = G01behaviour,
      intercept = intercept, euclidean = TRUE
    )
  )
  class(out) <- c("mobius_SvMF_prop4ii", "mobius_SvMF", "list")
  out
}

mobius_SvMF_prop4ii <- function(
    y, xs, xe = NULL, mean = NULL, k = NULL, a = NULL,
    G0 = NULL, G0reference = NULL, G01behaviour = "p1",
    doprelim = TRUE,
    det_constraint = c("orthogonal", "rotation"),
    algorithm = "NLOPT_LD_SLSQP", opts = list(), check = TRUE,
    intercept = TRUE, ...) {

  if (!.prop4ii_has_euc(xe)) {
    return(.mobius_SvMF_prop4ii_spherical_only(
      y = y, xs = xs, xe = NULL, mean = mean, k = k, a = a,
      G0 = G0, G0reference = G0reference,
      G01behaviour = G01behaviour, doprelim = doprelim,
      det_constraint = det_constraint, algorithm = algorithm,
      opts = utils::modifyList(opts, list(...)), check = check
    ))
  }

  det_constraint <- match.arg(det_constraint)
  G01behaviour <- match.arg(G01behaviour, c("p1", "free", "fixed"))
  dots <- list(...)
  if (length(dots)) opts <- utils::modifyList(opts, dots)

  if (doprelim) {
    preest <- mobius_SvMF_partransport_prelim_prop4ii(
      y = y, xs = xs, xe = xe, mean = mean, G0 = G0,
      G01behaviour = G01behaviour,
      det_constraint = det_constraint, check = check,
      intercept = intercept, algorithm = algorithm, opts = opts
    )
  } else {
    if (is.null(mean) || is.null(k) || is.null(a) || is.null(G0)) {
      stop("When doprelim=FALSE, mean, k, a, and G0 must all be supplied.")
    }
    mean <- as_mobius_link_prop4ii(mean)
    if (is.null(mean$P)) stop("mean must contain the Euclidean component.")
    preest <- list(
      mean = mean, Rtilde0 = mean$Rtilde0,
      P = mean$P, Qe_star = mean$Qe_star, Be = mean$Be,
      k = k, a = a, G0 = G0
    )
  }

  finalest <- mobius_SvMF_joint_fit_prop4ii(
    y = y, xs = xs, xe = xe,
    mean = preest$mean,
    k = if (is.null(k)) preest$k else k,
    a = if (is.null(a)) preest$a else a,
    G0 = preest$G0,
    G0reference = if (is.null(G0reference)) preest$G0 else G0reference,
    G01behaviour = G01behaviour,
    det_constraint = det_constraint,
    algorithm = algorithm,
    opts = opts,
    check = check,
    intercept = intercept
  )
  finalest$preest <- preest
  finalest
}

print.mobius_vMF_prop4ii <- function(x, ...) {
  cat("vMF regression with Proposition 4(ii) spherical link")
  if (isTRUE(x$linktype$euclidean)) cat(" and LinEuc covariates")
  cat("\n")
  cat("  det(Rtilde0) =", format(det(x$Rtilde0), digits = 6), "\n")
  if (!is.null(x$Be)) cat("  Be =", paste(format(x$Be, digits = 5), collapse = ", "), "\n")
  cat("  k =", format(x$k, digits = 6), "\n")
  cat("  log-likelihood =", format(x$lLik, digits = 7), "\n")
  cat("  AIC =", format(x$AIC, digits = 7), "\n")
  invisible(x)
}

print.mobius_SvMF_prop4ii <- function(x, ...) {
  cat("SvMF regression with Proposition 4(ii) spherical link")
  if (isTRUE(x$linktype$euclidean)) cat(" and LinEuc covariates")
  cat("\n")
  cat("  G0 estimation: original parallel-transport/Cayley procedure\n")
  cat("  G01behaviour =", x$linktype$G01behaviour, "\n")
  cat("  det(Rtilde0) =", format(det(x$Rtilde0), digits = 6), "\n")
  if (!is.null(x$Be)) cat("  Be =", paste(format(x$Be, digits = 5), collapse = ", "), "\n")
  cat("  k =", format(x$k, digits = 6), "\n")
  cat("  a =", paste(format(x$a, digits = 5), collapse = ", "), "\n")
  cat("  log-likelihood =", format(x$lLik, digits = 7), "\n")
  cat("  AIC =", format(x$AIC, digits = 7), "\n")
  invisible(x)
}

predict.mobius_vMF_prop4ii <- function(object, newdata = NULL,
                                        xs = NULL, xe = NULL, ...) {
  if (is.list(newdata) && is.null(xs)) {
    xs <- newdata$xs
    if (is.null(xe)) xe <- newdata$xe
  } else if (!is.null(newdata) && is.null(xs)) {
    xs <- newdata
  }
  if (is.null(xs)) xs <- object$xs
  if (is.null(xe) && !is.null(object$mean$P)) xe <- object$xe
  mobius_link(xs = as.matrix(xs), xe = xe, param = object$mean)
}

predict.mobius_SvMF_prop4ii <- function(object, newdata = NULL,
                                         xs = NULL, xe = NULL, ...) {
  predict.mobius_vMF_prop4ii(object, newdata = newdata, xs = xs, xe = xe, ...)
}

# Final wrapper.  The user-facing command remains type = "Prop4ii".  Passing
# xe includes the optional LinEuc component; xe = NULL fits the original
# spherical-only Proposition 4(ii) model.
mobius_SvMF <- function(
    y, xs = NULL, xe = NULL, mean = NULL, k = NULL, a = NULL,
    G0 = NULL, G0reference = NULL, G01behaviour = "p1",
    type = "LinEuc", fix_qs1 = FALSE,
    fix_qe1 = (type == "LinEuc"), intercept = TRUE,
    doprelim = TRUE,
    det_constraint = c("orthogonal", "rotation"),
    prop4ii_algorithm = "NLOPT_LD_SLSQP",
    prop4ii_opts = list(), ...) {

  if (.is_prop4ii_type(type)) {
    return(mobius_SvMF_prop4ii(
      y = y, xs = xs, xe = xe, mean = mean, k = k, a = a,
      G0 = G0, G0reference = G0reference,
      G01behaviour = G01behaviour, doprelim = doprelim,
      det_constraint = det_constraint,
      algorithm = prop4ii_algorithm,
      opts = prop4ii_opts,
      intercept = intercept,
      ...
    ))
  }

  .mobius_SvMF_definition1(
    y = y, xs = xs, xe = xe, mean = mean, k = k, a = a,
    G0 = G0, G0reference = G0reference,
    G01behaviour = G01behaviour,
    type = type, fix_qs1 = fix_qs1, fix_qe1 = fix_qe1,
    intercept = intercept, doprelim = doprelim, ...
  )
}

# =============================================================================
# Proposition 4(i) extension: spherical-only isotropic scaled Möbius link
#
# Paper parameterisation (strict Proposition 4(i)):
#   qe = 0, p = qs, Bs = beta_s I_(p-1), 0 < beta_s <= 1,
#
#   mu_s(xs) = M_S(xs : Rtilde0, psi),
#   psi = phi r_s1,   phi = (1 - beta_s)/(1 + beta_s),
#
# where Rtilde0 is orthogonal and ||psi|| < 1.  beta_s = 1 gives
# Proposition 4(ii).  We optimise psi through an unconstrained coordinate u,
#   psi = u / sqrt(1 + ||u||^2),
# so that ||psi|| < 1 automatically.
#
# This block is appended after the Proposition 4(ii) implementation so the
# existing type="Prop4ii" behaviour is unchanged.  The new user-facing type is
# type="Prop4i" (or "Proposition4i").
# =============================================================================

.is_prop4i_type <- function(type) {
  if (length(type) != 1L || is.na(type)) return(FALSE)
  key <- tolower(gsub("[^[:alnum:]]", "", as.character(type)))
  key %in% c("prop4i", "proposition4i")
}

.prop4i_u_to_psi <- function(u) {
  u <- as.numeric(u)
  den <- sqrt(1 + sum(u^2))
  u / den
}

.prop4i_psi_to_u <- function(psi, tol = 1e-12) {
  psi <- as.numeric(psi)
  r2 <- sum(psi^2)
  if (!is.finite(r2) || r2 >= 1) {
    r2 <- min(r2, 1 - tol)
    psi <- psi * sqrt(r2 / max(sum(psi^2), tol))
  }
  psi / sqrt(max(1 - sum(psi^2), tol))
}

.prop4i_beta_from_psi <- function(psi) {
  phi <- sqrt(sum(as.numeric(psi)^2))
  (1 - phi) / (1 + phi)
}

.prop4i_reference_direction <- function(mean, tol = 1e-10) {
  mean <- as_mobius_link_prop4i(mean, check = FALSE)
  phi <- sqrt(sum(mean$psi^2))
  if (phi > tol) {
    rs1 <- mean$psi / phi
    drop(mean$Rtilde0 %*% rs1)
  } else {
    # At beta_s = 1, r_s1 is not identifiable.  Use the same representative
    # convention as the Proposition 4(ii) implementation.
    mean$Rtilde0[, 1]
  }
}

mobius_link_prop4i <- function(Rtilde0, psi = NULL,
                               beta_s = NULL, rs1 = NULL,
                               check = TRUE) {
  Rtilde0 <- as.matrix(Rtilde0)
  storage.mode(Rtilde0) <- "double"
  p <- nrow(Rtilde0)

  if (is.null(psi)) {
    if (is.null(beta_s)) beta_s <- 1
    beta_s <- as.numeric(beta_s)[1]
    if (!is.finite(beta_s) || beta_s <= 0 || beta_s > 1) {
      stop("beta_s must satisfy 0 < beta_s <= 1.")
    }
    phi <- (1 - beta_s) / (1 + beta_s)
    if (is.null(rs1)) {
      rs1 <- c(1, rep(0, p - 1L))
    }
    rs1 <- as.numeric(rs1)
    if (length(rs1) != p) stop("rs1 must have length p.")
    nr <- sqrt(sum(rs1^2))
    if (!is.finite(nr) || nr <= 0) stop("rs1 must be nonzero and finite.")
    rs1 <- rs1 / nr
    psi <- phi * rs1
  } else {
    psi <- as.numeric(psi)
  }

  phi <- sqrt(sum(psi^2))
  rs1_out <- if (phi > 1e-10) psi / phi else rep(NA_real_, p)
  obj <- list(
    Rtilde0 = Rtilde0,
    psi = psi,
    phi = phi,
    beta_s = (1 - phi) / (1 + phi),
    rs1 = rs1_out
  )
  class(obj) <- c("mobius_link_prop4i", "list")
  if (check) mobius_link_prop4i_check(obj)
  obj
}

mobius_link_prop4i_check <- function(obj,
                                     tol = 100 * sqrt(.Machine$double.eps)) {
  if (!inherits(obj, "mobius_link_prop4i")) {
    stop("obj must inherit from 'mobius_link_prop4i'.")
  }
  R <- obj$Rtilde0
  if (!is.matrix(R) || nrow(R) != ncol(R) || nrow(R) < 2L) {
    stop("Rtilde0 must be a square p by p matrix with p >= 2.")
  }
  if (any(!is.finite(R))) stop("Rtilde0 contains non-finite values.")
  err <- max(abs(crossprod(R) - diag(ncol(R))))
  if (err > tol) {
    stop(sprintf("Rtilde0 is not orthogonal: max|t(R)R-I| = %.3e.", err))
  }
  psi <- obj$psi
  if (!is.numeric(psi) || length(psi) != nrow(R) || any(!is.finite(psi))) {
    stop("psi must be a finite numeric vector of length p.")
  }
  if (sum(psi^2) >= 1 - 1e-12) {
    stop("Proposition 4(i) requires ||psi|| < 1 (equivalently beta_s > 0).")
  }
  invisible(NULL)
}

as_mobius_link_prop4i <- function(obj, check = TRUE) {
  if (inherits(obj, "mobius_link_prop4i")) {
    if (check) mobius_link_prop4i_check(obj)
    return(obj)
  }
  if (is.list(obj) && !is.null(obj$Rtilde0)) {
    if (!is.null(obj$psi)) {
      return(mobius_link_prop4i(obj$Rtilde0, psi = obj$psi, check = check))
    }
    if (!is.null(obj$beta_s)) {
      return(mobius_link_prop4i(obj$Rtilde0, beta_s = obj$beta_s,
                                rs1 = obj$rs1, check = check))
    }
  }
  stop("obj must be a mobius_link_prop4i object or a compatible list.")
}

print.mobius_link_prop4i <- function(x, ...) {
  cat("Proposition 4(i) spherical Möbius link\n")
  cat("  p = qs =", nrow(x$Rtilde0), ", qe = 0\n")
  cat("  beta_s =", format(x$beta_s, digits = 6), "\n")
  cat("  phi =", format(x$phi, digits = 6), "\n")
  cat("  det(Rtilde0) =", format(det(x$Rtilde0), digits = 6), "\n")
  invisible(x)
}

dim.mobius_link_prop4i <- function(x) {
  p <- nrow(x$Rtilde0)
  c(p = p, qs = p, qe = 0L)
}

.prop4i_map <- function(xs, Rtilde0, psi) {
  xs <- as.matrix(xs)
  psi <- as.numeric(psi)
  z <- sweep(xs, 2L, psi, "+")
  den <- rowSums(z^2)
  if (any(!is.finite(den)) || any(den <= .Machine$double.eps)) {
    stop("Numerical singularity in the Proposition 4(i) Möbius map.")
  }
  fac <- (1 - sum(psi^2)) / den
  inner <- sweep(z, 1L, fac, "*")
  inner <- sweep(inner, 2L, psi, "+")
  pred <- inner %*% t(Rtilde0)
  nr <- sqrt(rowSums(pred^2))
  sweep(pred, 1L, nr, "/")
}

# Save the fully working Proposition 4(ii)+optional-Euclidean dispatcher before
# adding the new Proposition 4(i) class.
.mobius_link_before_prop4i <- mobius_link
.mobius_vMF_before_prop4i <- mobius_vMF
.mobius_SvMF_before_prop4i <- mobius_SvMF

mobius_link <- function(xs = NULL, xe = NULL, param = NULL, check = TRUE) {
  if (inherits(param, "mobius_link_prop4i")) {
    if (check) mobius_link_prop4i_check(param)
    if (is.null(xs)) stop("The Proposition 4(i) link requires xs.")
    if (!is.null(xe) && ncol(as.matrix(xe)) > 0L) {
      stop("Strict Proposition 4(i) has qe = 0; xe must be NULL.")
    }
    xs <- as.matrix(xs)
    if (ncol(xs) != nrow(param$Rtilde0)) {
      stop("ncol(xs) must equal nrow(Rtilde0).")
    }
    return(.prop4i_map(xs, param$Rtilde0, param$psi))
  }
  .mobius_link_before_prop4i(xs = xs, xe = xe, param = param, check = check)
}

# -----------------------------------------------------------------------------
# vMF preliminary fit for Proposition 4(i)
# -----------------------------------------------------------------------------

.prop4i_mean_spec <- function(Rref, psi_start = NULL) {
  p <- nrow(Rref)
  dR <- .prop4ii_skew_dim(p)
  idxR <- seq_len(dR)
  idxPsi <- dR + seq_len(p)
  par0 <- numeric(dR + p)
  if (!is.null(psi_start)) par0[idxPsi] <- .prop4i_psi_to_u(psi_start)
  list(
    p = p, Rref = Rref,
    idxR = idxR, idxPsi = idxPsi,
    par0 = par0,
    lower = rep(-Inf, length(par0)),
    upper = rep( Inf, length(par0))
  )
}

.prop4i_unpack_mean <- function(theta, spec) {
  R <- .prop4ii_cayley(theta[spec$idxR], spec$p) %*% spec$Rref
  psi <- .prop4i_u_to_psi(theta[spec$idxPsi])
  mean <- mobius_link_prop4i(R, psi = psi, check = FALSE)
  list(Rtilde0 = R, psi = psi, mean = mean)
}

.prop4i_fit_vmf_component <- function(y, xs, spec,
                                      algorithm = "NLOPT_LD_SLSQP",
                                      opts = list()) {
  eval_f <- function(theta) {
    pars <- .prop4i_unpack_mean(theta, spec)
    pred <- try(.prop4i_map(xs, pars$Rtilde0, pars$psi), silent = TRUE)
    if (inherits(pred, "try-error") || any(!is.finite(pred))) return(1e100)
    -mean(rowSums(y * pred))
  }

  default_opts <- list(
    algorithm = algorithm,
    xtol_rel = 1e-9,
    maxeval = 10000,
    print_level = 0
  )
  combined_opts <- utils::modifyList(default_opts, opts)
  uses_gradient <- grepl("^NLOPT_LD_", combined_opts$algorithm)
  eval_nloptr <- if (uses_gradient) {
    function(theta) {
      list(
        objective = eval_f(theta),
        gradient = .prop4ii_fd_gradient(
          eval_f, theta, spec$lower, spec$upper
        )
      )
    }
  } else {
    eval_f
  }

  fit <- try(
    nloptr::nloptr(
      x0 = spec$par0,
      eval_f = eval_nloptr,
      lb = spec$lower,
      ub = spec$upper,
      opts = combined_opts
    ),
    silent = TRUE
  )

  # The starting point (Proposition 4(ii), psi=0 by default) is always valid.
  # Retain it if an optimiser fails unexpectedly.
  theta_best <- spec$par0
  f_best <- eval_f(theta_best)
  if (!inherits(fit, "try-error") && is.finite(fit$objective) &&
      fit$objective < f_best) {
    theta_best <- fit$solution
    f_best <- fit$objective
  }
  pars <- .prop4i_unpack_mean(theta_best, spec)
  pred <- .prop4i_map(xs, pars$Rtilde0, pars$psi)

  if (inherits(fit, "try-error")) {
    fit <- list(status = -999L, message = as.character(fit),
                objective = f_best, solution = theta_best)
  } else {
    fit$objective <- f_best
    fit$solution <- theta_best
  }

  list(nlopt = fit, pars = pars, pred = pred, objective = f_best)
}

mobius_vMF_prop4i <- function(
    y, xs, xe = NULL, start = NULL,
    det_constraint = c("orthogonal", "rotation"),
    algorithm = "NLOPT_LD_SLSQP", opts = list(), check = TRUE) {

  det_constraint <- match.arg(det_constraint)
  y <- as.matrix(y)
  xs <- as.matrix(xs)
  storage.mode(y) <- storage.mode(xs) <- "double"

  if (!is.null(xe) && ncol(as.matrix(xe)) > 0L) {
    stop("Strict Proposition 4(i) has qe = 0; xe must be NULL.")
  }
  if (nrow(y) != nrow(xs) || ncol(y) != ncol(xs)) {
    stop("For Proposition 4(i), y and xs must have identical n by p dimensions.")
  }
  if (check) {
    .prop4ii_check_unit_rows(y, "y")
    .prop4ii_check_unit_rows(xs, "xs")
  }

  start_obj <- if (is.null(start)) NULL else as_mobius_link_prop4i(start, check = FALSE)
  start_R <- if (is.null(start_obj)) NULL else start_obj$Rtilde0
  start_psi <- if (is.null(start_obj)) rep(0, ncol(y)) else start_obj$psi

  target_signs <- if (det_constraint == "rotation") 1 else c(1, -1)
  fits <- vector("list", length(target_signs))
  for (j in seq_along(target_signs)) {
    sgn <- target_signs[j]
    Rref <- if (!is.null(start_R) &&
                (if (det(start_R) >= 0) 1 else -1) == sgn) {
      start_R
    } else {
      .prop4ii_procrustes(y, xs, det_sign = sgn)$R
    }
    spec <- .prop4i_mean_spec(Rref, psi_start = start_psi)
    fits[[j]] <- .prop4i_fit_vmf_component(
      y = y, xs = xs, spec = spec,
      algorithm = algorithm, opts = opts
    )
    fits[[j]]$det_sign <- sgn
  }

  values <- vapply(fits, function(z) z$objective, numeric(1))
  best <- fits[[which.min(values)]]
  pars <- best$pars
  pred <- best$pred
  dots <- pmax(-1, pmin(1, rowSums(y * pred)))
  kfit <- .prop4ii_estimate_k(mean(dots), ncol(y))
  n <- nrow(y)
  p <- ncol(y)
  dmean <- .prop4ii_skew_dim(p) + p
  DoF <- dmean + 1L
  lLik <- n * kfit$loglik_per_obs

  rr <- try(rotated_resid(y, pred, north_pole(p)), silent = TRUE)
  rresids <- if (inherits(rr, "try-error")) NULL else rr[, -1, drop = FALSE]
  if (!is.null(rresids)) {
    attr(rresids, "samehemisphere") <- attr(rr, "samehemisphere")
    colnames(rresids) <- paste0("r", seq_len(ncol(rresids)))
  }

  mean <- pars$mean
  out <- list(
    est = mean, mean = mean,
    Rtilde0 = mean$Rtilde0,
    psi = mean$psi, phi = mean$phi,
    beta_s = mean$beta_s, rs1 = mean$rs1,
    k = kfit$k,
    obj = best$objective,
    solution = best$nlopt$solution,
    nlopt = best$nlopt,
    y = y, xs = xs, xe = NULL,
    pred = pred, rresids = rresids,
    dists = acos(dots),
    DoF = DoF, AIC = 2 * DoF - 2 * lLik, lLik = lLik,
    start = start,
    component_fits = lapply(fits, function(z) {
      list(det_sign = z$det_sign, status = z$nlopt$status,
           message = z$nlopt$message, objective = z$objective)
    }),
    linktype = list(type = "Prop4i", det_constraint = det_constraint,
                    intercept = FALSE, euclidean = FALSE)
  )
  class(out) <- c("mobius_vMF_prop4i", "mobius_vMF", "list")
  out
}

# -----------------------------------------------------------------------------
# SvMF fit for Proposition 4(i)
# -----------------------------------------------------------------------------

mobius_SvMF_partransport_prelim_prop4i <- function(
    y, xs, xe = NULL, mean = NULL, G0 = NULL,
    G01behaviour = "p1",
    det_constraint = c("orthogonal", "rotation"),
    algorithm = "NLOPT_LD_SLSQP", opts = list(), check = TRUE, ...) {

  det_constraint <- match.arg(det_constraint)
  G01behaviour <- match.arg(G01behaviour, c("p1", "free", "fixed"))
  y <- as.matrix(y)
  xs <- as.matrix(xs)
  if (!is.null(xe) && ncol(as.matrix(xe)) > 0L) {
    stop("Strict Proposition 4(i) has qe = 0; xe must be NULL.")
  }
  if (check) {
    .prop4ii_check_unit_rows(y, "y")
    .prop4ii_check_unit_rows(xs, "xs")
  }

  prelim <- mobius_vMF_prop4i(
    y = y, xs = xs, xe = NULL, start = mean,
    det_constraint = det_constraint,
    algorithm = algorithm, opts = opts, check = FALSE
  )
  mean <- prelim$mean
  k <- prelim$k
  p1 <- .prop4i_reference_direction(mean)

  if (G01behaviour == "fixed" && is.null(G0)) {
    stop("At least the first column of G0 must be supplied when G01behaviour = 'fixed'.")
  }
  G01 <- switch(
    G01behaviour,
    p1 = p1,
    free = if (is.null(G0)) p1 else G0[, 1],
    fixed = G0[, 1]
  )

  rresid <- rotated_resid(y, prelim$pred, base = G01)
  if (!is.null(G0) && all(!is.na(G0))) {
    if (G01behaviour == "p1") {
      G0 <- cbind(G01, -jupp_Rmat(G0[, 1], G01) %*%
                        G0[, -1, drop = FALSE])
    }
    aremaining <- SvMF_prelim_scales(rresid, G0)
  } else {
    G0 <- SvMF_moment_axes(rresid, G01)
    aremaining <- SvMF_prelim_scales(rresid, G0)
  }

  list(
    mean = mean,
    Rtilde0 = mean$Rtilde0,
    psi = mean$psi,
    beta_s = mean$beta_s,
    k = k,
    a = c(1, aremaining),
    G0 = G0,
    nlopt = prelim$nlopt,
    pred = prelim$pred
  )
}

.prop4i_joint_loglik <- function(y, xs, mean, k, a, G0,
                                 approximate = TRUE) {
  if (!is.finite(k) || k <= 0 || any(!is.finite(a)) || any(a <= 0)) {
    return(-Inf)
  }
  pred <- try(mobius_link(xs = xs, xe = NULL, param = mean, check = FALSE),
              silent = TRUE)
  if (inherits(pred, "try-error") || any(!is.finite(pred))) return(-Inf)
  yrot <- try(undo_partransport(y, pred, G01 = G0[, 1]), silent = TRUE)
  if (inherits(yrot, "try-error") || any(!is.finite(yrot))) return(-Inf)

  ld <- NULL
  if (approximate && exists("ldSvMF_cann", mode = "function")) {
    z <- try(ldSvMF_cann(yrot, k = k, a = a, G = G0), silent = TRUE)
    if (!inherits(z, "try-error") && all(is.finite(z))) ld <- z
  }
  if (is.null(ld)) {
    z <- try(SvMF_log_lik_cann(yrot, SvMF_cann(k = k, a = a, G = G0)),
             silent = TRUE)
    if (!inherits(z, "try-error") && all(is.finite(z))) ld <- z
  }
  if (is.null(ld)) return(-Inf)
  sum(ld)
}

.prop4i_build_joint_spec <- function(Rref, prelim, G0reference,
                                     G01behaviour) {
  p <- nrow(Rref)
  G01behaviour <- match.arg(G01behaviour, c("p1", "free", "fixed"))
  dR <- .prop4ii_skew_dim(p)
  dPsi <- p
  dA <- max(p - 2L, 0L)
  dG <- if (G01behaviour == "free") {
    .prop4ii_skew_dim(p)
  } else {
    .prop4ii_skew_dim(p - 1L)
  }

  pos <- 0L
  take <- function(n) {
    if (n == 0L) return(integer(0))
    out <- pos + seq_len(n)
    pos <<- pos + n
    out
  }
  idxR <- take(dR)
  idxPsi <- take(dPsi)
  idxK <- take(1L)
  idxA <- take(dA)
  idxG <- take(dG)

  start_mean <- as_mobius_link_prop4i(prelim$mean, check = FALSE)
  Gstart <- .prop4ii_nearest_orthogonal(prelim$G0, det_sign = 1)
  Gcoord <- if (is.null(G0reference)) {
    Gstart
  } else {
    .prop4ii_nearest_orthogonal(G0reference, det_sign = 1)
  }

  mean_ref <- mobius_link_prop4i(Rref, psi = start_mean$psi, check = FALSE)
  target_base <- if (G01behaviour == "p1") {
    .prop4i_reference_direction(mean_ref)
  } else {
    Gstart[, 1]
  }

  if (G01behaviour == "free") {
    Gref <- Gcoord
    Q0 <- Gstart %*% t(Gref)
    thetaG0 <- .prop4ii_inverse_cayley(Q0)
    Ganchor <- NULL
    fixed_base <- NULL
  } else {
    Ganchor <- .prop4ii_rebase_G0(Gcoord, target_base)
    if (is.null(Ganchor)) Ganchor <- .prop4ii_rebase_G0(Gstart, target_base)
    if (is.null(Ganchor)) stop("Could not construct the G0 reference frame.")

    Gstart2 <- .prop4ii_rebase_G0(Gstart, target_base)
    if (is.null(Gstart2)) stop("Could not rebase the starting G0 frame.")
    H0 <- crossprod(Ganchor[, -1, drop = FALSE],
                    Gstart2[, -1, drop = FALSE])
    H0 <- .prop4ii_nearest_orthogonal(H0, det_sign = 1)
    thetaG0 <- .prop4ii_inverse_cayley(H0)
    Gref <- NULL
    fixed_base <- if (G01behaviour == "fixed") target_base else NULL
  }

  par0 <- numeric(pos)
  par0[idxR] <- 0
  par0[idxPsi] <- .prop4i_psi_to_u(start_mean$psi)
  par0[idxK] <- log(max(prelim$k, 1e-8))
  if (dA > 0L) par0[idxA] <- .prop4ii_scales_to_eta(prelim$a)
  if (dG > 0L) par0[idxG] <- thetaG0

  lower <- rep(-Inf, pos)
  upper <- rep( Inf, pos)
  lower[idxK] <- log(1e-8)
  upper[idxK] <- log(1e5)
  if (dA > 0L) {
    lower[idxA] <- -10
    upper[idxA] <- 10
  }

  list(
    p = p, Rref = Rref, a1 = prelim$a[1],
    G01behaviour = G01behaviour,
    Gref = Gref, Ganchor = Ganchor, fixed_base = fixed_base,
    idxR = idxR, idxPsi = idxPsi, idxK = idxK,
    idxA = idxA, idxG = idxG,
    par0 = par0, lower = lower, upper = upper
  )
}

.prop4i_unpack_joint <- function(theta, spec) {
  p <- spec$p
  R <- .prop4ii_cayley(theta[spec$idxR], p) %*% spec$Rref
  psi <- .prop4i_u_to_psi(theta[spec$idxPsi])
  mean <- mobius_link_prop4i(R, psi = psi, check = FALSE)
  k <- exp(theta[spec$idxK])
  eta <- if (length(spec$idxA)) theta[spec$idxA] else numeric(0)
  a <- .prop4ii_eta_to_scales(eta, spec$a1, p)

  if (spec$G01behaviour == "free") {
    G0 <- .prop4ii_cayley(theta[spec$idxG], p) %*% spec$Gref
  } else {
    Gbase <- spec$Ganchor
    if (spec$G01behaviour == "p1") {
      target <- .prop4i_reference_direction(mean)
      Gbase <- .prop4ii_rebase_G0(spec$Ganchor, target)
      if (is.null(Gbase)) return(NULL)
    }
    H <- if (p == 2L) {
      matrix(1, 1, 1)
    } else {
      .prop4ii_cayley(theta[spec$idxG], p - 1L)
    }
    B <- diag(p)
    B[-1, -1] <- H
    G0 <- Gbase %*% B
  }

  list(mean = mean, Rtilde0 = R, psi = psi, k = k, a = a, G0 = G0)
}

.prop4i_fit_component <- function(y, xs, Rref, prelim, G0reference,
                                  G01behaviour, algorithm, opts,
                                  approximate = TRUE) {
  spec <- .prop4i_build_joint_spec(
    Rref = Rref, prelim = prelim,
    G0reference = G0reference,
    G01behaviour = G01behaviour
  )

  eval_f <- function(theta) {
    pars <- .prop4i_unpack_joint(theta, spec)
    if (is.null(pars)) return(1e100)
    ll <- .prop4i_joint_loglik(
      y, xs, pars$mean, pars$k, pars$a, pars$G0,
      approximate = approximate
    )
    if (!is.finite(ll)) return(1e100)
    -ll / nrow(y)
  }

  default_opts <- list(
    algorithm = algorithm,
    xtol_rel = 1e-9,
    maxeval = 10000,
    print_level = 0
  )
  combined_opts <- utils::modifyList(default_opts, opts)
  uses_gradient <- grepl("^NLOPT_LD_", combined_opts$algorithm)
  eval_nloptr <- if (uses_gradient) {
    function(theta) {
      list(
        objective = eval_f(theta),
        gradient = .prop4ii_fd_gradient(
          eval_f, theta, spec$lower, spec$upper
        )
      )
    }
  } else {
    eval_f
  }

  fit <- try(
    nloptr::nloptr(
      x0 = spec$par0,
      eval_f = eval_nloptr,
      lb = spec$lower,
      ub = spec$upper,
      opts = combined_opts
    ),
    silent = TRUE
  )

  # Keep the valid preliminary point if the numerical optimiser fails.
  theta_best <- spec$par0
  f_best <- eval_f(theta_best)
  if (!inherits(fit, "try-error") && is.finite(fit$objective) &&
      fit$objective < f_best) {
    theta_best <- fit$solution
    f_best <- fit$objective
  }
  pars <- .prop4i_unpack_joint(theta_best, spec)
  ll <- if (is.null(pars)) -Inf else .prop4i_joint_loglik(
    y, xs, pars$mean, pars$k, pars$a, pars$G0,
    approximate = approximate
  )

  if (inherits(fit, "try-error")) {
    fit <- list(status = -999L, message = as.character(fit),
                objective = f_best, solution = theta_best)
  } else {
    fit$objective <- f_best
    fit$solution <- theta_best
  }
  list(nlopt = fit, pars = pars, lLik_approx = ll, spec = spec)
}

mobius_SvMF_joint_fit_prop4i <- function(
    y, xs, xe = NULL, mean, k, a, G0,
    G0reference = NULL, G01behaviour = "p1",
    det_constraint = c("orthogonal", "rotation"),
    algorithm = "NLOPT_LD_SLSQP", opts = list(), check = TRUE) {

  det_constraint <- match.arg(det_constraint)
  G01behaviour <- match.arg(G01behaviour, c("p1", "free", "fixed"))
  y <- as.matrix(y)
  xs <- as.matrix(xs)
  storage.mode(y) <- storage.mode(xs) <- "double"

  if (!is.null(xe) && ncol(as.matrix(xe)) > 0L) {
    stop("Strict Proposition 4(i) has qe = 0; xe must be NULL.")
  }
  if (nrow(y) != nrow(xs) || ncol(y) != ncol(xs)) {
    stop("For Proposition 4(i), y and xs must have identical dimensions.")
  }
  p <- ncol(y)
  if (check) {
    .prop4ii_check_unit_rows(y, "y")
    .prop4ii_check_unit_rows(xs, "xs")
  }

  start_mean <- as_mobius_link_prop4i(mean, check = FALSE)
  start_R <- start_mean$Rtilde0
  start_a <- .prop4ii_normalise_scales(a, p, a1 = a[1])
  start_G0 <- .prop4ii_nearest_orthogonal(G0, det_sign = 1)
  initial <- list(
    mean = start_mean, Rtilde0 = start_R, psi = start_mean$psi,
    k = k, a = start_a, G0 = start_G0
  )

  target_signs <- if (det_constraint == "rotation") 1 else c(1, -1)
  fits <- vector("list", length(target_signs))
  for (j in seq_along(target_signs)) {
    sgn <- target_signs[j]
    Rref <- if ((if (det(start_R) >= 0) 1 else -1) == sgn) {
      start_R
    } else {
      .prop4ii_procrustes(y, xs, det_sign = sgn)$R
    }
    pre_j <- initial
    pre_j$mean <- mobius_link_prop4i(
      Rref, psi = start_mean$psi, check = FALSE
    )
    pre_j$Rtilde0 <- Rref
    pre_j$psi <- start_mean$psi
    if (G01behaviour == "p1") {
      newbase <- .prop4i_reference_direction(pre_j$mean)
      rebased <- .prop4ii_rebase_G0(pre_j$G0, newbase)
      if (is.null(rebased)) stop("Could not rebase G0 for Proposition 4(i).")
      pre_j$G0 <- rebased
    }
    fits[[j]] <- .prop4i_fit_component(
      y = y, xs = xs, Rref = Rref, prelim = pre_j,
      G0reference = if (is.null(G0reference)) pre_j$G0 else G0reference,
      G01behaviour = G01behaviour,
      algorithm = algorithm, opts = opts,
      approximate = TRUE
    )
    fits[[j]]$det_sign <- sgn
  }

  ll_approx <- vapply(fits, function(z) z$lLik_approx, numeric(1))
  if (all(!is.finite(ll_approx))) {
    stop("All Proposition 4(i) optimisations failed.")
  }
  best <- fits[[which.max(ll_approx)]]
  pars <- best$pars

  # Same residual-axis identifiability convention as the existing models.
  aord <- order(pars$a[-1], decreasing = TRUE)
  pars$a[-1] <- pars$a[-1][aord]
  pars$G0[, -1] <- pars$G0[, -1, drop = FALSE][, aord, drop = FALSE]
  if (exists("standardise_col_signs", mode = "function")) {
    pars$G0[, -1] <- standardise_col_signs(pars$G0[, -1, drop = FALSE])
  }
  if (det(pars$G0) < 0) pars$G0[, p] <- -pars$G0[, p]

  pred <- mobius_link(xs = xs, xe = NULL, param = pars$mean, check = FALSE)
  if (p != 3L) {
    kres <- mobius_SvMF_konly(y = y, ymean = pred,
                              a = pars$a, G0 = pars$G0)
    pars$k <- kres$k
    lLik <- kres$lLik
  } else {
    yrot <- undo_partransport(y, pred, G01 = pars$G0[, 1])
    lLik <- sum(SvMF_log_lik_cann(
      yrot, SvMF_cann(k = pars$k, a = pars$a, G = pars$G0)
    ))
  }
  SvMF_cann_check(SvMF_cann(k = pars$k, a = pars$a, G = pars$G0))

  dots <- pmax(-1, pmin(1, rowSums(y * pred)))
  rrI <- rotated_resid(y, pred, north_pole(p))
  rresids_I <- rrI[, -1, drop = FALSE]
  attr(rresids_I, "samehemisphere") <- attr(rrI, "samehemisphere")
  colnames(rresids_I) <- paste0("r", seq_len(ncol(rresids_I)))
  rresids_G0 <- resid_SvMF_partransport(y, pred, G0 = pars$G0, scale = FALSE)
  rresids_std <- resid_SvMF_partransport(
    y, pred, pars$k, pars$a, pars$G0, scale = TRUE
  )

  dmean <- .prop4ii_skew_dim(p) + p
  dscale <- max(p - 2L, 0L)
  dG0 <- if (G01behaviour == "free") {
    .prop4ii_skew_dim(p)
  } else {
    .prop4ii_skew_dim(p - 1L)
  }
  DoF <- dmean + 1L + dscale + dG0
  mean <- pars$mean

  out <- list(
    mean = mean, est = mean,
    Rtilde0 = mean$Rtilde0,
    psi = mean$psi, phi = mean$phi,
    beta_s = mean$beta_s, rs1 = mean$rs1,
    k = pars$k, a = pars$a, G0 = pars$G0,
    obj = -lLik / nrow(y), nlopt = best$nlopt,
    y = y, xs = xs, xe = NULL, pred = pred,
    rresids_I = rresids_I,
    rresids_G0 = rresids_G0,
    rresids_std = rresids_std,
    dists = acos(dots),
    DoF = DoF, AIC = 2 * DoF - 2 * lLik, lLik = lLik,
    initial = initial,
    component_fits = lapply(fits, function(z) {
      list(det_sign = z$det_sign, status = z$nlopt$status,
           message = z$nlopt$message, lLik_approx = z$lLik_approx)
    }),
    linktype = list(
      type = "Prop4i", det_constraint = det_constraint,
      G01behaviour = G01behaviour,
      intercept = FALSE, euclidean = FALSE
    )
  )
  class(out) <- c("mobius_SvMF_prop4i", "mobius_SvMF", "list")
  out
}

mobius_SvMF_prop4i <- function(
    y, xs, xe = NULL, mean = NULL, k = NULL, a = NULL,
    G0 = NULL, G0reference = NULL, G01behaviour = "p1",
    doprelim = TRUE,
    det_constraint = c("orthogonal", "rotation"),
    algorithm = "NLOPT_LD_SLSQP", opts = list(), check = TRUE, ...) {

  if (!is.null(xe) && ncol(as.matrix(xe)) > 0L) {
    stop("Strict Proposition 4(i) has qe = 0; xe must be NULL.")
  }
  det_constraint <- match.arg(det_constraint)
  G01behaviour <- match.arg(G01behaviour, c("p1", "free", "fixed"))
  dots <- list(...)
  if (length(dots)) opts <- utils::modifyList(opts, dots)

  if (doprelim) {
    preest <- mobius_SvMF_partransport_prelim_prop4i(
      y = y, xs = xs, xe = NULL, mean = mean, G0 = G0,
      G01behaviour = G01behaviour,
      det_constraint = det_constraint,
      algorithm = algorithm, opts = opts, check = check
    )
  } else {
    if (is.null(mean) || is.null(k) || is.null(a) || is.null(G0)) {
      stop("When doprelim = FALSE, mean, k, a, and G0 must all be supplied.")
    }
    mean <- as_mobius_link_prop4i(mean)
    preest <- list(
      mean = mean, Rtilde0 = mean$Rtilde0, psi = mean$psi,
      beta_s = mean$beta_s, k = k, a = a, G0 = G0
    )
  }

  finalest <- mobius_SvMF_joint_fit_prop4i(
    y = y, xs = xs, xe = NULL,
    mean = preest$mean,
    k = if (is.null(k)) preest$k else k,
    a = if (is.null(a)) preest$a else a,
    G0 = preest$G0,
    G0reference = if (is.null(G0reference)) preest$G0 else G0reference,
    G01behaviour = G01behaviour,
    det_constraint = det_constraint,
    algorithm = algorithm, opts = opts, check = check
  )
  finalest$preest <- preest
  finalest
}

print.mobius_vMF_prop4i <- function(x, ...) {
  cat("vMF regression with Proposition 4(i) spherical Möbius link\n")
  cat("  beta_s =", format(x$beta_s, digits = 6), "\n")
  cat("  phi =", format(x$phi, digits = 6), "\n")
  cat("  det(Rtilde0) =", format(det(x$Rtilde0), digits = 6), "\n")
  cat("  k =", format(x$k, digits = 6), "\n")
  cat("  log-likelihood =", format(x$lLik, digits = 7), "\n")
  cat("  AIC =", format(x$AIC, digits = 7), "\n")
  invisible(x)
}

print.mobius_SvMF_prop4i <- function(x, ...) {
  cat("SvMF regression with Proposition 4(i) spherical Möbius link\n")
  cat("  G0 estimation: parallel-transport/Cayley procedure\n")
  cat("  G01behaviour =", x$linktype$G01behaviour, "\n")
  cat("  beta_s =", format(x$beta_s, digits = 6), "\n")
  cat("  phi =", format(x$phi, digits = 6), "\n")
  cat("  det(Rtilde0) =", format(det(x$Rtilde0), digits = 6), "\n")
  cat("  k =", format(x$k, digits = 6), "\n")
  cat("  a =", paste(format(x$a, digits = 5), collapse = ", "), "\n")
  cat("  log-likelihood =", format(x$lLik, digits = 7), "\n")
  cat("  AIC =", format(x$AIC, digits = 7), "\n")
  invisible(x)
}

predict.mobius_vMF_prop4i <- function(object, newdata = NULL,
                                       xs = newdata, ...) {
  if (is.null(xs)) xs <- object$xs
  mobius_link(xs = as.matrix(xs), xe = NULL, param = object$mean)
}

predict.mobius_SvMF_prop4i <- function(object, newdata = NULL,
                                        xs = newdata, ...) {
  predict.mobius_vMF_prop4i(object, newdata = newdata, xs = xs, ...)
}

# Final user-facing dispatchers.  Existing Proposition 4(ii) and Definition 1
# behaviour is delegated unchanged to the saved working wrappers.
mobius_vMF <- function(
    y, xs = NULL, xe = NULL, start = NULL,
    type = "SpEuc", fix_qs1 = FALSE,
    fix_qe1 = (type == "LinEuc"), intercept = TRUE,
    lb = NULL, ub = NULL,
    det_constraint = c("orthogonal", "rotation"),
    prop4i_algorithm = "NLOPT_LD_SLSQP",
    prop4i_opts = list(),
    prop4ii_algorithm = "NLOPT_LD_SLSQP",
    prop4ii_opts = list(), ...) {

  if (.is_prop4i_type(type)) {
    dots <- list(...)
    if (length(dots)) prop4i_opts <- utils::modifyList(prop4i_opts, dots)
    return(mobius_vMF_prop4i(
      y = y, xs = xs, xe = xe, start = start,
      det_constraint = det_constraint,
      algorithm = prop4i_algorithm, opts = prop4i_opts
    ))
  }

  .mobius_vMF_before_prop4i(
    y = y, xs = xs, xe = xe, start = start,
    type = type, fix_qs1 = fix_qs1, fix_qe1 = fix_qe1,
    intercept = intercept, lb = lb, ub = ub,
    det_constraint = det_constraint,
    prop4ii_algorithm = prop4ii_algorithm,
    prop4ii_opts = prop4ii_opts, ...
  )
}

mobius_SvMF <- function(
    y, xs = NULL, xe = NULL, mean = NULL, k = NULL, a = NULL,
    G0 = NULL, G0reference = NULL, G01behaviour = "p1",
    type = "LinEuc", fix_qs1 = FALSE,
    fix_qe1 = (type == "LinEuc"), intercept = TRUE,
    doprelim = TRUE,
    det_constraint = c("orthogonal", "rotation"),
    prop4i_algorithm = "NLOPT_LD_SLSQP",
    prop4i_opts = list(),
    prop4ii_algorithm = "NLOPT_LD_SLSQP",
    prop4ii_opts = list(), ...) {

  if (.is_prop4i_type(type)) {
    return(mobius_SvMF_prop4i(
      y = y, xs = xs, xe = xe, mean = mean, k = k, a = a,
      G0 = G0, G0reference = G0reference,
      G01behaviour = G01behaviour, doprelim = doprelim,
      det_constraint = det_constraint,
      algorithm = prop4i_algorithm,
      opts = prop4i_opts,
      ...
    ))
  }

  .mobius_SvMF_before_prop4i(
    y = y, xs = xs, xe = xe, mean = mean, k = k, a = a,
    G0 = G0, G0reference = G0reference,
    G01behaviour = G01behaviour,
    type = type, fix_qs1 = fix_qs1, fix_qe1 = fix_qe1,
    intercept = intercept, doprelim = doprelim,
    det_constraint = det_constraint,
    prop4ii_algorithm = prop4ii_algorithm,
    prop4ii_opts = prop4ii_opts,
    ...
  )
}

# =============================================================================
# Proposition 4(i)-constrained spherical component + optional Euclidean terms
#
# IMPORTANT:
# Proposition 4(i) in the paper itself is stated for qe = 0.  The extension
# below keeps the Proposition 4(i) restriction on the SPHERICAL component,
#
#     p = qs,   Bs = beta_s I_(p-1),   0 < beta_s <= 1,
#
# while adding the same linear Euclidean term used by the existing LinEuc
# implementation:
#
#   mu(xs,xe) = P S^{-1}{
#                 beta_s S(Qs^T xs) + Be Qe_star^T x_e,model
#               },
#
# with Rtilde0 = P Qs^T.  Hence beta_s = 1 reduces exactly to the existing
# Proposition 4(ii) + Euclidean model.
#
# User-facing usage:
#   type = "Prop4i", xe = NULL     : strict Proposition 4(i)
#   type = "Prop4i", xe != NULL    : Proposition 4(i)-constrained spherical
#                                    component + LinEuc covariates
# =============================================================================

# Save the working dispatchers before adding the optional-Euclidean Prop. 4(i)
# extension.  The spherical-only Prop. 4(i), Prop. 4(ii), and Definition 1
# behaviours remain unchanged.
.mobius_link_before_prop4i_euc <- mobius_link
.mobius_vMF_before_prop4i_euc <- mobius_vMF
.mobius_SvMF_before_prop4i_euc <- mobius_SvMF


# -----------------------------------------------------------------------------
# Parameter object and mean link
# -----------------------------------------------------------------------------

mobius_link_prop4i_euc <- function(
    Rtilde0, P, beta_s, Qe_star, Be,
    xe_center, xe_scale, intercept = TRUE, check = TRUE) {

  Rtilde0 <- as.matrix(Rtilde0)
  P <- as.matrix(P)
  Qe_star <- as.matrix(Qe_star)
  storage.mode(Rtilde0) <- "double"
  storage.mode(P) <- "double"
  storage.mode(Qe_star) <- "double"

  beta_s <- as.numeric(beta_s)[1L]
  Be <- as.numeric(Be)
  xe_center <- as.numeric(xe_center)
  xe_scale <- as.numeric(xe_scale)

  # Rtilde0 = P Qs^T, hence Qs = Rtilde0^T P.
  Qs <- t(Rtilde0) %*% P
  rs1 <- Qs[, 1L]
  phi <- (1 - beta_s) / (1 + beta_s)
  psi <- phi * rs1

  obj <- list(
    Rtilde0 = Rtilde0,
    P = P,
    beta_s = beta_s,
    phi = phi,
    psi = psi,
    rs1 = rs1,
    Qe_star = Qe_star,
    Be = Be,
    xe_center = xe_center,
    xe_scale = xe_scale,
    intercept = isTRUE(intercept)
  )
  class(obj) <- c("mobius_link_prop4i_euc", "list")
  if (check) mobius_link_prop4i_euc_check(obj)
  obj
}


mobius_link_prop4i_euc_check <- function(
    obj, tol = 100 * sqrt(.Machine$double.eps)) {

  if (!inherits(obj, "mobius_link_prop4i_euc")) {
    stop("obj must inherit from 'mobius_link_prop4i_euc'.")
  }

  R <- obj$Rtilde0
  P <- obj$P
  if (!is.matrix(R) || nrow(R) != ncol(R) || nrow(R) < 2L) {
    stop("Rtilde0 must be a square p by p matrix with p >= 2.")
  }
  p <- nrow(R)
  k <- p - 1L
  if (any(!is.finite(R)) ||
      max(abs(crossprod(R) - diag(p))) > tol) {
    stop("Rtilde0 must be orthogonal.")
  }

  if (!all(dim(P) == c(p, p)) || any(!is.finite(P)) ||
      max(abs(crossprod(P) - diag(p))) > tol || det(P) < 0) {
    stop("P must be a proper orthogonal p by p matrix.")
  }

  if (!is.finite(obj$beta_s) || obj$beta_s <= 0 || obj$beta_s > 1) {
    stop("beta_s must satisfy 0 < beta_s <= 1.")
  }

  V <- obj$Qe_star
  if (!is.matrix(V) || ncol(V) != k || nrow(V) < k ||
      any(!is.finite(V))) {
    stop("Qe_star must be m by (p-1), with m >= p-1.")
  }
  if (max(abs(crossprod(V) - diag(k))) > tol) {
    stop("Qe_star must have orthonormal columns.")
  }

  b <- obj$Be
  if (length(b) != k || any(!is.finite(b)) ||
      any(b <= 0) || any(b >= 1)) {
    stop("Be must contain p-1 values strictly between zero and one.")
  }

  q <- nrow(V) - as.integer(isTRUE(obj$intercept))
  if (q < 0L || length(obj$xe_center) != q ||
      length(obj$xe_scale) != q) {
    stop("Euclidean centering/scaling information has the wrong dimension.")
  }
  if (any(!is.finite(obj$xe_center)) ||
      any(!is.finite(obj$xe_scale)) ||
      any(obj$xe_scale <= 0)) {
    stop("Euclidean centering/scaling information is invalid.")
  }

  invisible(NULL)
}


as_mobius_link_prop4i_euc <- function(obj, check = TRUE) {
  if (inherits(obj, "mobius_link_prop4i_euc")) {
    if (check) mobius_link_prop4i_euc_check(obj)
    return(obj)
  }
  if (is.list(obj) &&
      all(c("Rtilde0", "P", "beta_s", "Qe_star", "Be",
            "xe_center", "xe_scale") %in% names(obj))) {
    return(mobius_link_prop4i_euc(
      Rtilde0 = obj$Rtilde0,
      P = obj$P,
      beta_s = obj$beta_s,
      Qe_star = obj$Qe_star,
      Be = obj$Be,
      xe_center = obj$xe_center,
      xe_scale = obj$xe_scale,
      intercept = if (is.null(obj$intercept)) TRUE else obj$intercept,
      check = check
    ))
  }
  stop("obj is not a compatible Proposition 4(i)+Euclidean parameter object.")
}


print.mobius_link_prop4i_euc <- function(x, ...) {
  cat("Proposition 4(i)-constrained spherical link with LinEuc covariates\n")
  cat("  p = qs =", nrow(x$Rtilde0), "\n")
  cat("  qe =", length(x$xe_center), "\n")
  cat("  beta_s =", format(x$beta_s, digits = 6), "\n")
  cat("  phi =", format(x$phi, digits = 6), "\n")
  cat("  det(Rtilde0) =", format(det(x$Rtilde0), digits = 6), "\n")
  cat("  Be =", paste(format(x$Be, digits = 5), collapse = ", "), "\n")
  invisible(x)
}


dim.mobius_link_prop4i_euc <- function(x) {
  p <- nrow(x$Rtilde0)
  c(p = p, qs = p, qe = length(x$xe_center))
}


# Mean-link dispatcher.  beta_s = 1 gives exactly the existing Prop. 4(ii)
# optional-Euclidean link.
mobius_link <- function(xs = NULL, xe = NULL, param = NULL, check = TRUE) {
  if (!inherits(param, "mobius_link_prop4i_euc")) {
    return(.mobius_link_before_prop4i_euc(
      xs = xs, xe = xe, param = param, check = check
    ))
  }

  if (check) mobius_link_prop4i_euc_check(param)
  if (is.null(xs)) stop("The Proposition 4(i)+Euclidean link requires xs.")
  if (is.null(xe)) stop("xe must be supplied for this fitted model.")

  xs <- as.matrix(xs)
  if (ncol(xs) != nrow(param$Rtilde0)) {
    stop("ncol(xs) must equal nrow(Rtilde0).")
  }

  ep <- .prop4ii_prepare_euc(
    xe,
    intercept = param$intercept,
    center = param$xe_center,
    scale = param$xe_scale,
    fitting = FALSE
  )
  if (nrow(xs) != nrow(ep$model)) {
    stop("xs and xe must have the same number of rows.")
  }

  Qs <- t(param$Rtilde0) %*% param$P
  spherical_part <- param$beta_s * Sp(xs %*% Qs)
  euclidean_part <- sweep(
    ep$model %*% param$Qe_star, 2L, param$Be, "*"
  )

  iSp(spherical_part + euclidean_part) %*% t(param$P)
}


# -----------------------------------------------------------------------------
# Mean-parameter optimisation
# -----------------------------------------------------------------------------

.prop4i_euc_mean_spec <- function(
    Rref, Pref, Vref, bstart, beta_start,
    xe_center, xe_scale, intercept) {

  base <- .prop4ii_mean_spec_euc(
    Rref = Rref,
    Pref = Pref,
    Vref = Vref,
    bstart = bstart,
    xe_center = xe_center,
    xe_scale = xe_scale,
    intercept = intercept
  )

  beta_start <- min(max(as.numeric(beta_start)[1L], 1e-5), 1 - 1e-5)
  idxBeta <- length(base$par0) + 1L

  list(
    base = base,
    base_length = length(base$par0),
    idxBeta = idxBeta,
    par0 = c(base$par0, qlogis(beta_start)),
    lower = c(base$lower, -12),
    upper = c(base$upper,  12)
  )
}


.prop4i_euc_unpack_mean <- function(theta, spec) {
  mp <- .prop4ii_unpack_mean_euc(
    theta[seq_len(spec$base_length)], spec$base
  )
  beta_s <- plogis(theta[spec$idxBeta])

  mean <- mobius_link_prop4i_euc(
    Rtilde0 = mp$Rtilde0,
    P = mp$P,
    beta_s = beta_s,
    Qe_star = mp$Qe_star,
    Be = mp$Be,
    xe_center = spec$base$xe_center,
    xe_scale = spec$base$xe_scale,
    intercept = spec$base$intercept,
    check = FALSE
  )

  list(
    mean = mean,
    Rtilde0 = mp$Rtilde0,
    P = mp$P,
    beta_s = beta_s,
    Qe_star = mp$Qe_star,
    Be = mp$Be
  )
}


.prop4i_euc_fit_vmf_component <- function(
    y, xs, xe, spec,
    algorithm = "NLOPT_LD_SLSQP", opts = list()) {

  eval_f <- function(theta) {
    pars <- .prop4i_euc_unpack_mean(theta, spec)
    pred <- try(
      mobius_link(xs = xs, xe = xe, param = pars$mean, check = FALSE),
      silent = TRUE
    )
    if (inherits(pred, "try-error") || any(!is.finite(pred))) {
      return(1e100)
    }
    -mean(rowSums(y * pred))
  }

  default_opts <- list(
    algorithm = algorithm,
    xtol_rel = 1e-9,
    maxeval = 10000,
    print_level = 0
  )
  combined_opts <- utils::modifyList(default_opts, opts)
  uses_gradient <- grepl("^NLOPT_LD_", combined_opts$algorithm)

  eval_nloptr <- if (uses_gradient) {
    function(theta) {
      list(
        objective = eval_f(theta),
        gradient = .prop4ii_fd_gradient(
          eval_f, theta, spec$lower, spec$upper
        )
      )
    }
  } else {
    eval_f
  }

  fit <- try(
    nloptr::nloptr(
      x0 = spec$par0,
      eval_f = eval_nloptr,
      lb = spec$lower,
      ub = spec$upper,
      opts = combined_opts
    ),
    silent = TRUE
  )

  theta_best <- spec$par0
  f_best <- eval_f(theta_best)
  if (!inherits(fit, "try-error") &&
      is.finite(fit$objective) && fit$objective < f_best) {
    theta_best <- fit$solution
    f_best <- fit$objective
  }

  pars <- .prop4i_euc_unpack_mean(theta_best, spec)
  pred <- mobius_link(
    xs = xs, xe = xe, param = pars$mean, check = FALSE
  )

  if (inherits(fit, "try-error")) {
    fit <- list(
      status = -999L,
      message = as.character(fit),
      objective = f_best,
      solution = theta_best
    )
  } else {
    fit$objective <- f_best
    fit$solution <- theta_best
  }

  list(
    nlopt = fit,
    pars = pars,
    pred = pred,
    objective = f_best
  )
}


# Construct a proper P starting value from a strict Proposition 4(i) mean.
# This is used only for numerical initialisation; the fitted P is free.
.prop4i_euc_P_from_strict <- function(mean) {
  mean <- as_mobius_link_prop4i(mean, check = FALSE)
  R <- mean$Rtilde0
  p <- nrow(R)

  rs1 <- mean$rs1
  if (length(rs1) != p || any(!is.finite(rs1))) {
    rs1 <- c(1, rep(0, p - 1L))
  }
  rs1 <- rs1 / sqrt(sum(rs1^2))

  # Obtain an orthonormal tangent basis perpendicular to rs1.
  vv <- svd(matrix(rs1, nrow = 1L), nu = 0, nv = p)$v
  tangent <- vv[, -1L, drop = FALSE]
  Qs <- cbind(rs1, tangent)

  # det(P)=det(Rtilde0)*det(Qs) must be +1.
  target_qs_sign <- if (det(R) >= 0) 1 else -1
  if ((if (det(Qs) >= 0) 1 else -1) != target_qs_sign) {
    Qs[, p] <- -Qs[, p]
  }

  P <- R %*% Qs
  if (det(P) < 0) {
    Qs[, p] <- -Qs[, p]
    P <- R %*% Qs
  }
  .prop4ii_nearest_orthogonal(P, det_sign = 1)
}


mobius_vMF_prop4i_euc <- function(
    y, xs, xe, start = NULL,
    det_constraint = c("orthogonal", "rotation"),
    intercept = TRUE,
    algorithm = "NLOPT_LD_SLSQP", opts = list(), check = TRUE) {

  det_constraint <- match.arg(det_constraint)
  y <- as.matrix(y)
  xs <- as.matrix(xs)
  xe <- as.matrix(xe)
  storage.mode(y) <- storage.mode(xs) <- storage.mode(xe) <- "double"

  if (nrow(y) != nrow(xs) || nrow(y) != nrow(xe) ||
      ncol(y) != ncol(xs)) {
    stop("For Proposition 4(i)+Euclidean, y and xs must be n by p and xe must have n rows.")
  }

  p <- ncol(y)
  kmean <- p - 1L
  if (check) {
    .prop4ii_check_unit_rows(y, "y")
    .prop4ii_check_unit_rows(xs, "xs")
  }

  ep <- .prop4ii_prepare_euc(xe, intercept = intercept, fitting = TRUE)
  if (ep$m < kmean) {
    stop("With Euclidean covariates, ncol(xe) + intercept must be at least p-1.")
  }

  if (!is.null(start) && is.list(start) && !is.null(start$mean)) {
    start <- start$mean
  }

  start_euc <- NULL
  start_strict <- NULL
  if (!is.null(start)) {
    if (inherits(start, "mobius_link_prop4i_euc")) {
      start_euc <- as_mobius_link_prop4i_euc(start, check = FALSE)
    } else if (inherits(start, "mobius_link_prop4i")) {
      start_strict <- as_mobius_link_prop4i(start, check = FALSE)
    }
  }

  start_R <- if (!is.null(start_euc)) {
    start_euc$Rtilde0
  } else if (!is.null(start_strict)) {
    start_strict$Rtilde0
  } else {
    NULL
  }

  start_P <- if (!is.null(start_euc)) {
    start_euc$P
  } else if (!is.null(start_strict)) {
    .prop4i_euc_P_from_strict(start_strict)
  } else {
    NULL
  }

  start_V <- if (!is.null(start_euc) &&
                 all(dim(start_euc$Qe_star) == c(ep$m, kmean))) {
    start_euc$Qe_star
  } else {
    diag(ep$m)[, seq_len(kmean), drop = FALSE]
  }

  start_Be <- if (!is.null(start_euc) &&
                  length(start_euc$Be) == kmean) {
    start_euc$Be
  } else {
    rep(0.05, kmean)
  }

  start_beta <- if (!is.null(start_euc)) {
    start_euc$beta_s
  } else if (!is.null(start_strict)) {
    start_strict$beta_s
  } else {
    0.8
  }

  target_signs <- if (det_constraint == "rotation") 1 else c(1, -1)
  fits <- vector("list", length(target_signs))

  for (j in seq_along(target_signs)) {
    sgn <- target_signs[j]
    Rref <- if (!is.null(start_R) &&
                (if (det(start_R) >= 0) 1 else -1) == sgn) {
      start_R
    } else {
      .prop4ii_procrustes(y, xs, det_sign = sgn)$R
    }

    Pref <- if (!is.null(start_P)) {
      .prop4ii_nearest_orthogonal(start_P, det_sign = 1)
    } else {
      .prop4ii_nearest_orthogonal(Rref, det_sign = 1)
    }

    spec <- .prop4i_euc_mean_spec(
      Rref = Rref,
      Pref = Pref,
      Vref = start_V,
      bstart = start_Be,
      beta_start = start_beta,
      xe_center = ep$center,
      xe_scale = ep$scale,
      intercept = intercept
    )

    fits[[j]] <- .prop4i_euc_fit_vmf_component(
      y = y, xs = xs, xe = xe, spec = spec,
      algorithm = algorithm, opts = opts
    )
    fits[[j]]$det_sign <- sgn
  }

  values <- vapply(fits, function(z) z$objective, numeric(1))
  best <- fits[[which.min(values)]]
  pars <- best$pars
  pred <- best$pred

  dots <- pmax(-1, pmin(1, rowSums(y * pred)))
  kfit <- .prop4ii_estimate_k(mean(dots), p)
  n <- nrow(y)

  dR <- .prop4ii_skew_dim(p)
  dP <- .prop4ii_skew_dim(p)
  dV <- .prop4ii_stiefel_dim(ep$m, kmean)
  dB <- kmean
  dBeta <- 1L
  dmean <- dR + dP + dV + dB + dBeta
  DoF <- dmean + 1L
  lLik <- n * kfit$loglik_per_obs

  rr <- try(rotated_resid(y, pred, north_pole(p)), silent = TRUE)
  rresids <- if (inherits(rr, "try-error")) NULL else rr[, -1, drop = FALSE]
  if (!is.null(rresids)) {
    attr(rresids, "samehemisphere") <- attr(rr, "samehemisphere")
    colnames(rresids) <- paste0("r", seq_len(ncol(rresids)))
  }

  mean_obj <- pars$mean
  out <- list(
    est = mean_obj, mean = mean_obj,
    Rtilde0 = mean_obj$Rtilde0,
    P = mean_obj$P,
    beta_s = mean_obj$beta_s,
    phi = mean_obj$phi,
    psi = mean_obj$psi,
    rs1 = mean_obj$rs1,
    Qe_star = mean_obj$Qe_star,
    Be = mean_obj$Be,
    k = kfit$k,
    obj = best$objective,
    solution = best$nlopt$solution,
    nlopt = best$nlopt,
    y = y, xs = xs, xe = xe,
    pred = pred,
    rresids = rresids,
    dists = acos(dots),
    DoF = DoF,
    AIC = 2 * DoF - 2 * lLik,
    lLik = lLik,
    start = start,
    component_fits = lapply(fits, function(z) {
      list(
        det_sign = z$det_sign,
        status = z$nlopt$status,
        message = z$nlopt$message,
        objective = z$objective
      )
    }),
    linktype = list(
      type = "Prop4i",
      det_constraint = det_constraint,
      intercept = intercept,
      euclidean = TRUE,
      spherical_constraint = "Bs = beta_s I"
    )
  )
  class(out) <- c(
    "mobius_vMF_prop4i_euc", "mobius_vMF_prop4i", "mobius_vMF", "list"
  )
  out
}


# -----------------------------------------------------------------------------
# SvMF fit
# -----------------------------------------------------------------------------

mobius_SvMF_partransport_prelim_prop4i_euc <- function(
    y, xs, xe, mean = NULL, G0 = NULL,
    G01behaviour = "p1",
    det_constraint = c("orthogonal", "rotation"),
    intercept = TRUE,
    algorithm = "NLOPT_LD_SLSQP", opts = list(),
    check = TRUE, ...) {

  det_constraint <- match.arg(det_constraint)
  G01behaviour <- match.arg(
    G01behaviour, c("p1", "free", "fixed")
  )

  prelim <- mobius_vMF_prop4i_euc(
    y = y, xs = xs, xe = xe, start = mean,
    det_constraint = det_constraint,
    intercept = intercept,
    algorithm = algorithm, opts = opts, check = check
  )

  p1 <- prelim$mean$P[, 1L]
  if (G01behaviour == "fixed" && is.null(G0)) {
    stop("At least the first column of G0 must be supplied when G01behaviour='fixed'.")
  }

  G01 <- switch(
    G01behaviour,
    p1 = p1,
    free = if (is.null(G0)) p1 else G0[, 1L],
    fixed = G0[, 1L]
  )

  rresid <- rotated_resid(as.matrix(y), prelim$pred, base = G01)
  if (!is.null(G0) && all(!is.na(G0))) {
    if (G01behaviour == "p1") {
      G0 <- cbind(
        G01,
        -jupp_Rmat(G0[, 1L], G01) %*%
          G0[, -1L, drop = FALSE]
      )
    }
    aremaining <- SvMF_prelim_scales(rresid, G0)
  } else {
    G0 <- SvMF_moment_axes(rresid, G01)
    aremaining <- SvMF_prelim_scales(rresid, G0)
  }

  list(
    mean = prelim$mean,
    Rtilde0 = prelim$Rtilde0,
    P = prelim$P,
    beta_s = prelim$beta_s,
    Qe_star = prelim$Qe_star,
    Be = prelim$Be,
    k = prelim$k,
    a = c(1, aremaining),
    G0 = G0,
    nlopt = prelim$nlopt,
    pred = prelim$pred
  )
}


.prop4i_euc_joint_spec <- function(
    mean_ref, prelim, G0reference, G01behaviour) {

  mean_ref <- as_mobius_link_prop4i_euc(mean_ref, check = FALSE)
  p <- nrow(mean_ref$Rtilde0)

  mspec <- .prop4i_euc_mean_spec(
    Rref = mean_ref$Rtilde0,
    Pref = mean_ref$P,
    Vref = mean_ref$Qe_star,
    bstart = mean_ref$Be,
    beta_start = mean_ref$beta_s,
    xe_center = mean_ref$xe_center,
    xe_scale = mean_ref$xe_scale,
    intercept = mean_ref$intercept
  )

  dA <- max(p - 2L, 0L)
  dG <- if (G01behaviour == "free") {
    .prop4ii_skew_dim(p)
  } else {
    .prop4ii_skew_dim(p - 1L)
  }

  pos <- length(mspec$par0)
  take <- function(n) {
    if (n == 0L) return(integer(0))
    z <- pos + seq_len(n)
    pos <<- pos + n
    z
  }

  idxK <- take(1L)
  idxA <- take(dA)
  idxG <- take(dG)

  Gstart <- .prop4ii_nearest_orthogonal(prelim$G0, det_sign = 1)
  Gcoord <- if (is.null(G0reference)) {
    Gstart
  } else {
    .prop4ii_nearest_orthogonal(G0reference, det_sign = 1)
  }

  if (G01behaviour == "free") {
    Gref <- Gcoord
    thetaG0 <- .prop4ii_inverse_cayley(Gstart %*% t(Gref))
    Ganchor <- NULL
  } else {
    target <- if (G01behaviour == "p1") {
      mean_ref$P[, 1L]
    } else {
      Gstart[, 1L]
    }

    Ganchor <- .prop4ii_rebase_G0(Gcoord, target)
    if (is.null(Ganchor)) {
      Ganchor <- .prop4ii_rebase_G0(Gstart, target)
    }
    if (is.null(Ganchor)) {
      stop("Could not construct the G0 reference frame.")
    }

    Gstart2 <- .prop4ii_rebase_G0(Gstart, target)
    H0 <- crossprod(
      Ganchor[, -1L, drop = FALSE],
      Gstart2[, -1L, drop = FALSE]
    )
    H0 <- .prop4ii_nearest_orthogonal(H0, det_sign = 1)
    thetaG0 <- .prop4ii_inverse_cayley(H0)
    Gref <- NULL
  }

  par0 <- c(mspec$par0, rep(0, pos - length(mspec$par0)))
  par0[idxK] <- log(max(prelim$k, 1e-8))
  if (dA > 0L) {
    par0[idxA] <- .prop4ii_scales_to_eta(prelim$a)
  }
  if (dG > 0L) {
    par0[idxG] <- thetaG0
  }

  lower <- c(mspec$lower, rep(-Inf, pos - length(mspec$lower)))
  upper <- c(mspec$upper, rep( Inf, pos - length(mspec$upper)))
  lower[idxK] <- log(1e-8)
  upper[idxK] <- log(1e5)
  if (dA > 0L) {
    lower[idxA] <- -10
    upper[idxA] <- 10
  }

  list(
    mean_spec = mspec,
    mean_length = length(mspec$par0),
    p = p,
    a1 = prelim$a[1L],
    G01behaviour = G01behaviour,
    Gref = Gref,
    Ganchor = Ganchor,
    idxK = idxK,
    idxA = idxA,
    idxG = idxG,
    par0 = par0,
    lower = lower,
    upper = upper
  )
}


.prop4i_euc_unpack_joint <- function(theta, spec) {
  mp <- .prop4i_euc_unpack_mean(
    theta[seq_len(spec$mean_length)], spec$mean_spec
  )

  kappa <- exp(theta[spec$idxK])
  eta <- if (length(spec$idxA)) theta[spec$idxA] else numeric(0)
  a <- .prop4ii_eta_to_scales(eta, spec$a1, spec$p)

  if (spec$G01behaviour == "free") {
    G0 <- .prop4ii_cayley(theta[spec$idxG], spec$p) %*% spec$Gref
  } else {
    Gbase <- spec$Ganchor
    if (spec$G01behaviour == "p1") {
      Gbase <- .prop4ii_rebase_G0(
        spec$Ganchor, mp$P[, 1L]
      )
      if (is.null(Gbase)) return(NULL)
    }

    H <- if (spec$p == 2L) {
      matrix(1, 1, 1)
    } else {
      .prop4ii_cayley(theta[spec$idxG], spec$p - 1L)
    }

    B <- diag(spec$p)
    B[-1L, -1L] <- H
    G0 <- Gbase %*% B
  }

  c(mp, list(k = kappa, a = a, G0 = G0))
}


.prop4i_euc_joint_loglik <- function(
    y, xs, xe, mean, k, a, G0, approximate = TRUE) {

  if (!is.finite(k) || k <= 0 ||
      any(!is.finite(a)) || any(a <= 0)) {
    return(-Inf)
  }

  pred <- try(
    mobius_link(xs = xs, xe = xe, param = mean, check = FALSE),
    silent = TRUE
  )
  if (inherits(pred, "try-error") || any(!is.finite(pred))) {
    return(-Inf)
  }

  yrot <- try(
    undo_partransport(y, pred, G01 = G0[, 1L]),
    silent = TRUE
  )
  if (inherits(yrot, "try-error") || any(!is.finite(yrot))) {
    return(-Inf)
  }

  ld <- NULL
  if (approximate && exists("ldSvMF_cann", mode = "function")) {
    z <- try(
      ldSvMF_cann(yrot, k = k, a = a, G = G0),
      silent = TRUE
    )
    if (!inherits(z, "try-error") && all(is.finite(z))) {
      ld <- z
    }
  }

  if (is.null(ld)) {
    z <- try(
      SvMF_log_lik_cann(
        yrot, SvMF_cann(k = k, a = a, G = G0)
      ),
      silent = TRUE
    )
    if (!inherits(z, "try-error") && all(is.finite(z))) {
      ld <- z
    }
  }

  if (is.null(ld)) return(-Inf)
  sum(ld)
}


.prop4i_euc_fit_joint_component <- function(
    y, xs, xe, mean_ref, prelim, G0reference,
    G01behaviour, algorithm, opts) {

  spec <- .prop4i_euc_joint_spec(
    mean_ref = mean_ref,
    prelim = prelim,
    G0reference = G0reference,
    G01behaviour = G01behaviour
  )

  eval_f <- function(theta) {
    pars <- .prop4i_euc_unpack_joint(theta, spec)
    if (is.null(pars)) return(1e100)

    ll <- .prop4i_euc_joint_loglik(
      y, xs, xe,
      pars$mean, pars$k, pars$a, pars$G0,
      approximate = TRUE
    )
    if (!is.finite(ll)) return(1e100)
    -ll / nrow(y)
  }

  default_opts <- list(
    algorithm = algorithm,
    xtol_rel = 1e-9,
    maxeval = 10000,
    print_level = 0
  )
  combined_opts <- utils::modifyList(default_opts, opts)
  uses_gradient <- grepl("^NLOPT_LD_", combined_opts$algorithm)

  eval_nloptr <- if (uses_gradient) {
    function(theta) {
      list(
        objective = eval_f(theta),
        gradient = .prop4ii_fd_gradient(
          eval_f, theta, spec$lower, spec$upper
        )
      )
    }
  } else {
    eval_f
  }

  fit <- try(
    nloptr::nloptr(
      x0 = spec$par0,
      eval_f = eval_nloptr,
      lb = spec$lower,
      ub = spec$upper,
      opts = combined_opts
    ),
    silent = TRUE
  )

  # Preserve the valid preliminary point if the optimiser fails.
  theta_best <- spec$par0
  f_best <- eval_f(theta_best)
  if (!inherits(fit, "try-error") &&
      is.finite(fit$objective) && fit$objective < f_best) {
    theta_best <- fit$solution
    f_best <- fit$objective
  }

  pars <- .prop4i_euc_unpack_joint(theta_best, spec)
  ll <- if (is.null(pars)) {
    -Inf
  } else {
    .prop4i_euc_joint_loglik(
      y, xs, xe,
      pars$mean, pars$k, pars$a, pars$G0,
      approximate = TRUE
    )
  }

  if (inherits(fit, "try-error")) {
    fit <- list(
      status = -999L,
      message = as.character(fit),
      objective = f_best,
      solution = theta_best
    )
  } else {
    fit$objective <- f_best
    fit$solution <- theta_best
  }

  list(
    nlopt = fit,
    pars = pars,
    lLik_approx = ll,
    spec = spec
  )
}


mobius_SvMF_joint_fit_prop4i_euc <- function(
    y, xs, xe, mean, k, a, G0,
    G0reference = NULL, G01behaviour = "p1",
    det_constraint = c("orthogonal", "rotation"),
    intercept = TRUE,
    algorithm = "NLOPT_LD_SLSQP", opts = list(),
    check = TRUE) {

  det_constraint <- match.arg(det_constraint)
  G01behaviour <- match.arg(
    G01behaviour, c("p1", "free", "fixed")
  )

  y <- as.matrix(y)
  xs <- as.matrix(xs)
  xe <- as.matrix(xe)
  p <- ncol(y)

  if (nrow(y) != nrow(xs) || nrow(y) != nrow(xe) ||
      ncol(xs) != p) {
    stop("For Proposition 4(i)+Euclidean, y and xs must be n by p and xe must have n rows.")
  }
  if (check) {
    .prop4ii_check_unit_rows(y, "y")
    .prop4ii_check_unit_rows(xs, "xs")
  }

  start_mean <- as_mobius_link_prop4i_euc(mean, check = FALSE)
  start_a <- .prop4ii_normalise_scales(a, p, a1 = a[1L])
  start_G0 <- .prop4ii_nearest_orthogonal(G0, det_sign = 1)

  initial <- list(
    mean = start_mean,
    Rtilde0 = start_mean$Rtilde0,
    P = start_mean$P,
    beta_s = start_mean$beta_s,
    Qe_star = start_mean$Qe_star,
    Be = start_mean$Be,
    k = k,
    a = start_a,
    G0 = start_G0
  )

  target_signs <- if (det_constraint == "rotation") 1 else c(1, -1)
  fits <- vector("list", length(target_signs))

  for (j in seq_along(target_signs)) {
    sgn <- target_signs[j]

    Rref <- if ((if (det(start_mean$Rtilde0) >= 0) 1 else -1) == sgn) {
      start_mean$Rtilde0
    } else {
      .prop4ii_procrustes(y, xs, det_sign = sgn)$R
    }

    mean_ref <- mobius_link_prop4i_euc(
      Rtilde0 = Rref,
      P = start_mean$P,
      beta_s = start_mean$beta_s,
      Qe_star = start_mean$Qe_star,
      Be = start_mean$Be,
      xe_center = start_mean$xe_center,
      xe_scale = start_mean$xe_scale,
      intercept = start_mean$intercept,
      check = FALSE
    )

    pre_j <- initial
    pre_j$mean <- mean_ref
    pre_j$Rtilde0 <- Rref

    if (G01behaviour == "p1") {
      rebased <- .prop4ii_rebase_G0(
        pre_j$G0, mean_ref$P[, 1L]
      )
      if (is.null(rebased)) {
        stop("Could not rebase G0 to P[,1].")
      }
      pre_j$G0 <- rebased
    }

    fits[[j]] <- .prop4i_euc_fit_joint_component(
      y = y, xs = xs, xe = xe,
      mean_ref = mean_ref,
      prelim = pre_j,
      G0reference = if (is.null(G0reference)) pre_j$G0 else G0reference,
      G01behaviour = G01behaviour,
      algorithm = algorithm,
      opts = opts
    )
    fits[[j]]$det_sign <- sgn
  }

  lla <- vapply(fits, function(z) z$lLik_approx, numeric(1))
  if (all(!is.finite(lla))) {
    stop("All Proposition 4(i)+Euclidean optimisations failed.")
  }

  best <- fits[[which.max(lla)]]
  pars <- best$pars

  # Same residual-axis identifiability convention as the existing models.
  aord <- order(pars$a[-1L], decreasing = TRUE)
  pars$a[-1L] <- pars$a[-1L][aord]
  pars$G0[, -1L] <- pars$G0[, -1L, drop = FALSE][, aord, drop = FALSE]
  if (exists("standardise_col_signs", mode = "function")) {
    pars$G0[, -1L] <- standardise_col_signs(
      pars$G0[, -1L, drop = FALSE]
    )
  }
  if (det(pars$G0) < 0) {
    pars$G0[, p] <- -pars$G0[, p]
  }

  pred <- mobius_link(
    xs = xs, xe = xe, param = pars$mean, check = FALSE
  )

  if (p != 3L) {
    kres <- mobius_SvMF_konly(
      y = y, ymean = pred,
      a = pars$a, G0 = pars$G0
    )
    pars$k <- kres$k
    lLik <- kres$lLik
  } else {
    yrot <- undo_partransport(
      y, pred, G01 = pars$G0[, 1L]
    )
    lLik <- sum(
      SvMF_log_lik_cann(
        yrot,
        SvMF_cann(k = pars$k, a = pars$a, G = pars$G0)
      )
    )
  }

  SvMF_cann_check(
    SvMF_cann(k = pars$k, a = pars$a, G = pars$G0)
  )

  dots <- pmax(-1, pmin(1, rowSums(y * pred)))
  rrI <- rotated_resid(y, pred, north_pole(p))
  rresids_I <- rrI[, -1L, drop = FALSE]
  attr(rresids_I, "samehemisphere") <- attr(rrI, "samehemisphere")
  colnames(rresids_I) <- paste0("r", seq_len(ncol(rresids_I)))

  rresids_G0 <- resid_SvMF_partransport(
    y, pred, G0 = pars$G0, scale = FALSE
  )
  rresids_std <- resid_SvMF_partransport(
    y, pred, pars$k, pars$a, pars$G0, scale = TRUE
  )

  ep_m <- nrow(pars$mean$Qe_star)
  kmean <- p - 1L

  dmean <- .prop4ii_skew_dim(p) +       # Rtilde0
    .prop4ii_skew_dim(p) +              # P
    .prop4ii_stiefel_dim(ep_m, kmean) + # Qe_star
    kmean +                             # Be
    1L                                  # beta_s

  dscale <- max(p - 2L, 0L)
  dG0 <- if (G01behaviour == "free") {
    .prop4ii_skew_dim(p)
  } else {
    .prop4ii_skew_dim(p - 1L)
  }
  DoF <- dmean + 1L + dscale + dG0

  mean_obj <- pars$mean
  out <- list(
    mean = mean_obj, est = mean_obj,
    Rtilde0 = mean_obj$Rtilde0,
    P = mean_obj$P,
    beta_s = mean_obj$beta_s,
    phi = mean_obj$phi,
    psi = mean_obj$psi,
    rs1 = mean_obj$rs1,
    Qe_star = mean_obj$Qe_star,
    Be = mean_obj$Be,
    k = pars$k,
    a = pars$a,
    G0 = pars$G0,
    obj = -lLik / nrow(y),
    nlopt = best$nlopt,
    y = y, xs = xs, xe = xe,
    pred = pred,
    rresids_I = rresids_I,
    rresids_G0 = rresids_G0,
    rresids_std = rresids_std,
    dists = acos(dots),
    DoF = DoF,
    AIC = 2 * DoF - 2 * lLik,
    lLik = lLik,
    initial = initial,
    component_fits = lapply(fits, function(z) {
      list(
        det_sign = z$det_sign,
        status = z$nlopt$status,
        message = z$nlopt$message,
        lLik_approx = z$lLik_approx
      )
    }),
    linktype = list(
      type = "Prop4i",
      det_constraint = det_constraint,
      G01behaviour = G01behaviour,
      intercept = intercept,
      euclidean = TRUE,
      spherical_constraint = "Bs = beta_s I"
    )
  )

  class(out) <- c(
    "mobius_SvMF_prop4i_euc",
    "mobius_SvMF_prop4i",
    "mobius_SvMF",
    "list"
  )
  out
}


mobius_SvMF_prop4i_euc <- function(
    y, xs, xe, mean = NULL, k = NULL, a = NULL,
    G0 = NULL, G0reference = NULL, G01behaviour = "p1",
    doprelim = TRUE,
    det_constraint = c("orthogonal", "rotation"),
    intercept = TRUE,
    algorithm = "NLOPT_LD_SLSQP", opts = list(),
    check = TRUE, ...) {

  det_constraint <- match.arg(det_constraint)
  G01behaviour <- match.arg(
    G01behaviour, c("p1", "free", "fixed")
  )

  dots <- list(...)
  if (length(dots)) {
    opts <- utils::modifyList(opts, dots)
  }

  if (doprelim) {
    preest <- mobius_SvMF_partransport_prelim_prop4i_euc(
      y = y, xs = xs, xe = xe,
      mean = mean, G0 = G0,
      G01behaviour = G01behaviour,
      det_constraint = det_constraint,
      intercept = intercept,
      algorithm = algorithm,
      opts = opts,
      check = check
    )
  } else {
    if (is.null(mean) || is.null(k) ||
        is.null(a) || is.null(G0)) {
      stop("When doprelim=FALSE, mean, k, a, and G0 must all be supplied.")
    }

    mean <- as_mobius_link_prop4i_euc(mean, check = FALSE)
    preest <- list(
      mean = mean,
      Rtilde0 = mean$Rtilde0,
      P = mean$P,
      beta_s = mean$beta_s,
      Qe_star = mean$Qe_star,
      Be = mean$Be,
      k = k,
      a = a,
      G0 = G0
    )
  }

  finalest <- mobius_SvMF_joint_fit_prop4i_euc(
    y = y, xs = xs, xe = xe,
    mean = preest$mean,
    k = if (is.null(k)) preest$k else k,
    a = if (is.null(a)) preest$a else a,
    G0 = preest$G0,
    G0reference = if (is.null(G0reference)) preest$G0 else G0reference,
    G01behaviour = G01behaviour,
    det_constraint = det_constraint,
    intercept = intercept,
    algorithm = algorithm,
    opts = opts,
    check = check
  )

  finalest$preest <- preest
  finalest
}


print.mobius_vMF_prop4i_euc <- function(x, ...) {
  cat("vMF regression with Proposition 4(i)-constrained spherical link and LinEuc covariates\n")
  cat("  beta_s =", format(x$beta_s, digits = 6), "\n")
  cat("  phi =", format(x$phi, digits = 6), "\n")
  cat("  det(Rtilde0) =", format(det(x$Rtilde0), digits = 6), "\n")
  cat("  Be =", paste(format(x$Be, digits = 5), collapse = ", "), "\n")
  cat("  k =", format(x$k, digits = 6), "\n")
  cat("  log-likelihood =", format(x$lLik, digits = 7), "\n")
  cat("  AIC =", format(x$AIC, digits = 7), "\n")
  invisible(x)
}


print.mobius_SvMF_prop4i_euc <- function(x, ...) {
  cat("SvMF regression with Proposition 4(i)-constrained spherical link and LinEuc covariates\n")
  cat("  G0 estimation: parallel-transport/Cayley procedure\n")
  cat("  G01behaviour =", x$linktype$G01behaviour, "\n")
  cat("  beta_s =", format(x$beta_s, digits = 6), "\n")
  cat("  phi =", format(x$phi, digits = 6), "\n")
  cat("  det(Rtilde0) =", format(det(x$Rtilde0), digits = 6), "\n")
  cat("  Be =", paste(format(x$Be, digits = 5), collapse = ", "), "\n")
  cat("  k =", format(x$k, digits = 6), "\n")
  cat("  a =", paste(format(x$a, digits = 5), collapse = ", "), "\n")
  cat("  log-likelihood =", format(x$lLik, digits = 7), "\n")
  cat("  AIC =", format(x$AIC, digits = 7), "\n")
  invisible(x)
}


predict.mobius_vMF_prop4i_euc <- function(
    object, newdata = NULL, xs = NULL, xe = NULL, ...) {

  if (is.null(xs)) xs <- if (is.null(newdata)) object$xs else newdata
  if (is.null(xe)) xe <- object$xe

  mobius_link(
    xs = as.matrix(xs),
    xe = as.matrix(xe),
    param = object$mean
  )
}


predict.mobius_SvMF_prop4i_euc <- function(
    object, newdata = NULL, xs = NULL, xe = NULL, ...) {

  predict.mobius_vMF_prop4i_euc(
    object, newdata = newdata, xs = xs, xe = xe, ...
  )
}


# -----------------------------------------------------------------------------
# Final user-facing dispatchers
# -----------------------------------------------------------------------------

mobius_vMF <- function(
    y, xs = NULL, xe = NULL, start = NULL,
    type = "SpEuc", fix_qs1 = FALSE,
    fix_qe1 = (type == "LinEuc"), intercept = TRUE,
    lb = NULL, ub = NULL,
    det_constraint = c("orthogonal", "rotation"),
    prop4i_algorithm = "NLOPT_LD_SLSQP",
    prop4i_opts = list(),
    prop4ii_algorithm = "NLOPT_LD_SLSQP",
    prop4ii_opts = list(), ...) {

  if (.is_prop4i_type(type) && .prop4ii_has_euc(xe)) {
    dots <- list(...)
    if (length(dots)) {
      prop4i_opts <- utils::modifyList(prop4i_opts, dots)
    }
    return(mobius_vMF_prop4i_euc(
      y = y, xs = xs, xe = xe, start = start,
      det_constraint = det_constraint,
      intercept = intercept,
      algorithm = prop4i_algorithm,
      opts = prop4i_opts
    ))
  }

  .mobius_vMF_before_prop4i_euc(
    y = y, xs = xs, xe = xe, start = start,
    type = type,
    fix_qs1 = fix_qs1,
    fix_qe1 = fix_qe1,
    intercept = intercept,
    lb = lb, ub = ub,
    det_constraint = det_constraint,
    prop4i_algorithm = prop4i_algorithm,
    prop4i_opts = prop4i_opts,
    prop4ii_algorithm = prop4ii_algorithm,
    prop4ii_opts = prop4ii_opts,
    ...
  )
}


mobius_SvMF <- function(
    y, xs = NULL, xe = NULL, mean = NULL, k = NULL, a = NULL,
    G0 = NULL, G0reference = NULL, G01behaviour = "p1",
    type = "LinEuc", fix_qs1 = FALSE,
    fix_qe1 = (type == "LinEuc"), intercept = TRUE,
    doprelim = TRUE,
    det_constraint = c("orthogonal", "rotation"),
    prop4i_algorithm = "NLOPT_LD_SLSQP",
    prop4i_opts = list(),
    prop4ii_algorithm = "NLOPT_LD_SLSQP",
    prop4ii_opts = list(), ...) {

  if (.is_prop4i_type(type) && .prop4ii_has_euc(xe)) {
    return(mobius_SvMF_prop4i_euc(
      y = y, xs = xs, xe = xe,
      mean = mean, k = k, a = a,
      G0 = G0, G0reference = G0reference,
      G01behaviour = G01behaviour,
      doprelim = doprelim,
      det_constraint = det_constraint,
      intercept = intercept,
      algorithm = prop4i_algorithm,
      opts = prop4i_opts,
      ...
    ))
  }

  .mobius_SvMF_before_prop4i_euc(
    y = y, xs = xs, xe = xe,
    mean = mean, k = k, a = a,
    G0 = G0, G0reference = G0reference,
    G01behaviour = G01behaviour,
    type = type,
    fix_qs1 = fix_qs1,
    fix_qe1 = fix_qe1,
    intercept = intercept,
    doprelim = doprelim,
    det_constraint = det_constraint,
    prop4i_algorithm = prop4i_algorithm,
    prop4i_opts = prop4i_opts,
    prop4ii_algorithm = prop4ii_algorithm,
    prop4ii_opts = prop4ii_opts,
    ...
  )
}

# =============================================================================
# End Proposition 4(i) + Euclidean extension
# =============================================================================
