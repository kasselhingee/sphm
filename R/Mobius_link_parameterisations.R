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
#' @param ce A single real value. The additive offset \eqn{c_e} in the denominator below Euclidean covariates in the 'H' link. `NULL` if no Euc covariates.
#' @details
#' # Cannonical Parameterisation
#' The parameters here are for a mean link defined as
#' \deqn{\mu_{SpEuc}(x) = P\mathcal{S}^{[-1]}\left(B_s \mathcal{S}(Q_s^\top x_s)  +  \frac{B_e(Q_e[,-1]^\top x_e}{Q_e[,1]^\top x_e + c_e}\right)}
#' where \eqn{\mu(x)} is a p-length unit vector,
#' \eqn{x_s} is a `qs`-length unit vector,
#' \eqn{x_e} is a `qe`-length vector,
#' `P` is a p x p rotation matrix,
#' `Qs` is a `qs x p` orthonormal matrix (`t(Qs) %*% Qs == diag(p)),
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

#' @describeIn mnlink_params Represent the parameters of the Mobius mean link [`mnlink()`] using the canonical form.
#' @export
mnlink_cann <- function(P, Bs = NULL, Qs = NULL, Be = NULL, Qe = NULL, ce = NULL, check = TRUE){
  stopifnot(is.matrix(P))
  obj <- list(P = P, Bs = Bs, Qs = Qs, Be = Be, Qe = Qe, ce = ce)
  obj <- lapply(obj, function(element){if (length(element) == 0){return(NULL)}else{return(element)}})

  # export
  class(obj) <- c("mnlink_cann", class(obj))
  if (check){mnlink_cann_check(obj)}
  return(obj)
}

#' @describeIn mnlink_params Represent the parameters of the Mobius mean link [`mnlink()`] using the canonical form.
#' @export
as_mnlink_cann <- function(obj){
  if (inherits(obj, "mnlink_cann")){return(obj)}
  if (inherits(obj, "mnlink_Omega")){return(Omega2cann(obj, check = FALSE))}
  if (!inherits(obj, "list")){stop("obj isn't a mnlink_cann or mnlink_Omega class, or a list.")}
  if ("P" %in% names(obj)){return(do.call(mnlink_cann, c(obj, list(check = FALSE))))}
  if ("p1" %in% names(obj)){return(Omega2cann(do.call(mnlink_Omega, c(obj, list(check = FALSE))), check = FALSE))}
  return(obj)
}

mnlink_cann_vec <- function(obj){
  stopifnot(inherits(obj, "mnlink_cann"))
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




#' @name mnlink_params
#' @param p1 First column of the P matrix (vector of length `p`)
#' @param qe1 First column of the Qe matrix (vector of length `qe`). `NULL` if no Euc covariates.
#' @param qs1 First column of the Qs matrix (vector of length `qs`). `NULL` if no Sph covariates.
#' @param Omega A `p` by `qe + qs` matrix representing 
#' \deqn{\Omega = [\Omega_s \Omega_e] = [P^* B_s {Q_s^*}^T \,   P^* B_e {Q_e^*}^T]}
#' @param ce The value of \eqn{c_e}.
NULL

#' @describeIn mnlink_params Represent the parameters of the Mobius mean link [`mnlink()`] using the Omega form.
#' @export
mnlink_Omega <- function(p1, qs1 = vector("numeric", 0), Omega, qe1 = vector("numeric", 0), ce = vector("numeric", 0), check = TRUE){
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
  class(obj) <- c("mnlink_Omega", class(obj))
  if (check) {mnlink_Omega_check(obj)}
  return(obj)
}

#' @describeIn mnlink_params Represent the parameters of the Mobius mean link [`mnlink()`] using the Omega form.
#' @export
as_mnlink_Omega <- function(obj){
  if (inherits(obj, "mnlink_cann")){return(cann2Omega(obj, check = FALSE))}
  if (inherits(obj, "mnlink_Omega")){return(obj)}
  if (!inherits(obj, "list")){stop("obj must be either a mnlink_cann or mnlink_Omega class, or a list.")}
  if ("P" %in% names(obj)){return(cann2Omega(do.call(mnlink_cann, c(obj, list(check = FALSE))), check = FALSE))}
  if ("p1" %in% names(obj)){return(do.call(mnlink_Omega, c(obj, list(check = FALSE))))}
  return(obj)
}

#' @noRd
#' @description Converts Omega parameterisation to a vector.
#' @param obj An `mnlink_Omega` parameter object.
mnlink_Omega_vec <- function(obj){
  stopifnot(inherits(obj, "mnlink_Omega"))
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
  class(out) <- "mnlink_Omega_vec"
  return(out)
}

#' @noRd
#' @description Inverse of `mnlink_Omega_vec()`
#' @param vec is a vector created by `mnlink_Omega_vec()`
#' @param p The dimension of the response (The dimension of covariates will be infered from `p`).
#' @param qe Number of Euclidean covariates
mnlink_Omega_unvec <- function(vec, p, qe = 0, check = TRUE){
  # length of vec = p + qs + qe + p*(qs + qe) + (qe>0)
  # l - 1*(qe>0) - qe - p*qe = qs + p*qs
  # (l - (qe>0) - qe - p*qe)/(1+p) = qs
  qs <- (length(vec) - p - (qe>0) - qe - p*qe)/(1+p)
  stopifnot(qs == as.integer(qs))
  names(vec) <- NULL
  
  mnlink_Omega(p1 = vec[1:p],
           qs1 = vec[p + seq.int(1, length.out = qs)],
           qe1 = vec[p + qs + seq.int(1, length.out = qe)],
           Omega = matrix(vec[p + qs + qe + seq.int(1, length.out = p*(qe + qs))], nrow = p, ncol = qs + qe, byrow = FALSE),
           ce = vec[p + qs + qe + p*(qs + qe) + seq.int(1, length.out = (qe > 0))],
           check = check)
}

#' @noRd
#' @describeIn mnlink_params For converting between parameterisations of the link function.
#' Sign of columns of P and Q are lost by the `cann2Omega()` transformation.
#' @export
cann2Omega <- function(obj, check = TRUE){
  if (check){mnlink_cann_check(obj)}
  p1 <- obj$P[, 1]
  qs1 <- obj$Qs[, 1]
  qe1 <- obj$Qe[, 1]
  Omega_s <- Omega_e <- ce <- NULL
  if (!is.null(qs1)){
    Omega_s <- obj$P[,-1] %*% obj$Bs %*% t(obj$Qs[, -1])
  }
  if (!is.null(qe1)){
    Omega_e <- obj$P[,-1] %*% obj$Be %*% t(obj$Qe[, -1])
    ce <- obj$ce
  }
  Omega <- cbind(Omega_s, Omega_e)
  return(mnlink_Omega(p1, qs1 = qs1, Omega = Omega, qe1 = qe1, ce = ce, check = FALSE))
}

#' @noRd
#' @description Converts an Omega parameterisation to the cannonical parameterisation.
#' Apart from p1 and q1, sign of columns of P and Q cannot be recovered from Omega.
Omega2cann <- function(obj, check = TRUE){
  if (check){mnlink_Omega_check(obj)}
  svdres <- svd(obj$Omega, nu = nrow(obj$Omega) - 1, nv = nrow(obj$Omega) - 1)

  P <- cbind(obj$p1, svdres$u)
  
  # much of the rest uses the SVD of Omega as written in the Euclidean link document (not ce)
  Qs <- Qe <- Bs <- Be <- ce <- NULL
  if (length(obj$qs1) > 0){
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
  
  out <- mnlink_cann(P, Bs = Bs, Qs = Qs, Be = Be, Qe = Qe, ce = ce, check = check)
  if (det(out$P) < 0){#if negative determinant, flip an axis to get positive determinant of P (i.e. a rotation matrix)
    out <- P_signswitch(out, ncol(out$P))
  }
  return(out)
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

mnlink_Omega_check <- function(obj){
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
  
  vals <- mnlink_Omega_check_numerical(obj)
  good <- (vals < sqrt(.Machine$double.eps))
  if (!all(good)){
    stop(paste("The following checks failed.", 
               paste0(names(vals)[!good], ": ", format(sqrt(vals[!good]), digits = 2), collapse = ", ") #sqrt here converts squared sizes to actual sizes
    ))
  }
  return(NULL)
}
mnlink_Omega_check_numerical <- function(obj){ #uses squared values for smoothness
  stopifnot(inherits(obj, "mnlink_Omega"))
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
  if ((qs > 0) & (qe > 0)) {# for commutivity check
    OmOm <- obj$Omega %*% t(obj$Omega) 
    Is_tilde <- diag(1, qs + qe, qs) # for commutivity check
    OmpartOmpart <- obj$Omega %*% (Is_tilde %*% t(Is_tilde)) %*% t(obj$Omega) # for commutivity check
    Omega_comm = (sum((OmOm %*% OmpartOmpart - OmpartOmpart %*% OmOm)^2)^2) # for commutivity check - Frobenius norm of 0
    # above second ^2 is a hack to make Omega_comm be on a similar scale to all the other tests. See the C++ implementation of this constraint for a vector that is a more refined representation of the constraint
    checkvals <- c(checkvals, Omega_comm = Omega_comm)
  }
  return(checkvals)
}

# if the values are close to satisfying the constraints, it might make sense to project and scale p1 and q1 to satisfy the constraints
# will use canonical parameterisation to do this because orthogonality of the columns will make for easier projections
Omega_proj <- function(obj){
  stopifnot(inherits(obj, "mnlink_Omega"))
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

# The Euc part of the Omega parameterisations is invariant to sign because ce can be anything (if it was fixed to +1 then there would be no invariance). This function switches the sign
Euc_signswitch <- function(obj){
  if (inherits(obj, "mnlink_Omega")){
    obj$qe1 <- -1 * obj$qe1
    obj$ce <- -1 * obj$ce
    obj$Omega[, length(obj$qs1) + (1:length(obj$qe1))] <- -1 * obj$Omega[, length(obj$qs1) + (1:length(obj$qe1))]
    return(obj)
  }
  if (inherits(obj, "mnlink_cann")){
    obj$ce <- -1 * obj$ce
    obj$Qe <- -1 * obj$Qe
    return(obj)
  }
}
P_signswitch <- function(obj, cols){
  stopifnot(inherits(obj, "mnlink_cann"))
  if (is.logical(cols)){cols <- which(cols)}
  if (1 %in% cols){stop("Sign of first column of P cannot be easily swapped")}
  obj$P[, cols] <- -1 * obj$P[,cols]
  if (!is.null(obj$Qs)){obj$Qs[,cols] <- -1 * obj$Qs[,cols]}
  if (!is.null(obj$Qe)){
    obj$Qe[,cols] <- -1 * obj$Qe[,cols]
  }
  return(obj)
}

#' @title Randomly generate mean link parameters
#' @description Uniformly generates orthogonal matrices `Qs`, `Qe` and `P`.
#' The diagonal elements of `Bs` and `Be` are drawn uniformly from between `0` and `1`.
#' The parameter `ce` is generated uniformly from between `0` and `1`.
#' @param preseed Before the generation of each parameter matrix, the `seed` is incremented by `1`, starting from `preseed`.
#' @param p The length of response vectors
#' @param qs The length of spherical covariate vectors
#' @param qe The length of Euclidean covariate vectors
#' @return A `mnlink_cann` object.
#' @export
rmnlink_cann <- function(p = 3, qs = 5, qe = 4, preseed = 0){
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
  paramobj <- mnlink_cann(P, Bs = Bs, Qs = Qs, Be = Be, Qe = Qe, ce = ce, check = FALSE)
  return(paramobj)
}

#' @noRd
#' @description To simplify testing places result of `rmnlink_cann()` into the calling environment.
#' @details
#' Places all the parameters in the environment by running list2env(obj, envir = environment())
rmnlink_cann__place_in_env <- function(p = 3, qs = 5, qe = 4, preseed = 0){
  paramobj <- rmnlink_cann(p = p, qs = qs, qe = qe, preseed = preseed)
  target_env <- if (is.null(environment(-1))) .GlobalEnv else environment(-1)
  list2env(c(paramobj, list(paramobj = paramobj, qs = qs, qe = qe, p = p)), envir = target_env)
  return(NULL)
}

is_LinEuc <- function(obj, tol = sqrt(.Machine$double.eps)){
  if (inherits(obj, "mnlink_Omega")){
    if (length(obj$qe1) > 0){
      checks <- c((obj$qe1 - c(1, rep(0, length(obj$qe1) - 1)))^2, (obj$ce - 1)^2)
      if (all(checks < tol)){return(TRUE)}else{return(FALSE)}
    }
  }
  if (inherits(obj, "mnlink_cann")){
    if (!is.null(obj$Qe)){
      checks <-  c((obj$Qe[,1] - c(1, rep(0, nrow(obj$Qe) - 1)))^2, (obj$ce[1] - 1)^2)
      if (all(checks < tol)){return(TRUE)}else{return(FALSE)}
    }
  }
  warning("obj doesn't have Euclidean component")
  return(NULL)
}

#' @title Obtain dimensions corresponding to a mean link parameter set
#' @param x A mean link parameter object of class `mnlink_Omega` or `mnlink_cann`.
#' @return An integer vector of p (length of response unit vectors), qs (length of spherical covariate unit vectors) and qe (length of Euclidean covariate vectors).
#' @name dim-mnlink_param
NULL

#' @rdname dim-mnlink_param
#' @export
#' @method dim mnlink_cann
dim.mnlink_cann <- function(x){
  c(p = ncol(x$P), 
    qs = switch(1 + is.null(x$Qs), nrow(x$Qs), 0),
    qe = switch(1 + is.null(x$Qe), nrow(x$Qe), 0))
}

#' @rdname dim-mnlink_param
#' @export
#' @method dim mnlink_Omega
dim.mnlink_Omega <- function(x){
  c(p = length(x$p1), 
    qs = length(x$qs1),
    qe = length(x$qe1))
}
