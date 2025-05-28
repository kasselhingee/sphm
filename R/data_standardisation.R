# data standardisations that can be incorporated into the mean link

#' Standardise Spherical Data to North Pole
#' @description Rotation according to the eigenvectors of the second moment of the data
#' projected perpendicular to the mean direction. Standardised data have mean of c(1, 0, 0,...). See Scealy and Wood 2019 for details.
#' @param y Data on the sphere. Each row is a data point in Cartesian coordinates.
#' @param G Axes of the second moment matrix of the data, projected so that the first column is the global mean.
#' @details Each returned data point is `t(G) %*% y` of the original data point y, where `G` is computed by `standardise_mat()`.
#' @export
standardise_sph <- function(y, tG = t(standardise_mat(y))){
  ystd <- y %*% t(tG)
  ystd <- unname(ystd)
  attr(ystd, "std_rotation") <- tG
  return(ystd)
}

destandardise_sph <- function(y, tG){
  ydestd <- y %*% tG
  attr(ydestd, "std_rotation") <- NULL
  return(ydestd)
}

#' @describeIn standardise_sph
#' @export
standardise_mat <- function(y){
  p <- ncol(y)
  mn <- colMeans(y)
  mn <- mn/sqrt(sum(mn^2))
  mnproj <- diag(1, p) - mn %*% t(mn)

  mom2 <- t(y) %*% y #quickly calculates the sum of projection matrices of rows of y
  projmom2 <- mnproj %*% mom2 %*% mnproj

  # use the non-degenerate directions in Ghat, then an arbitrary basis for the remaining except for the mean!!
  # Avoiding the mean direction is tricky.
  # one direction (corresponding to the mean) will ALWAYS be degenerate
  eig_projmom2 <- eigen(projmom2)
  degeneratedirections <- abs(eig_projmom2$values) <= sqrt(.Machine$double.eps)
  if (sum(degeneratedirections) > 1){
    # making the mom2 a little ridgier, then project again to remove the mean and continue
    eig_projmom2$values[degeneratedirections] <- eig_projmom2$values[degeneratedirections] + min(eig_projmom2$values[!degeneratedirections])/seq.int(2, length.out = sum(degeneratedirections))
    projmom2 <- mnproj %*% eig_projmom2$vectors %*% diag(eig_projmom2$values) %*% t(eig_projmom2$vectors) %*% mnproj
    eig_projmom2 <- eigen(projmom2)
  }
  Ghat <- cbind(mn, eigen(projmom2)$vectors[, 1:(p-1)])
  # make sure that Ghat is a rotation matrix, by making sure determinant is 1
  if (det(Ghat) < 0){Ghat[, p] <- -Ghat[, p]}
  return(Ghat)
}

# For covariates that are all a const != 0, leave unchanged
# Link isn't equivariant to scaling so dont do scaling
standardise_Euc <- function(xe){
  xe_names <- colnames(xe)
  constcovars <- apply(xe, 2, sd) < sqrt(.Machine$double.eps)
  # do shift and rotation of the other covariates
  xe_pcares <- princomp(xe[, !constcovars], cor = FALSE, scores = TRUE)
  # xe_pcares$scores is equivalent to xe[, !constcovars] %*% xe_pcares$loadings with centering
  
  # keep contsant covariats at the left
  xe <- cbind(xe[, constcovars, drop = FALSE], xe_pcares$scores)
  
  # update loadings and center for any constcovars
  xe_centers <- c(rep(0, sum(constcovars)), xe_pcares$center)
  names(xe_centers) <- xe_names
  xe_loadings <- diag(ncol(xe))
  xe_loadings[seq.int(1+sum(constcovars), length.out = ncol(xe_pcares$loadings)),
              seq.int(1+sum(constcovars), length.out = ncol(xe_pcares$loadings))] <-
    xe_pcares$loadings
  rownames(xe_loadings) <- xe_names
  colnames(xe_loadings)[constcovars] <- xe_names[constcovars]
  colnames(xe_loadings)[seq.int(1+sum(constcovars), length.out = ncol(xe_pcares$loadings))] <- colnames(xe_pcares$loadings)
  
  attr(xe, "std_rotation") <- t(xe_loadings)
  attr(xe, "std_center") <- xe_centers
  return(xe)
}

destandardise_Euc <- function(xe, center, rotation){
  stopifnot(is.vector(center))
  stopifnot(all(abs(t(rotation) %*% rotation - diag(ncol(rotation))) < sqrt(.Machine$double.eps)))
  xe <- xe %*% rotation #because xe is row vectors, destandardisation uses non-transpose of the forward rotation
  xe <- t(t(xe) + center)
  return(xe)
}

#' @noRd
#' @title Change parameters of mean link to match data standardisation
#' @description
#' Spherical covariates and response can be standardised by a rotation.
#' Euclidean covariates can be standardised by a shift and then rotation.
#' After the standardisation, the relationship between covariates and response
#' has the same form, but the parameters are different.
#' This function computes those parameters.
#' @details
#' Suppose spherical \eqn{xs} and Euclidean covariate \eqn{xe} are related to a spherical 
#' response \eqn{y} by [`meanlink()`].
#' Then `xsrot %*% xs` and `xerot %*% (xe - xecenter)` are related to
#' \eqn{`yrot` \times y} according to `meanlink()` with parameters given by
#' [`recoordinate_Omega()`].
#' 
#' Reversing this coordinate change can be performed by [`undo_recoordinate_Omega()`].
#' 
#' The standardisation of the data `xs`, `xe` and `y` can be performed by 
#' [`standardise_sph()`] and [`standardise_Euc()`].
#' 
#' I'm only sure that the following works when there is no shift of the 1s covariate and when param$qe1[onescovaridx] = 0
#' @param onescovaridx Gives the index in `xe` of the covariate that is identically 1 - needed whenever xecenter is non-zero.
recoordinate_Omega <- function(param, yrot = diag(length(param$p1)), 
                               xsrot = diag(length(param$qs1)),
                               xerot = diag(length(param$qe1)), 
                               xecenter = rep(0, length(param$qe1)),
                               onescovaridx = 1){
  stopifnot(inherits(param, "mnlink_Omega"))
  # in case arguments passed are NULL set to default
  if (is.null(yrot)){yrot <- diag(length(param$p1))}
  if (is.null(xsrot)){xsrot <- diag(length(param$qs1))}
  if (is.null(xerot)){xerot <- diag(length(param$qe1))}
  if (is.null(xecenter)){xecenter <- rep(0, length(param$qe1))}
  qs <- length(param$qs1)
  qe <- length(param$qe1)
  omstd <- om <- param
  
  # if xecenter non-zero, check onescovaridx
  if (any(abs(xecenter) > .Machine$double.eps)){
    # xecenter shouldn't shift the 1s covariate
    stopifnot(abs(xecenter[onescovaridx]) < .Machine$double.eps)
    # if qe1 is non-zero for the 1s covariate, then I dont know what it would mean
    if (abs(om$qe1[onescovaridx]) > .Machine$double.eps){
      warning(sprintf("qe1[onescovaridx]=%f is non-zero and incorporating a shift in the other covariates may not work", om$qe1[onescovaridx]))
    }
    # if the column of Omega related to the 1s covariate is non-zero
    # then the parameters already incorporate a shift
    if (any(abs(om$Omega[, qs+onescovaridx]) > .Machine$double.eps)){
      warning(paste("A shift of Euclidean covariates is already included in the parameters. Omega * shift =",  paste(om$Omega[, qs+onescovaridx], collapse = ", ")))
    }
  }
  
  # update for centering first
  omstd$ce <- om$ce + drop(omstd$qe1 %*% xecenter) 
  omstd$Omega[, qs + onescovaridx] <- om$Omega[, qs + onescovaridx, drop = FALSE] + om$Omega[,seq.int(qs + 1, length.out = qe)] %*% xecenter
  # Here unclear how to update qe1 based on centering if qe1[onescovaridx] is non-zero
  
  # update based on rotations xsrot and xerot
  omstd$qs1 <- drop(xsrot %*% omstd$qs1)
  omstd$qe1 <- drop(xerot %*% omstd$qe1)
  omstd$Omega[, seq.int(1, length.out = qs)] <- omstd$Omega[, seq.int(1, length.out = qs)] %*% t(xsrot)
  omstd$Omega[, seq.int(qs + 1, length.out = qe)] <- omstd$Omega[, seq.int(qs + 1, length.out = qe)] %*% t(xerot)
  
  # add yrot change
  omstd$Omega <- yrot %*% omstd$Omega
  omstd$p1 <- drop(yrot %*% omstd$p1)
  return(omstd)
}

undo_recoordinate_Omega <- function(param, yrot = diag(length(param$p1)), 
                                    xsrot = diag(length(param$qs1)),
                                    xerot = diag(length(param$qe1)), 
                                    xecenter = rep(0, length(param$qe1)),
                                    onescovaridx = 1){
  stopifnot(inherits(param, "mnlink_Omega"))
  # in case arguments passed are NULL set to default
  if (is.null(yrot)){yrot <- diag(length(param$p1))}
  if (is.null(xsrot)){xsrot <- diag(length(param$qs1))}
  if (is.null(xerot)){xerot <- diag(length(param$qe1))}
  if (is.null(xecenter)){xecenter <- rep(0, length(param$qe1))}
  qs <- length(param$qs1)
  qe <- length(param$qe1)
  om <- omstd <- param
  
  # if xecenter non-zero, check onescovaridx
  if (any(abs(xecenter) > .Machine$double.eps)){
    # xecenter shouldn't shift the 1s covariate
    stopifnot(abs(xecenter[onescovaridx]) < .Machine$double.eps)
    # if qe1 is non-zero for the 1s covariate, then I dont know what it would mean
    if (abs(om$qe1[onescovaridx]) > .Machine$double.eps){
      stop(sprintf("qe1[onescovaridx]=%f is non-zero and incorporating a shift in the other covariates may not work", om$qe1[onescovaridx]))
    }
  }
  
  # First undo effect of xsrot, xerot and yrot
  om$qs1 <- drop(t(xsrot) %*% omstd$qs1)
  om$qe1 <- drop(t(xerot) %*% omstd$qe1)
  om$Omega[, seq.int(1, length.out = qs)] <- omstd$Omega[, seq.int(1, length.out = qs)] %*% xsrot
  om$Omega[, seq.int(qs + 1, length.out = qe)] <- omstd$Omega[, seq.int(qs + 1, length.out = qe)] %*% xerot
  om$p1 <- drop(t(yrot) %*% omstd$p1)
  om$Omega <- t(yrot) %*% om$Omega
  
  # use the above calculated qe1 and Omega to get Omega 
  om$ce <- om$ce - drop(om$qe1 %*% xecenter) 
  om$Omega[, qs + onescovaridx] <- om$Omega[, qs + onescovaridx, drop = FALSE] - om$Omega[,seq.int(qs + 1, length.out = qe)] %*% xecenter
  
  return(om)
}
