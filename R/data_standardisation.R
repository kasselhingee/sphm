# data standardisations that can be incorporated into the mean link

#' Standardise Spherical Data to North Pole
#' @description Rotation according to the eigenvectors of the second moment of the data
#' projected perpendicular to the mean direction. Standardised data have mean of c(1, 0, 0,...). See Scealy and Wood 2019 for details.
#' @param y Data on the sphere. Each row is a data point in Cartesian coordinates.
#' @param G Axes of the second moment matrix of the data, projected so that the first column is the global mean.
#' @details Each returned data point is `t(G) %*% y` of the original data point y, where `G` is computed by `standardise_mat()`.
#' @export
standardise <- function(y, tG = t(standardise_mat(y))){
  ystd <- y %*% t(tG)
  ystd <- unname(ystd)
  attr(ystd, "std_rotation") <- tG
  return(ystd)
}

destandardise <- function(y, tG){
  ydestd <- y %*% tG
  attr(ydestd, "std_rotation") <- NULL
  return(ydestd)
}

#' @describeIn standardise
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

# For covariates that are all a const != 0, standardise to 1
standardise_Euc <- function(xe){
  xe_names <- colnames(xe)
  constcovars <- apply(xe, 2, sd) < sqrt(.Machine$double.eps)
  # do shift and rotation of the other covariates
  xe_pcares <- princomp(xe[, !constcovars], cor = FALSE, scores = TRUE)
  # xe_pcares$scores is equivalent to xe[, !constcovars] %*% xe_pcares$loadings with centering
  
  # keep contant covariats at the left
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
  xe <- xe %*% rotation #because xe is row vectors, destandardisation uses non-transpose of the forward rotation
  xe <- t(t(xe) + center)
  return(xe)
}

recoordinate_cann <- function(param, yrot = diag(nrow(param$P)), 
                              xsrot = diag(nrow(param$Qs)),
                              xerot = diag(nrow(param$Qe)), 
                              xecenter = rep(0, nrow(param$Qe))){
  stopifnot(inherits(param, "mnlink_cann"))
  paramstd <- param
  paramstd$Qs <- xsrot %*% param$Qs
  paramstd$Qe <- xerot %*% param$Qe
  paramstd$ce <- drop(t(param$Qe) %*% xecenter + param$ce)
  paramstd$P <-  yrot %*% param$P
  return(paramstd)
}

recoordinate_Omega <- function(param, yrot = diag(length(param$p1)), 
                               xsrot = diag(length(param$qs1)),
                               xerot = diag(length(param$qe1)), 
                               xecenter = rep(0, length(param$qe1))){
  stopifnot(inherits(param, "mnlink_Omega"))
  qs <- length(param$qs1)
  qe <- length(param$qe1)
  omstd <- om <- param
  omstd$qs1 <- drop(xsrot %*% om$qs1)
  omstd$qe1 <- drop(xerot %*% om$qe1)
  omstd$ce1 <- drop(t(om$qe1) %*% xecenter) + om$ce1
  omstd$PBce <- drop(om$Omega[, seq.int(qs + 1, length.out = qe)] %*% xecenter + om$PBce)
  omstd$Omega[, seq.int(1, length.out = qs)] <- om$Omega[, seq.int(1, length.out = qs)] %*% t(xsrot)
  omstd$Omega[, seq.int(qs + 1, length.out = qe)] <- om$Omega[, seq.int(qs + 1, length.out = qe)] %*% t(xerot)
  # add std of P
  omstd$PBce <- drop(yrot %*% omstd$PBce)
  omstd$Omega <- yrot %*% omstd$Omega
  omstd$p1 <- drop(yrot %*% omstd$p1)
  return(omstd)
}