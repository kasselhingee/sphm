# data standardisations that can be incorporated into the mean link

#' @title Standardise Spherical Data
#' @description Rotate spherical data according to the eigenvectors of the second moment of the data
#' projected perpendicular to the mean direction. Standardised data have mean of (1, 0, 0,...). See Scealy and Wood (2019, Section 4.1.3) for details.
#' @param y Matrix of row vectors representing data on the sphere as unit vectors in Cartesian coordinates.
#' @param rotation The rotation matrix to apply to `y`, as returned by [`second_moment_mat()`]
#'   (or its transpose). Defaults to `t(second_moment_mat(y))`, which computes the rotation
#'   from `y` itself. Supply this to apply the same rotation to a second dataset, e.g.
#'   `standardise_sph(xs, rotation = attr(standardise_sph(y), "std_rotation"))`.
#' @details Each returned data point is `rotation %*% y`. The mean direction of this returned
#'   data will be (1, 0, 0, ...) and [`second_moment_mat()`] of this standardised data will
#'   return the identity matrix. The rotation used is stored as the `"std_rotation"` attribute
#'   of the returned matrix.
#' @return A matrix of the same dimension as `y`.
#' @export
standardise_sph <- function(y, rotation = t(second_moment_mat(y))){
  ystd <- y %*% t(rotation)
  ystd <- unname(ystd)
  attr(ystd, "std_rotation") <- rotation
  return(ystd)
}

# Reverse the rotation applied by standardise_sph(). rotation must be the std_rotation
# attribute returned by standardise_sph() (i.e. the rotation matrix, not its transpose).
destandardise_sph <- function(y, rotation){
  ydestd <- y %*% rotation
  attr(ydestd, "std_rotation") <- NULL
  return(ydestd)
}

#' @describeIn standardise_sph Returns a matrix with columns that are the mean direction and then eigenvectors of \eqn{\sum_{i=1}^n y_i y_i^\top} projected orthogonal to the mean.
#' @export
second_moment_mat <- function(y){
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
  # standardise Ghat by first making rows positive
  Ghat[,-1] <- standardise_col_signs(Ghat[,-1])
  # make sure that Ghat is a rotation matrix, by making sure determinant is 1
  if (det(Ghat) < 0){Ghat[, p] <- -Ghat[, p]}
  return(Ghat)
}

# Centre and rotate Euclidean covariates by PCA (principal components). Constant
# covariates are left unchanged. Stores the PCA centre and loadings as attributes
# "std_center" and "std_rotation" so that undo_recoordinate_Omega() can reverse the
# effect on the mean link parameters. Scaling by SD is NOT applied because the
# link function is not equivariant to scaling.
standardise_Euc <- function(xe){
  xe_names <- colnames(xe)
  constcovars <- apply(xe, 2, sd) < sqrt(.Machine$double.eps)
  # do shift and rotation of the other covariates
  xe_pcares <- princomp(xe[, !constcovars], cor = FALSE, scores = TRUE)
  # xe_pcares$scores is equivalent to xe[, !constcovars] %*% xe_pcares$loadings with centering
  
  # update only the non-constant covariates
  xe[,!constcovars] <- xe_pcares$scores
  
  # update loadings and center for any constcovars
  xe_centers <- rep(0, ncol(xe))
  xe_centers[!constcovars] <- xe_pcares$center
  names(xe_centers) <- xe_names
  xe_loadings <- diag(ncol(xe))
  xe_loadings[!constcovars, !constcovars] <- xe_pcares$loadings
  rownames(xe_loadings) <- xe_names
  colnames(xe_loadings)[constcovars] <- xe_names[constcovars]
  colnames(xe_loadings)[!constcovars] <- colnames(xe_pcares$loadings)
  
  attr(xe, "std_rotation") <- t(xe_loadings)
  attr(xe, "std_center") <- xe_centers
  return(xe)
}

# Reverse the PCA centering and rotation applied by standardise_Euc().
# rotation must be the std_rotation attribute (the PCA loadings matrix).
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
#' response \eqn{y} by [`mobius_link()`].
#' Then `xsrot %*% xs` and `xerot %*% (xe - xecenter)` are related to
#' \eqn{`yrot` \times y} according to `mobius_link()` with parameters given by
#' `recoordinate_Omega()`.
#'
#' Reversing this coordinate change can be performed by `undo_recoordinate_Omega()`.
#'
#' The standardisation of the data `xs`, `xe` and `y` can be performed by
#' [`standardise_sph()`] and `standardise_Euc()`.
#' 
#' Note: only verified correct when there is no shift of the 1s covariate
#' and when `param$qe1[onescovaridx] = 0`.
#' `onescovaridx` gives the index in `xe` of the covariate that is identically 1.
recoordinate_Omega <- function(param, yrot = diag(length(param$p1)), 
                               xsrot = diag(length(param$qs1)),
                               xerot = diag(length(param$qe1)), 
                               xecenter = rep(0, length(param$qe1)),
                               onescovaridx = 0){
  stopifnot(inherits(param, "mobius_link_Omega"))
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
    if ((length(onescovaridx) == 0) || !is.finite(onescovaridx)){
      stop("Please specify the index corresponding to the 1s covariate.")
    }
    if ((onescovaridx < 1) || (onescovaridx > qe)){stop("onescovardix is outside the range possible Euc covariate indexes")}
    # xecenter shouldn't shift the 1s covariate
    stopifnot(abs(xecenter[onescovaridx]) < .Machine$double.eps)
    # if qe1 is non-zero for the 1s covariate, then I dont know what it would mean
    if (abs(om$qe1[onescovaridx]) > .Machine$double.eps){
      warning(sprintf("qe1[onescovaridx]=%f is non-zero and incorporating a shift in the other covariates may not work", om$qe1[onescovaridx]))
    }
  }
  
  # update for centering first
  omstd$ce <- om$ce + drop(omstd$qe1 %*% xecenter) 
  if ((onescovaridx > 0) && (onescovaridx < qe)){
    omstd$Omega[, qs + onescovaridx] <- om$Omega[, qs + onescovaridx, drop = FALSE] + om$Omega[,seq.int(qs + 1, length.out = qe)] %*% xecenter
  }
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

# Reverse the parameter transformation applied by recoordinate_Omega().
# Given parameters in standardised coordinates, returns parameters in original coordinates.
undo_recoordinate_Omega <- function(param, yrot = diag(length(param$p1)),
                                    xsrot = diag(length(param$qs1)),
                                    xerot = diag(length(param$qe1)), 
                                    xecenter = rep(0, length(param$qe1)),
                                    onescovaridx = 0){
  stopifnot(inherits(param, "mobius_link_Omega"))
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
    if ((onescovaridx < 1) || (onescovaridx > qe)){stop("onescovardix is outside the range possible Euc covariate indexes")}
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
  if ((onescovaridx > 0) && (onescovaridx < qe)){
    om$Omega[, qs + onescovaridx] <- om$Omega[, qs + onescovaridx, drop = FALSE] - om$Omega[,seq.int(qs + 1, length.out = qe)] %*% xecenter
  }
  
  return(om)
}

# Standardise all data in preplist (y, xs, xe) and update the starting link parameters
# to match the standardised coordinate system. If intercept = FALSE, xe is not standardised.
standardise_data <- function(preplist, intercept){
  # standardise inputs
  preplist$y <- standardise_sph(preplist$y)
  if (!is.null(preplist$xs)){preplist$xs <- standardise_sph(preplist$xs)}
  if (!is.null(preplist$xe) && intercept){preplist$xe <- standardise_Euc(preplist$xe)}
  
  # update starting coordinates if provided
  if (!is.null(preplist$start)){
    preplist$start <- recoordinate_Omega(as_mobius_link_Omega(preplist$start), 
                                yrot = attr(preplist$y, "std_rotation"), 
                                xsrot = attr(preplist$xs, "std_rotation"), #if xs/xe is NULL then attr(xs/xe, ..) is NULL too
                                xerot = attr(preplist$xe, "std_rotation"), 
                                xecenter = attr(preplist$xe, "std_center"),
                                onescovaridx = preplist$onescovaridx)
  }
  return(preplist)
}

# I dont think I use this function at all
destandardise_data <- function(preplist, intercept){
  preplist$y = destandardise_sph(preplist$y, attr(preplist$y, "std_rotation"))
  preplist$xs = if (!is.null(preplist$xs)){destandardise_sph(preplist$xs, attr(preplist$xs, "std_rotation"))}
  preplist$xe = if (!is.null(preplist$xe)){destandardise_Euc(preplist$xe, attr(preplist$xe, "std_center"), attr(preplist$xe, "std_rotation"))}
  preplist$start <- undo_recoordinate_Omega(as_mobius_link_Omega(preplist$start), 
                          yrot = attr(preplist$y, "std_rotation"), 
                          xsrot = attr(preplist$xs, "std_rotation"), #if xs/xe is NULL then attr(xs/xe, ..) is NULL too
                          xerot = attr(preplist$xe, "std_rotation"), 
                          xecenter = attr(preplist$xe, "std_center"),
                          onescovaridx = preplist$onescovaridx)
  return(preplist)
}



# Add dummy Euclidean covariates as required by the link type:
# - For "LinEuc" links, prepend a column of zeros ("dummyzero") so the Euclidean
#   part is treated as a linear predictor (not stereographic).
# - If intercept = TRUE and no constant covariate exists, append a column of ones.
# Stores onescovaridx (the column index of the constant-1 covariate) in preplist.
addEuccovars <- function(preplist, type, intercept){
  # if LinEuc add a zeros covariate
  if ((type == "LinEuc") && (!is.null(preplist$xe))){
    if (!is.null(preplist$xe)){
      if (any(preplist$xe[,1]^2 > .Machine$double.eps)){
        preplist$xe <- cbind("dummyzero" = 0, preplist$xe)
      }
    }
    if (!is.null(preplist$start)){stopifnot(is_LinEuc(preplist$start))}
  }
 
  # if intercept==TRUE add a ones 
  onescovaridx <- 0
  if (intercept && is.null(preplist$xe)){warning("Intercept==TRUE will be ignored because no Euclidean covariates")}
  if (intercept && !is.null(preplist$xe)){
    # search for 1s covariate, otherwise add it to end
    constxe <- (apply(preplist$xe, 2, sd) < sqrt(.Machine$double.eps))
    onescovaridx <- which(constxe)[which(abs(colMeans(preplist$xe[,constxe, drop = FALSE]) - 1) < .Machine$double.eps)]
    if (length(onescovaridx) == 0){
      preplist$xe <- cbind(preplist$xe, "ones" = 1)
      onescovaridx <- ncol(preplist$xe)
      if (!is.null(preplist$start)){
        preplist$start <- as_mobius_link_cann(preplist$start)
        # if start is wrong dimension add a row of zeros
        if (dim(preplist$start)["qe"] != ncol(preplist$xe)){
          stopifnot(dim(preplist$start)["qe"] == ncol(preplist$xe) - 1)
          preplist$start$Qe <- rbind(preplist$start$Qe, ones = 0)
        }
      }
    }
    if (length(onescovaridx) > 1){onescovaridx <- onescovaridx[1]}
  }
  return(c(
    preplist,
    list(onescovaridx = onescovaridx)
  ))
}

# Generate default starting parameters for the mean link optimisation when no start is supplied.
# Assumes data have been standardised (so identity-like starting values are sensible).
# Sets P = I, Bs = 0.9 I, Be = 0.9 I, Qs and Qe to the first p columns of identity.
defaultstart <- function(preplist, type){
  if (is.null(preplist$start)){
    p <- ncol(preplist$y)
    if (!is.null(preplist$xe)){stopifnot(ncol(preplist$xe) >= p)}
    if (!is.null(preplist$xs)){stopifnot(ncol(preplist$xs) >= p)}
    preplist$start <- mobius_link_cann(
                P = diag(p),
                Bs = if (!is.null(preplist$xs)){diag(0.9, p-1)},
                Qs = if (!is.null(preplist$xs)){diag(1, ncol(preplist$xs), p)},
                Be = if (!is.null(preplist$xe)){diag(0.9, p-1)},
                Qe = if (!is.null(preplist$xe)){diag(1, ncol(preplist$xe), p)},
                ce = if (!is.null(preplist$xe)){1}
    )
    if ((type == "SpEuc") && !is.null(preplist$xe)){ #default to be larger an all values of -xe
      preplist$start$ce[1] <- max(-preplist$xe)  +  0.1*IQR(preplist$xe)
    }
  }
  return(preplist)
}
