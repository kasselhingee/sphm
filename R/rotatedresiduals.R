#' @title Jupp's Rotated Residuals
#' @description Compute rotated (parallel-transported) residuals on the sphere.
#' For each observation, the crude residual (projection of `y` onto the tangent
#' plane at `ypred`) is parallel transported along the geodesic to a common
#' `base` location, making residuals comparable across observations.
#' @param y Matrix of observed unit vectors (one per row).
#' @param ypred Matrix of predicted unit vectors (one per row, same size as `y`).
#' @param base Unit vector giving the base location to transport residuals to.
#'   Typically `north_pole(p)` (see [`north_pole()`]).
#' @param path Transport method: `"geo"` or `"Amaral"` use [`parallel_transport_mat()`];
#'   `"Jupp"` uses `JuppRmat()` (internal; gives the negative); `"Absil"` uses a third formula.
#' @return A matrix of rotated residuals (one per row), with attribute
#'   `"samehemisphere"` indicating observations on the same side as `ypred`.
#' @seealso [parallel_transport_mat()]
#' @family residuals
#' @export
rotatedresid <- function(y, ypred, base, path = "geo"){
  #crude residuals
  cresids <- cruderesid(y, ypred)

  # rotate them
  transportmat <- switch(path, Jupp = JuppRmat, geo = parallel_transport_mat, Amaral = parallel_transport_mat, Absil = partransportmat)
  rresids <- t(sapply(1:nrow(y), function(i){
    tmat <- transportmat(ypred[i, ], base)
    if (!all(is.finite(tmat))){warning("A predicted mean is at the antipode of the rotated residual base location")}
    tmat %*%  cresids[i, ]
  }))
  attr(rresids, "samehemisphere") <- attr(cresids, "samehemisphere")
  return(rresids)
}

# Compute crude residuals: for each row, project y onto the tangent plane at ypred.
# Returns a matrix with attribute "samehemisphere" (TRUE if y and ypred are in the same hemisphere).
cruderesid <- function(y, ypred){
  out <- t(sapply(1:nrow(y), function(i){
    projmat <- ypred[i,] %*% t(ypred[i, ])
    y[i, ] - projmat %*% y[i, ]
  }))
  samehemisphere <- rowSums(y * ypred) >= 0
  attr(out, "samehemisphere") <- samehemisphere
  return(out)
}

# Jupp (1988) rotation matrix for parallel transport from `y` to `base`.
# Equivalent to the negative of parallel_transport_mat() — i.e. JuppRmat(a,b) = -parallel_transport_mat(a,b).
# Used internally by rotatedresid() when path = "Jupp".
JuppRmat <- function(y, base){
  (y+base) %*% t(y+base)/(1+drop(y%*%base)) - diag(1, length(y))
}

#' @title Parallel transport matrix
#' @description Rotation matrix that performs parallel transport along the minimum-distance geodesic
#' on the sphere from `start` to `end`. When applied to a tangent vector at `start`, it moves
#' that vector to the corresponding tangent vector at `end`. When applied to `start` itself,
#' it returns `end`.
#'
#' This is the matrix from Amaral et al. (2007, Lemma 2), which is equivalent to the negative
#' of Jupp's (1988) rotation matrix (internal function `JuppRmat()`). Both perform the same
#' parallel transport; this version uses a rotation-based formula.
#' @param start Starting location as a unit vector.
#' @param end End location as a unit vector.
#' @return A square orthogonal matrix.
#' @details If `end = -start` (antipodal points), there is no unique geodesic and a reflection
#' matrix in the `end` direction is returned.
#' @seealso [rotatedresid()]
#' @family residuals
#' @export
parallel_transport_mat <- function(start, end){ #assumes a and b are unit vectors
  ab <- (end %*% start)[[1]]
  if (ab > 1){ab <- min(ab, 1)}
  if (ab < -1){ab <- max(ab, -1)}
  alpha <- acos(ab)
  cvec <- start - end*ab
  if ((cvec %*% cvec)[[1]] > 0){
    cvec <- cvec/sqrt(cvec %*% cvec)[[1]] #if statement here is useful when start==end
  }
  A <- end%o%cvec - cvec%o%end
  Q = diag(length(end)) + sin(alpha)*A + (cos(alpha) - 1)*(end%o%end + cvec%o%cvec)
  return(Q)
}

# Parallel transport matrix from Absil, Mahony and Sepulchre (2008) equation (8.4).
# Alternative formula to parallel_transport_mat() and JuppRmat(); used when path = "Absil".
partransportmat <- function(start, end){
  alpha <- drop(acos(start %*% end))
  u <- (end - cos(alpha) * start)/sin(alpha)
  out <- diag(length(start)) - sin(alpha) * start %*% t(u) + (cos(alpha) - 1) * u %*% t(u)
}

# Compute SvMF-standardised rotated residuals. Parallel transports crude residuals
# to the base location G0[,1], then projects onto G0 axes and optionally scales
# by the SvMF concentration k and scale vector a to give approximately unit-variance residuals.
resid_SvMF_partransport <- function(y, ypred, k = NULL, a = NULL, G0, scale = TRUE){
  rresids_std <- rresids_tmp <- rotatedresid(y, ypred, G0[,1])
  rresids_std <- rresids_std %*% G0
  rresids_std <- rresids_std[, -1]
  if (scale){rresids_std <- sqrt(k) * rresids_std %*% diag(a[1]/a[-1])}
  attr(rresids_std, "samehemisphere") <-  attr(rresids_tmp, "samehemisphere")
  colnames(rresids_std) <- paste0("r", 1:ncol(rresids_std))
  return(rresids_std)
}
