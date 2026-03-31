#' Jupp's Rotated Residuals

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

cruderesid <- function(y, ypred){
  out <- t(sapply(1:nrow(y), function(i){
    projmat <- ypred[i,] %*% t(ypred[i, ])
    y[i, ] - projmat %*% y[i, ]
  }))
  samehemisphere <- rowSums(y * ypred) >= 0
  attr(out, "samehemisphere") <- samehemisphere
  return(out)
}

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
#' of Jupp's (1988) rotation matrix [`JuppRmat()`]. Both perform the same parallel transport;
#' this version uses a rotation-based formula.
#' @param start Starting location as a unit vector.
#' @param end End location as a unit vector.
#' @return A square orthogonal matrix.
#' @details If `end = -start` (antipodal points), there is no unique geodesic and a reflection
#' matrix in the `end` direction is returned.
#' @seealso [rotatedresid()], [JuppRmat()]
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

# Parallel transport matrix from Absil, Mahony and Sepulchre equation (8.4)
partransportmat <- function(start, end){
  alpha <- drop(acos(start %*% end))
  u <- (end - cos(alpha) * start)/sin(alpha)
  out <- diag(length(start)) - sin(alpha) * start %*% t(u) + (cos(alpha) - 1) * u %*% t(u)
}

resid_SvMF_partransport <- function(y, ypred, k = NULL, a = NULL, G0, scale = TRUE){
  rresids_std <- rresids_tmp <- rotatedresid(y, ypred, G0[,1])
  rresids_std <- rresids_std %*% G0
  rresids_std <- rresids_std[, -1]
  if (scale){rresids_std <- sqrt(k) * rresids_std %*% diag(a[1]/a[-1])}
  attr(rresids_std, "samehemisphere") <-  attr(rresids_tmp, "samehemisphere")
  colnames(rresids_std) <- paste0("r", 1:ncol(rresids_std))
  return(rresids_std)
}
