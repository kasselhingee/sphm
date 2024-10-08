vnorm2=function(x) sum(x^2)
vnorm=function(x) sqrt(vnorm2(x))


#' Cayley Transform Reparameterisation
#' @description Map the upper triangle of a skew-symmetric matrix to an
#' orthogonal matrix with (non-negative?) determinant.
#' @examples
#' x <- 1:3
cayley <- function(x){
  # length(x) <- d * (d-1)/2
  d <- (1 + sqrt(8*length(x) + 1))/2
  p <- matrix(0., d, d)
  p[upper.tri(p, diag = FALSE)] <- x
  p[lower.tri(p, diag = FALSE)] <- -t(p)[lower.tri(p, diag = FALSE)]
  P=(diag(1,d)-p)%*%solve(diag(1,d)+p)
  return(P)
}

#' Standardise Data to North Pole
#' @description Rotation according to the eigenvectors of the second moment of the data
#' projected perpendicular to the mean direction. Standardised data have mean of c(1, 0, 0,...). See Scealy and Wood 2019 for details.
#' @param y Data on the sphere. Each row is a data point in Cartesian coordinates.
#' @param G Axes of the second moment matrix of the data, projected so that the first column is the global mean.
#' @details Each returned data point is `t(G)` of the original data point, where `G` is computed by `standardise_mat()`.
#' @export
standardise <- function(y, G = standardise_mat(y)){
  ystd <- y %*% G
  ystd <- unname(ystd)
  return(ystd)
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

  Ghat <- cbind(mn, eigen(projmom2)$vectors[, 1:(p-1)])
  return(Ghat)
}

nthpole <- function(p){c(1, rep(0, p-1))}
