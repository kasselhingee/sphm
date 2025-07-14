#' Simple estimators of SvMF Parameters
#' @param y iid observations of a SvMF. Each row an observation.
#' @param mu Mean
#' @param G The orientation paramter matrix Gamma, which has the mean as the first column.
#' @details
#' `SvMF_mom_axes` follows Scealy and Wood (2019) Section 4.1.1
SvMF_mom_axes <- function(y, mu){
  projmat <- diag(length(mu)) - mu %*% t(mu)
  projy <- y %*% t(projmat)
  mom2nd <- t(projy) %*% projy #t() %*% () quickly calculates the sum of projection matrices of rows of projy
  mon2ndeigen <- eigen(mom2nd, symmetric = TRUE)
  Gest <- cbind(mu, mon2ndeigen$vectors[,-length(mu)])
  #Put Gest into somewhat standard from
  Gest <- toBigPosElRot_keepfirst(Gest)
  return(Gest)
}

#' `SvMF_prelim_scales` follows Scealy and Wood (2019) Proposition 2
SvMF_prelim_scales <- function(y, G){
  # rotate y to star coordinates, where the columns of G are the basis coordinates
  # G^T . G = I. So G^T . y does this transformation
  # since y is stored as row vectors
  stary <- y %*% G
  staryL <- stary[,-1]
  # at high concentrations staryL is Gaussian with a diagonal matrix related to the scales a
  # such that the variance is (a_i/(k *a_1))^2
  aremaining <- apply(staryL, 2, sd) #this is actually a_i/(k *a_1)
  # since we know they must multiply to 1 we can avoid knowing k and a_1
  aremaining <- aremaining/prod(aremaining)^(1/length(aremaining))
  return(aremaining)
}