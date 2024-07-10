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

