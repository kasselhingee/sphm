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

nthpole <- function(p){c(1, rep(0, p-1))}

#' Standardise sign of columns of a matrix to have positive first element
topos1strow <- function(mat){
  mat <- t(t(mat) * sign(mat[1, ]))
  return(mat)
}

# gives the degree of freedom of an object with n rows and p columns (i.e. p orthonormal vectors in n space)
# This formula is from wikipedia and (2.1) of Edelman et al 1998.
DoF_Stiefel <- function(n, p){
  if (n == 0){return(0)}
  if (n < p){return(NA_integer_)} #it is not possible to have more orthogonal vectors (p) than the dimension of the ambient space (n)
  n*p - p * (p+1)/2
}