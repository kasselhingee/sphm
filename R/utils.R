#' @useDynLib sphm, .registration = TRUE
#' @importFrom Rcpp evalCpp

#' @title Euclidean norm
#' @description Returns the Euclidean norm (square root of the sum of squares) of a vector. `vnorm2()` returns the square of the Euclidean norm.
#' @param x a vector
#' @return A single numeric value.
#' @export
vnorm=function(x) sqrt(vnorm2(x))

#' @rdname vnorm
#' @export
vnorm2=function(x) sum(x^2)

#' @noRd
#' @title Stereographic projection
#' @param x is a matrix of row vectors
#' @noRd
Sp=function(x) {
  if (is.vector(x)){x <- matrix(x, nrow = 1)}
  # detect -e1 vectors, remembering the x may be in a disc
  is_me1 <- colSums(t(x) != c(1, rep(0, ncol(x) - 1))) == 0
  out <- x[, -1, drop = FALSE]
  out[is_me1, ] <- 1e+9
  out[!is_me1, ] <- out[!is_me1, , drop = FALSE]/(1+x[!is_me1, 1, drop = TRUE])
  if (nrow(out) == 1){return(as.vector(out))}
  else {return(out)}
}

#' @noRd
#' @title Inverse stereographic projection
#' @param y is a matrix of row vectors
iSp=function(y){
  if (is.vector(y)){y <- matrix(y, nrow = 1)}
  norms2 <- rowSums(y^2)
  out <- cbind(1-norms2, 2*y)/(1+norms2)
  if (nrow(out) == 1){return(as.vector(out))}
  else {return(out)}
}

#' @noRd
#' @title Cayley Transform
#' @description Map the upper triangle of a skew-symmetric matrix to an
#' orthogonal matrix with (non-negative) determinant.
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

#' @title 'North pole' vector
#' @description The vector \eqn{(1,0,0,...)^\top}{(1,0,0,...)} for given dimension.
#' @param p The dimension of the space/length of the vector.
#' @return A vector of length `p`.
#' @export
nthpole <- function(p){c(1, rep(0, p-1))}

#' @noRd
#' Standardise sign of columns of a matrix to have positive first element, or unchanged sign if 0 first element
topos1strow <- function(mat){
  neg <- mat[1, ] < -.Machine$double.eps #only negative values found
  sgn <- rep(1, length(neg))
  sgn[neg] <- -1
  mat <- t(t(mat) * sgn)
  return(mat)
}

#' @title Standardise signs of columns
#' @description
#' Standardise signs of columns so that largest (absolute value) element is positive in each column, or unchanged if all elements are 0.
#' This should not change the signs of diagonal elements in diagonal matrices.
#' @param mat
#' @export
toBigPosEl <- function(mat){
  maxidx <- apply(abs(mat), 2, which.max)
  neg <- mat[cbind(maxidx, 1:length(maxidx))] < -.Machine$double.eps #only negative values found
  mat[,neg] <- -mat[,neg]
  return(mat)
}

#' @title Convert orthogonal matrix to rotation matrix
#' @description Convert from an orthogonal matrix to rotation matrix
#' by switching sign of final column to make determinant positive
#' @param mat
#' @export
toRot <- function(mat){
  mat[,ncol(mat)] <- sign(det(mat)) * mat[,ncol(mat)] #make mat a rotation matrix
  return(mat)
}

#' @noRd
#' @description Applies toBigPosEl() and toRot, not doing any sign switchin on the first column
toBigPosElRot_keepfirst <- function(mat){
  mat[,-1] <- toBigPosEl(mat[,-1]) # make Gamma consistentish signs, and a rotation matrix
  mat <- toRot(mat)
  return(mat)
}

# gives the degree of freedom of an object with n rows and p columns (i.e. p orthonormal vectors in n space)
# This formula is from wikipedia and (2.1) of Edelman et al 1998.
DoF_Stiefel <- function(n, p){
  if (n == 0){return(0)}
  if (n < p){return(NA_integer_)} #it is not possible to have more orthogonal vectors (p) than the dimension of the ambient space (n)
  n*p - p * (p+1)/2
}
