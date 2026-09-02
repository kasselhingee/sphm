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

#' @title Stereographic projection from the south pole
#' @description `Sp()` is the stereographic projection from the south pole
#' \eqn{-e_1 = (-1, 0, \ldots, 0)} of \eqn{S^{p-1}} to \eqn{R^{p-1}}: for
#' \eqn{x = (x_1, \ldots, x_p)} on the sphere,
#' \deqn{Sp(x)_j = x_{j+1} / (1 + x_1), \quad j = 1, \ldots, p-1.}
#' `iSp()` is its inverse.
#'
#' This is the projection used by the scaled Mobius link [`mobius_link()`], which applies
#' `Sp()` to the spherical covariates, accumulates the scaled displacements linearly in
#' \eqn{R^{p-1}}, and maps the result back to the sphere with `iSp()`.
#' @details `Sp()` is undefined at the projection point itself, the south pole
#' \eqn{-e_1}. Rows of `x` equal to \eqn{-e_1} are returned as `1e9` rather than
#' `Inf`, so that downstream optimisation sees a large finite number.
#' @param x A matrix of row vectors on the sphere \eqn{S^{p-1}}, or a single such vector.
#' @param y A matrix of row vectors in \eqn{R^{p-1}}, or a single such vector.
#' @return A matrix with the same number of rows as the input: `p-1` columns for `Sp()`
#' and `p` columns for `iSp()`. If the input is a single vector (or a matrix with one
#' row) the result is returned as a plain vector rather than a one-row matrix.
#' @family link-function
#' @examples
#' x <- rbind(north_pole(3), c(0, 1, 0), c(0, 0, 1))
#' Sp(x)
#' iSp(Sp(x))  # recovers x
#'
#' Sp(north_pole(3))          # the north pole maps to the origin
#' Sp(c(-1, 0, 0))            # the south pole: undefined, returned as 1e9
#' @export
Sp=function(x) {
  if (is.vector(x)){x <- matrix(x, nrow = 1)}
  # detect south-pole vectors (x = -e1 = (-1,0,...,0)), where the projection is undefined
  is_me1 <- colSums(t(x) != c(-1, rep(0, ncol(x) - 1))) == 0
  out <- x[, -1, drop = FALSE]
  out[is_me1, ] <- 1e+9
  out[!is_me1, ] <- out[!is_me1, , drop = FALSE]/(1+x[!is_me1, 1, drop = TRUE])
  if (nrow(out) == 1){return(as.vector(out))}
  else {return(out)}
}

#' @rdname Sp
#' @description `iSp()` maps \eqn{y} in \eqn{R^{p-1}} back to \eqn{S^{p-1}} via
#' \deqn{iSp(y) = (1 - \|y\|^2, 2y) / (1 + \|y\|^2).}
#' It is defined everywhere on \eqn{R^{p-1}} and never returns the south pole.
#' @export
iSp=function(y){
  if (is.vector(y)){y <- matrix(y, nrow = 1)}
  norms2 <- rowSums(y^2)
  # inverse stereographic formula from the south pole: (1-||y||^2, 2y) / (1+||y||^2)
  out <- cbind(1-norms2, 2*y)/(1+norms2)
  if (nrow(out) == 1){return(as.vector(out))}
  else {return(out)}
}

#' @noRd
#' @title Cayley Transform
#' @description Maps the upper-triangle entries (in vector form) of a skew-symmetric
#' matrix A to the orthogonal matrix `(I - A)(I + A)^{-1}`. The Cayley transform is
#' guaranteed to produce an orthogonal matrix for any skew-symmetric A, making it useful
#' for parameterising rotation matrices without explicit orthogonality constraints.
#' The inverse operation is `inverse_cayley_transform()` (exported from C++).
#'
#' No longer called anywhere in the package: the Proposition 4 code was its last caller and
#' now goes through the C++ `cayley_transform()` instead, so that every skew coordinate in
#' the package uses the one strict-lower-triangle convention. Note the two differ by a
#' transpose -- `cayley()` computes \eqn{(I-A)(I+A)^{-1}}, `cayley_transform()` computes
#' \eqn{(I-A)^{-1}(I+A)} -- and fill opposite triangles. Kept for reference and use at the
#' console.
#' @param x numeric vector of length `d*(d-1)/2`: the upper-triangle entries of a
#'   `d x d` skew-symmetric matrix, filled column-by-column.
#' @examples
#' x <- 1:3
cayley <- function(x){
  # length(x) = d*(d-1)/2; recover d from the quadratic formula
  d <- (1 + sqrt(8*length(x) + 1))/2
  p <- matrix(0., d, d)
  p[upper.tri(p, diag = FALSE)] <- x
  p[lower.tri(p, diag = FALSE)] <- -t(p)[lower.tri(p, diag = FALSE)]
  # Cayley transform: (I - A)(I + A)^{-1}, guaranteed orthogonal for skew-symmetric A
  P=(diag(1,d)-p)%*%solve(diag(1,d)+p)
  return(P)
}

#' @title North pole vector
#' @description Returns the north pole unit vector \eqn{(1,0,0,...)^\top}{(1,0,0,...)} of length `p`.
#' This is the default reference direction used for standardising spherical data.
#' @param p The dimension of the space/length of the vector.
#' @return A numeric vector of length `p`.
#' @export
north_pole <- function(p){c(1, rep(0, p-1))}

#' @noRd
#' @description Standardise sign of columns of a matrix to have positive first element, or unchanged sign if 0 first element
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
#' @param mat A matrix
#' @export
standardise_col_signs <- function(mat){
  maxidx <- apply(abs(mat), 2, which.max)
  neg <- mat[cbind(maxidx, 1:length(maxidx))] < -.Machine$double.eps #only negative values found
  mat[,neg] <- -mat[,neg]
  return(mat)
}

#' @title Convert orthogonal matrix to rotation matrix
#' @description Convert from an orthogonal matrix to rotation matrix
#' by switching sign of final column to make determinant positive
#' @param mat an orthogonal matrix
#' @export
as_rotation_mat <- function(mat){
  mat[,ncol(mat)] <- sign(det(mat)) * mat[,ncol(mat)] #make mat a rotation matrix
  return(mat)
}

#' @noRd
#' @description Applies standardise_col_signs() and as_rotation_mat, not doing any sign switchin on the first column
toBigPosElRot_keepfirst <- function(mat){
  mat[,-1] <- standardise_col_signs(mat[,-1]) # make Gamma consistentish signs, and a rotation matrix
  mat <- as_rotation_mat(mat)
  return(mat)
}

# gives the degree of freedom of an object with n rows and p columns (i.e. p orthonormal vectors in n space)
# This formula is from wikipedia and (2.1) of Edelman et al 1998.
DoF_Stiefel <- function(n, p){
  if (n == 0){return(0)}
  if (n < p){return(NA_integer_)} #it is not possible to have more orthogonal vectors (p) than the dimension of the ambient space (n)
  n*p - p * (p+1)/2
}
