#' Preliminary objective function for S2S Link
#' @details Assumes that the distribution is isotropic around the mean with constant concentration, thus maximising
#' \deqn{\sum_i=1^n y_i^T \mu(x_i).}
#' @param y A set of unit vectors in embedded coordinates, each row corresponds to a single unit vector.
#' @param x A set of covariate vectors (also unit vectors), each row corresponds to the same row in `y`.
#' @param paramobj A set of link parameters. See [`cannS2S()`] and [`OmegaS2S()`].
#' @export
pobjS2S <- function(y, x, paramobj){
  predictedmeans <- meanlinkS2S(x = x, paramobj = paramobj)
  stopifnot(nrow(y) == nrow(predictedmeans))
  stopifnot(ncol(y) == ncol(predictedmeans))
  return(-sum(rowSums(y * predictedmeans)))
}

optim_pobjS2S <- function(y, x, paramobj0){ #paramobj0 is the starting parameter object
  p <- ncol(y)
  om0 <- as_OmegaS2S(paramobj0)
  globopt <- nloptr(
    x0 = OmegaS2S_vec(om0),
    eval_f = function(theta, y, x){
      pobjS2S(y, x, OmegaS2S_unvec(theta, ncol(y)))
    },
    eval_g_eq = function(theta, y, x){
      om <- OmegaS2S_unvec(theta, p)
      sum(OmegaS2S_check_internal(om))
    },
    lb = om0 *0 - 10, #10 is just a guess here. Since everything is related to spheres, I suspect most values are well below 1.
    ub = om0 *0 + 10,
    opts = list(algorithm = "NLOPT_GN_ISRES"), #the only algorthim that natively handles non-linear equality constraints - all the others have to use augmented Lagrangian ideas.
    y,
    x
  )
  
  # re do with a local optimisation to polish (NLopt docs)
}

#' Preliminary objective function for S2S Link with p=q
#' @details Uses Cayley transform to parameterise P and Q. 
#' + Could be more accurate the closer P and Q are to the identity (check notes with Andy).
#' + Might be missing sign stuff to get negative determinants for Cayley transform
#' @export
pre_est3_mod=function(y,x,theta){
  
  b1=theta[7]
  b2=theta[8]
  
  P=cayley(theta[1:3])
  Q=cayley(theta[4:6])
  B=b1*diag(c(1,b2))
  
  means <- meanlinkS2S(t(x), P = P, Q = Q, B = B, check = FALSE)
  return(-sum(rowSums(t(y) * means)))
}


