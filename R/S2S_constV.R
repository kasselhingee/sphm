#' SvMF Regression with Parallel Transported Axes
#' @description
#' This function is designed to eventually be an internal function.
#' No standardisation is performed.
#' 
#' @details
#' The mean is assumed to follow the usual mean link.
#' The concentration and scaling in the SvMF is assumed constant across observations.
#' The scaling axes of the SvMF at location \eqn{\mu} are assumed to be the parallel transport along the geodesic of axes at the first column of the matrix `P` from the mean link. These axes specified at first column of the matrix `P` are to be estimated and constant with respect to covariates (and \eqn{\mu})
#' __Warning: C++ function still uses Jupp's transport rather than Amaral matrix and p!=3 only handled approximately__
#' __Warning: resulting relevant matrices are only very close to orthogonal - something do with with p1 not being precisely orthogonal to singular vectors despite the projection of Omega (I could do a projection again!?)__
#' __Warning: estimation could be improved by using moment estimators for Gstar and a etc if possible to start the optimisation__
#' 
#' The first element of each column of Gstar will have positive value.
#' @param y Response data on a sphere
#' @param xs Covariate data on a sphere
#' @param xe Covariate data in Euclidean Space
#' @param k Starting concentration. I suspect lower means less chance of finding a local minimum.
#' @param a The scaling vector `a`. `a[1]` is a fixed tuning parameter and the remainining is used as a starting guess.
#' @param mean Parameters for the mean link, used as a starting guess.
#' @param G0 A `p x p` orthonormal matrix specifying the starting guess of the axes of the SvMF distribution. G0 should have positive determinant because in the estimatino routine G0 or parts of G0 are representented using Cayley transforms.
#' @param G0reference A `p x p` rotation matrix specifying a set of coordinates to represent G0 in for the estimation. Ideally the columns of `G0reference` will be close to the best `G0` because the Cayley transform representation has the best performance when applied to matrices close to the identity.
#' @param G01behaviour "p1" identifies `G0[,1]` with `p1`. "fixed" fixes `G0[,1]` to its initial value. "free" allows `G0[,1]` to be estimated freely.
#' @param ... Named optional arguments passed as a list to the `opts` argument of [`nloptr::nloptr()`].
#' @details
#' From the starting parameters, optimises everything. For p != 3, the concentration is approximated.
#' No standardisation is performed.
#' @export
optim_constV <- function(y, xs, xe, mean, k, a, G0, G0reference = diag(p), G01behaviour = "p1", ...){
  p <- ncol(y)
  qs <- ncol(xs)
  qe <- ncol(xe)
  # checks
  om0 <- as_mnlink_Omega(mean)
  mnlink_Omega_check(om0)
  stopifnot(length(a) == p)
  a1 = a[1]
  aremaining = a[-1]
  stopifnot(isTRUE(all.equal(prod(aremaining), 1)))
  
  initial <-  list(
    mean = om0,
    k = k, 
    a = a, 
    G0 = G0
  )
  
  # estimation prep
  dims_in <- c(p, length(om0$qe1))
  omvec0 <- mnlink_Omega_vec(om0)
  ulltape <- tape_ull_S2S_constV_nota1(omvec = omvec0, k = k,
                                       a1 = a1, 
                                       aremaining = aremaining,
                                       G0 = G0,
                                       p = p, qe = qe, 
                                       yx = cbind(y, xs, xe),
                                       referencecoords = G0reference,
                                       G01behaviour = G01behaviour)
  ulltape <- scorematchingad::avgrange(ulltape)
  constraint_tape <- tape_namedfun("Omega_constraints_wrap", omvec0, vector(mode = "numeric"), dims_in, matrix(nrow = 0, ncol = 0), check_for_nan = FALSE)
  ineqconstraint_tape <- tape_namedfun("Omega_ineqconstraints", omvec0, vector(mode = "numeric"), dims_in, matrix(nrow = 0, ncol = 0), check_for_nan = FALSE)
  x0 <- ulltape$xtape
  
  # check Jacobians of constraints when constraints satisfied
  Jac_eq <- matrix(constraint_tape$Jacobian(omvec0), byrow = TRUE, ncol = length(omvec0))
  stopifnot(all(abs(svd(Jac_eq)$d) > sqrt(.Machine$double.eps))) 
  Jac_ineq <- matrix(ineqconstraint_tape$Jacobian(omvec0), byrow = TRUE, ncol = length(omvec0))
  stopifnot(all(abs(svd(Jac_ineq)$d) > sqrt(.Machine$double.eps)))

  # prepare nloptr options
  default_opts <- list(algorithm = "NLOPT_LD_SLSQP",
                       xtol_rel = 1E-10, #1E-04,
                       tol_constraints_eq = rep(1E-1, constraint_tape$range),
                       # check_derivatives = TRUE, check_derivatives_print = 'errors', check_derivatives_tol = 1E-3,
                       # print_level = 3,
                       maxeval = 1E4)
  ellipsis_args <- list(...)
  combined_opts <- utils::modifyList(default_opts, ellipsis_args)
  
  # lower bound for concentration k
  lb <- rep(-Inf, ulltape$domain)
  lb[length(omvec0) + 1] <- 0
  
  # Optimisation
  est <- nloptr::nloptr(
    x0 = x0,
    eval_f = function(theta){-ulltape$eval(theta, a1)},
    eval_grad_f = function(theta){-ulltape$Jac(theta, a1)},
    lb = lb,
    eval_g_eq =  function(theta){constraint_tape$forward(0, theta[1:length(omvec0)])},
    eval_jac_g_eq =  function(theta){
      out <- cbind(matrix(constraint_tape$Jac(theta[1:length(omvec0)], vector(mode = "numeric")), byrow = TRUE, ncol = length(omvec0)),
            matrix(0, nrow = constraint_tape$range, ncol = length(theta) - length(omvec0)))
      if (any(is.nan(out))){browser()}
      return(out)
    },
    opts = combined_opts
  )
  
  if (!(est$status %in% c(0, 1, 2, 3, 4))){warning("Optimistation did not finish properly.")}
  estparamlist <- S2S_constV_nota1_fromvecparamsR(est$solution, p, qs, qe, 
                                                  referencecoords = G0reference,
                                                  G01behaviour = G01behaviour,
                                                  G01 = initial$G0[,1])
  
  # project Omega to satisfy orthogonality constraint
  est_om <- Omega_proj(mnlink_Omega_unvec(estparamlist$omvec, p, qe = qe, check = FALSE))
  
  #make first element of each vector positive
  estparamlist$G0[,-1] <- topos1strow(estparamlist$G0[,-1])
  
  # make sure aremaining is in decreasing order
  aord <- order(estparamlist$aremaining, decreasing = TRUE)
  estparamlist$aremaining <- estparamlist$aremaining[aord]
  estparamlist$G0[,-1] <- estparamlist$G0[,-1][, aord]
  
  
  outsolution <- list(
    mean = est_om,
    k = estparamlist$k,
    a = c(a1, estparamlist$aremaining),
    G0 = estparamlist$G0
  )
  
  return(list(
    solution = outsolution,
    nlopt = est,
    initial = initial
  ))
}

#' Function for simulating data given mean link and SvMF parameters
rS2S_constV <- function(xs, xe, mnparam, k, a, G0){
  ymean <- mnlink(xs = xs, xe = xe, param = mnparam)
  
  # simulate noise
  y_ld <- t(apply(ymean, 1, function(mn){
    G <- cbind(mn, -JuppRmat(G0[,1], mn) %*% G0[,-1])
    obs <- rSvMF(1, SvMFcann(k, a, G))
    ld <- uldSvMF_cann(obs, k = k, a = a, G = G)
    return(c(obs, ld))
  }))
  return(y_ld)
}
