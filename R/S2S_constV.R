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
optim_constV <- function(y, xs, xe, mean, k, a, G0, G0reference = diag(p), G01behaviour = "p1", fix_qs1 = FALSE, fix_qe1 = FALSE, ssqOmbuffer = 2, ...){
  om0 <- as_mnlink_Omega(mean)
  # check inputs:
  try(mnlink_Omega_check(om0))
  p <- ncol(y)
  stopifnot(p == length(om0$p1))
  if (!is.null(xs)){
    stopifnot(ncol(xs) == length(om0$qs1))
  } else {
    stopifnot(length(om0$qs1) == 0)
  }
  # check Euc info if Euc is part of the link
  if (!is.null(xe)){
    stopifnot(ncol(xe) == length(om0$qe1))
  } else {
    stopifnot(length(om0$qe1) == 0)
  }
  stopifnot(length(a) == p)
  a1 = a[1]
  aremaining = a[-1]
  stopifnot(isTRUE(all.equal(prod(aremaining), 1)))

  qs <- length(om0$qs1)
  qe <- length(om0$qe1)

  
  initial <-  list(
    mean = om0,
    k = k, 
    a = a, 
    G0 = G0
  )
  
  # Prepare objective tapes
  dims_in <- c(p, length(om0$qe1))
  om0vec <- mnlink_Omega_vec(om0)
  objtape <- tape_ull_S2S_constV_nota1(omvec = om0vec, k = k,
                                       a1 = a1, 
                                       aremaining = aremaining,
                                       G0 = G0,
                                       p = p, qe = qe, 
                                       yx = cbind(y, xs, xe),
                                       referencecoords = G0reference,
                                       G01behaviour = G01behaviour)
  objtape <- scorematchingad::avgrange(objtape) #objtape initially returns a value for each measurement. Average here to get average over all data.
  constraint_tape <- tape_namedfun("Omega_constraints_wrap", om0vec, vector(mode = "numeric"), dims_in, matrix(nrow = 0, ncol = 0), check_for_nan = FALSE)
  ineqconstraint_tape <- tape_namedfun("Omega_ineqconstraints", om0vec, vector(mode = "numeric"), dims_in, matrix(nrow = 0, ncol = 0), check_for_nan = FALSE)

  # fix mean link parameters depending on arguments
  # use the starting parameters om0 to detect whether we have xs and xe as their form is more predictable due to the Omega class
  omfixed <- lapply(om0, function(x) x * 0)
  if (fix_qe1 && (length(om0$qe1) > 0)){
    # if shogo and Euc, fix some elements
    omfixed$qe1 <- omfixed$qe1 + 1
    omfixed$ce <- omfixed$ce + 1
  }
  
  # fix qs1 to the starting values if qs1 exists
  # use the starting parameters to check as their form is more constrained by the Omega class
  if (fix_qs1 && (length(om0$qs1) > 0)){
    omfixed$qs1 <- omfixed$qs1 + 1
  }
  
  # update constraint tapes based on omfixed
  isfixed <- mnlink_Omega_vec(as_mnlink_Omega(omfixed)) > 0.5
  constraint_tape <- scorematchingad::fixindependent(constraint_tape, om0vec, isfixed)
  ineqconstraint_tape <- scorematchingad::fixindependent(ineqconstraint_tape, om0vec, isfixed)
  # drop constraint returns that are constant:
  keep <- which(vapply((1:constraint_tape$range)-1, function(i){!constraint_tape$parameter(i)}, FUN.VALUE = FALSE))
  constraint_tape <- scorematchingad::keeprange(constraint_tape, keep)
  
  # update objtape based on fixed values
  objtape <- scorematchingad::fixindependent(objtape, objtape$xtape, c(isfixed, rep(0, length(objtape$xtape) - length(isfixed))))
  # using tapeing values as starting parameters
  x0 <- objtape$xtape

  # check Jacobians of constraints are non-singular for the starting parameters.
  # For pathological params (e.g. the default starting params of no rotations), it can be zero.
  # If it is singular, perturb start very slightly
  Jac_eq <- matrix(constraint_tape$Jacobian(x0[1:constraint_tape$domain]), byrow = TRUE, ncol = constraint_tape$domain)
  if (any(abs(svd(Jac_eq)$d) < sqrt(.Machine$double.eps))){x0[1:constraint_tape$domain] <- x0[1:constraint_tape$domain] - 1E-4}
  Jac_eq <- matrix(constraint_tape$Jacobian(x0[1:constraint_tape$domain]), byrow = TRUE, ncol = constraint_tape$domain)
  stopifnot(all(abs(svd(Jac_eq)$d) > sqrt(.Machine$double.eps))) 

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
  lb <- rep(-Inf, objtape$domain)
  lb[length(om0vec) + 1] <- 0
  
  # Optimisation
  # current dynamic parameter values of tapes will be used
  locopt <- nloptr::nloptr(
    x0 = x0,
    eval_f = function(theta){-objtape$forward(0, theta)},
    eval_grad_f = function(theta){-objtape$Jacobian(theta)},
    lb = lb,
    eval_g_eq =  function(theta){constraint_tape$forward(0, theta[1:length(om0vec)])},
    eval_jac_g_eq =  function(theta){
      Jac <- cbind(matrix(constraint_tape$Jacobian(theta[1:length(om0vec)]), byrow = TRUE, ncol = length(om0vec)),
             matrix(0, nrow = constraint_tape$range, ncol = length(theta) - length(om0vec)))
      # colnames(Jac) <- names(om0vec)
      # print(round(Jac, 3))
      # print(apply(Jac, 1, function(x)max(abs(x))))
      Jac
    },
    opts = combined_opts
  )
  if (!(locopt$status %in% 1:4)){warning(locopt$message)}

  #output some diagnostics - vector names would be nice here
  locopt$solution_grad_f <- -objtape$Jacobian(locopt$solution)
  locopt$solution_jac_g_eq <- matrix(constraint_tape$Jacobian(locopt$solution),
                                     byrow = TRUE, ncol = length(locopt$solution))
  locopt$solution_Hes_f <- matrix(-objtape$Hessian0(locopt$solution),
         nrow = objtape$domain,
         byrow = TRUE)

  # remove the tapes from the return to save on memory
  locopt$eval_f <- locopt$eval_g_eq <- locopt$eval_g_ineq <- locopt$nloptr_environment <- NULL

  # insert any fixed values of mean parameters
  meanpars <- locopt$solution[1:(sum(!is.na(isfixed)))]
  meanpars <- scorematchingad:::t_sfi2u(locopt$solution, om0vec, isfixed)
  fullparam <- c(meanpars, locopt$solution[-(1:(sum(!is.na(isfixed))))])

  
  estparamlist <- S2S_constV_nota1_fromvecparamsR(fullparam, p, qs, qe, 
                                                  referencecoords = G0reference,
                                                  G01behaviour = G01behaviour,
                                                  G01 = initial$G0[,1])
  
  #project mean pars to have correct orthogonality
  projectedom <- Omega_proj(mnlink_Omega_unvec(estparamlist$omvec, p, length(om0$qe1), check = FALSE))
  try({mnlink_Omega_check(projectedom)})

  # make sure aremaining is in decreasing order
  aord <- order(estparamlist$aremaining, decreasing = TRUE)
  aremaining <- estparamlist$aremaining[aord]
  estparamlist$G0[,-1] <- estparamlist$G0[,-1][, aord]

  #For axes G0 standardise the return by
  # (1) make first element of each vector positive (except the first column)
  G0 <- estparamlist$G0
  G0[,-1] <- topos1strow(G0[,-1])
  # (2) make rotation matrix by flipping final column according to determinant
  if (det(G0) < 0){G0[,p] <- -G0[,p]}
  
  outsolution <- list(
    mean = est_om,
    k = estparamlist$k,
    a = c(a1, aremaining),
    G0 = G0
  )
  

  return(list(
    solution = outsolution,
    nlopt = locopt,
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
