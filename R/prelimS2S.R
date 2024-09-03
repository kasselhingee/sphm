#' Preliminary objective function for S2S Link
#' @details Assumes that the distribution is isotropic around the mean with constant concentration, thus maximising
#' \deqn{\sum_i=1^n y_i^T \mu(x_i).}
#' @param y A set of unit vectors in embedded coordinates, each row corresponds to a single unit vector.
#' @param x A set of covariate vectors (also unit vectors), each row corresponds to the same row in `y`.
#' @param paramobj A set of link parameters. See [`cannS2S()`] and [`OmegaS2S()`].
#' @export
pobjS2S <- function(y, x, paramobj){
  predictedmeans <- meanlinkS2S(x = x, paramobj = paramobj, check = FALSE)
  stopifnot(nrow(y) == nrow(predictedmeans))
  stopifnot(ncol(y) == ncol(predictedmeans))
  return(-mean(rowSums(y * predictedmeans)))
}


optim_pobjS2S_f <- function(theta, y, x){
  pobjS2S(y, x, OmegaS2S_proj(OmegaS2S_unvec(theta, ncol(y), check = FALSE), method = "Omega"))
}
optim_pobjS2S_g_eq <- function(theta, y, x){
  om <- OmegaS2S_unvec(theta, ncol(y), check = FALSE)
  OmegaS2S_check_internal(om)
}

#' Optimisation of the Preliminary Objectivf Function for S2S Link
#' @details Uses `nloptr`. First a global optimisation with +/-10 of `paramobj0` using algorthim `NLOPT_GN_ISRES` (the only algorthim that natively handles non-linear equality constraints - all the others have to use augmented Lagrangian ideas). The a local optimisation using algorithm `NLOPT_LN_COBYLA`.
#' @param paramobj0 is a starting parameter object.
#' @param global If `TRUE` will do a global search.
#' @param local If `TRUE` will do a local search
optim_pobjS2S_pureR <- function(y, x, paramobj0, global = TRUE, local = TRUE){ #paramobj0 is the starting parameter object
  p <- ncol(y)
  om0 <- as_OmegaS2S(paramobj0)
  om0_local <- om0
  globopt <- locopt <- NULL
  if (global){
    globopt <- nloptr::nloptr(
      x0 = OmegaS2S_vec(om0),
      eval_f = optim_pobjS2S_f,
      eval_g_eq = optim_pobjS2S_g_eq,
      lb = OmegaS2S_vec(om0) *0 - 10, #10 is just a guess here. Since everything is related to spheres, I suspect most values are well below 1.
      ub = OmegaS2S_vec(om0) *0 + 10,
      opts = list(algorithm = "NLOPT_GN_ISRES",
                  xtol_rel = 1E-04,
                  maxeval = 1E4), #the only algorthim that natively handles non-linear equality constraints - all the others have to use augmented Lagrangian ideas.
      y = y,
      x = x
    )
    om0_local <- OmegaS2S_unvec(globopt$solution, p, check = FALSE)
  }
 
  if (local){
    # re do with a local optimisation to polish (NLopt docs suggest this)
    locopt <- nloptr::nloptr(
      x0 = OmegaS2S_vec(om0_local),
      eval_f = optim_pobjS2S_f,
      eval_g_eq = optim_pobjS2S_g_eq,
      lb = OmegaS2S_vec(om0) *0 - 10, #10 is just a guess here. Since everything is related to spheres, I suspect most values are well below 1.
      ub = OmegaS2S_vec(om0) *0 + 10,
      opts = list(algorithm = "NLOPT_LN_COBYLA",
                  xtol_rel = 1E-04,
                  maxeval = 1E3),
      y = y,
      x = x
    )
  } 
  
  solutionvec <- switch(1 + global + 2*local,
                        stop("At least one of global or local must be TRUE"),
                        globopt$solution,
                        locopt$solution,
                        locopt$solution
                        )
  
  return(list(
    solution = OmegaS2S_proj(OmegaS2S_unvec(solutionvec, p, check = FALSE), method = "Omega"),
    glob_nloptr = globopt,
    loc_nloptr = locopt
  ))
}

#' Optimisation of the Preliminary Objectivf Function for S2S Link
#' @details Uses `nloptr`. `NLopt` doesn't have any algorithms for global optimisation with non-linear equality constraints that use provided gradients. So `_parttape` only does local optimisation and uses `NLOPT_LD_SLSQP` which is the only algorithm that takes advantage of derivatives and can handle non-linear equality constraints.
#' @param paramobj0 is a starting parameter object.
#' @param ... Passed as options to [`nloptr()`]. Default is
optim_pobjS2S_parttape <- function(y, x, paramobj0, ...){ #paramobj0 is the starting parameter object
  p <- ncol(y)
  om0 <- as_OmegaS2S(paramobj0)
  
  obj_tape <- tape_namedfun("pobjS2Scpp", OmegaS2S_vec(om0), vector(mode = "numeric"), p, cbind(y,x))
  constraint_tape <- tape_namedfun("wrap_OmegaS2S_constraints", OmegaS2S_vec(om0), vector(mode = "numeric"), p, matrix(nrow = 0, ncol = 0))

  # prepare nloptr options
  default_opts <- list(algorithm = "NLOPT_LD_SLSQP",
                xtol_rel = 1E-10, #1E-04,
                tol_constraints_eq = rep(1E-1, 2),
                maxeval = 1E4)
  ellipsis_args <- list(...)
  combined_opts <- utils::modifyList(default_opts, ellipsis_args)
  
  locopt <- nloptr::nloptr(
    x0 = OmegaS2S_vec(om0),
    eval_f = function(theta){scorematchingad:::pForward0(obj_tape, theta, vector(mode = "numeric"))},
    eval_grad_f = function(theta){scorematchingad:::pJacobian(obj_tape, theta, vector(mode = "numeric"))},
    eval_g_eq =  function(theta){scorematchingad:::pForward0(constraint_tape, theta, vector(mode = "numeric"))[1:2]},
    eval_jac_g_eq =  function(theta){matrix(scorematchingad:::pJacobian(constraint_tape, theta, vector(mode = "numeric")), byrow = TRUE, ncol = length(theta))},
    opts = combined_opts
  )
  
  return(list(
    solution = OmegaS2S_proj(OmegaS2S_unvec(locopt$solution, p, check = FALSE), method = "Omega"),
    loc_nloptr = locopt
  ))
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


