# preliminary estimator using purely R (not automatic differentiation)

optim_pobjS2S_f <- function(theta, y, x){
  prelimobj(y, x, Omega_proj(mnlink_Omega_unvec(theta, ncol(y), check = FALSE)))
}
optim_pobjS2S_g_eq <- function(theta, y, x){
  om <- mnlink_Omega_unvec(theta, ncol(y), check = FALSE)
  mnlink_Omega_check_numerical(om)
}

#' Optimisation of the Preliminary Objectivf Function for S2S Link
#' @details Uses `nloptr`. First a global optimisation with +/-10 of `paramobj0` using algorthim `NLOPT_GN_ISRES` (the only algorthim that natively handles non-linear equality constraints - all the others have to use augmented Lagrangian ideas). The a local optimisation using algorithm `NLOPT_LN_COBYLA`.
#' @param paramobj0 is a starting parameter object.
#' @param global If `TRUE` will do a global search.
#' @param local If `TRUE` will do a local search
optim_pobjS2S_pureR <- function(y, x, paramobj0, global = TRUE, local = TRUE){ #paramobj0 is the starting parameter object
  p <- ncol(y)
  om0 <- as_mnlink_Omega(paramobj0)
  om0_local <- om0
  globopt <- locopt <- NULL
  if (global){
    globopt <- nloptr::nloptr(
      x0 = mnlink_Omega_vec(om0),
      eval_f = optim_pobjS2S_f,
      eval_g_eq = optim_pobjS2S_g_eq,
      lb = mnlink_Omega_vec(om0) *0 - 10, #10 is just a guess here. Since everything is related to spheres, I suspect most values are well below 1.
      ub = mnlink_Omega_vec(om0) *0 + 10,
      opts = list(algorithm = "NLOPT_GN_ISRES",
                  xtol_rel = 1E-04,
                  maxeval = 1E4), #the only algorthim that natively handles non-linear equality constraints - all the others have to use augmented Lagrangian ideas.
      y = y,
      x = x
    )
    om0_local <- mnlink_Omega_unvec(globopt$solution, p, check = FALSE)
  }
 
  if (local){
    # re do with a local optimisation to polish (NLopt docs suggest this)
    locopt <- nloptr::nloptr(
      x0 = mnlink_Omega_vec(om0_local),
      eval_f = optim_pobjS2S_f,
      eval_g_eq = optim_pobjS2S_g_eq,
      lb = mnlink_Omega_vec(om0) *0 - 10, #10 is just a guess here. Since everything is related to spheres, I suspect most values are well below 1.
      ub = mnlink_Omega_vec(om0) *0 + 10,
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
    solution = Omega_proj(mnlink_Omega_unvec(solutionvec, p, check = FALSE)),
    glob_nloptr = globopt,
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
  
  means <- mnlink(xs = x, param = mnlink_cann(P = P, Qs = Q, Bs = B, check = FALSE))
  return(-sum(rowSums(y * means)))
}
