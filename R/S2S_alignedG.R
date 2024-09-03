#' MLE for the SvMF with S2S Mobius Link and alignedG_

#' @param y Response data on a sphere
#' @param x Covariate data on a sphere
#' @param a1 The first element of the vector a, which is tuning parameter.
#' @param aremaining The remaining vector a, used as a starting guess.
#' @param param_mean Parameters for the mean link, used as a starting guess.
optim_alignedG <- function(y, x, a1, param_mean, k, aremaining, xtol_rel = 1E-2, ...){ #all the parameters are used as starting guesses, except a[1] that is a tuning parameter
  p <- ncol(y)
  om0 <- as_OmegaS2S(param_mean)
  P <- Omega2cann(om0)$P
  
  # preliminary estimate of the mean
  prelim <- optim_pobjS2S_parttape(y, x, om0)
  
  # generate tapes of ll that can be reused
  ll_mean <- tape_namedfun("ll_SvMF_S2S_alignedG_mean", 
                           OmegaS2S_vec(om0),
                           c(k, a1, aremaining, as.vector(P)),
                           p,
                           cbind(y, x))
  ll_mean_constraint <- tape_namedfun("wrap_OmegaS2S_constraints", OmegaS2S_vec(om0), vector(mode = "numeric"), p, matrix(nrow = 0, ncol = 0))
  ll_k <- tape_namedfun("ll_SvMF_S2S_alignedG_k",
                        k,
                        c(OmegaS2S_vec(om0), a1, aremaining, as.vector(P)),
                        p,
                        cbind(y, x))
  #skipped here: the tape for a must be recreated for each new version of Omega
  
  # prepare nloptr options
  default_opts <- list(xtol_rel = xtol_rel, #1E-04,
                       tol_constraints_eq = rep(1E-1, 2),
                       maxeval = 1E4)
  ellipsis_args <- list(...)
  combined_opts <- utils::modifyList(default_opts, ellipsis_args)
  
  # optimise iteratively. Start with k, then do a, then mean and repeat
  est0 <- list(
    mean = OmegaS2S_vec(prelim$solution),
    k = k,
    aremaining = aremaining
  )
  est <- est0 #iteratively update est
  diff <- abs(unlist(est)) * xtol_rel * 3
  browser()
  while (max(abs(diff)/abs(unlist(est))) > xtol_rel){ #tol is a relative finish
    estprev <- est
    
    # update k
    P <- Omega2cann(OmegaS2S_unvec(est$mean, p))$P
    newk <- nloptr::nloptr(
      x0 = est$k,
      eval_f = function(k){mean(scorematchingad:::pForward0(ll_k, k, c(est$mean, a1, est$aremaining, as.vector(P))))},
      eval_grad_f = function(k){mean(scorematchingad:::pJacobian(ll_k, k, c(est$mean, a1, est$aremaining, as.vector(P))))},
      opts = c(list(algorithm = "NLOPT_LD_SLSQP"), combined_opts[names(combined_opts) != "tol_constraints_eq"])
    )
    
    
    diff <- unlist(est) - unlist(estprev)
  }
  newk <- nloptr::nloptr(
    x0 = est$k,
    eval_f = function(theta){scorematchingad:::pForward0(ll_k, theta, vector(mode = "numeric"))},
    eval_grad_f = function(theta){scorematchingad:::pJacobian(obj_tape, theta, vector(mode = "numeric"))},
    eval_g_eq =  function(theta){scorematchingad:::pForward0(constraint_tape, theta, vector(mode = "numeric"))[1:2]},
    eval_jac_g_eq =  function(theta){matrix(scorematchingad:::pJacobian(constraint_tape, theta, vector(mode = "numeric")), byrow = TRUE, ncol = length(theta))},
    opts = combined_opts
  )
  
  

  
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