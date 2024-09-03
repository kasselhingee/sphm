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
                       maxeval = 1E4,
                       check_derivatives = TRUE)
  ellipsis_args <- list(...)
  combined_opts <- utils::modifyList(default_opts, ellipsis_args)
  
  # optimise iteratively. Start with k, then do a, then mean and repeat
  est0 <- list(
    mean = OmegaS2S_vec(om0),
    k = k,
    aremaining = aremaining
  )
  est <- est0 #iteratively update est
  diff <- abs(unlist(est)) * xtol_rel * 3
  iter <- 0
  while ( (iter < combined_opts$maxeval)  & (max(abs(diff)/abs(unlist(est))) > xtol_rel)){
    estprev <- est
    
    # update k
    P <- Omega2cann(OmegaS2S_unvec(est$mean, p))$P
    newk <- nloptr::nloptr(
      x0 = est$k,
      eval_f = function(k){-sum(scorematchingad:::pForward0(ll_k, k, c(est$mean, a1, est$aremaining, as.vector(P))))},
      eval_grad_f = function(k){-sum(scorematchingad:::pJacobian(ll_k, k, c(est$mean, a1, est$aremaining, as.vector(P))))},
      opts = c(list(algorithm = "NLOPT_LD_SLSQP"), combined_opts[names(combined_opts) != "tol_constraints_eq"])
    )
    est$k <- newk$solution
    
    # update aremaining
    ll_aremaining <- tape_namedfun("ll_SvMF_S2S_alignedG_a",
                          est$aremaining,
                          c(est$k, a1),
                          c(p, est$mean),
                          cbind(y, x))
    newaremaining <- nloptr::nloptr(
      x0 = est$aremaining,
      eval_f = function(aremaining){-sum(scorematchingad:::pForward0(ll_aremaining, aremaining, c(est$k, a1)))},
      eval_grad_f = function(aremaining){-colSums(matrix(scorematchingad:::pJacobian(ll_aremaining, aremaining, c(est$k, a1)), byrow = TRUE, ncol = length(aremaining)))},
      opts = c(list(algorithm = "NLOPT_LD_SLSQP"), combined_opts[names(combined_opts) != "tol_constraints_eq"])
    )
    est$aremaining <- newaremaining$solution
    
    #update mean link
    P <- Omega2cann(OmegaS2S_unvec(est$mean, p))$P
    newmean <- nloptr::nloptr(
      x0 = est$mean,
      eval_f = function(theta){-sum(scorematchingad:::pForward0(ll_mean, theta, c(est$k, a1, est$aremaining, as.vector(P))))},
      eval_grad_f = function(theta){-colSums(matrix(scorematchingad:::pJacobian(ll_mean, theta, c(est$k, a1, est$aremaining, as.vector(P))), byrow = TRUE, ncol = length(theta)))},
      opts = c(list(algorithm = "NLOPT_LD_SLSQP"), combined_opts[names(combined_opts) != "tol_constraints_eq"])
    )
    est$mean <- newmean$solution
    
    iter <- iter + 1
    diff <- unlist(est) - unlist(estprev)
    print(est)
  }
  
  return(list(
    solution = est,
    iter = iter
  ))
}