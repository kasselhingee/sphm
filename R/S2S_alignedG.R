#' MLE for the SvMF with S2S Mobius Link and alignedG_

#' @param y Response data on a sphere
#' @param x Covariate data on a sphere
#' @param a1 The first element of the vector a, which is tuning parameter.
#' @param aremaining The remaining vector a, used as a starting guess.
#' @param param_mean Parameters for the mean link, used as a starting guess.
optim_alignedG <- function(y, x, a1, param_mean, k, aremaining, xtol_rel = 1E-5, ...){ #all the parameters are used as starting guesses, except a[1] that is a tuning parameter
  p <- ncol(y)
  om0 <- as_OmegaS2S(param_mean)
  OmegaS2S_check(om0)
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
                       maxeval = 1E4,
                       check_derivatives = FALSE)
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
  times <- data.frame(list(k = NA, aremaining = NA, mean = NA)) #record times
  ests <- list() #track estimates
  while ( (iter < combined_opts$maxeval)  & (max(abs(diff)/abs(unlist(est))) > xtol_rel)){
    estprev <- est
    iter <- iter + 1
    P <- Omega2cann(OmegaS2S_unvec(est$mean, p))$P
    
    # update k
    ktime <- system.time({
    newk <- nloptr::nloptr(
      x0 = est$k,
      eval_f = function(k){-sum(scorematchingad:::pForward0(ll_k, k, c(est$mean, a1, est$aremaining, as.vector(P))))},
      eval_grad_f = function(k){-sum(scorematchingad:::pJacobian(ll_k, k, c(est$mean, a1, est$aremaining, as.vector(P))))},
      opts = c(list(algorithm = "NLOPT_LD_SLSQP"), combined_opts[names(combined_opts) != "tol_constraints_eq"])
    )})
    est$k <- newk$solution
    times[iter, "k"] <- sum(ktime[c("user.self", "user.child")])
    
    # update aremaining
    atime_start <- proc.time() # for timing
    if (length(est$aremaining) > 1) {
      #if its length one then the value must be exactly 1
      #otherwise here aremaining[1] will be derived from all the others
      ll_aremaining <- tape_namedfun("ll_SvMF_S2S_alignedG_a",
                          log(est$aremaining[-1]),
                          c(est$k, a1),
                          c(p, est$mean, as.vector(P)),
                          cbind(y, x))
      newlaremaining <- nloptr::nloptr( #searching for log(a) to avoid non-zero bound AND deriving the first value from all the others to avoid prod=1 requirement
        x0 = log(est$aremaining[-1]),
        eval_f = function(theta){
          -sum(scorematchingad:::pForward0(ll_aremaining, theta, c(est$k, a1)))},
        eval_grad_f = function(theta){
          -colSums(matrix(scorematchingad:::pJacobian(ll_aremaining, theta, c(est$k, a1)), byrow = TRUE, ncol = length(theta)))
          },
        opts = c(list(algorithm = "NLOPT_LD_SLSQP"), combined_opts)
      )
      est$aremaining <- exp(c(-sum(newlaremaining$solution), newlaremaining$solution))
    }
    # optimizing log a's in the following (optimising the aremaining would have first steps that went to below zero or super high)
    atime <- proc.time() - atime_start
    times[iter, "aremaining"] <- sum(atime[c("user.self", "user.child")])
    
    #update mean link
    mntime_start <- proc.time()
    newmean <- nloptr::nloptr(
      x0 = est$mean,
      eval_f = function(theta){-sum(scorematchingad:::pForward0(ll_mean, theta, c(est$k, a1, est$aremaining, as.vector(P))))},
      eval_grad_f = function(theta){-colSums(matrix(scorematchingad:::pJacobian(ll_mean, theta, c(est$k, a1, est$aremaining, as.vector(P))), byrow = TRUE, ncol = length(theta)))},
      eval_g_eq =  function(theta){scorematchingad:::pForward0(ll_mean_constraint, theta, vector(mode = "numeric"))[1:2]},
      eval_jac_g_eq =  function(theta){matrix(scorematchingad:::pJacobian(ll_mean_constraint, theta, vector(mode = "numeric")), byrow = TRUE, ncol = length(theta))},
      opts = c(list(algorithm = "NLOPT_LD_SLSQP", tol_constraints_eq = rep(1E-1, 2)), combined_opts)
    )
    est$mean <- OmegaS2S_vec(OmegaS2S_proj(OmegaS2S_unvec(newmean$solution, p, check = FALSE), method = "Omega"))
    mntime <- proc.time() - mntime_start
    times[iter, "mean"] <- sum(mntime[c("user.self", "user.child")])
    
    ests[[iter]] <- est
    diff <- unlist(est) - unlist(estprev)
    cat(".")
  }
  
  return(list(
    solution = est,
    iter = iter,
    ests = ests,
    times = times
  ))
}
