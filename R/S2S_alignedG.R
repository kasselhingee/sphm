#' MLE for the SvMF with S2S Mobius Link and alignedG_

#' @details The `NLOPT_LD_SLSQP` algorithm is a bit inconvenient as it seems to choose huge values of the parameters to try even though it is meant to local optimisation! I've had to put upper bounds of log(.Machine$double.xmax)/2-10 (~344) on aremaining
#' @param y Response data on a sphere
#' @param x Covariate data on a sphere
#' @param a1 The first element of the vector a, which is tuning parameter.
#' @param aremaining The remaining vector a, used as a starting guess.
#' @param param_mean Parameters for the mean link, used as a starting guess.
#' @param verbose `0` means no extra output. `1` means objective is printed each time `aremaining`, `k` and `param_mean` have all been updated. `2` further prints values of the parameters. (this roughly mirrors the behaviour of `print_level` for `nloptr()`).
#' @export
optim_alignedG <- function(y, x, a1, param_mean, k, aremaining, xtol_rel = 1E-5, verbose = 0, ...){ #all the parameters are used as starting guesses, except a[1] that is a tuning parameter
  p <- ncol(y)
  om0 <- as_OmegaS2S(param_mean)
  OmegaS2S_check(om0)
  P <- Omega2cann(om0)$P
  
  # generate tapes of ll that can be reused
  ll_mean <- tape_namedfun("ull_S2S_alignedG_mean", 
                           OmegaS2S_vec(om0),
                           c(k, a1, aremaining, as.vector(P)),
                           p,
                           cbind(y, x),
                           check_for_nan = FALSE)
  ll_mean_constraint <- tape_namedfun("wrap_OmegaS2S_constraints", OmegaS2S_vec(om0), vector(mode = "numeric"), p, matrix(nrow = 0, ncol = 0), check_for_nan = FALSE)
  ll_mean_ineqconstraint <- tape_namedfun("OmegaS2S_ineqconstaints", OmegaS2S_vec(om0), vector(mode = "numeric"), p, matrix(nrow = 0, ncol = 0), check_for_nan = FALSE)
  withCallingHandlers({
  ll_k <- tape_namedfun("ull_S2S_alignedG_k",
                        k,
                        c(OmegaS2S_vec(om0), a1, aremaining, as.vector(P)),
                        p,
                        cbind(y, x),
                        check_for_nan = FALSE)
  },
    warning = function(w){if (grepl("p=3", conditionMessage(w))){invokeRestart("muffleWarning")}}
  ) #this muffles a warning that the ull_S2S_alignedG_k function doesn't compute the normalising constant for p!=3.

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
  estinfo_k <- data.frame(list("status" = NA, "message" = NA, "iterations" = NA, "objective" = NA)) #record iterations etc of each
  estinfo_mn <- estinfo_a <- estinfo_k
  ests <- list() #track estimates
  while ( (iter < combined_opts$maxeval)  & (max(abs(diff)/abs(unlist(est))) > xtol_rel)){
    estprev <- est
    iter <- iter + 1
    
    # update k
    ktime <- system.time({
      if (p==3){newk <- optim_alignedG_k_3(est, ll_k, P, a1, combined_opts)}
      else {newk <- optim_alignedG_k_nograd(est, ll_k, P, a1, combined_opts)}
    })
    est$k <- newk$solution
    times[iter, "k"] <- sum(ktime[c("user.self", "user.child")])
    estinfo_k[iter,] <- newk[c("status", "message", "iterations", "objective")]
    
    # update aremaining
    atime_start <- proc.time() # for timing
    if (length(est$aremaining) > 1) {
      #if its length one then the value must be exactly 1
      #otherwise here aremaining[1] will be derived from all the others
      ll_aremaining <- tape_namedfun("ull_S2S_alignedG_a",
                          log(est$aremaining[-1]),
                          c(est$k, a1),
                          c(p, est$mean, as.vector(P)),
                          cbind(y, x),
                          check_for_nan = FALSE)
      newlaremaining <- nloptr::nloptr( #searching for log(a) to avoid non-zero bound AND deriving the first value from all the others to avoid prod=1 requirement
        x0 = log(est$aremaining[-1]),
        eval_f = function(theta){
          -sum(scorematchingad:::pForward0(ll_aremaining, theta, c(est$k, a1)))},
        eval_grad_f = function(theta){
          -colSums(matrix(scorematchingad:::pJacobian(ll_aremaining, theta, c(est$k, a1)), byrow = TRUE, ncol = length(theta)))
          },
        ub = rep(1, length(est$aremaining)-1) * log(.Machine$double.xmax), #this keeps the tapes evaluating to non infinite values. Note that R's (and I suspect C++'s generally) limit is log(.Machine$double.xmax), not sure why half is needed here. It isn't really needed if the tapes return the nan to R rather than erroring (check_for_nan = FALSE)
        opts = c(list(algorithm = "NLOPT_LD_SLSQP"), combined_opts)
      )
      est$aremaining <- exp(c(-sum(newlaremaining$solution), newlaremaining$solution))
    }
    atime <- proc.time() - atime_start
    times[iter, "aremaining"] <- sum(atime[c("user.self", "user.child")])
    if (length(est$aremaining) > 1) {
      estinfo_a[iter,] <- newlaremaining[c("status", "message", "iterations", "objective")]
    } else {
      estinfo_a[iter,] <- NA
    }
    
    #update mean link
    mntime_start <- proc.time()
    newmean <- nloptr::nloptr(
      x0 = est$mean,
      eval_f = function(theta){-sum(scorematchingad:::pForward0(ll_mean, theta, c(est$k, a1, est$aremaining, as.vector(P))))},
      eval_grad_f = function(theta){-colSums(matrix(scorematchingad:::pJacobian(ll_mean, theta, c(est$k, a1, est$aremaining, as.vector(P))), byrow = TRUE, ncol = length(theta)))},
      eval_g_eq =  function(theta){scorematchingad:::pForward0(ll_mean_constraint, theta, vector(mode = "numeric"))[1:2]},
      eval_jac_g_eq =  function(theta){matrix(scorematchingad:::pJacobian(ll_mean_constraint, theta, vector(mode = "numeric")), byrow = TRUE, ncol = length(theta))},
      eval_g_ineq =  function(theta){scorematchingad:::pForward0(ll_mean_ineqconstraint, theta, vector(mode = "numeric"))},
      eval_jac_g_ineq =  function(theta){matrix(scorematchingad:::pJacobian(ll_mean_ineqconstraint, theta, vector(mode = "numeric")), byrow = TRUE, ncol = length(theta))},
      opts = c(list(algorithm = "NLOPT_LD_SLSQP", tol_constraints_eq = rep(1E-1, 2)), combined_opts)
    )
    est$mean <- OmegaS2S_vec(OmegaS2S_proj(OmegaS2S_unvec(newmean$solution, p, check = FALSE), method = "Omega"))
    mntime <- proc.time() - mntime_start
    times[iter, "mean"] <- sum(mntime[c("user.self", "user.child")])
    estinfo_mn[iter,] <- newmean[c("status", "message", "iterations", "objective")]
    
    # Update P and check params
    Omegapar <- OmegaS2S_unvec(est$mean, p, check = FALSE)
    cannpar <- Omega2cann(Omegapar, check = FALSE)
    tryCatch({cannS2S_check(cannpar)}, error = function(e){warning(conditionMessage(e)); return(NULL)})
    P <- cannpar$P

    # get likelihood from newmean computations
    ll <- -newmean$objective
    if (p != 3){ll <- ll - nrow(y) * lvMFnormconst(est$k, length(est$aremaining) + 1)} #if p!=3 account for normalising constant outside the optimisation of the mean parameters
    
    ests[[iter]] <- est
    diff <- unlist(est) - unlist(estprev)
    if (verbose < 0.1){
      cat(".")
    } else {
      print(-ll)
      if (verbose > 1.1){
        print(cannpar)
        print(est$k)
        print(est$aremaining)
      } 
    }
  }
  
  colnames(times) <- paste0(colnames(times), ".usertime")
  estinfo <- cbind(times, k = estinfo_k, aremaining = estinfo_a, mean = estinfo_mn)
  
  return(list(
    solution = est,
    objective = -ll,
    iter = iter,
    ests = ests,
    estinfo = estinfo
  ))
}

#optimisation of k for p=3 using gradient
optim_alignedG_k_3 <- function(est, ll_k, P, a1, combined_opts){
    newk <- nloptr::nloptr(
      x0 = est$k,
      eval_f = function(k){-sum(scorematchingad:::pForward0(ll_k, k, c(est$mean, a1, est$aremaining, as.vector(P))))},
      eval_grad_f = function(k){-sum(scorematchingad:::pJacobian(ll_k, k, c(est$mean, a1, est$aremaining, as.vector(P))))},
      lb = 0, ub = .Machine$double.xmax,
      opts = c(list(algorithm = "NLOPT_LD_SLSQP"), combined_opts[names(combined_opts) != "tol_constraints_eq"])
    )
  return(newk)
}

# optimisation without using gradient
optim_alignedG_k_nograd <- function(est, ll_k, P, a1, combined_opts){
    newk <- nloptr::nloptr(
      x0 = est$k,
      eval_f = function(k){
        if (length(est$aremaining) + 1 == 3){ #p == 3
          lnormconst <- 0
        } else {
          lnormconst <- lvMFnormconst(k, length(est$aremaining) + 1)
        }
        -sum(-lnormconst + scorematchingad:::pForward0(ll_k, k, c(est$mean, a1, est$aremaining, as.vector(P))))},
      lb = 0, ub = .Machine$double.xmax,
      opts = c(list(algorithm = "NLOPT_LN_COBYLA"), combined_opts[names(combined_opts) != "tol_constraints_eq"])
    )
  return(newk)
}





