#' Homosckedastic SvMF Regression
#' @details
#' The mean is assumed to follow the usual mean link.
#' The concentration and scaling in the SvMF is assumed constant across observations.
#' The scaling axes of the SvMF at location \eqn{\mu} are assumed to be the parallel transport along the geodesic of axes at the first column of the matrix `P` from the mean link. These axes specified at first column of the matrix `P` are to be estimated and constant with respect to covariates (and \eqn{\mu})
#' __Warning: C++ function still uses Jupp's transport rather than Amaral matrix and p !=3 not handled yet__
#' @param y Response data on a sphere
#' @param x Covariate data on a sphere
#' @param a1 The first element of the vector a, which is tuning parameter.
#' @param aremaining The remaining vector a, used as a starting guess.
#' @param mean Parameters for the mean link, used as a starting guess.
#' @param Gstar starting guess of the axes at `p1`.
#' @export
optim_constV <- function(y, x, mean, k, a, Gstar, xtol_rel = 1E-5, verbose = 0, ...){
  p <- ncol(y)
  q <- ncol(x)
  # checks
  om0 <- as_OmegaS2S(mean)
  OmegaS2S_check(om0)
  if (!isTRUE(all.equal(cbind(om0$p1, Gstar) %*% t(cbind(om0$p1, Gstar)), diag(1, p), check.attributes = FALSE))){ # p1 orthogonal to Vstar
    stop("Gstar is not orthogonal to p1.")
  }
  stopifnot(length(a) == p)
  a1 = a[1]
  aremaining = a[-1]
  stopifnot(isTRUE(all.equal(prod(aremaining), 1)))
  
  # standardisation of data
  stdmat <- standardise_mat(y)
  ystd <- y %*% stdmat
  # apply same operation to initial parameters
  cann0 <- as_cannS2S(om0)
  om0std <- as_OmegaS2S(cannS2S(t(stdmat) %*% cann0$P, cann0$Q, cann0$B))
  stdGstar <- t(stdmat) %*% Gstar #Because stdmat performs a rigid transformation, it is really just a change in basis for the whole problem, so I think this is what we want for the axes too.
  stdKstar <- t(getHstar(om0std$p1)) %*% stdGstar
  stdKstar[, 1] <- det(stdKstar) * stdKstar[,1] #because Cayley transform only works on det of +1
  
  # preliminary estimate of mean link
  estprelim <- optim_pobjS2S_parttape(ystd, x, om0std)
  if (!(estprelim$loc_nloptr$status %in% c(0, 1, 2, 3, 4))){warning("Preliminary optimistation did not finish properly.")}
  om0prelim <- estprelim$solution
  
  # estimation any p prep
  omvec0 <- OmegaS2S_vec(om0prelim)
  ll_mean_constraint <- tape_namedfun("wrap_OmegaS2S_constraints", omvec0, vector(mode = "numeric"), p, matrix(nrow = 0, ncol = 0), check_for_nan = FALSE)
  # prepare nloptr options
  default_opts <- list(xtol_rel = xtol_rel, #1E-04,
                       maxeval = 1E4,
                       check_derivatives = FALSE)
  ellipsis_args <- list(...)
  combined_opts <- utils::modifyList(default_opts, ellipsis_args)
  
  #estimation for p=3 specific
  ulltape <- tape_ull_S2S_constV_nota1(omvec = omvec0, k = k,
                                       a1 = a1, aremaining = aremaining, Kstar = stdKstar,
                                       p = p, cbind(ystd, x))
  # lower bound for concentration k
  lb <- rep(-Inf, ulltape$domain)
  lb[p + q + p*q + 1] <- 0
  est <- nloptr::nloptr(
    x0 = S2S_constV_nota1_tovecparams(omvec = omvec0, k = k, aremaining = aremaining, Kstar = stdKstar),
    eval_f = function(theta){
      # print(unlist(S2S_constV_nota1_fromvecparamsR(theta, p, q)[c("k", "aremaining")]))
      -sum(ulltape$eval(theta, a[1]))
      },
    eval_grad_f = function(theta){-colSums(matrix(ulltape$Jac(theta, a[1]), byrow = TRUE, ncol = length(theta)))},
    lb = lb,
    eval_g_eq =  function(theta){ll_mean_constraint$eval(theta[1:length(omvec0)], vector(mode = "numeric"))},
    eval_jac_g_eq =  function(theta){
      cbind(matrix(ll_mean_constraint$Jac(theta[1:length(omvec0)], vector(mode = "numeric")), byrow = TRUE, ncol = length(omvec0)),
            matrix(0, nrow = 2, ncol = length(ulltape$xtape) - length(ll_mean_constraint$xtape)))
    },
    opts = c(list(algorithm = "NLOPT_LD_SLSQP", tol_constraints_eq = rep(1E-1, 2)), combined_opts)
  )
  
  # estimation for p!=3 (alternating between k and others) ## NOT COMPLETE
  
  if (!(est$status %in% c(0, 1, 2, 3, 4))){warning("Optimistation did not finish properly.")}
  estparamlist <- S2S_constV_nota1_fromvecparamsR(est$solution, p, q)
  
  # project Omega to satisfy orthogonality constraint
  est_om <- OmegaS2S_proj(OmegaS2S_unvec(estparamlist$omvec, p, check = FALSE))
  
  # calculate Gstar now (because getHstar is sensitive to changes of basis A * H(p1) != H(A*p1))
  Gstar <- getHstar(est_om$p1) %*% estparamlist$Kstar
  
  # undo standardisation coordinate change
  est_cann <- as_cannS2S(est_om)
  est_om <- as_OmegaS2S(cannS2S(stdmat %*% est_cann$P, est_cann$Q, est_cann$B))
  Gstar <- stdmat %*% Gstar
  
  
  outsolution <- list(
    mean = est_om,
    k = estparamlist$k,
    a = c(a1, estparamlist$aremaining),
    Gstar = Gstar
  )
  
  return(list(
    solution = outsolution,
    stdmat = stdmat,
    nlopt_prelim = estprelim,
    nlopt_final = est,
    initial = list(
      mean = mean,
      k = k, 
      a = a, 
      Gstar = Gstar
    )
  ))
}


