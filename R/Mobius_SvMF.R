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
#' @param doprelim When `FALSE` the preliminary von Mises-Fisher regression and subsequent moment estimation of `a` and `G0` is omitted. The provided parameters are used as the initial values for an optimisation of all parameters of the SvMF regression all together using [`nloptr::nloptr()`].
#' @param ... Named optional arguments passed as a list to the `opts` argument of [`nloptr::nloptr()`].
#' @details
#' From the starting parameters, optimises everything. For p != 3, the concentration is approximated.
#' No standardisation is performed.
#' @export
mobius_SvMF <- function(y, xs, xe, mean = NULL, k = NULL, a = NULL, G0 = NULL, G0reference = NULL, G01behaviour = "p1", type = "Shogo", fix_qs1 = FALSE, fix_qe1 = (type == "Shogo"), intercept = TRUE, doprelim = TRUE, ...){

  if (doprelim){
  preest <- mobius_SvMF_partransport_prelim(y, xs, xe, 
                                            mean = mean,
                                            G0 = G0, G01behaviour = G01behaviour, 
                                            type = type, fix_qs1 = fix_qs1, fix_qe1 = fix_qe1,
                                            intercept = intercept, ...)
  } else {
     stopifnot(!is.null(mean))
     stopifnot(!is.null(k))
     stopifnot(!is.null(a))
     stopifnot(!is.null(G0))
     preest <- list(
       mean = mean,
       k = k,
       a = a,
       G0 = G0)
  }

  finalest <- optim_constV(y, xs, xe, 
                           mean = preest$mean,
                           k = if(!is.null(k)){k}else{preest$k},
                           a = if(!is.null(a)){a}else{preest$a},
                           G0 = preest$G0, 
                           G0reference = if(!is.null(G0reference)){G0reference}else{preest$G0},
                           G01behaviour = G01behaviour, 
                           type = type, fix_qs1 = fix_qs1, fix_qe1 = fix_qe1,
                           intercept = intercept, ...)
  return(c(finalest, list(preest = preest)))
}

optim_constV <- function(y, xs, xe, mean, k, a, G0 = NULL, G0reference = NULL, G01behaviour = "p1", type = "Shogo", fix_qs1 = FALSE, fix_qe1 = (type == "Shogo"), intercept = TRUE, lb = NULL, ub = NULL, ...){
  initial <-  list(
    mean = mean,
    k = k, 
    a = a, 
    G0 = G0
  )
  SvMFcann_check(SvMFcann(k = initial$k, a = initial$a, G = initial$G0))

  p <- ncol(y)
  preplist <- list(y = y, xs = xs, xe = xe, start = mean)
  # if needed, add Euclidean covariates and update start accordingly
  preplist <- addEuccovars(preplist, type = type, intercept = intercept)
  # standardise y, xe and xe and update start accordingly. Dont standardise xe if intercept = FALSE
  preplist <- standardise_data(preplist, intercept)
  if (!is.null(G0)){preplist$G0 <- attr(preplist$y, "std_rotation") %*% G0} #update G0 too
  if (!is.null(G0reference)){preplist$G0reference <- attr(preplist$y, "std_rotation") %*% G0reference} #update G0 too
  # If start not supplied, choose start close to identities since data standardised
  preplist <- defaultstart(preplist, type)
  if (!is.null(a)){preplist$a <- a}
  if (!is.null(k)){preplist$k <- k}

  # Check Shogo link initiated properly
  if ((type == "Shogo") && (!is.null(preplist$xe))){
    stopifnot(is_Shogo(preplist$start))
    stopifnot(all(preplist$xe[, 1]^2 < sqrt(.Machine$double.eps)))
  }

  ### More detailed preparation ###
  om0 <- as_mnlink_Omega(preplist$start)
  # check inputs:
  check_meanlink(preplist$y, preplist$xs, preplist$xe, om0)
  stopifnot(length(preplist$a) == p)
  a1 = preplist$a[1]
  aremaining = preplist$a[-1]
  stopifnot(isTRUE(all.equal(prod(aremaining), 1)))
  if (is.null(preplist$G0reference)){
    G0reference <- preplist$G0
  } else {
    G0reference <- preplist$G0reference
  }
   
  qs <- length(om0$qs1)
  qe <- length(om0$qe1)


  # Prepare constraint tape
  conprep <- estprep_meanconstraints(om0, fix_qs1, fix_qe1)
  # below updates om0vec with x0 values according to isfixed
  om0vec <- scorematchingad:::t_sfi2u(conprep$x0, conprep$om0vec, conprep$isfixed)
  
  # Prepare objective tape
  objtape_ind <- tape_ull_S2S_constV_nota1(omvec = om0vec, k = preplist$k,
                                       a1 = a1, 
                                       aremaining = aremaining,
                                       G0 = preplist$G0,
                                       p = p, qe = qe, 
                                       yx = cbind(preplist$y, preplist$xs, preplist$xe),
                                       referencecoords = G0reference,
                                       G01behaviour = G01behaviour)
  # update objtape based on fixed values
  objtape_ind <- scorematchingad::fixindependent(objtape_ind, objtape_ind$xtape, c(conprep$isfixed, rep(0, length(objtape_ind$xtape) - length(conprep$isfixed))))
 
  #objtape initially returns a value for each measurement. Average here to get average over all data.
  objtape <- scorematchingad::avgrange(objtape_ind)
  # using tapeing values as starting parameters
  x0 <- objtape$xtape

  # prepare nloptr options
  default_opts <- list(algorithm = "NLOPT_LD_SLSQP",
                       xtol_rel = 1E-10, #1E-04,
                       tol_constraints_eq = rep(1E-1, conprep$constraint_tape$range),
                       # check_derivatives = TRUE, check_derivatives_print = 'errors', check_derivatives_tol = 1E-3,
                       # print_level = 3,
                       maxeval = 1E4)
  ellipsis_args <- list(...)
  combined_opts <- utils::modifyList(default_opts, ellipsis_args)
  
  # set lower bound for concentration k, unless lb passed manually
  if (is.null(lb)){
    lb <- rep(-Inf, objtape$domain)
    lb[length(conprep$x0) + 1] <- 0
  }

  # activate a progress bar
  pb <- progress::progress_bar$new(total = combined_opts$maxeval + 5, format = ":bar :percent :current :tick_rate elapsed::elapsedfull eta::eta")
  
  
  # Optimisation
  # current dynamic parameter values of tapes will be used
  nlopt <- nloptr::nloptr(
    x0 = x0,
    eval_f = function(theta){
      if (!pb$finished) pb$tick()
      list(objective = -objtape$forward(0, theta), gradient = -objtape$Jacobian(theta))
      },
    lb = lb,
    ub = ub,
    eval_g_eq =  function(theta){list(
      constraints = conprep$constraint_tape$forward(0, theta[1:conprep$constraint_tape$domain]),
      jacobian = cbind(matrix(conprep$constraint_tape$Jacobian(theta[1:conprep$constraint_tape$domain]), byrow = TRUE, ncol = conprep$constraint_tape$domain),
             matrix(0, nrow = conprep$constraint_tape$range, ncol = length(theta) - conprep$constraint_tape$domain))
    )},
    opts = combined_opts
  )
  if (!(nlopt$status %in% 1:4)){warning(nlopt$message)}

  #output some diagnostics - vector names would be nice here
  nlopt$solution_grad_f <- -objtape$Jacobian(nlopt$solution)
  nlopt$solution_g_eq <- conprep$constraint_tape$forward(0, nlopt$solution[1:conprep$constraint_tape$domain])
  nlopt$solution_jac_g_eq <- matrix(conprep$constraint_tape$Jacobian(nlopt$solution[1:conprep$constraint_tape$domain]),
                                     byrow = TRUE, ncol = length(nlopt$solution[1:conprep$constraint_tape$domain]))
  nlopt$solution_Hes_f <- matrix(-objtape$Hessian0(nlopt$solution),
         nrow = objtape$domain,
         byrow = TRUE)
  nlopt$ldens <- drop(objtape_ind$forward(0,nlopt$solution))
  lLik <- sum(nlopt$ldens)

  # remove the tapes from the return to save on memory
  nlopt$eval_f <- nlopt$eval_g_eq <- nlopt$eval_g_ineq <- nlopt$nloptr_environment <- NULL

  # insert any fixed values of mean parameters
  meanpars <- nlopt$solution[1:length(conprep$x0)]
  meanpars <- scorematchingad:::t_sfi2u(meanpars, om0vec, conprep$isfixed)
  fullparam <- c(meanpars, nlopt$solution[-(1:length(conprep$x0))])

  
  estparamlist <- S2S_constV_nota1_fromvecparamsR(fullparam, p, qs, qe, 
                                                  referencecoords = G0reference,
                                                  G01behaviour = G01behaviour,
                                                  G01 = preplist$G0[,1])
  
  #project mean pars to have correct orthogonality
  projectedom <- Omega_proj(mnlink_Omega_unvec(estparamlist$omvec, p, length(om0$qe1), check = FALSE))
  try({mnlink_Omega_check(projectedom)})

  # make sure aremaining is in decreasing order
  aord <- order(estparamlist$aremaining, decreasing = TRUE)
  aremaining <- estparamlist$aremaining[aord]
  estparamlist$G0[,-1] <- estparamlist$G0[,-1][, aord]
  
  # Polishing optimisation of just concentration if vMF normalising constant being approximated
  pred <- mnlink(xs = preplist$xs, xe = preplist$xe, param = projectedom)
  if (p!=3){
    result <- mobius_SvMF_konly(y = preplist$y, ymean = pred, a = c(a1, aremaining), G0 = estparamlist$G0)
    # update parameters
    estparamlist$k <- result$k
    # update log-likelood
    lLik <- result$lLik
  }
  
  ### Making nicer return objects ###
  # Aspects of the fit that are invariant to coordinates used
  # distances in response space
  
  dists <- acos(rowSums(pred * preplist$y))
  # get residuals as coordinates wrt G0. So under high concentration these residuals follow something multivariate normal.
  rresids_std <- resid_SvMF_partransport(preplist$y, pred, estparamlist$k, c(a1, estparamlist$aremaining), estparamlist$G0, scale = TRUE)
  rresids_G0 <- resid_SvMF_partransport(preplist$y, pred, G0 = estparamlist$G0, scale = FALSE) # par transport to G01, coords G0
  rresids_I_tmp <- rotatedresid(preplist$y, pred, nthpole(ncol(preplist$y))) # par transport to nthpole, coords cannonical
  rresids_I <- rresids_I_tmp[, -1]
  attr(rresids_I, "samehemisphere") <-  attr(rresids_I_tmp, "samehemisphere")
  colnames(rresids_I) <- paste0("r", 1:ncol(rresids_I))
  

  
  ### revert estimated parameters and pred to pre-standardisation coordinates ###
  est <- undo_recoordinate_Omega(projectedom, 
                          yrot = attr(preplist$y, "std_rotation"), 
                          xsrot = attr(preplist$xs, "std_rotation"), #if xs/xe is NULL then attr(xs/xe, ..) is NULL too
                          xerot = attr(preplist$xe, "std_rotation"), 
                          xecenter = attr(preplist$xe, "std_center"),
                          onescovaridx = preplist$onescovaridx)

  #put G0 into same coordinates as y
  G0 <- estparamlist$G0
  G0 <- t(attr(preplist$y, "std_rotation")) %*% G0
  #For axes G0 standardise the return by
  # (1) make first element of each vector positive (except the first column)
  G0[,-1] <- toBigPosEl(G0[,-1])
  # (2) make rotation matrix by flipping final column according to determinant
  if (det(G0) < 0){G0[,p] <- -G0[,p]}
  # DoF
  DoF <- mobius_DoF(p, length(est$qs1), length(est$qe1), fix_qs1 = fix_qs1, fix_qe1 = fix_qe1) + 
    1 + #concentration
    (p-1)-1 + #aremaining given that prod(aremaining) = 1
    if (G01behaviour == "free"){ #G0 freedom
      DoF_Stiefel(p,p)
    } else {
      DoF_Stiefel(p-1, p-1)
    }
  # AIC
  AIC = 2*DoF - 2 * lLik
  
  if (p==3) {
    if (estparamlist$k < 1E-15){
    	warning("Estimated concentration is very small and computation of the vMF normalising constant may be breaking down.")
    }
  }
  
  #Scealy and Wood (2019) Proposition 1 check for unimodality
  SvMFcann_check(SvMFcann(k = estparamlist$k, a = c(a1, aremaining), G = G0))
  
  niceout <- list(
    mean = est,
    k = estparamlist$k,
    a = c(a1, aremaining),
    G0 = G0,
    obj = nlopt$objective,
    nlopt = nlopt,
    y = y,
    xs = xs,
    xe = if (!is.null(xe)){if (intercept){destandardise_Euc(preplist$xe, attr(preplist$xe, "std_center"), attr(preplist$xe, "std_rotation"))} else {xe}}, #this recovers any added covariates too
    pred = destandardise_sph(pred, tG = attr(preplist$y, "std_rotation")),
    rresids = rresids_I,
    rresids_I = rresids_I,
    rresids_G0 = rresids_G0,
    rresids_std = rresids_std,
    dists = dists,
    DoF = DoF,
    AIC = AIC,
    lLik = lLik,
    initial = initial
  )
  return(niceout)
}

mobius_SvMF_partransport_prelim <- function(y, xs, xe, mean = NULL, G0 = NULL, G01behaviour = "p1", type = "Shogo", fix_qs1 = FALSE, fix_qe1 = (type == "Shogo"), intercept = TRUE, ...){
  prelim <- mobius_vMF(y = y, xs = xs, xe = xe, 
             start = mean, 
             type = type, fix_qs1 = fix_qs1, fix_qe1 = fix_qe1, intercept = intercept, ...)
  # update starting values accordingly
  mean <- as_mnlink_Omega(prelim$est)
  k <- prelim$k
  p <- ncol(y)
  
  # get/choose G01 depending on behaviour
  if (G01behaviour == "fixed" && is.null(G0)){stop("At least the first column of G0 must be supplied when G01behaviour = 'fixed'.")}
  G01 <- switch(G01behaviour,
         p1 = prelim$est$p1,
         free = if(is.null(G0)){prelim$est$p1}else{G0[,1]},
         fixed = G0[,1])
  # get rotated residuals
  rresid <- rotatedresid(y, prelim$pred, base = G01)
  if (!is.null(G0) && all(!is.na(G0))){
    # if G0 fully supplied just use rresid to approximate scales a using the high concentration approximation
    if (G01behaviour == "p1"){
      G0 <- cbind(G01, -JuppRmat(G0[,1], G01) %*% G0[,-1])
    }
    aremaining <- SvMF_prelim_scales(rresid, G0)
  } else {
    # axes:
    G0 <- SvMF_mom_axes(rresid, G01)
    # estimate the scales
    aremaining <- SvMF_prelim_scales(rresid, G0)
  }
  
  prelim <- list(
    mean = mean,
    k = k,
    a = c(1, aremaining),
    G0 = G0,
    nlopt = prelim$nlopt
  )
  return(prelim)
}

mobius_SvMF_konly <- function(y, ymean, a, G0){
  yrot <- undo_partransport(y = y, ymean = ymean, G01 = G0[,1])
  res <- optimise(function(k){
    sum(SvMF_ll_cann(yrot, SvMFcann(k = k, a = a, G = G0)))
  }, lower = 1E-8, upper = 1E5, maximum = TRUE)
  SvMFcann_check(SvMFcann(k = res$maximum, a = a, G = G0))
  if (res$maximum == 1E-8){warning("Concentration at numerical lower limit of 1E-8")}
  if (res$maximum == 1E5){warning("Concentration at numerical upper limit of 1E5")}
  return(list(
    k = res$maximum,
    a = a,
    G0 = G0,
    lLik = res$objective
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


# @param y matrix of observations
# @param ymean matrix of predicted means
# @param G01 first column of the G0 matrix
# Parallel transports and rotates y so that it has ymean = G01
undo_partransport <- function(y, ymean, G01){
  #rotate all observations to reverse the transport from G0[,1] to ymean
  yrot <- lapply(1:nrow(y), function(idx){
    drop(t(rotationmat_amaral(G01, ymean[idx, ])) %*% y[idx, ])
  })
  yrot <- do.call(rbind, yrot)
  return(yrot)
}

# log-density of each row of y according to a mobius_SvMF regression model
dS2S_constV <- function(y, xs, xe, mean, k, a, G0){
  ymean <- mnlink(xs = xs, xe = xe, param = mean)
  diff <- lvMFnormconst_approx(k, ncol(y)) - lvMFnormconst(k, ncol(y))
  
  #rotate all observations so that ymean --> G0[,1]
  yrot <- undo_partransport(y, ymean, G01 = G0[,1])
  ldCpp <- uldSvMF_cann(yrot, k = k, a = a, G = G0)
  ldR <- SvMF_ll_cann(yrot, SvMFcann(k = k, a = a, G = G0))
  ld <- cbind(Cpp = ldCpp, R = ldR)
  attr(ld, "error") <- diff
  return(ld)
}
