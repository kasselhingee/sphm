#' @title Objective function for vMF regression
#' @description 
#' For given parameters of [`mnlink()`], computes \deqn{\frac{1}{n}\sum_{i=1}^n y_i^T \mu(x_i)}
#' where \eqn{y_i} are observed unit vectors and \eqn{\mu(x_i)} is the corresponding predicted unit vector.
# @param y A matrix of unit vectors, each row corresponds to a single unit vector.
# @param xs A matrix of spherical covariate vectors (also unit vectors), each row corresponds to the same row in `y`.
# @param xe A matrix of Euclidean covariate vectors, each row corresponds to the same row in `y`.
# @param param A set of link parameters. See [`mnlink_cann()`] and [`mnlink_Omega()`].
#' @inheritParams mnlink
#' @export
prelimobj <- function(y, xs = NULL, xe = NULL, param){
  predictedmeans <- mnlink(xs = xs, xe = xe, param = param, check = FALSE)
  stopifnot(nrow(y) == nrow(predictedmeans))
  stopifnot(ncol(y) == ncol(predictedmeans))
  return(-mean(rowSums(y * predictedmeans)))
}

#' @title Regression with vMF Error and Scaled Mobius Link
#' @importClassesFrom scorematchingad Rcpp_ADFun
#' @description Regression of spherical response with vMF error and scaled Mobius link [`mnlink()`].
#' @details Optimisation of link parameters is by maximising 
#' \deqn{\sum_i=1^n y_i^T \mu(x_i)}
#' where \eqn{y_i} are observed unit vectors and \eqn{\mu(x_i)} is the corresponding predicted unit vector.
#' The concentration of the vMF is then estimated from the residuals.
#' 
#' Uses `nloptr`.
#' 
#' Before fitting, rotates y, xs and xe into standard form. If supplied, `start`, is updated accordingly.
#' 
#' If `type == "LinEuc"` a column of zeros called `'dummyzero'` is added to the front of `xe`.
#' 
#' @param y A matrix of unit vectors, each row corresponds to a single unit vector.
#' @param xs A matrix of spherical covariate vectors (also unit vectors), each row corresponds to the same row in `y`.
#' @param xe A matrix of Euclidean covariate vectors, each row corresponds to the same row in `y`.
#' @param start is a starting parameter object. For LinEuc mean link the Qe matrix must have an extra row and column that at the front/top, with 1 in the first entry (and zero elsewhere).
#' @param ... Passed as options to [`nloptr()`]. 
#' @param intercept `TRUE` to include a Euclidean intercept term using a covariate that is always `1`. This is needed for centering of Euclidean covariates, which is part of standardising the covariates. If `intercept = FALSE` then the Euclidean covariates will not be standardised.
#' @param lb Passed to [`nloptr()`]. Usually not used.
#' @param ub Passed to [`nloptr()`]. Usually not used.
#' @inheritParams mobius_SvMF
#' @export
mobius_vMF <- function(y, xs = NULL, xe = NULL, start = NULL, type = "SpEuc", fix_qs1 = FALSE, fix_qe1 = (type == "LinEuc"), intercept = TRUE, lb = NULL, ub = NULL, ...){
  p <- ncol(y)
  preplist <- list(y = y, xs = xs, xe = xe, start = start)
  # if needed, add Euclidean covariates and update start accordingly
  preplist <- addEuccovars(preplist, type = type, intercept = intercept)
  # standardise y, xe and xe and update start accordingly. Dont standardise xe if intercept = FALSE
  preplist <- standardise_data(preplist, intercept)
  # If start not supplied, choose start close to identities since data standardised
  preplist <- defaultstart(preplist, type)

  # Check LinEuc link initiated properly
  if ((type == "LinEuc") && (!is.null(preplist$xe))){
    stopifnot(is_LinEuc(preplist$start))
    stopifnot(all(preplist$xe[, 1]^2 < sqrt(.Machine$double.eps)))
  }

  ### More detailed preparation ###
  om0 <- as_mnlink_Omega(preplist$start)
  # check inputs:
  check_meanlink(preplist$y, preplist$xs, preplist$xe, om0)

  # Prepare constraint tape
  conprep <- estprep_meanconstraints(om0, fix_qs1, fix_qe1)
  # below updates om0vec with x0 values according to isfixed
  om0vec <- scorematchingad:::t_sfi2u(conprep$x0, conprep$om0vec, conprep$isfixed)

  # Prepare objective tape.
  objtape <- tape_namedfun("prelimobj_cpp", om0vec, vector(mode = "numeric"), c(p, length(om0$qe1)), cbind(preplist$y,preplist$xs,preplist$xe), check_for_nan = FALSE)
  objtape <- scorematchingad::avgrange(objtape) #objtape initially returns a value for each measurement. Average here to get average over all data.
  
  # update objtape based on fixed values
  objtape <- scorematchingad::fixindependent(objtape, objtape$xtape, conprep$isfixed)
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
    eval_g_eq =  function(theta){
      list(constraints = conprep$constraint_tape$forward(0, theta),
           jacobian = matrix(conprep$constraint_tape$Jacobian(theta), byrow = TRUE, ncol = length(theta)))},
    opts = combined_opts,
    lb = lb,
    ub = ub
  )
  if (!(nlopt$status %in% 1:4)){warning(nlopt$message)}
  
  # Estimate concentration
  # Note that the objective is average of y.ypred
  res <- optimise(function(k){
    -lvMFnormconst(k, p) + k * (-nlopt$objective) #full vMF log-likelihood (standardised by number of observations)
  }, lower = 1E-8, upper = 1E5, maximum = TRUE)
  k <- res$maximum
  
  #output some diagnostics - vector names would be nice here
  nlopt$solution_grad_f <- -objtape$Jacobian(nlopt$solution)
  nlopt$solution_g_eq <- conprep$constraint_tape$forward(0, nlopt$solution)
  nlopt$solution_jac_g_eq <- matrix(conprep$constraint_tape$Jacobian(nlopt$solution),
                                     byrow = TRUE, ncol = length(nlopt$solution))
  nlopt$solution_Hes_f <- matrix(-objtape$Hessian0(nlopt$solution),
         nrow = objtape$domain,
         byrow = TRUE)

  # remove the tapes from the return to save on memory
  nlopt$eval_f <- nlopt$eval_g_eq <- nlopt$nloptr_environment <- NULL
   
  # insert any fixed values of mean parameters
  fullparam <- scorematchingad:::t_sfi2u(nlopt$solution, om0vec, conprep$isfixed)
 
  #project mean pars to have correct orthogonality
  projectedom <- Omega_proj(mnlink_Omega_unvec(fullparam, p, length(om0$qe1), check = FALSE))
  try({mnlink_Omega_check(projectedom)})

  ### Making nicer return objects ###
  # Aspects of the fit that are invariant to coordinates used
  # distances in response space
  pred <- mnlink(xs = preplist$xs, xe = preplist$xe, param = projectedom)
  dists <- acos(rowSums(pred * preplist$y))
  rresids_tmp <- rotatedresid(preplist$y, pred, nthpole(ncol(preplist$y)))
  rresids <- rresids_tmp[, -1]
  attr(rresids, "samehemisphere") <-  attr(rresids_tmp, "samehemisphere")
  colnames(rresids) <- paste0("r", 1:ncol(rresids))
  
  ### revert estimated parameters and pred to pre-standardisation coordinates ###
  est <- undo_recoordinate_Omega(projectedom, 
                          yrot = attr(preplist$y, "std_rotation"), 
                          xsrot = attr(preplist$xs, "std_rotation"), #if xs/xe is NULL then attr(xs/xe, ..) is NULL too
                          xerot = attr(preplist$xe, "std_rotation"), 
                          xecenter = attr(preplist$xe, "std_center"),
                          onescovaridx = preplist$onescovaridx)
  ### Stabilise sign of Euc est based on ce ###
  if (isTRUE(est$ce < 0)){est <- Euc_signswitch(est)}
  
  # DoF
  DoF <- mobius_DoF(p, length(est$qs1), length(est$qe1), fix_qs1 = fix_qs1, fix_qe1 = fix_qe1) + 
    1 #concentration
  # AIC using result from concentration search
  AIC <- 2*DoF - 2 * res$objective * nrow(y)
  lLik <- res$objective * nrow(y)
  
  niceout <- list(
    est = est,
    mean = est,
    k = k,
    obj = nlopt$objective,
    solution = projectedom, #non-standardised solution
    nlopt = nlopt,
    y = y,
    xs = xs,
    xe = if (!is.null(xe)){if (intercept){destandardise_Euc(preplist$xe, attr(preplist$xe, "std_center"), attr(preplist$xe, "std_rotation"))} else {xe}}, #this recovers any added covariates too
    pred = destandardise_sph(pred, tG = attr(preplist$y, "std_rotation")),
    rresids = rresids,
    dists = dists,
    DoF = DoF,
    AIC = AIC,
    lLik = lLik,
    start = start,
    linktype = list(type = type, fix_qs1 = fix_qs1, fix_qe1 = fix_qe1, intercept = intercept)
  )
  return(niceout)
}

#' @title Degrees of freedom of the Mobius mean link function
#' @description The parameters of the Mobius mean link function [`mnlink()`] have a number of constraints.
#' This function incorporates these constraints to obtain the total degrees of freedom of the parameters.
#' @param p The length of response vectors
#' @param qs The length of spherical covariate vectors
#' @param qe The length of Euclidean covariate vectors
#' @param fix_qs1 Whether the first column of `Qs` is fixed (i.e. not estimated and not free).
#' @param fix_qe1 Whether `ce` and the first column of `Qe` is fixed (i.e. not estimated and not free).
#' @return An integer
#' @export
mobius_DoF <- function(p, qs = 0, qe = 0, fix_qs1 = FALSE, fix_qe1 = FALSE){
  if (qs == 0){fix_qs1 <- FALSE} #ignore fix_qs1
  if (qe == 0){fix_qe1 <- FALSE} #ignore fix_qs1
  DoF <- DoF_Stiefel(p, p) + #P
    DoF_Stiefel(qs-fix_qs1, p-fix_qs1) + #Qs
    DoF_Stiefel(qe-fix_qe1, p-fix_qe1) + #Qe
    (qs>0)*(p-1) + #Bs
    (qe>0)*(p-1) + #Be
    1*((qe>0) & (!fix_qe1)) #ce
}

check_meanlink <- function(y, xs, xe, om0){
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
}

estprep_meanconstraints <- function(om0, fix_qs1, fix_qe1){
  dims_in <- c(length(om0$p1), length(om0$qe1))
  om0vec <- mnlink_Omega_vec(om0)
 
  # generate tapes 
  constraint_tape <- tape_namedfun("Omega_constraints_wrap", om0vec, vector(mode = "numeric"), dims_in, matrix(nrow = 0, ncol = 0), check_for_nan = FALSE)

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
  # drop constraint returns that are constant:
  keep <- which(vapply((1:constraint_tape$range)-1, function(i){!constraint_tape$parameter(i)}, FUN.VALUE = FALSE))
  constraint_tape <- scorematchingad::keeprange(constraint_tape, keep)
  
  # check Jacobians of constraints are non-singular for the starting parameters.
  # For pathological params (e.g. the default starting params of no rotations), it can be zero.
  # If it is singular, perturb starting omega very slightly
  x0 <- om0vec[!isfixed]
  Jac_eq <- matrix(constraint_tape$Jacobian(x0), byrow = TRUE, ncol = constraint_tape$domain)
  if (any(abs(svd(Jac_eq)$d) < sqrt(.Machine$double.eps))){
    # modify Qs and Qe so that they dont individually form orthogonal vectors
    cann <- as_mnlink_cann(om0)
    scalemodifier <- max(abs(rbind(cann$Qs, cann$Qe)))/1E3
    cann$Qs[,-1] <- cann$Qs[,-1] + scalemodifier*savednoisemat[1:nrow(cann$Qs), 1:(ncol(cann$Qs)-1)]
    cann$Qe[,-1] <- cann$Qe[,-1] + scalemodifier*savednoisemat[1:nrow(cann$Qe), 1:(ncol(cann$Qe)-1)]
    om0new <- as_mnlink_Omega(cann)
    x0 <- mnlink_Omega_vec(om0new)[!isfixed]
    }
  Jac_eq <- matrix(constraint_tape$Jacobian(x0), byrow = TRUE, ncol = constraint_tape$domain)
  if (any(abs(svd(Jac_eq)$d) < sqrt(.Machine$double.eps))){
    warning("Initial parameters lead to a singular constraint Jacobian")
  }
  
  return(list(
    om0vec = om0vec,
    x0 = x0, #x0 may be perturbed to avoid singular Jac_eq 
    isfixed = isfixed,
    constraint_tape = constraint_tape
    ))
}

#' @title Try other initial starting parameters for a given regression
#' @description Given a vMF regression, repeat the optimisation from initial parameters randomly generated by [`rmnlink_cann()`].
#' @param mod_vMF Result of [`mobius_vMF()`]
#' @param preseed Passed to [`rmnlink_cann()`]
#' @details `fix_qe1` and `fix_qs1` of `mod_vMF` will be respected.
#' @return An vMF regression using the same data as `mod_vMF`.
#' @export
mobius_vMF_restart <- function(mod_vMF, preseed = 1){
  inparam <- as_mnlink_Omega(mod_vMF$est)
  dims <- dim.mnlink_Omega(inparam)
  start <- rmnlink_cann(p = dims["p"], 
                           qs = dims[["qs"]] - (dims[["qs"]] > 0 & mod_vMF$linktype$fix_qs1), 
                           qe = dims[["qe"]] - (dims[["qe"]] > 0 & mod_vMF$linktype$fix_qe1), 
                           preseed = preseed)
  # for sitatuations with fixed elements, treat simulated matrices as part of the full matrix
  if (dims[["qe"]] > 0 & mod_vMF$linktype$fix_qe1){
    start$ce <- inparam$ce
    # convert qe1 to nth pole
    rotmat <- rotationmat_amaral(inparam$qe1, nthpole(dims[["qe"]]))
    # build Qe from deciding it is random in the space orthogonal to qe1
    rotQe <- rbind(0, start$Qe)
    if (ncol(rotQe) == dims["p"]){
      rotQe <- cbind(0, rotQe[,-1, drop = FALSE])
    } else if (ncol(rotQe) == dims["p"] - 1) { #case when qe=p and fix_qe1=TRUE
      rotQe <- cbind(0, rotQe)
    }
    rotQe[1,1] <- 1
    start$Qe <- t(rotmat) %*% rotQe
    # put inparam$qe1 back in case inparam$qe1 = - nthpole(dims[["qe"]]) and rotmat is then a reflection not a rotation
    start$Qe[,1] <- inparam$qe1
  }
  if (dims[["qs"]] > 0 & mod_vMF$linktype$fix_qs1){
    # convert qs1 to nth pole
    rotmat <- rotationmat_amaral(inparam$qs1, nthpole(dims[["qs"]]))
    rotQs <- rbind(0, start$Qs)
    if (ncol(rotQs) == dims["p"]){
      rotQs <- cbind(0, rotQs[,-1, drop = FALSE])
    } else if (ncol(rotQs) == dims["p"] - 1) { #case when qs=p and fix_qs1=TRUE
      rotQs <- cbind(0, rotQs)
    }
    rotQs[1,1] <- 1
    start$Qs <- t(rotmat) %*% rotQs
    # put inparam$qs1 back in case inparam$qs1 = - nthpole(dims[["qs"]]) and rotmat is then a reflection not a rotation
    start$Qs[,1] <- inparam$qs1
  }
  start <- as_mnlink_Omega(start)
  stopifnot(all(dim.mnlink_Omega(start) == dims))
  mobius_vMF(y = mod_vMF$y,
             xs = mod_vMF$xs,
             xe = mod_vMF$xe,
             fix_qs1 = mod_vMF$linktype$fix_qs1,
             fix_qe1 = mod_vMF$linktype$fix_qe1,
             type = mod_vMF$linktype$type,
             intercept = mod_vMF$linktype$intercept,
             start = start)
}
