#' Optimise Mobius Link with vMF Error
#' @importClassesFrom scorematchingad Rcpp_ADFun
#' @details Assumes that the distribution is isotropic around the mean with constant concentration, thus maximising
#' \deqn{\sum_i=1^n y_i^T \mu(x_i).}
#' @param y A set of unit vectors in embedded coordinates, each row corresponds to a single unit vector.
#' @param x A set of covariate vectors (also unit vectors), each row corresponds to the same row in `y`.
#' @param paramobj A set of link parameters. See [`cannS2S()`] and [`OmegaS2S()`].
#' @inheritParams mnlink
#' @export
prelimobj <- function(y, xs = NULL, xe = NULL, param){
  predictedmeans <- mnlink(xs = xs, xe = xe, param = param, check = FALSE)
  stopifnot(nrow(y) == nrow(predictedmeans))
  stopifnot(ncol(y) == ncol(predictedmeans))
  return(-mean(rowSums(y * predictedmeans)))
}

#' Optimisation of the Preliminary Objective Function
#' @details Uses `nloptr`. `NLopt` doesn't have any algorithms for global optimisation with non-linear equality constraints that use provided gradients. So `_parttape` only does local optimisation and uses `NLOPT_LD_SLSQP` which is the only algorithm that takes advantage of derivatives and can handle non-linear equality constraints.
#' 
#' Before fitting, standardises y, xs and xe (*the latter needs implementing*). If supplied, `start`, is updated accordingly.
#' Note that if standardised y has a vMF distribution with the given means, the unstandardised y *does not* because of the second-moment standardisation (I would expect is to not be isotropic).
#' 
#' If `type == "Shogo"` a column of zeros called `'dummyzero'` is added to the front of `xe`.
#' 
#' Default scaling of 0.9 avoids being on the inequality boundary at the start of the search.
#' @param start is a starting parameter object. For Shogo mean link the Qe matrix must have an extra row and column that at the front/top, with 1 in the first entry (and zero elsewhere).
#' @param ... Passed as options to [`nloptr()`]. 
#' @param intercept `TRUE` to include a Euclidean intercept term using a covariate that is always `1`. This is needed for centering of Euclidean covariates, which is part of standardising the covariates. If `intercept = FALSE` then the Euclidean covariates will not be standardised.
#' @export
mobius_vMF <- function(y, xs = NULL, xe = NULL, start = NULL, type = "Kassel", fix_qs1 = FALSE, fix_qe1 = (type == "Shogo"), intercept = TRUE, ...){
  p <- ncol(y)
  preplist <- list(y = y, xs = xs, xe = xe, start = start)
  # if needed, add Euclidean covariates and update start accordingly
  preplist <- addEuccovars(preplist, type = type, intercept = intercept)
  # standardise y, xe and xe and update start accordingly. Dont standardise xe if intercept = FALSE
  preplist <- standardise_data(preplist, intercept)
  # If start not supplied, choose start close to identities since data standardised
  preplist <- defaultstart(preplist, type)

  # Check Shogo link initiated properly
  if ((type == "Shogo") && (!is.null(preplist$xe))){
    stopifnot(is_Shogo(preplist$start))
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
  
  # Optimisation
  # current dynamic parameter values of tapes will be used
  nlopt <- nloptr::nloptr(
    x0 = x0,
    eval_f = function(theta){-objtape$forward(0, theta)},
    eval_grad_f = function(theta){-objtape$Jacobian(theta)},
    eval_g_eq =  function(theta){conprep$constraint_tape$forward(0, theta)},
    eval_jac_g_eq =  function(theta){
      Jac <- matrix(conprep$constraint_tape$Jacobian(theta), byrow = TRUE, ncol = length(theta))
      # colnames(Jac) <- names(om0vec)
      # print(round(Jac, 3))
      # print(apply(Jac, 1, function(x)max(abs(x))))
      Jac
      },
    opts = combined_opts
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
  
  niceout <- list(
    est = est,
    k = k,
    obj = nlopt$objective,
    solution = projectedom, #non-standardised solution
    nlopt = nlopt,
    y = y,
    xs = xs,
    xe = if (!is.null(xe)){if (intercept){destandardise_Euc(preplist$xe, attr(preplist$xe, "std_center"), attr(preplist$xe, "std_rotation"))} else {xe}}, #this recovers any added covariates too
    pred = destandardise_sph(pred, tG = attr(preplist$y, "std_rotation")),
    rresids = rresids,
    dists = dists
  )
  return(niceout)
}

mobius_DoF <- function(p, qs = 0, qe = 0, fix_qs1 = FALSE, fix_qe1 = FALSE){
  DoF <- DoF_Stiefel(p, p) + #P
    DoF_Stiefel(qs, p) + #Qs
    DoF_Stiefel(qe, p) + #Qe
    (qs>0)*(p-1) + #Bs
    (qe>0)*(p-1) #Be
}

# Standard errors using the Fisher Information Matrix
# If k not supplied, estimates k
vMF_SE <- function(y, xs = NULL, xe = NULL, k = NULL, param, type = "Kassel"){
  p <- ncol(y)
  om <- as_mnlink_Omega(param)
  dims_in <- c(p, length(om$qe1))
  vec_om <- mnlink_Omega_vec(om)
  # Prepare objective tape.
  objtape_long <- tape_namedfun("prelimobj_cpp", vec_om, vector(mode = "numeric"), dims_in, cbind(y,xs,xe), check_for_nan = FALSE)
  
  if ((!is.null(xe)) && (type == "Shogo")){
    # if shogo and Euc, fix some elements
    omfixed <- lapply(om, function(x) x * 0)
    omfixed$qe1 <- omfixed$qe1 + 1
    omfixed$ce <- omfixed$ce + 1
    isfixed <- mnlink_Omega_vec(as_mnlink_Omega(omfixed)) > 0.5
    objtape_long <- scorematchingad::fixindependent(objtape_long, vec_om, isfixed)
    vec_om <- vec_om[!isfixed]
  }
  
  # Estimate concentration
  if (is.null(k)){
    mu_y <- mean(objtape_long$eval(vec_om, vector(mode = "numeric"))) #average of mu.y
    res <- optimise(function(k){
      -lvMFnormconst(k, p) + k * mu_y #full vMF log-likelihood (standardised by number of observations)
    }, lower = 1E-8, upper = 1E5, maximum = TRUE)
    k <- res$maximum
  }
  stop("SEs not implemented correctly")
  
  # Fisher Information Matrix is the variance of the gradient of the log-likelihood
  # And in MLE is equal to the expected double derivative of the log-likelihood
  # But I'm missing something because the covariance of the gradients gets scaled by k twice while the average hessian is scaled by k only once.
  # I think I need to project the gradients/hessians to the tangent of the parameter space: i.e. tangent to p1, qe1, qs1, and *somehow* Omega's special orthogonality constraints
  grads <- matrix(objtape_long$Jacobian(vec_om), byrow = TRUE, ncol = objtape_long$domain) * k #each row is the gradient at a data point
  FisherI <- stats::cov(grads) #also called the 'variablility matrix'
  
  # Sensitivity Matrix
  # E(d^2(ll)/dtheta^2) under certain regularity conditions (passing derivatives outside an intergal) will be equal to d^2(E[ll])/dtheta^2
  # Here, I dont assume the regularity conditions:
  jactape <- scorematchingad::tape_Jacobian(objtape_long) #rowwise fill of gradient of each data point
  allhess <- matrix(jactape$Jacobian(vec_om), byrow = TRUE, ncol = objtape_long$domain^2)
  sensitivitymat <- -matrix(colMeans(allhess), nrow = objtape_long$domain, ncol = objtape_long$domain) * k
  
  # method assuming that E(d^2(ll)/dtheta^2) = d^2(E[ll])/dtheta^2
  objtape <- scorematchingad::avgrange(objtape_long) #Average of mu.y.
  sensitivitymatb <- -matrix(objtape$Hessian0(vec_om), ncol = objtape$domain, nrow = objtape$domain) * k
  
  es <- eigen(sensitivitymat)
  
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
  # If it is singular, perturb start very slightly
  x0 <- constraint_tape$xtape
  Jac_eq <- matrix(constraint_tape$Jacobian(x0), byrow = TRUE, ncol = constraint_tape$domain)
  if (any(abs(svd(Jac_eq)$d) < sqrt(.Machine$double.eps))){x0 <- x0 - 1E-4}
  Jac_eq <- matrix(constraint_tape$Jacobian(x0), byrow = TRUE, ncol = constraint_tape$domain)
  stopifnot(all(abs(svd(Jac_eq)$d) > sqrt(.Machine$double.eps))) 
  
  return(list(
    om0vec = om0vec,
    x0 = x0, #x0 may be perturbed to avoid singular Jac_eq 
    isfixed = isfixed,
    constraint_tape = constraint_tape
    ))
}
