
#' @param ssqOmbuffer The sum of squared singular values of Omega is allowed to go `ssqOmbuffer` above the limit given by singular values of 1 (or 2 if there are both Euclidean and spherical coordinates).
prelim_ad <- function(y, xs = NULL, xe = NULL, paramobj0, type = "Kassel", globalfirst = FALSE, ssqOmbuffer = 2, ...){ #paramobj0 is the starting parameter object
  om0 <- as_mnlink_Omega(paramobj0)
  # check inputs:
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
    if (type == "Shogo"){
      stopifnot(is_Shogo(om0))
      stopifnot(all(xe[, 1]^2 < sqrt(.Machine$double.eps)))
    }
  } else {
    stopifnot(length(om0$qe1) == 0)
  }

  dims_in <- c(p, length(om0$qe1))
  vec_om0 <- mnlink_Omega_vec(om0)
  # Prepare objective tape.
  obj_tape <- tape_namedfun("prelimobj_cpp", vec_om0, vector(mode = "numeric"), dims_in, cbind(y,xs,xe), check_for_nan = FALSE)
  obj_tape <- scorematchingad::avgrange(obj_tape) #Average of mu.y.
  constraint_tape <- tape_namedfun("Omega_constraints_wrap", vec_om0, vector(mode = "numeric"), dims_in, matrix(nrow = 0, ncol = 0), check_for_nan = FALSE)
  ineqconstraint_tape <- tape_namedfun("Omega_ineqconstraints", vec_om0, vector(mode = "numeric"), dims_in, matrix(nrow = 0, ncol = 0), check_for_nan = FALSE)
  
  if ((!is.null(xe)) && (type == "Shogo")){
    # if shogo and Euc, fix some elements
    omfixed <- lapply(om0, function(x) x * 0)
    omfixed$qe1 <- omfixed$qe1 + 1
    omfixed$ce1 <- omfixed$ce1 + 1
    isfixed <- mnlink_Omega_vec(as_mnlink_Omega(omfixed)) > 0.5
    obj_tape <- scorematchingad::fixindependent(obj_tape, vec_om0, isfixed)
    constraint_tape <- scorematchingad::fixindependent(constraint_tape, vec_om0, isfixed)
    ineqconstraint_tape <- scorematchingad::fixindependent(ineqconstraint_tape, vec_om0, isfixed)
    vec_om0 <- vec_om0[!isfixed]
    
    # drop the qe constraint
    keep <- setdiff(1:constraint_tape$range, 2 + !is.null(xs))
    constraint_tape <- scorematchingad:::keeprange(constraint_tape, keep)
  }
  
  # check Jacobians of constraints are non-singular for the starting parameters.
  # For pathological params (e.g. the default starting params of no rotations), it can be zero.
  # If it is singular, perturb start very slightly
  Jac_eq <- matrix(constraint_tape$Jacobian(vec_om0), byrow = TRUE, ncol = length(vec_om0))
  if (any(abs(svd(Jac_eq)$d) < sqrt(.Machine$double.eps))){vec_om0 <- vec_om0+1E-4}
  Jac_eq <- matrix(constraint_tape$Jacobian(vec_om0), byrow = TRUE, ncol = length(vec_om0))
  stopifnot(all(abs(svd(Jac_eq)$d) > sqrt(.Machine$double.eps)))
  Jac_ineq <- matrix(ineqconstraint_tape$Jacobian(vec_om0), byrow = TRUE, ncol = length(vec_om0))
  stopifnot(all(abs(svd(Jac_ineq)$d) > sqrt(.Machine$double.eps)))
  

  if (globalfirst){ #do a quick global search first #nlopt recommends doing a local search afterwards
    default_opts <- list(algorithm = "NLOPT_GN_ISRES", #the only algorthim that natively handles non-linear equality constraints - all the others have to use augmented Lagrangian ideas.
	        xtol_rel = 1E-04,
	        tol_constraints_eq = rep(1E-1, constraint_tape$range),
	        maxeval = 1E2)
    ellipsis_args <- list(...)
    combined_opts <- utils::modifyList(default_opts, ellipsis_args)
    globopt <- nloptr::nloptr(
      x0 = vec_om0,
      eval_f = function(theta){-obj_tape$eval(theta, vector(mode = "numeric"))},
      eval_g_eq =  function(theta){constraint_tape$eval(theta, vector(mode = "numeric"))},
      eval_g_ineq =  function(theta){ineqconstraint_tape$eval(theta, vector(mode = "numeric")) - ssqOmbuffer},
      lb = vec_om0 * 0 - 10, #10 is just a guess here. For the spherical covariate stuff, I suspect most values are well below 1. *Euc will be different*
      ub = vec_om0 * 0 + 10,
      opts = combined_opts
    )
    warning("lb and ub not properly set")
    if (!(globopt$status %in% 1:4)){warning(globopt$message)}
    vec_om0 <- globopt$solution
  }
  
  # prepare nloptr options
  default_opts <- list(algorithm = "NLOPT_LD_SLSQP",
                xtol_rel = 1E-10, #1E-04,
                tol_constraints_eq = rep(1E-1, constraint_tape$range),
                # check_derivatives = TRUE, check_derivatives_print = 'errors', check_derivatives_tol = 1E-3,
                # print_level = 3,
                maxeval = 1E4)
  ellipsis_args <- list(...)
  combined_opts <- utils::modifyList(default_opts, ellipsis_args)
  
  locopt <- nloptr::nloptr(
    x0 = vec_om0,
    eval_f = function(theta){-obj_tape$eval(theta, vector(mode = "numeric"))},
    eval_grad_f = function(theta){-obj_tape$Jac(theta, vector(mode = "numeric"))},
    eval_g_eq =  function(theta){constraint_tape$eval(theta, vector(mode = "numeric"))},
    eval_jac_g_eq =  function(theta){
      Jac <- matrix(constraint_tape$Jacobian(theta), byrow = TRUE, ncol = length(theta))
      # colnames(Jac) <- names(vec_om0)
      # print(round(Jac, 3))
      # print(apply(Jac, 1, function(x)max(abs(x))))
      Jac
      },
    eval_g_ineq =  function(theta){ineqconstraint_tape$eval(theta, vector(mode = "numeric")) - ssqOmbuffer},
    eval_jac_g_ineq =  function(theta){
      Jac <- matrix(ineqconstraint_tape$Jacobian(theta), byrow = TRUE, ncol = length(theta))
      # print(apply(Jac, 1, function(x)max(abs(x))))
      Jac
      },
    opts = combined_opts
  )
  if (!(locopt$status %in% 1:4)){warning(locopt$message)}

  # Because locopt$solution isnt the full Omega parameterisation for 'Shogo' type, we need to rebuild it.
  if (type == "Shogo"){
    fullparam <- isfixed
    fullparam[!isfixed] <- locopt$solution
    fullparam[isfixed] <- mnlink_Omega_vec(om0)[isfixed]
  } else {
    fullparam <- locopt$solution
  }
 
  unprojresult <- mnlink_Omega_unvec(fullparam, p, length(om0$qe1), check = FALSE)
  # mnlink_Omega_check(unprojresult)
  projresult <- Omega_proj(unprojresult)
  # mnlink_Omega_check(projresult)
  
  # remove the tapes from the return to save on memory
  locopt$eval_f <- locopt$eval_g_eq <- locopt$eval_g_ineq <- locopt$nloptr_environment <- NULL
  
  return(list(
    solution = projresult,
    loc_nloptr = locopt
  ))
}

# Standard errors using the Fisher Information Matrix
# If k not supplied, estimates k
vMF_SE <- function(y, xs = NULL, xe = NULL, k = NULL, param, type = "Kassel"){
  p <- ncol(y)
  om <- as_mnlink_Omega(param)
  dims_in <- c(p, length(om$qe1))
  vec_om <- mnlink_Omega_vec(om)
  # Prepare objective tape.
  obj_tape_long <- tape_namedfun("prelimobj_cpp", vec_om, vector(mode = "numeric"), dims_in, cbind(y,xs,xe), check_for_nan = FALSE)
  
  if ((!is.null(xe)) && (type == "Shogo")){
    # if shogo and Euc, fix some elements
    omfixed <- lapply(om, function(x) x * 0)
    omfixed$qe1 <- omfixed$qe1 + 1
    omfixed$ce1 <- omfixed$ce1 + 1
    isfixed <- mnlink_Omega_vec(as_mnlink_Omega(omfixed)) > 0.5
    obj_tape_long <- scorematchingad::fixindependent(obj_tape_long, vec_om, isfixed)
    vec_om <- vec_om[!isfixed]
  }
  
  # Estimate concentration
  if (is.null(k)){
    mu_y <- mean(obj_tape_long$eval(vec_om, vector(mode = "numeric"))) #average of mu.y
    res <- optimise(function(k){
      -lvMFnormconst(k, p) + k * mu_y #full vMF log-likelihood (standardised by number of observations)
    }, lower = 1E-8, upper = 1E5, maximum = TRUE)
    k <- res$maximum
  }
  browser()
  
  # Fisher Information Matrix is the variance of the gradient of the log-likelihood
  # And in MLE is equal to the expected double derivative of the log-likelihood
  # But I'm missing something because the covariance of the gradients gets scaled by k twice while the average hessian is scaled by k only once.
  # I think I need to project the gradients/hessians to the tangent of the parameter space: i.e. tangent to p1, qe1, qs1, and *somehow* Omega's special orthogonality constraints
  grads <- matrix(obj_tape_long$Jacobian(vec_om), byrow = TRUE, ncol = obj_tape_long$domain) * k #each row is the gradient at a data point
  FisherI <- stats::cov(grads) #also called the 'variablility matrix'
  
  # Sensitivity Matrix
  # E(d^2(ll)/dtheta^2) under certain regularity conditions (passing derivatives outside an intergal) will be equal to d^2(E[ll])/dtheta^2
  # Here, I dont assume the regularity conditions:
  jactape <- scorematchingad::tape_Jacobian(obj_tape_long) #rowwise fill of gradient of each data point
  allhess <- matrix(jactape$Jacobian(vec_om), byrow = TRUE, ncol = obj_tape_long$domain^2)
  sensitivitymat <- -matrix(colMeans(allhess), nrow = obj_tape_long$domain, ncol = obj_tape_long$domain) * k
  
  # method assuming that E(d^2(ll)/dtheta^2) = d^2(E[ll])/dtheta^2
  obj_tape <- scorematchingad::avgrange(obj_tape_long) #Average of mu.y.
  sensitivitymatb <- -matrix(obj_tape$Hessian0(vec_om), ncol = obj_tape$domain, nrow = obj_tape$domain) * k
  
  es <- eigen(sensitivitymat)
  
}