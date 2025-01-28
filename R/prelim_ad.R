

prelim_ad <- function(y, xs = NULL, xe = NULL, paramobj0, ...){ #paramobj0 is the starting parameter object
  om0 <- as_mnlink_Omega(paramobj0)
  # check inputs:
  mnlink_Omega_check(om0)
  p <- ncol(y)
  stopifnot(p == length(om0$p1))
  if (!is.null(xs)){
    stopifnot(ncol(xs) == length(om0$qs1))
  } else {
    stopifnot(length(om0$qs1) == 0)
  }
  if (!is.null(xe)){
    stopifnot(ncol(xe) == length(om0$qe1))
  } else {
    stopifnot(length(om0$qe1) == 0)
  }

  dims_in <- c(p, length(om0$qe1))
  obj_tape <- tape_namedfun("prelimobj_cpp", mnlink_Omega_vec(om0), vector(mode = "numeric"), dims_in, cbind(y,xs,xe), check_for_nan = FALSE)
  constraint_tape <- tape_namedfun("Omega_constraints_wrap", mnlink_Omega_vec(om0), vector(mode = "numeric"), dims_in, matrix(nrow = 0, ncol = 0), check_for_nan = FALSE)
  ineqconstraint_tape <- tape_namedfun("Omega_constraints_wrap", mnlink_Omega_vec(om0), vector(mode = "numeric"), dims_in, matrix(nrow = 0, ncol = 0), check_for_nan = FALSE)

  # prepare nloptr options
  default_opts <- list(algorithm = "NLOPT_LD_SLSQP",
                xtol_rel = 1E-10, #1E-04,
                tol_constraints_eq = rep(1E-1, 1 + (length(om0$qs1) > 0) + (length(om0$qe1) > 0)),
                maxeval = 1E4)
  ellipsis_args <- list(...)
  combined_opts <- utils::modifyList(default_opts, ellipsis_args)
  
  locopt <- nloptr::nloptr(
    x0 = mnlink_Omega_vec(om0),
    eval_f = function(theta){obj_tape$eval(theta, vector(mode = "numeric"))},
    eval_grad_f = function(theta){obj_tape$Jac(theta, vector(mode = "numeric"))},
    eval_g_eq =  function(theta){constraint_tape$eval(theta, vector(mode = "numeric"))},
    eval_jac_g_eq =  function(theta){matrix(constraint_tape$Jacobian(theta), byrow = TRUE, ncol = length(theta))},
    eval_g_ineq =  function(theta){ineqconstraint_tape$eval(theta, vector(mode = "numeric")) - 2},
    eval_jac_g_ineq =  function(theta){matrix(ineqconstraint_tape$Jacobian(theta), byrow = TRUE, ncol = length(theta))},
    opts = combined_opts
  )
  
  return(list(
    solution = Omega_proj(mnlink_Omega_unvec(locopt$solution, p, length(om0$qe1), check = FALSE)),
    loc_nloptr = locopt
  ))
}


