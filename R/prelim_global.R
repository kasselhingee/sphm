prelim_global <- function(y, xs = NULL, xe = NULL, paramobj0, type = "Kassel", ...){ #paramobj0 is the starting parameter object
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
  obj_tape <- tape_namedfun("prelimobj_cpp", vec_om0, vector(mode = "numeric"), dims_in, cbind(y,xs,xe), check_for_nan = FALSE)
  constraint_tape <- tape_namedfun("Omega_constraints_wrap", vec_om0, vector(mode = "numeric"), dims_in, matrix(nrow = 0, ncol = 0), check_for_nan = FALSE)
  ineqconstraint_tape <- tape_namedfun("Omega_ineqconstraints", vec_om0, vector(mode = "numeric"), dims_in, matrix(nrow = 0, ncol = 0), check_for_nan = FALSE)
  
  if (type == "Shogo"){
    browser()
    omfixed <- lapply(om0, function(x) x * 0)
    omfixed$qe1 <- omfixed$qe1 + 1
    omfixed$ce1 <- omfixed$ce1 + 1
    isfixed <- mnlink_Omega_vec(as_mnlink_Omega(omfixed)) > 0.5
    obj_tape <- scorematchingad::fixindependent(obj_tape, vec_om0, isfixed)
    constraint_tape <- scorematchingad::fixindependent(constraint_tape, vec_om0, isfixed)
    ineqconstraint_tape <- scorematchingad::fixindependent(ineqconstraint_tape, vec_om0, isfixed)
    vec_om0 <- vec_om0[!isfixed]
  }

  # prepare nloptr options
  default_opts <- list(algorithm = "NLOPT_GN_ISRES", #the only algorthim that natively handles non-linear equality constraints - all the others have to use augmented Lagrangian ideas.
                xtol_rel = 1E-04,
                tol_constraints_eq = rep(1E-1, 1 + (length(om0$qs1) > 0) + (length(om0$qe1) > 0)),
                print_level = 3,
                maxeval = 1E4)
  ellipsis_args <- list(...)
  combined_opts <- utils::modifyList(default_opts, ellipsis_args)
  
  globopt <- nloptr::nloptr(
    x0 = vec_om0,
    eval_f = function(theta){obj_tape$eval(theta, vector(mode = "numeric"))},
    eval_g_eq =  function(theta){constraint_tape$eval(theta, vector(mode = "numeric"))},
    eval_g_ineq =  function(theta){ineqconstraint_tape$eval(theta, vector(mode = "numeric")) - 2},
    lb = vec_om0 * 0 - 10, #10 is just a guess here. For the spherical covariate stuff, I suspect most values are well below 1. *Euc will be different*
    ub = vec_om0 * 0 + 10,
    opts = combined_opts
  )
  warning("lb and ub not properly set")
  
  if (!(globopt$status %in% 1:4)){warning(globopt$message)}
  
  return(list(
    solution = Omega_proj(mnlink_Omega_unvec(globopt$solution, p, length(om0$qe1), check = FALSE)),
    glob_nloptr = globopt
  ))
}
