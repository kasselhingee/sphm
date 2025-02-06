

prelim_ad <- function(y, xs = NULL, xe = NULL, paramobj0, type = "Kassel", ...){ #paramobj0 is the starting parameter object
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
  default_opts <- list(algorithm = "NLOPT_LD_SLSQP",
                xtol_rel = 1E-10, #1E-04,
                tol_constraints_eq = rep(1E-1, 1 + (length(om0$qs1) > 0) + ((type != "Shogo") && (length(om0$qe1) > 0))),
                check_derivatives = TRUE, check_derivatives_print = 'errors',
                # print_level = 3,
                maxeval = 1E4)
  ellipsis_args <- list(...)
  combined_opts <- utils::modifyList(default_opts, ellipsis_args)
  
  browser()
  locopt <- nloptr::nloptr(
    x0 = vec_om0,
    eval_f = function(theta){out <- obj_tape$eval(theta, vector(mode = "numeric")); print(out); out},
    eval_grad_f = function(theta){out <- obj_tape$Jac(theta, vector(mode = "numeric")); print(out); out},
    eval_g_eq =  function(theta){out <- constraint_tape$eval(theta, vector(mode = "numeric"))[-3]; print(out); out},
    eval_jac_g_eq =  function(theta){out <- matrix(constraint_tape$Jacobian(theta), byrow = TRUE, ncol = length(theta))[-3, ]; print(out); out},
    eval_g_ineq =  function(theta){out <- ineqconstraint_tape$eval(theta, vector(mode = "numeric")) - 2; print(out); out},
    eval_jac_g_ineq =  function(theta){out <- matrix(ineqconstraint_tape$Jacobian(theta), byrow = TRUE, ncol = length(theta)); print(out); out},
    opts = combined_opts
  )
  
  if (!(locopt$status %in% 1:4)){warning(locopt$message)}
  
  browser()
  return(list(
    solution = Omega_proj(mnlink_Omega_unvec(locopt$solution, p, length(om0$qe1), check = FALSE)),
    loc_nloptr = locopt
  ))
}


