#' @title Objective function for vMF regression
#' @description
#' Internal objective function used by [`mobius_vMF()`].
#' For given parameters of [`mobius_link()`], computes \deqn{\frac{1}{n}\sum_{i=1}^n y_i^T \mu(x_i)}
#' where \eqn{y_i} are observed unit vectors and \eqn{\mu(x_i)} is the corresponding predicted unit vector.
#' Maximising this is equivalent to maximising the vMF log-likelihood over the mean link parameters.
#' @inheritParams mobius_link
#' @keywords internal
prelimobj <- function(y, xs = NULL, xe = NULL, param){
  predictedmeans <- mobius_link(xs = xs, xe = xe, param = param, check = FALSE)
  stopifnot(nrow(y) == nrow(predictedmeans))
  stopifnot(ncol(y) == ncol(predictedmeans))
  return(-mean(rowSums(y * predictedmeans)))
}

#' @title Regression with vMF Error and Scaled Mobius Link
#' @importClassesFrom scorematchingad Rcpp_ADFun
#' @description Regression of spherical response with vMF error and scaled Mobius link [`mobius_link()`].
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
#' @param type The mean link form. `"SpEuc"` and `"LinEuc"` are the general links (see
#'   [`mobius_link_params`]). `"Prop4ii"` constrains `p = qs` and `Bs = I_(p-1)`, and
#'   `"Prop4i"` constrains `Bs = beta_s I_(p-1)` with `0 < beta_s <= 1` --- the submodels of
#'   Proposition 4 of the manuscript, fitted by dedicated optimisers that parameterise the
#'   constrained manifolds directly (see [`mobius_link_prop4i`], [`mobius_link_prop4ii`]).
#'   Whether `xe` is supplied selects the variant with or without Euclidean covariates, and
#'   `fix_qs1`/`fix_qe1`/`lb`/`ub` are ignored for these two. Note that the Proposition 4
#'   fitters use central finite-difference gradients rather than the C++ automatic
#'   differentiation used by `"SpEuc"`/`"LinEuc"`, so they are slower.
#' @param det_constraint Only used when `type` is `"Prop4i"` or `"Prop4ii"`. Whether
#'   `Rtilde0` is constrained to be orthogonal (both determinant signs are tried and the
#'   better fit kept) or to be a rotation (determinant `+1`).
#' @param prop4_algorithm Only used when `type` is `"Prop4i"` or `"Prop4ii"`. The nloptr
#'   algorithm for the constrained fit.
#' @param prop4_opts Only used when `type` is `"Prop4i"` or `"Prop4ii"`. A list of nloptr
#'   options; anything passed through `...` is merged into it.
#' @inheritParams mobius_SvMF
#' @family regression
#' @export
mobius_vMF <- function(y, xs = NULL, xe = NULL, start = NULL, type = "SpEuc", fix_qs1 = FALSE, fix_qe1 = (type == "LinEuc"), intercept = TRUE, lb = NULL, ub = NULL,
                       det_constraint = c("orthogonal", "rotation"),
                       prop4_algorithm = "NLOPT_LD_SLSQP", prop4_opts = list(), ...){
  # The constrained links of Proposition 4 are fitted by their own optimisers, which
  # parameterise the constrained manifolds directly rather than through Omega.
  # Whether Euclidean covariates are present selects the variant.
  if (.is_prop4i_type(type) || .is_prop4ii_type(type)){
    dots <- list(...)
    if (length(dots)){prop4_opts <- utils::modifyList(prop4_opts, dots)}
    fitter <- if (.is_prop4i_type(type)){
      if (.is_xe_nonempty(xe)) mobius_vMF_prop4i_euc else mobius_vMF_prop4i
    } else {
      if (.is_xe_nonempty(xe)) mobius_vMF_prop4ii_euc else mobius_vMF_prop4ii
    }
    # Only pass what each fitter takes: the spherical-only fitters have no `intercept`
    # (there is no Euclidean term to centre against), and mobius_vMF_prop4ii has no
    # optimiser settings because it is solved in closed form by Procrustes rotation.
    args <- list(y = y, xs = xs, xe = xe, start = start,
                 det_constraint = det_constraint, intercept = intercept,
                 algorithm = prop4_algorithm, opts = prop4_opts)
    return(do.call(fitter, args[names(args) %in% names(formals(fitter))]))
  }
  mobius_vMF_general(y = y, xs = xs, xe = xe, start = start, type = type,
                     fix_qs1 = fix_qs1, fix_qe1 = fix_qe1, intercept = intercept,
                     lb = lb, ub = ub, ...)
}

#' @noRd
#' @title Unconstrained Mobius vMF regression
#' @description The `SpEuc`/`LinEuc` fit. [`mobius_vMF()`] dispatches here for every
#' `type` other than the constrained links of Proposition 4.
mobius_vMF_general <- function(y, xs = NULL, xe = NULL, start = NULL, type = "SpEuc", fix_qs1 = FALSE, fix_qe1 = (type == "LinEuc"), intercept = TRUE, lb = NULL, ub = NULL, ...){
  p <- ncol(y)
  preplist <- list(y = y, xs = xs, xe = xe, start = start)
  # For LinEuc: prepend a dummy-zero column to xe so the first covariate direction is always
  # fixed at the north pole (required by the LinEuc parameterisation); also append an
  # intercept column of ones if intercept = TRUE. Updates start to match the augmented xe.
  preplist <- addEuccovars(preplist, type = type, intercept = intercept)
  # Rotate y so its empirical mean direction aligns with the north pole; center and whiten xe
  # (if intercept = TRUE). Update start to match the new coordinates. This makes the default
  # starting point (near-identity matrices) a sensible initial guess.
  preplist <- standardise_data(preplist, intercept)
  # If start not supplied, choose start close to identities since data standardised
  preplist <- defaultstart(preplist, type)

  # For LinEuc: verify that qe1 equals the north pole (is_LinEuc) and that the dummy column
  # of xe is numerically zero, both required by the LinEuc parameterisation.
  if ((type == "LinEuc") && (!is.null(preplist$xe))){
    stopifnot(is_LinEuc(preplist$start))
    stopifnot(all(preplist$xe[, 1]^2 < sqrt(.Machine$double.eps)))
  }

  ### More detailed preparation ###
  om0 <- as_mobius_link_Omega(preplist$start)
  # check inputs:
  check_meanlink(preplist$y, preplist$xs, preplist$xe, om0)

  # Set up constraints: tape, fixed-parameter mask, and (possibly perturbed) starting point
  conprep <- estprep_meanconstraints(om0, fix_qs1, fix_qe1)
  # conprep$om0vec          - starting parameters as a vector: c(p1, qs1, qe1, Omega, ce)
  # conprep$isfixed         - logical mask: which positions in om0vec are held fixed
  # conprep$x0              - free-parameter starting values (om0vec[!isfixed]), possibly
  #                           perturbed to avoid a singular constraint Jacobian
  # conprep$constraint_tape - AD tape for orthonormality constraints on P, Qs, Qe
  # Rebuild om0vec incorporating any perturbation in conprep$x0. om0vec then serves as the
  # fixed-value carrier in the later t_sfi2u call: t_sfi2u(nlopt$solution, om0vec, isfixed)
  # reconstructs the full parameter vector after optimisation using om0vec's fixed positions.
  om0vec <- scorematchingad:::t_sfi2u(conprep$x0, conprep$om0vec, conprep$isfixed)

  # Prepare objective tape.
  # tape_namedfun records prelimobj_cpp under CppAD, producing a tape that evaluates the
  # objective and its gradient at any parameter value via $forward(0, x) and $Jacobian(x).
  # om0vec is the independent (differentiable) variable; the data (y, xs, xe) and dimensions
  # (p, qe) are baked in as constants.
  objtape <- tape_namedfun("prelimobj_cpp", om0vec, vector(mode = "numeric"), c(p, length(om0$qe1)), cbind(preplist$y,preplist$xs,preplist$xe), check_for_nan = FALSE)
  objtape <- scorematchingad::avgrange(objtape) #objtape initially returns a value for each measurement. Average here to get average over all data.
  
  # Re-tape the objective with fixed parameters baked in as constants, reducing the tape's
  # domain from length(om0vec) to sum(!isfixed) (free parameters only).
  objtape <- scorematchingad::fixindependent(objtape, objtape$xtape, conprep$isfixed)
  # Use the taping values as the starting point for optimisation. These are the (possibly
  # perturbed) starting values baked into the tape by estprep_meanconstraints.
  x0 <- objtape$xtape
  
  # prepare nloptr options
  # NLOPT_LD_SLSQP: Sequential Least Squares Programming; gradient-based constrained
  # optimiser suited here because the AD tapes supply exact gradients at low cost.
  # tol_constraints_eq = 1E-1: constraints (orthonormality) are satisfied only approximately
  # during optimisation; exact orthonormality is restored afterwards by Omega_proj.
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
  # nloptr minimises; the objective sum_i y_i . mu(x_i) is to be maximised, hence the negation.
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
  
  # Estimate concentration k by a separate 1D search.
  # Given the optimised mean link, the per-observation vMF log-likelihood is
  #   -log_norm_const(k, p) + k * mean(y . ypred).
  # The first term depends only on k and p; the second is k times the already-maximised
  # objective (-nlopt$objective = mean(y . ypred)). So k can be found by a scalar optimise
  # over k, treating mean(y . ypred) as fixed.
  res <- optimise(function(k){
    -vMF_log_norm_const_exact(k, p) + k * (-nlopt$objective) #full vMF log-likelihood (standardised by number of observations)
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
 
  # SLSQP satisfies orthonormality constraints only up to tol_constraints_eq = 1E-1;
  # Omega_proj projects P, Qs, Qe back onto the exact Stiefel manifold.
  projectedom <- Omega_proj(mobius_link_Omega_unvec(fullparam, p, length(om0$qe1), check = FALSE))
  try({mobius_link_Omega_check(projectedom)})

  ### Making nicer return objects ###
  # Residuals and distances are computed in the standardised coordinate system (before
  # reverting), since they are invariant to the choice of coordinates.
  pred <- mobius_link(xs = preplist$xs, xe = preplist$xe, param = projectedom)
  dists <- acos(rowSums(pred * preplist$y))
  # rotated_resid parallel-transports each observation to the north pole. The returned matrix
  # has p columns: the first is the (normalised) response direction (always ~1 for small
  # residuals). [, -1] drops it, keeping only the p-1 tangent components.
  rresids_tmp <- rotated_resid(preplist$y, pred, north_pole(ncol(preplist$y)))
  rresids <- rresids_tmp[, -1]
  attr(rresids, "samehemisphere") <-  attr(rresids_tmp, "samehemisphere")
  colnames(rresids) <- paste0("r", 1:ncol(rresids))
  
  # Invert the rotations and centering applied by standardise_data, expressing the estimated
  # parameters in the original user-supplied coordinates.
  est <- undo_recoordinate_Omega(projectedom,
                          yrot = attr(preplist$y, "std_rotation"),
                          xsrot = attr(preplist$xs, "std_rotation"), #if xs/xe is NULL then attr(xs/xe, ..) is NULL too
                          xerot = attr(preplist$xe, "std_rotation"),
                          xecenter = attr(preplist$xe, "std_center"),
                          onescovaridx = preplist$onescovaridx)
  # The Euclidean link has an inherent sign ambiguity: replacing (ce, Qe, Be) with
  # (-ce, -Qe, -Be) leaves the mean link unchanged. Convention is ce >= 0; Euc_signswitch
  # flips all three if ce < 0.
  if (isTRUE(est$ce < 0)){est <- Euc_signswitch(est)}
  
  # DoF
  DoF <- mobius_dof(p, length(est$qs1), length(est$qe1), fix_qs1 = fix_qs1, fix_qe1 = fix_qe1) + 
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
    # destandardise_Euc inverts the centering and whitening from standardise_data; using
    # preplist$xe also recovers any columns added by addEuccovars (dummy-zero and intercept).
    xe = if (!is.null(xe)){if (intercept){destandardise_Euc(preplist$xe, attr(preplist$xe, "std_center"), attr(preplist$xe, "std_rotation"))} else {xe}},
    pred = destandardise_sph(pred, rotation = attr(preplist$y, "std_rotation")),
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
#' @description The parameters of the Mobius mean link function [`mobius_link()`] have a number of constraints.
#' This function incorporates these constraints to obtain the total degrees of freedom of the parameters.
#' @param p The length of response vectors
#' @param qs The length of spherical covariate vectors
#' @param qe The length of Euclidean covariate vectors
#' @param fix_qs1 Whether the first column of `Qs` is fixed (i.e. not estimated and not free).
#' @param fix_qe1 Whether `ce` and the first column of `Qe` is fixed (i.e. not estimated and not free).
#' @return An integer
#' @family regression
#' @export
mobius_dof <- function(p, qs = 0, qe = 0, fix_qs1 = FALSE, fix_qe1 = FALSE){
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
  try(mobius_link_Omega_check(om0))
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

#' @title Prepare mean-link constraint setup for optimisation
#' @description
#' Prepares the constraint tape, fixed-parameter mask, and starting point for constrained
#' optimisation of the Mobius mean link parameters in [`mobius_vMF()`].
#' Specifically, the function:
#' 1. Vectorises `om0` to `om0vec`, the starting parameters as `c(p1, qs1, qe1, Omega, ce)`.
#' 2. Creates an AD constraint tape evaluating the orthonormality constraints on `P`, `Qs`,
#'    `Qe` and their Jacobians.
#' 3. Determines `isfixed` — a logical mask of which positions in `om0vec` are held constant,
#'    based on `fix_qs1` and `fix_qe1`.
#' 4. Reduces the tape via `fixindependent` and `keeprange`, dropping constraint outputs that
#'    become constant once fixed parameters are removed.
#' 5. Computes `x0 = om0vec[!isfixed]`, the free-parameter starting vector.
#' 6. Checks the constraint Jacobian at `x0` for singularity (near-zero singular values).
#' 7. If singular, perturbs `x0` by adding small noise (from a pregenerated matrix in
#'    `R/sysdata.rda`) to the off-diagonal columns of `Qs` and `Qe`, then recomputes `x0`.
#' @param om0 Starting [`mobius_link_Omega`] object.
#' @param fix_qs1 Logical; whether to hold the first column of `Qs` fixed.
#' @param fix_qe1 Logical; whether to hold `ce` and the first column of `Qe` fixed.
#' @return A list:
#' \describe{
#'   \item{om0vec}{Starting parameters as a flat vector: `c(p1, qs1, qe1, Omega, ce)`.}
#'   \item{isfixed}{Logical vector marking fixed positions in `om0vec`.}
#'   \item{x0}{Starting values for the free (non-fixed) parameters, possibly perturbed to
#'     avoid a singular constraint Jacobian.}
#'   \item{constraint_tape}{AD tape evaluating orthonormality constraints on `P`, `Qs`, `Qe`
#'     and their Jacobians, with fixed-parameter dimensions already reduced.}
#' }
#' @keywords internal
estprep_meanconstraints <- function(om0, fix_qs1, fix_qe1){
  dims_in <- c(length(om0$p1), length(om0$qe1))
  om0vec <- mobius_link_Omega_vec(om0)
 
  # generate tapes
  # tape_namedfun records Omega_constraints_wrap under CppAD: evaluates the orthonormality
  # constraints and their Jacobian at any om0vec via $forward(0, x) and $Jacobian(x).
  constraint_tape <- tape_namedfun("Omega_constraints_wrap", om0vec, vector(mode = "numeric"), dims_in, matrix(nrow = 0, ncol = 0), check_for_nan = FALSE)

  # fix mean link parameters depending on arguments
  # use the starting parameters om0 to detect whether we have xs and xe as their form is more predictable due to the Omega class
  omfixed <- lapply(om0, function(x) x * 0)
  if (fix_qe1 && (length(om0$qe1) > 0)){
    # if qe1 is fixed and qe1 is a non-empty vector, then record the corresponding omfixed elements as 1
    omfixed$qe1 <- omfixed$qe1 + 1
    omfixed$ce <- omfixed$ce + 1
  }
  
  # fix qs1 to the starting values if qs1 exists
  # use the starting parameters to check as their form is more constrained by the Omega class
  if (fix_qs1 && (length(om0$qs1) > 0)){
    omfixed$qs1 <- omfixed$qs1 + 1
  }
  
  # Vectorise the omfixed marker through the Omega class so the logical mask aligns exactly
  # with the positions in om0vec.
  isfixed <- mobius_link_Omega_vec(as_mobius_link_Omega(omfixed)) > 0.5
  # Bake the fixed-parameter values into the constraint tape, reducing its domain to free
  # parameters only (so the constraint Jacobian is over free parameters only).
  constraint_tape <- scorematchingad::fixindependent(constraint_tape, om0vec, isfixed)
  # After fixing parameters, some constraint equations become identically zero (e.g. the
  # orthonormality constraint on a fixed column is pre-satisfied). Drop those to keep the
  # constraint Jacobian full-rank.
  keep <- which(vapply((1:constraint_tape$range)-1, function(i){!constraint_tape$parameter(i)}, FUN.VALUE = FALSE))
  constraint_tape <- scorematchingad::keeprange(constraint_tape, keep)
  
  # check Jacobians of constraints are non-singular for the starting parameters.
  # For pathological params (e.g. the default starting params of no rotations), it can be zero.
  # If it is singular, perturb starting omega very slightly
  x0 <- om0vec[!isfixed]
  Jac_eq <- matrix(constraint_tape$Jacobian(x0), byrow = TRUE, ncol = constraint_tape$domain)
  if (any(abs(svd(Jac_eq)$d) < sqrt(.Machine$double.eps))){
    # modify Qs and Qe so that they dont individually form orthogonal vectors by adding a little bit of noise
    cann <- as_mobius_link_cann(om0)
    scalemodifier <- max(abs(rbind(cann$Qs, cann$Qe)))/1E3
    cann$Qs[,-1] <- cann$Qs[,-1] + scalemodifier*savednoisemat[1:nrow(cann$Qs), 1:(ncol(cann$Qs)-1)] #use pregenerated 100 x 100 noise matrix in R/sysdata.rda, built by data-raw/generate-savednoisemat.R
    cann$Qe[,-1] <- cann$Qe[,-1] + scalemodifier*savednoisemat[1:nrow(cann$Qe), 1:(ncol(cann$Qe)-1)] #use pregenerated 100 x 100 noise matrix in R/sysdata.rda, built by data-raw/generate-savednoisemat.R
    om0new <- as_mobius_link_Omega(cann)
    x0 <- mobius_link_Omega_vec(om0new)[!isfixed]
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
#' @description Given a vMF regression, repeat the optimisation from initial parameters randomly generated by [`rand_mobius_link_cann()`].
#' @param mod_vMF Result of [`mobius_vMF()`]
#' @param preseed Passed to [`rand_mobius_link_cann()`]
#' @details `fix_qe1` and `fix_qs1` of `mod_vMF` will be respected.
#' @return An vMF regression using the same data as `mod_vMF`.
#' @family regression
#' @export
mobius_vMF_refit <- function(mod_vMF, preseed = 1){
  .prop4_reject(mod_vMF, "mobius_vMF_refit")
  inparam <- as_mobius_link_Omega(mod_vMF$est)
  dims <- dim.mobius_link_Omega(inparam)
  start <- rand_mobius_link_cann(p = dims["p"], 
                           qs = dims[["qs"]] - (dims[["qs"]] > 0 & mod_vMF$linktype$fix_qs1), 
                           qe = dims[["qe"]] - (dims[["qe"]] > 0 & mod_vMF$linktype$fix_qe1), 
                           preseed = preseed)
  # Goal: for each fixed first column (qe1 or qs1), construct a random starting matrix whose
  # first column equals the fixed value and whose remaining columns are random in the
  # orthogonal complement. Strategy:
  #   1. rand_mobius_link_cann already generated a random matrix for the free columns only.
  #   2. Prepend a zero placeholder column to form the full p x qe (or p x qs) matrix.
  #   3. Apply the inverse parallel transport from the north pole to qe1 (or qs1), mapping
  #      the placeholder first column to qe1 and rotating the free columns accordingly.
  #   4. Overwrite the first column with the exact fixed value (guards against the edge case
  #      where qe1 = -north_pole and parallel transport is a reflection, not a rotation).
  if (dims[["qe"]] > 0 & mod_vMF$linktype$fix_qe1){
    start$ce <- inparam$ce
    rotmat <- parallel_transport_mat(inparam$qe1, north_pole(dims[["qe"]]))
    rotQe <- rbind(0, start$Qe)
    # Two cases depending on whether rand_mobius_link_cann returned a square or tall matrix
    if (ncol(rotQe) == dims["p"]){
      rotQe <- cbind(0, rotQe[,-1, drop = FALSE])
    } else if (ncol(rotQe) == dims["p"] - 1) { #case when qe=p and fix_qe1=TRUE
      rotQe <- cbind(0, rotQe)
    }
    rotQe[1,1] <- 1
    start$Qe <- t(rotmat) %*% rotQe
    # Overwrite first column: parallel transport may be a reflection when qe1 = -north_pole
    start$Qe[,1] <- inparam$qe1
  }
  if (dims[["qs"]] > 0 & mod_vMF$linktype$fix_qs1){
    rotmat <- parallel_transport_mat(inparam$qs1, north_pole(dims[["qs"]]))
    rotQs <- rbind(0, start$Qs)
    if (ncol(rotQs) == dims["p"]){
      rotQs <- cbind(0, rotQs[,-1, drop = FALSE])
    } else if (ncol(rotQs) == dims["p"] - 1) { #case when qs=p and fix_qs1=TRUE
      rotQs <- cbind(0, rotQs)
    }
    rotQs[1,1] <- 1
    start$Qs <- t(rotmat) %*% rotQs
    # Overwrite first column: parallel transport may be a reflection when qs1 = -north_pole
    start$Qs[,1] <- inparam$qs1
  }
  start <- as_mobius_link_Omega(start)
  stopifnot(all(dim.mobius_link_Omega(start) == dims))
  mobius_vMF(y = mod_vMF$y,
             xs = mod_vMF$xs,
             xe = mod_vMF$xe,
             fix_qs1 = mod_vMF$linktype$fix_qs1,
             fix_qe1 = mod_vMF$linktype$fix_qe1,
             type = mod_vMF$linktype$type,
             intercept = mod_vMF$linktype$intercept,
             start = start)
}
