#' SvMF Regression with Scaled Mobius Link and Parallel Transported Axes
#' @description
#' Performs regression for spherical response data with a scaled Mobius mean link and scale von Mises-Fisher error distribution (SvMF).
#' @details
#' The mean is assumed to follow [`mobius_link()`].
#' The concentration and scaling in the SvMF is assumed constant across observations.
#' The symmetry axes of the SvMF at location \eqn{\mu} are assumed to be the result of parallel transport along the geodesic of axes `G0[,-1]` from  `G0[,1]` to \eqn{\mu}. The axis at this base location `G0[,1]` are estimated.
#'
#' Unless requested by `doprelim = FALSE`, a preliminary regression using a vMF error distribution is performed.
#'
#' When p!=3, the optimisation first approximates the normalising constant of the vMF distribution, then the concentration is optimised seperately.
#' @param y Response data on a sphere
#' @param xs Covariate data on a sphere
#' @param xe Covariate data in Euclidean Space
#' @param k Starting concentration. If omitted, then concentration of the preliminary vMF regression is used.
#' @param a The scaling vector `a`. `a[1]` is a fixed tuning parameter and the remainining is used as an initial value for the optimisation. If omitted, then a moment estimate using residuals of from the preliminary vMF regression is used.
#' @param mean Parameters for the mean link, used as a starting guess for the preliminary vMF regression if `doprelim = TRUE`.
#' @param G0 A `p x p` rotation matrix specifying the starting guess of the axes of the SvMF distribution. G0 should have positive determinant because in the estimation routine G0 or parts of G0 are representented using Cayley transforms. If omitted, then a moment estimate using residuals of from the preliminary vMF regression is used.
#' @param G0reference A `p x p` rotation matrix specifying a set of coordinates with which to represent G0 for the estimation. This is because internally `G0` or `G0[,-1]` is represented using a Cayley transform, which performs better the closer the matrix is to the identity. Ideally the columns of `G0reference` will be close to the best `G0`.
#' @param G01behaviour "p1" identifies `G0[,1]` with `p1`. "fixed" fixes `G0[,1]` to its initial value. "free" allows `G0[,1]` to be estimated freely.
#' @param doprelim When `FALSE` the preliminary von Mises-Fisher regression and subsequent moment estimation of `a` and `G0` is omitted. The provided parameters are used as the initial values for an optimisation of all parameters of the SvMF regression all together using [`nloptr::nloptr()`].
#' @param ... Named optional arguments passed as a list to the `opts` argument of [`nloptr::nloptr()`].
#' @param type Specify the link type ("SpEuc" or "LinEuc")
#' @param fix_qs1 Fix qs1 to the starting values.
#' @param fix_qe1 Fix qe1 to the starting values - typically used in when "LinEuc" link used.
#' @param intercept Include an intercept (Euclidean covariate that is always `1`).
#' @return A list:
#' \describe{
#'   \item{mean}{Estimated mean parameters}
#'   \item{k}{Estimated concentration parameter k}
#'   \item{a}{Estimated vector of scales `a` (except `a[1]` which is fixed)}
#'   \item{G0}{Estimated base location and symmetry axes.}
#'   \item{obj}{Final objective value from NLopt}
#'   \item{nlopt}{NLopt optimization result object}
#'   \item{y}{Response data}
#'   \item{xs}{Spherical covariate values}
#'   \item{xe}{Euclidean covariate values (including any automatically added columns like dummy zeros and an intercept)}
#'   \item{pred}{Predicted values on sphere}
#'   \item{rresids_I}{Residuals rotated (by parallel transport) to the north pole}
#'   \item{rresids_G0, rresids_std}{Residuals rotated to `G0[,1]` and expressed according to the axes given by `G0[,-1]`. `rresids_std` additionaly scales the rotated residuals by \eqn{\sqrt{k} a[1]/a[j]}.}
#'   \item{dists}{Residual geodesic distances (geodesic distance between predicted value and observed value)}
#'   \item{DoF}{Degrees of freedom of the optimisation}
#'   \item{AIC}{Akaike Information Criterion}
#'   \item{lLik}{Log-likelihood}
#'   \item{initial}{Initial values}
#'   \item{preest}{Results of the preliminary vMF regression}
#'   \item{linktype}{The mean link settings used for the fit: `type`, `fix_qs1`, `fix_qe1` and `intercept`. Refit helpers such as [`mobius_SvMF_signflip_refit()`] read these back so they refit the same model.}
#'   \item{G01behaviour}{The `G01behaviour` used for the fit}
#' }
#' @inheritParams mobius_vMF
#' @family regression
#' @export
mobius_SvMF <- function(y, xs = NULL, xe = NULL, mean = NULL, k = NULL, a = NULL, G0 = NULL, G0reference = NULL, G01behaviour = "p1", type = "LinEuc", fix_qs1 = FALSE, fix_qe1 = (type == "LinEuc"), intercept = TRUE, doprelim = TRUE,
                        det_constraint = c("orthogonal", "rotation"),
                        prop4_algorithm = "NLOPT_LD_SLSQP", prop4_opts = list(), ...){
  # See mobius_vMF() for the rationale; the constrained links of Proposition 4 have their
  # own two-phase estimators in Mobius_prop4i.R / Mobius_prop4ii.R.
  if (.is_prop4i_type(type) || .is_prop4ii_type(type)){
    fitter <- if (.is_prop4i_type(type)){
      if (.prop4_has_euc(xe)) mobius_SvMF_prop4i_euc else mobius_SvMF_prop4i
    } else {
      if (.prop4_has_euc(xe)) mobius_SvMF_prop4ii_euc else mobius_SvMF_prop4ii
    }
    # The spherical-only fitters take no `intercept` (there is no Euclidean term).
    args <- list(y = y, xs = xs, xe = xe, mean = mean, k = k, a = a,
                 G0 = G0, G0reference = G0reference, G01behaviour = G01behaviour,
                 doprelim = doprelim, det_constraint = det_constraint,
                 intercept = intercept, algorithm = prop4_algorithm,
                 opts = prop4_opts)
    args <- args[names(args) %in% names(formals(fitter))]
    return(do.call(fitter, c(args, list(...))))
  }
  mobius_SvMF_general(y = y, xs = xs, xe = xe, mean = mean, k = k, a = a, G0 = G0,
                      G0reference = G0reference, G01behaviour = G01behaviour, type = type,
                      fix_qs1 = fix_qs1, fix_qe1 = fix_qe1, intercept = intercept,
                      doprelim = doprelim, ...)
}

#' @noRd
#' @title Unconstrained Mobius SvMF regression
#' @description The `SpEuc`/`LinEuc` fit. [`mobius_SvMF()`] dispatches here for every
#' `type` other than the constrained links of Proposition 4.
mobius_SvMF_general <- function(y, xs = NULL, xe = NULL, mean = NULL, k = NULL, a = NULL, G0 = NULL, G0reference = NULL, G01behaviour = "p1", type = "LinEuc", fix_qs1 = FALSE, fix_qe1 = (type == "LinEuc"), intercept = TRUE, doprelim = TRUE, ...){

  # Two-phase estimation:
  # Phase 1 (prelim): vMF regression for the mean link, followed by moment estimation of G0
  # and a. This provides a warm start for the joint optimisation.
  # Phase 2 (finalest): joint maximum likelihood over all parameters (mean, k, a, G0).
  if (doprelim){
    preest <- mobius_SvMF_partransport_prelim(y, xs, xe,
                                              mean = mean,
                                              G0 = G0, G01behaviour = G01behaviour,
                                              type = type, fix_qs1 = fix_qs1, fix_qe1 = fix_qe1,
                                              intercept = intercept, ...)
  } else {
    # doprelim = FALSE: caller must supply all starting values directly.
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

  # Explicit caller-supplied k, a, G0reference take precedence over the prelim estimates.
  finalest <- mobius_SvMF_joint_fit(y, xs, xe,
                           mean = preest$mean,
                           k = if(!is.null(k)){k}else{preest$k},
                           a = if(!is.null(a)){a}else{preest$a},
                           G0 = preest$G0,
                           G0reference = if(!is.null(G0reference)){G0reference}else{preest$G0},
                           G01behaviour = G01behaviour,
                           type = type, fix_qs1 = fix_qs1, fix_qe1 = fix_qe1,
                           intercept = intercept, ...)
  # Record how the model was fitted alongside the estimates, mirroring mobius_vMF(). Refit
  # helpers read these back rather than requiring the caller to restate them, which would let
  # a refit silently optimise a different model than the one it was handed. G01behaviour is a
  # sibling of linktype rather than a member because it describes the error distribution
  # rather than the mean link.
  return(c(finalest, list(
    preest = preest,
    linktype = list(type = type, fix_qs1 = fix_qs1, fix_qe1 = fix_qe1, intercept = intercept),
    G01behaviour = G01behaviour)))
}

#' @title Sign-flip refit for a SvMF regression
#' @description Re-runs the SvMF optimisation from starting points obtained by flipping
#' the signs of the three "pole" directions (p1, qs1, qe1) in the fitted mean-link
#' parameters, analogous to \code{\link{mobius_vMF_signflip_refit}}.
#' Starting points are generated by \code{signflip_starts()}: sign-flips of p1, qs1, qe1,
#' plus an alternative \code{ce} value from the opposite size regime when xe is present and
#' \code{ce} is free (\code{fix_qe1 = FALSE}).
#' Each restart calls \code{\link{mobius_SvMF}} with \code{doprelim = FALSE}, so the
#' sign-flipped mean and the alternative ce are used directly as the joint-optimiser starting
#' values rather than being discarded by a vMF preliminary step.
#' The link type, \code{fix_qs1}, \code{fix_qe1}, \code{intercept} and \code{G01behaviour} are
#' taken from \code{mod_SvMF} rather than supplied as arguments, so the restarts always
#' optimise the same model as the fit being refitted.
#' @param mod_SvMF Result of \code{\link{mobius_SvMF}}. Must carry the \code{linktype} and
#'   \code{G01behaviour} elements that \code{\link{mobius_SvMF}} records.
#' @param xtol_rel Relative convergence tolerance passed to nloptr (default \code{1e-4},
#'   coarser than the default \code{1e-10} so enumeration is fast).
#' @param maxeval Maximum number of objective evaluations per restart (default 500).
#' @param ... Additional arguments passed to \code{\link{mobius_SvMF}}.
#' @return The SvMF fit with the lowest \code{obj} value across all restarts, or
#'   \code{NULL} if every restart errored.
#' @family regression
#' @export
mobius_SvMF_signflip_refit <- function(mod_SvMF, xtol_rel = 1e-4, maxeval = 500, ...) {
  .prop4_reject(mod_SvMF, "mobius_SvMF_signflip_refit")
  # Settings are read off the fit rather than taken as arguments, so a restart always
  # optimises the same model as mod_SvMF. Fail loudly rather than let NULL settings fall
  # through to mobius_SvMF()'s defaults: results of mobius_SvMF_joint_fit(), and fits cached
  # before these fields were recorded, have no linktype.
  if (is.null(mod_SvMF$linktype)) {
    stop("mod_SvMF has no 'linktype' element: it was not produced by mobius_SvMF() (or predates the recording of fit settings). Refit with mobius_SvMF() first.")
  }
  lt     <- mod_SvMF$linktype
  cann0  <- as_mobius_link_cann(mod_SvMF$mean)
  starts <- signflip_starts(cann0, lt$fix_qs1, lt$fix_qe1)

  best_fit <- NULL
  best_obj <- Inf
  for (cann_trial in starts) {
    # suppressWarnings: warnings from the sign-flip fits are not informative to the caller.
    # tryCatch: a sign-flipped start can be numerically ill-conditioned; skip those silently.
    # doprelim = FALSE: bypass vMF prelim so the modified starting ce (and sign-flipped poles)
    # feed directly into the joint SvMF optimizer. With doprelim = TRUE the vMF prelim would
    # re-estimate ce from data, negating the ce variation in signflip_starts.
    fit <- tryCatch(
      suppressWarnings(mobius_SvMF(
        y            = mod_SvMF$y,
        xs           = mod_SvMF$xs,
        xe           = mod_SvMF$xe,
        mean         = cann_trial,
        k            = mod_SvMF$k,
        a            = mod_SvMF$a,
        G0           = mod_SvMF$G0,
        G0reference  = mod_SvMF$G0,
        G01behaviour = mod_SvMF$G01behaviour,
        type         = lt$type,
        fix_qs1      = lt$fix_qs1,
        fix_qe1      = lt$fix_qe1,
        intercept    = lt$intercept,
        doprelim     = FALSE,
        xtol_rel     = xtol_rel,
        maxeval      = maxeval,
        ...
      )),
      error = function(e) NULL
    )
    # nloptr minimises -log-lik, so a lower obj = better fit.
    if (!is.null(fit) && fit$obj < best_obj) {
      best_obj <- fit$obj
      best_fit <- fit
    }
  }
  return(best_fit)
}

#' @title SvMF regression with sign-flip multistart exploration
#' @description Fits the SvMF regression in three phases:
#' (1) \code{\link{mobius_vMF_multistart}} finds a reliable basin for the mean link
#'     then a coarse SvMF fit is run from that starting point;
#' (2) sign-flip restarts (see \code{\link{mobius_SvMF_signflip_refit}}) at the same
#'     coarse tolerance to explore other basins again;
#' (3) a tight final SvMF optimisation from the best starting point found.
#' @param y,xs,xe,G0,G0reference,G01behaviour,type,fix_qs1,fix_qe1,intercept,lb,ub
#'   As in \code{\link{mobius_SvMF}}.
#' @param xtol_rel Relative tolerance for phases 1 and 2 (default \code{1e-4}).
#' @param maxeval Maximum evaluations per fit in phases 1 and 2 (default \code{500}).
#' @param ... Additional arguments forwarded to every internal \code{\link{mobius_SvMF}} and
#'   \code{\link{mobius_vMF_multistart}} call.
#' @return Same structure as \code{\link{mobius_SvMF}}.
#' @family regression
#' @export
mobius_SvMF_multistart <- function(y, xs = NULL, xe = NULL,
                                    G0 = NULL, G0reference = NULL,
                                    G01behaviour = "p1",
                                    type = "LinEuc", fix_qs1 = FALSE,
                                    fix_qe1 = (type == "LinEuc"),
                                    intercept = TRUE, lb = NULL, ub = NULL,
                                    xtol_rel = 1e-4, maxeval = 500, ...) {
  .prop4_reject(type, "mobius_SvMF_multistart")
  # Phase 1: vMF multistart to find a reliable starting basin for the mean link,
  # then coarse SvMF from that starting point.
  # A plain coarse SvMF from a default start can settle in a low-k local optimum;
  # vMF multistart navigates a simpler landscape (mean only, fewer params) and is
  # more reliable at finding the globally best basin for the mean link.
  vMF_mod <- suppressWarnings(mobius_vMF_multistart(y, xs = xs, xe = xe,
                                                     type = type,
                                                     fix_qs1 = fix_qs1, fix_qe1 = fix_qe1,
                                                     intercept = intercept,
                                                     xtol_rel = xtol_rel, maxeval = maxeval, ...))
  mod0 <- suppressWarnings(mobius_SvMF(y, xs = xs, xe = xe,
                                        mean = vMF_mod$est,
                                        G0 = G0, G0reference = G0reference,
                                        G01behaviour = G01behaviour,
                                        type = type, fix_qs1 = fix_qs1, fix_qe1 = fix_qe1,
                                        intercept = intercept, lb = lb, ub = ub,
                                        xtol_rel = xtol_rel, maxeval = maxeval, ...))

  # Phase 2: sign-flip sweep from the coarse fit's landing point.
  mod1 <- mobius_SvMF_signflip_refit(mod0, xtol_rel = xtol_rel, maxeval = maxeval, ...)

  best <- if (!is.null(mod1) && mod1$obj < mod0$obj) mod1 else mod0

  # Phase 3: tight final SvMF from the best basin found in phases 1–2.
  # OPEN QUESTION (was an interactive browser() here, removed because it makes
  # testthat::test_local() and any non-interactive run halt): does best$G0[,1] ever differ
  # from best$mean$p1 when G01behaviour == "p1"? To investigate, re-insert
  #   if (G01behaviour == "p1") browser()
  # here, or log max(abs(best$G0[,1] - best$mean$p1)).
  # Slight re-orthogonalise G0; for G01behaviour == "p1" pin column 1 to p1 exactly
  # to satisfy the constraint that G0[,1] must equal p1 exactly.
  first_col <- if (G01behaviour == "p1") best$mean$p1 else {
    v <- best$G0[, 1L]; v / sqrt(sum(v^2))
  }
  G0_rest <- best$G0[, -1L, drop = FALSE]
  G0_rest <- G0_rest - first_col %*% (t(first_col) %*% G0_rest)
  best_G0 <- as_rotation_mat(cbind(first_col, qr.Q(qr(G0_rest))))
  mobius_SvMF(y, xs = xs, xe = xe,
              mean = best$mean, k = best$k, a = best$a, G0 = best_G0,
              G0reference = if (!is.null(G0reference)) G0reference else best_G0,
              G01behaviour = G01behaviour,
              type = type, fix_qs1 = fix_qs1, fix_qe1 = fix_qe1,
              intercept = intercept, lb = lb, ub = ub,
              doprelim = FALSE, ...)
}

# Joint maximum likelihood optimisation of all SvMF regression parameters for fixed
# scale structure (a and G0 are estimated, but the constraint prod(a[-1]) = 1 is enforced).
# Called by mobius_SvMF() after the preliminary vMF estimate. Uses automatic differentiation
# (CppAD via scorematchingad) and a gradient-based solver (nloptr/SLSQP).
mobius_SvMF_joint_fit <- function(y, xs, xe, mean, k, a, G0 = NULL, G0reference = NULL, G01behaviour = "p1", type = "LinEuc", fix_qs1 = FALSE, fix_qe1 = (type == "LinEuc"), intercept = TRUE, lb = NULL, ub = NULL, ...){
  initial <-  list(
    mean = mean,
    k = k, 
    a = a, 
    G0 = G0
  )
  SvMF_cann_check(SvMF_cann(k = initial$k, a = initial$a, G = initial$G0))

  p <- ncol(y)
  preplist <- list(y = y, xs = xs, xe = xe, start = mean)
  # if needed, add Euclidean covariates and update start accordingly
  preplist <- addEuccovars(preplist, type = type, intercept = intercept)
  # standardise y, xe and xe and update start accordingly. Dont standardise xe if intercept = FALSE
  preplist <- standardise_data(preplist, intercept)
  # G0 lives in the same space as y, so it must be rotated by the same std_rotation used on y.
  if (!is.null(G0)){preplist$G0 <- attr(preplist$y, "std_rotation") %*% G0}
  if (!is.null(G0reference)){preplist$G0reference <- attr(preplist$y, "std_rotation") %*% G0reference}
  # If start not supplied, choose start close to identities since data standardised
  preplist <- defaultstart(preplist, type)
  # Caller-supplied a and k take precedence over any defaults set by defaultstart.
  if (!is.null(a)){preplist$a <- a}
  if (!is.null(k)){preplist$k <- k}

  # Check LinEuc link initiated properly
  if ((type == "LinEuc") && (!is.null(preplist$xe))){
    stopifnot(is_LinEuc(preplist$start))
    stopifnot(all(preplist$xe[, 1]^2 < sqrt(.Machine$double.eps)))
  }

  ### More detailed preparation ###
  om0 <- as_mobius_link_Omega(preplist$start)
  # check inputs:
  check_meanlink(preplist$y, preplist$xs, preplist$xe, om0)
  stopifnot(length(preplist$a) == p)
  a1 <- preplist$a[1]
  aremaining <- preplist$a[-1]
  stopifnot(isTRUE(all.equal(prod(aremaining), 1)))
  if (is.null(preplist$G0reference)){
    G0reference <- preplist$G0
  } else {
    G0reference <- preplist$G0reference
  }
   
  qs <- length(om0$qs1)
  qe <- length(om0$qe1)


  # Set up constraints: tape, fixed-parameter mask, and (possibly perturbed) starting point
  conprep <- estprep_meanconstraints(om0, fix_qs1, fix_qe1)
  # conprep$om0vec          - starting parameters as a vector: c(p1, qs1, qe1, Omega, ce)
  # conprep$isfixed         - logical mask: which positions in om0vec are held fixed
  # conprep$x0              - free-parameter starting values (om0vec[!isfixed]), possibly
  #                           perturbed to avoid a singular constraint Jacobian
  # conprep$constraint_tape - AD tape for orthonormality constraints on P, Qs, Qe
  # Rebuild om0vec incorporating any perturbation in conprep$x0. om0vec then serves as the
  # fixed-value carrier in the later t_sfi2u call: t_sfi2u(meanpars, om0vec, isfixed)
  # reconstructs the full mean-link parameter vector after optimisation.
  om0vec <- scorematchingad:::t_sfi2u(conprep$x0, conprep$om0vec, conprep$isfixed)
  
  # Prepare objective tape.
  # tape_ld_mobius_SvMF_partransport_nota1 records the full SvMF log-density under CppAD.
  # Independent variables (differentiated w.r.t.): free mean-link parameters from om0vec,
  # followed by k, aremaining, and the Cayley-transform coordinates of G0 — all in one flat
  # vector. The Cayley transform maps the rotation manifold to an unconstrained space, so G0
  # requires no extra constraints in the optimiser.
  # Constants baked in: data (y, xs, xe), dimensions, G0reference, G01behaviour.
  objtape_ind <- tape_ld_mobius_SvMF_partransport_nota1(omvec = om0vec, k = preplist$k,
                                       a1 = a1,
                                       aremaining = aremaining,
                                       G0 = preplist$G0,
                                       p = p, qe = qe,
                                       yx = cbind(preplist$y, preplist$xs, preplist$xe),
                                       referencecoords = G0reference,
                                       G01behaviour = G01behaviour)
  # Re-tape with fixed mean-link parameters baked in. The isfixed vector is padded with zeros
  # for the non-mean-link components (k, aremaining, Cayley G0) since those are always free.
  objtape_ind <- scorematchingad::fixindependent(objtape_ind, objtape_ind$xtape, c(conprep$isfixed, rep(0, length(objtape_ind$xtape) - length(conprep$isfixed))))

  # objtape initially returns a value for each measurement. Average here to get average over all data.
  objtape <- scorematchingad::avgrange(objtape_ind)
  # Use the taping values as the starting point for optimisation.
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

  # The free mean-link parameters occupy positions 1:length(conprep$x0) in the flat vector;
  # position length(conprep$x0) + 1 is k. Enforce k > 0 as a lower bound.
  if (is.null(lb)){
    lb <- rep(-Inf, objtape$domain)
    lb[length(conprep$x0) + 1] <- 0
  }

  # activate a progress bar
  pb <- progress::progress_bar$new(total = combined_opts$maxeval + 5, format = ":bar :percent :current :tick_rate elapsed::elapsedfull eta::eta")
  
  
  # Optimisation
  # current dynamic parameter values of tapes will be used
  # nloptr minimises; the SvMF log-likelihood is to be maximised, hence the negation.
  nlopt <- nloptr::nloptr(
    x0 = x0,
    eval_f = function(theta){
      if (!pb$finished) pb$tick()
      list(objective = -objtape$forward(0, theta), gradient = -objtape$Jacobian(theta))
      },
    lb = lb,
    ub = ub,
    eval_g_eq =  function(theta){list(
      # Orthonormality constraints apply only to the mean-link parameters (first
      # constraint_tape$domain elements of theta). The zero padding gives zero Jacobian
      # for k, aremaining, and the Cayley G0 components, which are unconstrained.
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

  # The solution vector is partitioned as [free mean-link | k | aremaining | Cayley G0].
  # Reinsert fixed mean-link values (t_sfi2u uses om0vec's fixed positions as the template).
  meanpars <- nlopt$solution[1:length(conprep$x0)]
  meanpars <- scorematchingad:::t_sfi2u(meanpars, om0vec, conprep$isfixed)
  fullparam <- c(meanpars, nlopt$solution[-(1:length(conprep$x0))])

  # Reconstruct the full typed parameter list (om, k, aremaining, G0) from the flat vector,
  # applying the inverse Cayley transform to recover G0 from its unconstrained coordinates.
  estparamlist <- mobius_SvMF_partransport_nota1_fromvecparams_forR(fullparam, p, qs, qe,
                                                  referencecoords = G0reference,
                                                  G01behaviour = G01behaviour,
                                                  G01 = preplist$G0[,1])
  
  #project mean pars to have correct orthogonality
  projectedom <- Omega_proj(mobius_link_Omega_unvec(estparamlist$omvec, p, length(om0$qe1), check = FALSE))
  try({mobius_link_Omega_check(projectedom)})

  # The model is identifiable only up to permutation of (a_j, G0[:,j]) pairs. By convention
  # the scales are returned in decreasing order with the G0 columns reordered to match.
  aord <- order(estparamlist$aremaining, decreasing = TRUE)
  aremaining <- estparamlist$aremaining[aord]
  estparamlist$G0[,-1] <- estparamlist$G0[,-1][, aord]

  pred <- mobius_link(xs = preplist$xs, xe = preplist$xe, param = projectedom)
  if (p!=3){
    # When p != 3 the joint AD optimisation uses a numerical approximation to the vMF
    # normalising constant (required for AD differentiability). After convergence,
    # mobius_SvMF_konly refines k by an exact univariate search using the more accurate
    # Bessel function implementation in base R, correcting for the approximation error in k.
    result <- mobius_SvMF_konly(y = preplist$y, ymean = pred, a = c(a1, aremaining), G0 = estparamlist$G0)
    estparamlist$k <- result$k
    lLik <- result$lLik
  }
  
  ### Making nicer return objects ###
  # Residuals and distances are computed in the standardised coordinate system (before
  # reverting), since they are invariant to the choice of coordinates.
  dists <- acos(rowSums(pred * preplist$y))
  # rresids_std: parallel-transported to G0[:,1] and scaled by sqrt(k)*a[1]/a[j] —
  # under high concentration these are approximately i.i.d. standard normal.
  rresids_std <- resid_SvMF_partransport(preplist$y, pred, estparamlist$k, c(a1, estparamlist$aremaining), estparamlist$G0, scale = TRUE)
  # rresids_G0: same parallel transport to G0[:,1] but unscaled.
  rresids_G0 <- resid_SvMF_partransport(preplist$y, pred, G0 = estparamlist$G0, scale = FALSE)
  # rresids_I: parallel-transported to the north pole; [, -1] drops the response component
  # (always ~1 for small residuals), keeping only the p-1 tangent components.
  rresids_I_tmp <- rotated_resid(preplist$y, pred, north_pole(ncol(preplist$y)))
  rresids_I <- rresids_I_tmp[, -1]
  attr(rresids_I, "samehemisphere") <-  attr(rresids_I_tmp, "samehemisphere")
  colnames(rresids_I) <- paste0("r", 1:ncol(rresids_I))

  # Invert the rotations and centering applied by standardise_data, expressing the estimated
  # mean-link parameters in the original user-supplied coordinates.
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

  # G0 is in standardised y coordinates; t(std_rotation) reverts it to the original space.
  G0 <- estparamlist$G0
  G0 <- t(attr(preplist$y, "std_rotation")) %*% G0
  # Sign convention for G0 axes:
  # (1) make first element of each non-first column positive, removing the per-axis sign ambiguity
  G0[,-1] <- standardise_col_signs(G0[,-1])
  # (2) ensure G0 is a proper rotation matrix (det = +1) by flipping the final column if needed
  if (det(G0) < 0){G0[,p] <- -G0[,p]}
  # DoF: mean link + concentration + scales + G0 axes
  # (p-1)-1: aremaining has p-1 components but prod(aremaining) = 1 removes one degree of freedom
  # G0: DoF_Stiefel(p,p) when G0[,1] is free; DoF_Stiefel(p-1, p-1) when G0[,1] is
  #     identified with p1 or fixed (one fewer free column)
  DoF <- mobius_dof(p, length(est$qs1), length(est$qe1), fix_qs1 = fix_qs1, fix_qe1 = fix_qe1) +
    1 + #concentration
    (p-1)-1 + #aremaining given that prod(aremaining) = 1
    if (G01behaviour == "free"){ #G0 freedom
      DoF_Stiefel(p,p)
    } else {
      DoF_Stiefel(p-1, p-1)
    }
  # AIC
  AIC <- 2*DoF - 2 * lLik
  
  if (p==3) {
    if (estparamlist$k < 1E-15){
    	warning("Estimated concentration is very small and computation of the vMF normalising constant may be breaking down.")
    }
  }
  
  #Scealy and Wood (2019) Proposition 1 check for unimodality
  SvMF_cann_check(SvMF_cann(k = estparamlist$k, a = c(a1, aremaining), G = G0))
  
  niceout <- list(
    mean = est,
    k = estparamlist$k,
    a = c(a1, aremaining),
    G0 = G0,
    obj = nlopt$objective,
    nlopt = nlopt,
    y = y,
    xs = xs,
    # destandardise_Euc inverts the centering and whitening from standardise_data; using
    # preplist$xe also recovers any columns added by addEuccovars (dummy-zero and intercept).
    xe = if (!is.null(xe)){if (intercept){destandardise_Euc(preplist$xe, attr(preplist$xe, "std_center"), attr(preplist$xe, "std_rotation"))} else {xe}},
    pred = destandardise_sph(pred, rotation = attr(preplist$y, "std_rotation")),
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

# Preliminary estimation step for mobius_SvMF(): runs a vMF regression to get a good
# starting mean link, then estimates G0 axes (parallel transport base) and shape
# parameters a using method-of-moments on the rotated residuals (Scealy & Wood 2019, Sec 4.1).
mobius_SvMF_partransport_prelim <- function(y, xs, xe, mean = NULL, G0 = NULL, G01behaviour = "p1", type = "LinEuc", fix_qs1 = FALSE, fix_qe1 = (type == "LinEuc"), intercept = TRUE, ...){
  prelim <- mobius_vMF(y = y, xs = xs, xe = xe, 
             start = mean, 
             type = type, fix_qs1 = fix_qs1, fix_qe1 = fix_qe1, intercept = intercept, ...)
  # update starting values accordingly
  mean <- as_mobius_link_Omega(prelim$est)
  k <- prelim$k
  p <- ncol(y)
  
  # get/choose G01 depending on behaviour
  if (G01behaviour == "fixed" && is.null(G0)){stop("At least the first column of G0 must be supplied when G01behaviour = 'fixed'.")}
  G01 <- switch(G01behaviour,
         p1 = prelim$est$p1,
         free = if(is.null(G0)){prelim$est$p1}else{G0[,1]},
         fixed = G0[,1])
  # get rotated residuals
  rresid <- rotated_resid(y, prelim$pred, base = G01)
  if (!is.null(G0) && all(!is.na(G0))){
    # G0 is supplied: re-express its axes relative to the new G01 (prelim$est$p1) by parallel-
    # transporting the columns of G0[,-1] from G0[,1] to G01. The minus sign is because
    # -jupp_Rmat(a,b) equals parallel_transport_mat(a,b). Then estimate scales from the rotated residuals.
    if (G01behaviour == "p1"){
      G0 <- cbind(G01, -jupp_Rmat(G0[,1], G01) %*% G0[,-1])
    }
    aremaining <- SvMF_prelim_scales(rresid, G0)
  } else {
    # G0 not supplied: estimate principal-axis directions from the sample second-moment
    # matrix of the rotated residuals (Scealy & Wood 2019, §4.1).
    G0 <- SvMF_moment_axes(rresid, G01)
    # Estimate scale parameters a using the high-concentration approximation.
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

# Estimate only the concentration parameter k by univariate MLE, holding a and G0 fixed.
# Used after the main joint optimisation converges when p != 3 (because the normalising
# constant can be computed more accurately for a scalar search than during joint AD optimisation).
mobius_SvMF_konly <- function(y, ymean, a, G0){
  yrot <- undo_partransport(y = y, ymean = ymean, G01 = G0[,1])
  res <- optimise(function(k){
    sum(SvMF_log_lik_cann(yrot, SvMF_cann(k = k, a = a, G = G0)))
  }, lower = 1E-8, upper = 1E5, maximum = TRUE)
  SvMF_cann_check(SvMF_cann(k = res$maximum, a = a, G = G0))
  if (res$maximum == 1E-8){warning("Concentration at numerical lower limit of 1E-8")}
  if (res$maximum == 1E5){warning("Concentration at numerical upper limit of 1E5")}
  return(list(
    k = res$maximum,
    a = a,
    G0 = G0,
    lLik = res$objective
  ))
}

#' @title Function for simulating data given mean link and SvMF parameters
#' @inheritParams mobius_SvMF_log_lik
#' @return A matrix of `p+1` columns and the same number of rows as `xs` or `xe`. The final column is the log-density of the simulated response.
#' @family regression
#' @export
rmobius_SvMF <- function(xs, xe, mean, k, a, G0){
  ymean <- mobius_link(xs = xs, xe = xe, param = mean)
  
  # simulate noise
  y_ld <- t(apply(ymean, 1, function(mn){
    # Construct local axes at the predicted mean mn by parallel-transporting the columns of
    # G0[,-1] from G0[,1] to mn. (-jupp_Rmat(a,b) equals parallel_transport_mat(a,b).)
    G <- cbind(mn, -jupp_Rmat(G0[,1], mn) %*% G0[,-1])
    obs <- rSvMF(1, SvMF_cann(k, a, G))
    ld <- ldSvMF_cann(obs, k = k, a = a, G = G)
    return(c(obs, ld))
  }))
  return(y_ld)
}


# @param y matrix of observations
# @param ymean matrix of predicted means
# @param G01 first column of the G0 matrix
# Parallel-transports each observation y[i,] back from ymean[i,] to G01, bringing all
# observations into a common reference frame (as if they all had mean direction G01).
# This is required to evaluate the pooled SvMF log-likelihood via SvMF_log_lik_cann /
# ldSvMF_cann, which assume the mean direction is G01.
undo_partransport <- function(y, ymean, G01){
  #rotate all observations to reverse the transport from G0[,1] to ymean
  yrot <- lapply(1:nrow(y), function(idx){
    drop(t(parallel_transport_mat(G01, ymean[idx, ])) %*% y[idx, ])
  })
  yrot <- do.call(rbind, yrot)
  return(yrot)
}

#' @title log-density of data for a given SvMF regression
#' @param mean Parameter object (see [`mobius_link_params`]) specifying mean link.
#' @param k Concentration of the SvMF error distribution
#' @param a Scales of the SvMF error distribution
#' @param G0 The base location of parallel transport along with axes \eqn{\gamma_{0j}}, \eqn{j=2,...,p}.
#' @inheritParams mobius_SvMF
#' @description The log-density of each row of `y` for a given SvMF regression. Two methods are used. The approximate method used in the optimisation of all parameters (labelled `Cpp`) and an exact log-density using highly accurate Bessel function implementations from base `R` (labelled `R`).
#' @return A matrix with two columns and the same number of rows as `y`.
#' @family regression
#' @export
mobius_SvMF_log_lik <- function(y, xs, xe, mean, k, a, G0){
  ymean <- mobius_link(xs = xs, xe = xe, param = mean)
  # diff: the discrepancy between the approximate normalising constant used in the C++ AD tape
  # and the exact value; stored as an attribute for diagnostic use.
  diff <- vMF_log_norm_const(k, ncol(y)) - vMF_log_norm_const_exact(k, ncol(y))

  yrot <- undo_partransport(y, ymean, G01 = G0[,1])
  # ldCpp: uses the same approximation as the optimiser (fast, suitable for AD).
  ldCpp <- ldSvMF_cann(yrot, k = k, a = a, G = G0)
  # ldR: uses exact Bessel functions from base R (accurate, not AD-differentiable).
  ldR <- SvMF_log_lik_cann(yrot, SvMF_cann(k = k, a = a, G = G0))
  ld <- cbind(Cpp = ldCpp, R = ldR)
  attr(ld, "error") <- diff
  return(ld)
}
