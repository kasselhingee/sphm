# Proposition 4(i): the scaled Moebius link with an isotropic spherical scale,
#     p = qs,  Bs = beta_s I_(p-1),  0 < beta_s <= 1.
#
# Spherical covariate only (the form Proposition 4(i) is stated for, qe = 0):
#     mu_s(xs) = M_S(xs : Rtilde0, psi),   psi = phi r_s1,  phi = (1 - beta_s)/(1 + beta_s)
# with Rtilde0 orthogonal and ||psi|| < 1. beta_s = 1 gives Proposition 4(ii). The
# optimiser works in an unconstrained coordinate u with psi = u / sqrt(1 + ||u||^2), so
# ||psi|| < 1 holds automatically.
#
# Proposition 4(i) itself is stated for qe = 0. The Euclidean variant here keeps the
# restriction on the SPHERICAL component and adds the same linear Euclidean term the
# LinEuc implementation uses:
#     mu(xs, xe) = P iSp{ beta_s Sp(Qs^T xs) + Be Qe*^T xe_tilde },   Rtilde0 = P Qs^T
# so beta_s = 1 reduces exactly to the Proposition 4(ii)+Euclidean model.
#
# Reached through mobius_vMF(type = "Prop4i") / mobius_SvMF(type = "Prop4i"); the
# dispatchers in Mobius_vMF.R and Mobius_SvMF.R pick the _euc variant when xe is supplied.
# Shared internals live in Mobius_prop4_helpers.R.
#
# Future unification candidate: as in Mobius_prop4ii.R, the spherical-only and _euc
# joint-fit engines below are structurally parallel but parameterise different manifolds.

.prop4i_u_to_psi <- function(u) {
  u <- as.numeric(u)
  den <- sqrt(1 + sum(u^2))
  u / den
}

.prop4i_psi_to_u <- function(psi, tol = 1e-12) {
  psi <- as.numeric(psi)
  r2 <- sum(psi^2)
  if (!is.finite(r2) || r2 >= 1) {
    r2 <- min(r2, 1 - tol)
    psi <- psi * sqrt(r2 / max(sum(psi^2), tol))
  }
  psi / sqrt(max(1 - sum(psi^2), tol))
}

.prop4i_beta_from_psi <- function(psi) {
  phi <- sqrt(sum(as.numeric(psi)^2))
  (1 - phi) / (1 + phi)
}

.prop4i_reference_direction <- function(mean, tol = 1e-10) {
  mean <- as_mobius_link_prop4i(mean, check = FALSE)
  phi <- sqrt(sum(mean$psi^2))
  if (phi > tol) {
    rs1 <- mean$psi / phi
    drop(mean$Rtilde0 %*% rs1)
  } else {
    # At beta_s = 1, r_s1 is not identifiable.  Use the same representative
    # convention as the Proposition 4(ii) implementation.
    mean$Rtilde0[, 1]
  }
}

.prop4i_map <- function(xs, Rtilde0, psi) {
  xs <- as.matrix(xs)
  psi <- as.numeric(psi)
  z <- sweep(xs, 2L, psi, "+")
  den <- rowSums(z^2)
  if (any(!is.finite(den)) || any(den <= .Machine$double.eps)) {
    stop("Numerical singularity in the Proposition 4(i) M\u00f6bius map.")
  }
  fac <- (1 - sum(psi^2)) / den
  inner <- sweep(z, 1L, fac, "*")
  inner <- sweep(inner, 2L, psi, "+")
  pred <- inner %*% t(Rtilde0)
  nr <- sqrt(rowSums(pred^2))
  sweep(pred, 1L, nr, "/")
}

# The parameter class carries the Euclidean fields optionally, matching
# mobius_link_prop4ii. beta_s = 1 (equivalently psi = 0) recovers Proposition 4(ii).
#' @name mobius_link_prop4i
#' @title Proposition 4(i) Mean Link Parameters
#' @description Parameters of the scaled Mobius mean link constrained by \eqn{p = q_s} and
#' an isotropic spherical scale \eqn{B_s = \beta_s I_{p-1}}, \eqn{0 < \beta_s \le 1} ---
#' Proposition 4(i) of the manuscript. Without Euclidean covariates the link is the
#' isotropic Mobius map \eqn{\mu_s(x_s) = M_S(x_s : \tilde R_0, \psi)} with
#' \eqn{\psi = \phi r_{s1}} and \eqn{\phi = (1-\beta_s)/(1+\beta_s)}.
#' \eqn{\beta_s = 1} gives [`mobius_link_prop4ii`].
#' @param Rtilde0 A `p` x `p` orthonormal matrix.
#' @param psi The Mobius shift, a vector of length `p` with `||psi|| < 1`. Supply either
#'   `psi`, or `beta_s` and `rs1`.
#' @param beta_s The isotropic spherical scale, in `(0, 1]`.
#' @param rs1 The reference direction, a unit vector of length `p`. Not identified when
#'   `beta_s == 1`.
#' @param P Final rotation on the response sphere: a `p` x `p` rotation matrix. Supplied
#'   only when the fit has Euclidean covariates, in which case \eqn{\tilde R_0 = P Q_s^\top}
#'   and `rs1`/`psi` are derived from it rather than taken from the arguments.
#' @param Qe_star An `m` x `(p-1)` matrix with orthonormal columns acting on the
#'   standardised Euclidean covariates, where `m` counts them plus any intercept.
#' @param Be The `p-1` diagonal entries of the Euclidean scaling matrix, each in `(0, 1)`.
#' @param xe_center,xe_scale The centre and scale used to standardise the Euclidean
#'   covariates during fitting, so that predictions reproduce them.
#' @param intercept `TRUE` if an intercept column was appended to the Euclidean covariates.
#' @param check Whether to validate the result.
#' @param obj An object to coerce or validate.
#' @return `mobius_link_prop4i()` and `as_mobius_link_prop4i()` return an object of class
#'   "mobius_link_prop4i". `dim()` returns the named vector `c(p, qs, qe)`.
#' @details
#' Proposition 4(i) itself is stated for \eqn{q_e = 0}. Supplying `P`, `Qe_star` and `Be`
#' (together, or not at all) keeps the restriction on the spherical component and adds the
#' same linear Euclidean term the `LinEuc` link uses, so \eqn{\beta_s = 1} reduces exactly
#' to the Proposition 4(ii)+Euclidean model.
#'
#' Use [`as_mobius_link_cann()`] to convert to the canonical parameterisation for
#' interpretation; see [`mobius_link_params`].
#' @seealso [`mobius_link_prop4ii`], [`mobius_link()`], [`mobius_SvMF()`].
#' @family link-function
#' @export
mobius_link_prop4i <- function(Rtilde0, psi = NULL,
                               beta_s = NULL, rs1 = NULL,
                               P = NULL, Qe_star = NULL, Be = NULL,
                               xe_center = NULL, xe_scale = NULL,
                               intercept = TRUE, check = TRUE) {
  Rtilde0 <- as.matrix(Rtilde0)
  storage.mode(Rtilde0) <- "double"
  p <- nrow(Rtilde0)
  has_euc <- !is.null(P) || !is.null(Qe_star) || !is.null(Be)

  if (has_euc) {
    if (is.null(P) || is.null(Qe_star) || is.null(Be)) {
      stop("P, Qe_star, and Be must either all be supplied or all be NULL.")
    }
    if (is.null(beta_s)) {
      stop("beta_s must be supplied for the Proposition 4(i)+Euclidean form.")
    }
    P <- as.matrix(P)
    Qe_star <- as.matrix(Qe_star)
    storage.mode(P) <- storage.mode(Qe_star) <- "double"
    beta_s <- as.numeric(beta_s)[1L]
    if (!is.finite(beta_s) || beta_s <= 0 || beta_s > 1) {
      stop("beta_s must satisfy 0 < beta_s <= 1.")
    }
    # Rtilde0 = P Qs^T, hence Qs = Rtilde0^T P and r_s1 = Qs[,1]: with Euclidean
    # covariates the reference direction is determined by P, so psi/rs1 arguments
    # are derived here rather than taken from the caller.
    Qs <- t(Rtilde0) %*% P
    rs1 <- Qs[, 1L]
    phi <- (1 - beta_s) / (1 + beta_s)
    psi <- phi * rs1
  } else {
    if (is.null(psi)) {
      if (is.null(beta_s)) beta_s <- 1
      beta_s <- as.numeric(beta_s)[1]
      if (!is.finite(beta_s) || beta_s <= 0 || beta_s > 1) {
        stop("beta_s must satisfy 0 < beta_s <= 1.")
      }
      phi <- (1 - beta_s) / (1 + beta_s)
      if (is.null(rs1)) {
        rs1 <- c(1, rep(0, p - 1L))
      }
      rs1 <- as.numeric(rs1)
      if (length(rs1) != p) stop("rs1 must have length p.")
      nr <- sqrt(sum(rs1^2))
      if (!is.finite(nr) || nr <= 0) stop("rs1 must be nonzero and finite.")
      rs1 <- rs1 / nr
      psi <- phi * rs1
    } else {
      psi <- as.numeric(psi)
    }
    phi <- sqrt(sum(psi^2))
    beta_s <- (1 - phi) / (1 + phi)
    # at beta_s = 1 the reference direction is not identifiable
    rs1 <- if (phi > 1e-10) psi / phi else rep(NA_real_, p)
  }

  obj <- list(
    Rtilde0 = Rtilde0,
    psi = psi,
    phi = phi,
    beta_s = beta_s,
    rs1 = rs1,
    P = if (has_euc) P else NULL,
    Qe_star = if (has_euc) Qe_star else NULL,
    Be = if (has_euc) as.numeric(Be) else NULL,
    xe_center = if (has_euc) as.numeric(xe_center) else NULL,
    xe_scale = if (has_euc) as.numeric(xe_scale) else NULL,
    intercept = if (has_euc) isTRUE(intercept) else FALSE
  )
  class(obj) <- c("mobius_link_prop4i", "list")
  if (check) mobius_link_prop4i_check(obj)
  obj
}

mobius_link_prop4i_check <- function(obj,
                                     tol = 100 * sqrt(.Machine$double.eps)) {
  if (!inherits(obj, "mobius_link_prop4i")) {
    stop("obj must inherit from 'mobius_link_prop4i'.")
  }
  R <- obj$Rtilde0
  if (!is.matrix(R) || nrow(R) != ncol(R) || nrow(R) < 2L) {
    stop("Rtilde0 must be a square p by p matrix with p >= 2.")
  }
  if (any(!is.finite(R))) stop("Rtilde0 contains non-finite values.")
  err <- max(abs(crossprod(R) - diag(ncol(R))))
  if (err > tol) {
    stop(sprintf("Rtilde0 is not orthogonal: max|t(R)R-I| = %.3e.", err))
  }
  psi <- obj$psi
  if (!is.numeric(psi) || length(psi) != nrow(R) || any(!is.finite(psi))) {
    stop("psi must be a finite numeric vector of length p.")
  }
  if (sum(psi^2) >= 1 - 1e-12) {
    stop("Proposition 4(i) requires ||psi|| < 1 (equivalently beta_s > 0).")
  }

  has_euc <- !is.null(obj$P)
  if (!has_euc) {
    if (!is.null(obj$Qe_star) || !is.null(obj$Be)) {
      stop("P, Qe_star, and Be must either all be supplied or all be NULL.")
    }
    return(invisible(NULL))
  }

  p <- nrow(R)
  k <- p - 1L
  P <- obj$P
  if (!all(dim(P) == c(p, p)) || any(!is.finite(P)) ||
      max(abs(crossprod(P) - diag(p))) > tol || det(P) < 0) {
    stop("P must be a proper orthogonal p by p matrix.")
  }
  if (!is.finite(obj$beta_s) || obj$beta_s <= 0 || obj$beta_s > 1) {
    stop("beta_s must satisfy 0 < beta_s <= 1.")
  }
  V <- obj$Qe_star
  if (!is.matrix(V) || ncol(V) != k || nrow(V) < k || any(!is.finite(V))) {
    stop("Qe_star must be m by (p-1), with m >= p-1.")
  }
  if (max(abs(crossprod(V) - diag(k))) > tol) {
    stop("Qe_star must have orthonormal columns.")
  }
  b <- obj$Be
  if (length(b) != k || any(!is.finite(b)) || any(b <= 0) || any(b >= 1)) {
    stop("Be must contain p-1 values strictly between zero and one.")
  }
  q <- nrow(V) - as.integer(isTRUE(obj$intercept))
  if (q < 0L || length(obj$xe_center) != q || length(obj$xe_scale) != q) {
    stop("Euclidean centering/scaling information has the wrong dimension.")
  }
  if (any(!is.finite(obj$xe_center)) || any(!is.finite(obj$xe_scale)) ||
      any(obj$xe_scale <= 0)) {
    stop("Euclidean centering/scaling information is invalid.")
  }
  invisible(NULL)
}

#' @describeIn mobius_link_prop4i Coerce a compatible object to the Proposition 4(i)
#'   parameterisation.
#' @export
as_mobius_link_prop4i <- function(obj, check = TRUE) {
  if (inherits(obj, "mobius_link_prop4i")) {
    if (check) mobius_link_prop4i_check(obj)
    return(obj)
  }
  if (is.list(obj) &&
      all(c("Rtilde0", "P", "beta_s", "Qe_star", "Be",
            "xe_center", "xe_scale") %in% names(obj))) {
    return(mobius_link_prop4i(
      Rtilde0 = obj$Rtilde0, P = obj$P, beta_s = obj$beta_s,
      Qe_star = obj$Qe_star, Be = obj$Be,
      xe_center = obj$xe_center, xe_scale = obj$xe_scale,
      intercept = if (is.null(obj$intercept)) TRUE else obj$intercept,
      check = check
    ))
  }
  if (is.list(obj) && !is.null(obj$Rtilde0)) {
    if (!is.null(obj$psi)) {
      return(mobius_link_prop4i(obj$Rtilde0, psi = obj$psi, check = check))
    }
    if (!is.null(obj$beta_s)) {
      return(mobius_link_prop4i(obj$Rtilde0, beta_s = obj$beta_s,
                                rs1 = obj$rs1, check = check))
    }
  }
  stop("obj must be a mobius_link_prop4i object or a compatible list.")
}

#' @export
print.mobius_link_prop4i <- function(x, ...) {
  if (is.null(x$P)) {
    cat("Proposition 4(i) spherical M\u00f6bius link\n")
    cat("  p = qs =", nrow(x$Rtilde0), ", qe = 0\n")
  } else {
    cat("Proposition 4(i)-constrained spherical link with LinEuc covariates\n")
    cat("  p = qs =", nrow(x$Rtilde0), "\n")
    cat("  qe =", length(x$xe_center), "\n")
  }
  cat("  beta_s =", format(x$beta_s, digits = 6), "\n")
  cat("  phi =", format(x$phi, digits = 6), "\n")
  cat("  det(Rtilde0) =", format(det(x$Rtilde0), digits = 6), "\n")
  if (!is.null(x$Be)) {
    cat("  Be =", paste(format(x$Be, digits = 5), collapse = ", "), "\n")
  }
  invisible(x)
}

#' @export
dim.mobius_link_prop4i <- function(x) {
  p <- nrow(x$Rtilde0)
  c(p = p, qs = p, qe = length(x$xe_center))
}

#' @noRd
#' @description Evaluate the Proposition 4(i) mean link. Called by [`mobius_link()`].
mobius_link_pred_prop4i <- function(xs = NULL, xe = NULL, param = NULL, check = TRUE) {
  if (check) mobius_link_prop4i_check(param)
  if (is.null(xs)) stop("The Proposition 4(i) link requires xs.")
  xs <- as.matrix(xs)
  if (ncol(xs) != nrow(param$Rtilde0)) {
    stop("ncol(xs) must equal nrow(Rtilde0).")
  }

  if (is.null(param$P)) {
    # strict Proposition 4(i): qe = 0
    if (.is_xe_nonempty(xe)) {
      stop("This fitted Proposition 4(i) link has no Euclidean component.")
    }
    return(.prop4i_map(xs, param$Rtilde0, param$psi))
  }

  if (is.null(xe)) stop("xe must be supplied for this fitted model.")
  ep <- .prop4_prepare_euc(
    xe, intercept = param$intercept,
    center = param$xe_center, scale = param$xe_scale,
    fitting = FALSE
  )
  if (nrow(xs) != nrow(ep$model)) {
    stop("xs and xe must have the same number of rows.")
  }
  # Rtilde0 = P Qs^T, hence Qs = Rtilde0^T P.
  Qs <- t(param$Rtilde0) %*% param$P
  spherical_part <- param$beta_s * Sp(xs %*% Qs)
  euclidean_part <- sweep(ep$model %*% param$Qe_star, 2L, param$Be, "*")
  iSp(spherical_part + euclidean_part) %*% t(param$P)
}

.prop4i_mean_spec <- function(Rref, psi_start = NULL) {
  p <- nrow(Rref)
  dR <- .prop4_skew_dim(p)
  idxR <- seq_len(dR)
  idxPsi <- dR + seq_len(p)
  par0 <- numeric(dR + p)
  if (!is.null(psi_start)) par0[idxPsi] <- .prop4i_psi_to_u(psi_start)
  list(
    p = p, Rref = Rref,
    idxR = idxR, idxPsi = idxPsi,
    par0 = par0,
    lower = rep(-Inf, length(par0)),
    upper = rep( Inf, length(par0))
  )
}

.prop4i_unpack_mean <- function(theta, spec) {
  R <- .prop4_cayley(theta[spec$idxR], spec$p) %*% spec$Rref
  psi <- .prop4i_u_to_psi(theta[spec$idxPsi])
  mean <- mobius_link_prop4i(R, psi = psi, check = FALSE)
  list(Rtilde0 = R, psi = psi, mean = mean)
}

.prop4i_fit_vmf_component <- function(y, xs, spec,
                                      algorithm = "NLOPT_LD_SLSQP",
                                      opts = list()) {
  eval_f <- function(theta) {
    pars <- .prop4i_unpack_mean(theta, spec)
    pred <- try(.prop4i_map(xs, pars$Rtilde0, pars$psi), silent = TRUE)
    if (inherits(pred, "try-error") || any(!is.finite(pred))) return(1e100)
    -mean(rowSums(y * pred))
  }

  default_opts <- list(
    algorithm = algorithm,
    xtol_rel = 1e-9,
    maxeval = 10000,
    print_level = 0
  )
  combined_opts <- utils::modifyList(default_opts, opts)
  uses_gradient <- grepl("^NLOPT_LD_", combined_opts$algorithm)
  eval_nloptr <- if (uses_gradient) {
    function(theta) {
      list(
        objective = eval_f(theta),
        gradient = .prop4_fd_gradient(
          eval_f, theta, spec$lower, spec$upper
        )
      )
    }
  } else {
    eval_f
  }

  fit <- try(
    nloptr::nloptr(
      x0 = spec$par0,
      eval_f = eval_nloptr,
      lb = spec$lower,
      ub = spec$upper,
      opts = combined_opts
    ),
    silent = TRUE
  )

  # The starting point (Proposition 4(ii), psi=0 by default) is always valid.
  # Retain it if an optimiser fails unexpectedly.
  theta_best <- spec$par0
  f_best <- eval_f(theta_best)
  if (!inherits(fit, "try-error") && is.finite(fit$objective) &&
      fit$objective < f_best) {
    theta_best <- fit$solution
    f_best <- fit$objective
  }
  pars <- .prop4i_unpack_mean(theta_best, spec)
  pred <- .prop4i_map(xs, pars$Rtilde0, pars$psi)

  if (inherits(fit, "try-error")) {
    fit <- list(status = -999L, message = as.character(fit),
                objective = f_best, solution = theta_best)
  } else {
    fit$objective <- f_best
    fit$solution <- theta_best
  }

  list(nlopt = fit, pars = pars, pred = pred, objective = f_best)
}

mobius_vMF_prop4i <- function(
    y, xs, xe = NULL, start = NULL,
    det_constraint = c("orthogonal", "rotation"),
    algorithm = "NLOPT_LD_SLSQP", opts = list(), check = TRUE) {

  det_constraint <- match.arg(det_constraint)
  y <- as.matrix(y)
  xs <- as.matrix(xs)
  storage.mode(y) <- storage.mode(xs) <- "double"

  if (!is.null(xe) && ncol(as.matrix(xe)) > 0L) {
    stop("Strict Proposition 4(i) has qe = 0; xe must be NULL.")
  }
  if (nrow(y) != nrow(xs) || ncol(y) != ncol(xs)) {
    stop("For Proposition 4(i), y and xs must have identical n by p dimensions.")
  }
  if (check) {
    .prop4_check_unit_rows(y, "y")
    .prop4_check_unit_rows(xs, "xs")
  }

  start_obj <- if (is.null(start)) NULL else as_mobius_link_prop4i(start, check = FALSE)
  start_R <- if (is.null(start_obj)) NULL else start_obj$Rtilde0
  start_psi <- if (is.null(start_obj)) rep(0, ncol(y)) else start_obj$psi

  target_signs <- if (det_constraint == "rotation") 1 else c(1, -1)
  fits <- vector("list", length(target_signs))
  for (j in seq_along(target_signs)) {
    sgn <- target_signs[j]
    Rref <- if (!is.null(start_R) &&
                (if (det(start_R) >= 0) 1 else -1) == sgn) {
      start_R
    } else {
      .prop4_procrustes(y, xs, det_sign = sgn)$R
    }
    spec <- .prop4i_mean_spec(Rref, psi_start = start_psi)
    fits[[j]] <- .prop4i_fit_vmf_component(
      y = y, xs = xs, spec = spec,
      algorithm = algorithm, opts = opts
    )
    fits[[j]]$det_sign <- sgn
  }

  values <- vapply(fits, function(z) z$objective, numeric(1))
  best <- fits[[which.min(values)]]
  pars <- best$pars
  pred <- best$pred
  dots <- pmax(-1, pmin(1, rowSums(y * pred)))
  kfit <- .prop4_estimate_k(mean(dots), ncol(y))
  n <- nrow(y)
  p <- ncol(y)
  dmean <- .prop4_skew_dim(p) + p
  DoF <- dmean + 1L
  lLik <- n * kfit$loglik_per_obs

  rr <- try(rotated_resid(y, pred, north_pole(p)), silent = TRUE)
  rresids <- if (inherits(rr, "try-error")) NULL else rr[, -1, drop = FALSE]
  if (!is.null(rresids)) {
    attr(rresids, "samehemisphere") <- attr(rr, "samehemisphere")
    colnames(rresids) <- paste0("r", seq_len(ncol(rresids)))
  }

  mean <- pars$mean
  out <- list(
    est = mean, mean = mean,
    Rtilde0 = mean$Rtilde0,
    psi = mean$psi, phi = mean$phi,
    beta_s = mean$beta_s, rs1 = mean$rs1,
    k = kfit$k,
    obj = best$objective,
    solution = best$nlopt$solution,
    nlopt = best$nlopt,
    y = y, xs = xs, xe = NULL,
    pred = pred, rresids = rresids,
    dists = acos(dots),
    DoF = DoF, AIC = 2 * DoF - 2 * lLik, lLik = lLik,
    start = start,
    component_fits = lapply(fits, function(z) {
      list(det_sign = z$det_sign, status = z$nlopt$status,
           message = z$nlopt$message, objective = z$objective)
    }),
    linktype = list(type = "Prop4i", det_constraint = det_constraint,
                    intercept = FALSE, euclidean = FALSE)
  )
  class(out) <- c("mobius_vMF_prop4i", "mobius_vMF", "list")
  out
}

.prop4i_euc_mean_spec <- function(
    Rref, Pref, Vref, bstart, beta_start,
    xe_center, xe_scale, intercept) {

  base <- .prop4ii_mean_spec_euc(
    Rref = Rref,
    Pref = Pref,
    Vref = Vref,
    bstart = bstart,
    xe_center = xe_center,
    xe_scale = xe_scale,
    intercept = intercept
  )

  beta_start <- min(max(as.numeric(beta_start)[1L], 1e-5), 1 - 1e-5)
  idxBeta <- length(base$par0) + 1L

  list(
    base = base,
    base_length = length(base$par0),
    idxBeta = idxBeta,
    par0 = c(base$par0, qlogis(beta_start)),
    lower = c(base$lower, -12),
    upper = c(base$upper,  12)
  )
}

.prop4i_euc_unpack_mean <- function(theta, spec) {
  mp <- .prop4ii_unpack_mean_euc(
    theta[seq_len(spec$base_length)], spec$base
  )
  beta_s <- plogis(theta[spec$idxBeta])

  mean <- mobius_link_prop4i(
    Rtilde0 = mp$Rtilde0,
    P = mp$P,
    beta_s = beta_s,
    Qe_star = mp$Qe_star,
    Be = mp$Be,
    xe_center = spec$base$xe_center,
    xe_scale = spec$base$xe_scale,
    intercept = spec$base$intercept,
    check = FALSE
  )

  list(
    mean = mean,
    Rtilde0 = mp$Rtilde0,
    P = mp$P,
    beta_s = beta_s,
    Qe_star = mp$Qe_star,
    Be = mp$Be
  )
}

.prop4i_euc_fit_vmf_component <- function(
    y, xs, xe, spec,
    algorithm = "NLOPT_LD_SLSQP", opts = list()) {

  eval_f <- function(theta) {
    pars <- .prop4i_euc_unpack_mean(theta, spec)
    pred <- try(
      mobius_link(xs = xs, xe = xe, param = pars$mean, check = FALSE),
      silent = TRUE
    )
    if (inherits(pred, "try-error") || any(!is.finite(pred))) {
      return(1e100)
    }
    -mean(rowSums(y * pred))
  }

  default_opts <- list(
    algorithm = algorithm,
    xtol_rel = 1e-9,
    maxeval = 10000,
    print_level = 0
  )
  combined_opts <- utils::modifyList(default_opts, opts)
  uses_gradient <- grepl("^NLOPT_LD_", combined_opts$algorithm)

  eval_nloptr <- if (uses_gradient) {
    function(theta) {
      list(
        objective = eval_f(theta),
        gradient = .prop4_fd_gradient(
          eval_f, theta, spec$lower, spec$upper
        )
      )
    }
  } else {
    eval_f
  }

  fit <- try(
    nloptr::nloptr(
      x0 = spec$par0,
      eval_f = eval_nloptr,
      lb = spec$lower,
      ub = spec$upper,
      opts = combined_opts
    ),
    silent = TRUE
  )

  theta_best <- spec$par0
  f_best <- eval_f(theta_best)
  if (!inherits(fit, "try-error") &&
      is.finite(fit$objective) && fit$objective < f_best) {
    theta_best <- fit$solution
    f_best <- fit$objective
  }

  pars <- .prop4i_euc_unpack_mean(theta_best, spec)
  pred <- mobius_link(
    xs = xs, xe = xe, param = pars$mean, check = FALSE
  )

  if (inherits(fit, "try-error")) {
    fit <- list(
      status = -999L,
      message = as.character(fit),
      objective = f_best,
      solution = theta_best
    )
  } else {
    fit$objective <- f_best
    fit$solution <- theta_best
  }

  list(
    nlopt = fit,
    pars = pars,
    pred = pred,
    objective = f_best
  )
}

.prop4i_euc_P_from_strict <- function(mean) {
  mean <- as_mobius_link_prop4i(mean, check = FALSE)
  R <- mean$Rtilde0
  p <- nrow(R)

  rs1 <- mean$rs1
  if (length(rs1) != p || any(!is.finite(rs1))) {
    rs1 <- c(1, rep(0, p - 1L))
  }
  rs1 <- rs1 / sqrt(sum(rs1^2))

  # Obtain an orthonormal tangent basis perpendicular to rs1.
  vv <- svd(matrix(rs1, nrow = 1L), nu = 0, nv = p)$v
  tangent <- vv[, -1L, drop = FALSE]
  Qs <- cbind(rs1, tangent)

  # det(P)=det(Rtilde0)*det(Qs) must be +1.
  target_qs_sign <- if (det(R) >= 0) 1 else -1
  if ((if (det(Qs) >= 0) 1 else -1) != target_qs_sign) {
    Qs[, p] <- -Qs[, p]
  }

  P <- R %*% Qs
  if (det(P) < 0) {
    Qs[, p] <- -Qs[, p]
    P <- R %*% Qs
  }
  .prop4_nearest_orthogonal(P, det_sign = 1)
}

mobius_vMF_prop4i_euc <- function(
    y, xs, xe, start = NULL,
    det_constraint = c("orthogonal", "rotation"),
    intercept = TRUE,
    algorithm = "NLOPT_LD_SLSQP", opts = list(), check = TRUE) {

  det_constraint <- match.arg(det_constraint)
  y <- as.matrix(y)
  xs <- as.matrix(xs)
  xe <- as.matrix(xe)
  storage.mode(y) <- storage.mode(xs) <- storage.mode(xe) <- "double"

  if (nrow(y) != nrow(xs) || nrow(y) != nrow(xe) ||
      ncol(y) != ncol(xs)) {
    stop("For Proposition 4(i)+Euclidean, y and xs must be n by p and xe must have n rows.")
  }

  p <- ncol(y)
  kmean <- p - 1L
  if (check) {
    .prop4_check_unit_rows(y, "y")
    .prop4_check_unit_rows(xs, "xs")
  }

  ep <- .prop4_prepare_euc(xe, intercept = intercept, fitting = TRUE)
  if (ep$m < kmean) {
    stop("With Euclidean covariates, ncol(xe) + intercept must be at least p-1.")
  }

  if (!is.null(start) && is.list(start) && !is.null(start$mean)) {
    start <- start$mean
  }

  start_euc <- NULL
  start_strict <- NULL
  if (!is.null(start)) {
    # Since the two Proposition 4(i) parameter forms now share one class, the presence of
    # the Euclidean block -- not the class -- says which form was supplied. A strict start
    # contributes only Rtilde0; P is then derived by .prop4i_euc_P_from_strict().
    start <- as_mobius_link_prop4i(start, check = FALSE)
    if (is.null(start$P)) {
      start_strict <- start
    } else {
      start_euc <- start
    }
  }

  start_R <- if (!is.null(start_euc)) {
    start_euc$Rtilde0
  } else if (!is.null(start_strict)) {
    start_strict$Rtilde0
  } else {
    NULL
  }

  start_P <- if (!is.null(start_euc)) {
    start_euc$P
  } else if (!is.null(start_strict)) {
    .prop4i_euc_P_from_strict(start_strict)
  } else {
    NULL
  }

  start_V <- if (!is.null(start_euc) &&
                 all(dim(start_euc$Qe_star) == c(ep$m, kmean))) {
    start_euc$Qe_star
  } else {
    diag(ep$m)[, seq_len(kmean), drop = FALSE]
  }

  start_Be <- if (!is.null(start_euc) &&
                  length(start_euc$Be) == kmean) {
    start_euc$Be
  } else {
    rep(0.05, kmean)
  }

  start_beta <- if (!is.null(start_euc)) {
    start_euc$beta_s
  } else if (!is.null(start_strict)) {
    start_strict$beta_s
  } else {
    0.8
  }

  target_signs <- if (det_constraint == "rotation") 1 else c(1, -1)
  fits <- vector("list", length(target_signs))

  for (j in seq_along(target_signs)) {
    sgn <- target_signs[j]
    Rref <- if (!is.null(start_R) &&
                (if (det(start_R) >= 0) 1 else -1) == sgn) {
      start_R
    } else {
      .prop4_procrustes(y, xs, det_sign = sgn)$R
    }

    Pref <- if (!is.null(start_P)) {
      .prop4_nearest_orthogonal(start_P, det_sign = 1)
    } else {
      .prop4_nearest_orthogonal(Rref, det_sign = 1)
    }

    spec <- .prop4i_euc_mean_spec(
      Rref = Rref,
      Pref = Pref,
      Vref = start_V,
      bstart = start_Be,
      beta_start = start_beta,
      xe_center = ep$center,
      xe_scale = ep$scale,
      intercept = intercept
    )

    fits[[j]] <- .prop4i_euc_fit_vmf_component(
      y = y, xs = xs, xe = xe, spec = spec,
      algorithm = algorithm, opts = opts
    )
    fits[[j]]$det_sign <- sgn
  }

  values <- vapply(fits, function(z) z$objective, numeric(1))
  best <- fits[[which.min(values)]]
  pars <- best$pars
  pred <- best$pred

  dots <- pmax(-1, pmin(1, rowSums(y * pred)))
  kfit <- .prop4_estimate_k(mean(dots), p)
  n <- nrow(y)

  dR <- .prop4_skew_dim(p)
  dP <- .prop4_skew_dim(p)
  dV <- DoF_Stiefel(ep$m, kmean)
  dB <- kmean
  dBeta <- 1L
  dmean <- dR + dP + dV + dB + dBeta
  DoF <- dmean + 1L
  lLik <- n * kfit$loglik_per_obs

  rr <- try(rotated_resid(y, pred, north_pole(p)), silent = TRUE)
  rresids <- if (inherits(rr, "try-error")) NULL else rr[, -1, drop = FALSE]
  if (!is.null(rresids)) {
    attr(rresids, "samehemisphere") <- attr(rr, "samehemisphere")
    colnames(rresids) <- paste0("r", seq_len(ncol(rresids)))
  }

  mean_obj <- pars$mean
  out <- list(
    est = mean_obj, mean = mean_obj,
    Rtilde0 = mean_obj$Rtilde0,
    P = mean_obj$P,
    beta_s = mean_obj$beta_s,
    phi = mean_obj$phi,
    psi = mean_obj$psi,
    rs1 = mean_obj$rs1,
    Qe_star = mean_obj$Qe_star,
    Be = mean_obj$Be,
    k = kfit$k,
    obj = best$objective,
    solution = best$nlopt$solution,
    nlopt = best$nlopt,
    y = y, xs = xs, xe = xe,
    pred = pred,
    rresids = rresids,
    dists = acos(dots),
    DoF = DoF,
    AIC = 2 * DoF - 2 * lLik,
    lLik = lLik,
    start = start,
    component_fits = lapply(fits, function(z) {
      list(
        det_sign = z$det_sign,
        status = z$nlopt$status,
        message = z$nlopt$message,
        objective = z$objective
      )
    }),
    linktype = list(
      type = "Prop4i",
      det_constraint = det_constraint,
      intercept = intercept,
      euclidean = TRUE,
      spherical_constraint = "Bs = beta_s I"
    )
  )
  class(out) <- c(
    "mobius_vMF_prop4i_euc", "mobius_vMF_prop4i", "mobius_vMF", "list"
  )
  out
}

mobius_SvMF_partransport_prelim_prop4i <- function(
    y, xs, xe = NULL, mean = NULL, G0 = NULL,
    G01behaviour = "p1",
    det_constraint = c("orthogonal", "rotation"),
    algorithm = "NLOPT_LD_SLSQP", opts = list(), check = TRUE, ...) {

  det_constraint <- match.arg(det_constraint)
  G01behaviour <- match.arg(G01behaviour, c("p1", "free", "fixed"))
  y <- as.matrix(y)
  xs <- as.matrix(xs)
  if (!is.null(xe) && ncol(as.matrix(xe)) > 0L) {
    stop("Strict Proposition 4(i) has qe = 0; xe must be NULL.")
  }
  if (check) {
    .prop4_check_unit_rows(y, "y")
    .prop4_check_unit_rows(xs, "xs")
  }

  prelim <- mobius_vMF_prop4i(
    y = y, xs = xs, xe = NULL, start = mean,
    det_constraint = det_constraint,
    algorithm = algorithm, opts = opts, check = FALSE
  )
  mean <- prelim$mean
  k <- prelim$k
  p1 <- .prop4i_reference_direction(mean)

  if (G01behaviour == "fixed" && is.null(G0)) {
    stop("At least the first column of G0 must be supplied when G01behaviour = 'fixed'.")
  }
  G01 <- switch(
    G01behaviour,
    p1 = p1,
    free = if (is.null(G0)) p1 else G0[, 1],
    fixed = G0[, 1]
  )

  rresid <- rotated_resid(y, prelim$pred, base = G01)
  if (!is.null(G0) && all(!is.na(G0))) {
    if (G01behaviour == "p1") {
      G0 <- cbind(G01, -jupp_Rmat(G0[, 1], G01) %*%
                        G0[, -1, drop = FALSE])
    }
    aremaining <- SvMF_prelim_scales(rresid, G0)
  } else if (G01behaviour == "free") {
    # A base far from the true G0[,1] leaves the moment estimates near isotropic and strands
    # the joint optimiser, so search for a better one -- see .prop4_choose_G0_start().
    # Only for "free": "p1" pins the base to the identifiable representative and "fixed"
    # takes it from the caller, so in neither is the base ours to choose.
    Gstart <- .prop4_choose_G0_start(as.matrix(y), prelim$pred, prelim$k, G01)
    G0 <- Gstart$G0
    aremaining <- Gstart$a[-1]
  } else {
    G0 <- SvMF_moment_axes(rresid, G01)
    aremaining <- SvMF_prelim_scales(rresid, G0)
  }

  list(
    mean = mean,
    Rtilde0 = mean$Rtilde0,
    psi = mean$psi,
    beta_s = mean$beta_s,
    k = k,
    a = c(1, aremaining),
    G0 = G0,
    nlopt = prelim$nlopt,
    pred = prelim$pred
  )
}

mobius_SvMF_partransport_prelim_prop4i_euc <- function(
    y, xs, xe, mean = NULL, G0 = NULL,
    G01behaviour = "p1",
    det_constraint = c("orthogonal", "rotation"),
    intercept = TRUE,
    algorithm = "NLOPT_LD_SLSQP", opts = list(),
    check = TRUE, ...) {

  det_constraint <- match.arg(det_constraint)
  G01behaviour <- match.arg(
    G01behaviour, c("p1", "free", "fixed")
  )

  prelim <- mobius_vMF_prop4i_euc(
    y = y, xs = xs, xe = xe, start = mean,
    det_constraint = det_constraint,
    intercept = intercept,
    algorithm = algorithm, opts = opts, check = check
  )

  p1 <- prelim$mean$P[, 1L]
  if (G01behaviour == "fixed" && is.null(G0)) {
    stop("At least the first column of G0 must be supplied when G01behaviour='fixed'.")
  }

  G01 <- switch(
    G01behaviour,
    p1 = p1,
    free = if (is.null(G0)) p1 else G0[, 1L],
    fixed = G0[, 1L]
  )

  rresid <- rotated_resid(as.matrix(y), prelim$pred, base = G01)
  if (!is.null(G0) && all(!is.na(G0))) {
    if (G01behaviour == "p1") {
      G0 <- cbind(
        G01,
        -jupp_Rmat(G0[, 1L], G01) %*%
          G0[, -1L, drop = FALSE]
      )
    }
    aremaining <- SvMF_prelim_scales(rresid, G0)
  } else if (G01behaviour == "free") {
    # A base far from the true G0[,1] leaves the moment estimates near isotropic and strands
    # the joint optimiser, so search for a better one -- see .prop4_choose_G0_start().
    # Only for "free": "p1" pins the base to the identifiable representative and "fixed"
    # takes it from the caller, so in neither is the base ours to choose.
    Gstart <- .prop4_choose_G0_start(as.matrix(y), prelim$pred, prelim$k, G01)
    G0 <- Gstart$G0
    aremaining <- Gstart$a[-1]
  } else {
    G0 <- SvMF_moment_axes(rresid, G01)
    aremaining <- SvMF_prelim_scales(rresid, G0)
  }

  list(
    mean = prelim$mean,
    Rtilde0 = prelim$Rtilde0,
    P = prelim$P,
    beta_s = prelim$beta_s,
    Qe_star = prelim$Qe_star,
    Be = prelim$Be,
    k = prelim$k,
    a = c(1, aremaining),
    G0 = G0,
    nlopt = prelim$nlopt,
    pred = prelim$pred
  )
}

.prop4i_joint_loglik <- function(y, xs, mean, k, a, G0,
                                 approximate = TRUE) {
  if (!is.finite(k) || k <= 0 || any(!is.finite(a)) || any(a <= 0)) {
    return(-Inf)
  }
  pred <- try(mobius_link(xs = xs, xe = NULL, param = mean, check = FALSE),
              silent = TRUE)
  if (inherits(pred, "try-error") || any(!is.finite(pred))) return(-Inf)
  yrot <- try(undo_partransport(y, pred, G01 = G0[, 1]), silent = TRUE)
  if (inherits(yrot, "try-error") || any(!is.finite(yrot))) return(-Inf)

  ld <- NULL
  if (approximate) {
    z <- try(ldSvMF_cann(yrot, k = k, a = a, G = G0), silent = TRUE)
    if (!inherits(z, "try-error") && all(is.finite(z))) ld <- z
  }
  if (is.null(ld)) {
    z <- try(SvMF_log_lik_cann(yrot, SvMF_cann(k = k, a = a, G = G0)),
             silent = TRUE)
    if (!inherits(z, "try-error") && all(is.finite(z))) ld <- z
  }
  if (is.null(ld)) return(-Inf)
  sum(ld)
}

.prop4i_build_joint_spec <- function(Rref, prelim, G0reference,
                                     G01behaviour) {
  p <- nrow(Rref)
  G01behaviour <- match.arg(G01behaviour, c("p1", "free", "fixed"))
  dR <- .prop4_skew_dim(p)
  dPsi <- p
  dA <- max(p - 2L, 0L)
  dG <- if (G01behaviour == "free") {
    .prop4_skew_dim(p)
  } else {
    .prop4_skew_dim(p - 1L)
  }

  pos <- 0L
  take <- function(n) {
    if (n == 0L) return(integer(0))
    out <- pos + seq_len(n)
    pos <<- pos + n
    out
  }
  idxR <- take(dR)
  idxPsi <- take(dPsi)
  idxK <- take(1L)
  idxA <- take(dA)
  idxG <- take(dG)

  start_mean <- as_mobius_link_prop4i(prelim$mean, check = FALSE)
  Gstart <- .prop4_nearest_orthogonal(prelim$G0, det_sign = 1)
  Gcoord <- if (is.null(G0reference)) {
    Gstart
  } else {
    .prop4_nearest_orthogonal(G0reference, det_sign = 1)
  }

  mean_ref <- mobius_link_prop4i(Rref, psi = start_mean$psi, check = FALSE)
  target_base <- if (G01behaviour == "p1") {
    .prop4i_reference_direction(mean_ref)
  } else {
    Gstart[, 1]
  }

  if (G01behaviour == "free") {
    Gref <- Gcoord
    Q0 <- Gstart %*% t(Gref)
    thetaG0 <- .prop4_inverse_cayley(Q0)
    Ganchor <- NULL
    fixed_base <- NULL
  } else {
    Ganchor <- .prop4_rebase_G0(Gcoord, target_base)
    if (is.null(Ganchor)) Ganchor <- .prop4_rebase_G0(Gstart, target_base)
    if (is.null(Ganchor)) stop("Could not construct the G0 reference frame.")

    Gstart2 <- .prop4_rebase_G0(Gstart, target_base)
    if (is.null(Gstart2)) stop("Could not rebase the starting G0 frame.")
    H0 <- crossprod(Ganchor[, -1, drop = FALSE],
                    Gstart2[, -1, drop = FALSE])
    H0 <- .prop4_nearest_orthogonal(H0, det_sign = 1)
    thetaG0 <- .prop4_inverse_cayley(H0)
    Gref <- NULL
    fixed_base <- if (G01behaviour == "fixed") target_base else NULL
  }

  par0 <- numeric(pos)
  par0[idxR] <- 0
  par0[idxPsi] <- .prop4i_psi_to_u(start_mean$psi)
  par0[idxK] <- log(max(prelim$k, 1e-8))
  if (dA > 0L) par0[idxA] <- .prop4_scales_to_eta(prelim$a)
  if (dG > 0L) par0[idxG] <- thetaG0

  lower <- rep(-Inf, pos)
  upper <- rep( Inf, pos)
  lower[idxK] <- log(1e-8)
  upper[idxK] <- log(1e5)
  if (dA > 0L) {
    lower[idxA] <- -10
    upper[idxA] <- 10
  }

  list(
    p = p, Rref = Rref, a1 = prelim$a[1],
    G01behaviour = G01behaviour,
    Gref = Gref, Ganchor = Ganchor, fixed_base = fixed_base,
    idxR = idxR, idxPsi = idxPsi, idxK = idxK,
    idxA = idxA, idxG = idxG,
    par0 = par0, lower = lower, upper = upper
  )
}

.prop4i_unpack_joint <- function(theta, spec) {
  p <- spec$p
  R <- .prop4_cayley(theta[spec$idxR], p) %*% spec$Rref
  psi <- .prop4i_u_to_psi(theta[spec$idxPsi])
  mean <- mobius_link_prop4i(R, psi = psi, check = FALSE)
  k <- exp(theta[spec$idxK])
  eta <- if (length(spec$idxA)) theta[spec$idxA] else numeric(0)
  a <- .prop4_eta_to_scales(eta, spec$a1, p)

  if (spec$G01behaviour == "free") {
    G0 <- .prop4_cayley(theta[spec$idxG], p) %*% spec$Gref
  } else {
    Gbase <- spec$Ganchor
    if (spec$G01behaviour == "p1") {
      target <- .prop4i_reference_direction(mean)
      Gbase <- .prop4_rebase_G0(spec$Ganchor, target)
      if (is.null(Gbase)) return(NULL)
    }
    H <- if (p == 2L) {
      matrix(1, 1, 1)
    } else {
      .prop4_cayley(theta[spec$idxG], p - 1L)
    }
    B <- diag(p)
    B[-1, -1] <- H
    G0 <- Gbase %*% B
  }

  list(mean = mean, Rtilde0 = R, psi = psi, k = k, a = a, G0 = G0)
}

.prop4i_fit_component <- function(y, xs, Rref, prelim, G0reference,
                                  G01behaviour, algorithm, opts,
                                  approximate = TRUE) {
  spec <- .prop4i_build_joint_spec(
    Rref = Rref, prelim = prelim,
    G0reference = G0reference,
    G01behaviour = G01behaviour
  )

  eval_f <- function(theta) {
    pars <- .prop4i_unpack_joint(theta, spec)
    if (is.null(pars)) return(1e100)
    ll <- .prop4i_joint_loglik(
      y, xs, pars$mean, pars$k, pars$a, pars$G0,
      approximate = approximate
    )
    if (!is.finite(ll)) return(1e100)
    -ll / nrow(y)
  }

  default_opts <- list(
    algorithm = algorithm,
    xtol_rel = 1e-9,
    maxeval = 10000,
    print_level = 0
  )
  combined_opts <- utils::modifyList(default_opts, opts)
  uses_gradient <- grepl("^NLOPT_LD_", combined_opts$algorithm)
  eval_nloptr <- if (uses_gradient) {
    function(theta) {
      list(
        objective = eval_f(theta),
        gradient = .prop4_fd_gradient(
          eval_f, theta, spec$lower, spec$upper
        )
      )
    }
  } else {
    eval_f
  }

  fit <- try(
    nloptr::nloptr(
      x0 = spec$par0,
      eval_f = eval_nloptr,
      lb = spec$lower,
      ub = spec$upper,
      opts = combined_opts
    ),
    silent = TRUE
  )

  # Keep the valid preliminary point if the numerical optimiser fails.
  theta_best <- spec$par0
  f_best <- eval_f(theta_best)
  if (!inherits(fit, "try-error") && is.finite(fit$objective) &&
      fit$objective < f_best) {
    theta_best <- fit$solution
    f_best <- fit$objective
  }
  pars <- .prop4i_unpack_joint(theta_best, spec)
  ll <- if (is.null(pars)) -Inf else .prop4i_joint_loglik(
    y, xs, pars$mean, pars$k, pars$a, pars$G0,
    approximate = approximate
  )

  if (inherits(fit, "try-error")) {
    fit <- list(status = -999L, message = as.character(fit),
                objective = f_best, solution = theta_best)
  } else {
    fit$objective <- f_best
    fit$solution <- theta_best
  }
  list(nlopt = fit, pars = pars, lLik_approx = ll, spec = spec)
}

mobius_SvMF_joint_fit_prop4i <- function(
    y, xs, xe = NULL, mean, k, a, G0,
    G0reference = NULL, G01behaviour = "p1",
    det_constraint = c("orthogonal", "rotation"),
    algorithm = "NLOPT_LD_SLSQP", opts = list(), check = TRUE) {

  det_constraint <- match.arg(det_constraint)
  G01behaviour <- match.arg(G01behaviour, c("p1", "free", "fixed"))
  y <- as.matrix(y)
  xs <- as.matrix(xs)
  storage.mode(y) <- storage.mode(xs) <- "double"

  if (!is.null(xe) && ncol(as.matrix(xe)) > 0L) {
    stop("Strict Proposition 4(i) has qe = 0; xe must be NULL.")
  }
  if (nrow(y) != nrow(xs) || ncol(y) != ncol(xs)) {
    stop("For Proposition 4(i), y and xs must have identical dimensions.")
  }
  p <- ncol(y)
  if (check) {
    .prop4_check_unit_rows(y, "y")
    .prop4_check_unit_rows(xs, "xs")
  }

  start_mean <- as_mobius_link_prop4i(mean, check = FALSE)
  start_R <- start_mean$Rtilde0
  start_a <- .prop4_normalise_scales(a, p, a1 = a[1])
  start_G0 <- .prop4_nearest_orthogonal(G0, det_sign = 1)
  initial <- list(
    mean = start_mean, Rtilde0 = start_R, psi = start_mean$psi,
    k = k, a = start_a, G0 = start_G0
  )

  target_signs <- if (det_constraint == "rotation") 1 else c(1, -1)
  fits <- vector("list", length(target_signs))
  for (j in seq_along(target_signs)) {
    sgn <- target_signs[j]
    Rref <- if ((if (det(start_R) >= 0) 1 else -1) == sgn) {
      start_R
    } else {
      .prop4_procrustes(y, xs, det_sign = sgn)$R
    }
    pre_j <- initial
    pre_j$mean <- mobius_link_prop4i(
      Rref, psi = start_mean$psi, check = FALSE
    )
    pre_j$Rtilde0 <- Rref
    pre_j$psi <- start_mean$psi
    if (G01behaviour == "p1") {
      newbase <- .prop4i_reference_direction(pre_j$mean)
      rebased <- .prop4_rebase_G0(pre_j$G0, newbase)
      if (is.null(rebased)) stop("Could not rebase G0 for Proposition 4(i).")
      pre_j$G0 <- rebased
    }
    fits[[j]] <- .prop4i_fit_component(
      y = y, xs = xs, Rref = Rref, prelim = pre_j,
      G0reference = if (is.null(G0reference)) pre_j$G0 else G0reference,
      G01behaviour = G01behaviour,
      algorithm = algorithm, opts = opts,
      approximate = TRUE
    )
    fits[[j]]$det_sign <- sgn
  }

  ll_approx <- vapply(fits, function(z) z$lLik_approx, numeric(1))
  if (all(!is.finite(ll_approx))) {
    stop("All Proposition 4(i) optimisations failed.")
  }
  best <- fits[[which.max(ll_approx)]]
  pars <- best$pars

  # Same residual-axis identifiability convention as the existing models.
  aord <- order(pars$a[-1], decreasing = TRUE)
  pars$a[-1] <- pars$a[-1][aord]
  pars$G0[, -1] <- pars$G0[, -1, drop = FALSE][, aord, drop = FALSE]
  pars$G0[, -1] <- standardise_col_signs(pars$G0[, -1, drop = FALSE])
  if (det(pars$G0) < 0) pars$G0[, p] <- -pars$G0[, p]

  pred <- mobius_link(xs = xs, xe = NULL, param = pars$mean, check = FALSE)
  if (p != 3L) {
    kres <- mobius_SvMF_konly(y = y, ymean = pred,
                              a = pars$a, G0 = pars$G0)
    pars$k <- kres$k
    lLik <- kres$lLik
  } else {
    yrot <- undo_partransport(y, pred, G01 = pars$G0[, 1])
    lLik <- sum(SvMF_log_lik_cann(
      yrot, SvMF_cann(k = pars$k, a = pars$a, G = pars$G0)
    ))
  }
  SvMF_cann_check(SvMF_cann(k = pars$k, a = pars$a, G = pars$G0))

  dots <- pmax(-1, pmin(1, rowSums(y * pred)))
  rrI <- rotated_resid(y, pred, north_pole(p))
  rresids_I <- rrI[, -1, drop = FALSE]
  attr(rresids_I, "samehemisphere") <- attr(rrI, "samehemisphere")
  colnames(rresids_I) <- paste0("r", seq_len(ncol(rresids_I)))
  rresids_G0 <- resid_SvMF_partransport(y, pred, G0 = pars$G0, scale = FALSE)
  rresids_std <- resid_SvMF_partransport(
    y, pred, pars$k, pars$a, pars$G0, scale = TRUE
  )

  dmean <- .prop4_skew_dim(p) + p
  dscale <- max(p - 2L, 0L)
  dG0 <- if (G01behaviour == "free") {
    .prop4_skew_dim(p)
  } else {
    .prop4_skew_dim(p - 1L)
  }
  DoF <- dmean + 1L + dscale + dG0
  mean <- pars$mean

  out <- list(
    mean = mean, est = mean,
    Rtilde0 = mean$Rtilde0,
    psi = mean$psi, phi = mean$phi,
    beta_s = mean$beta_s, rs1 = mean$rs1,
    k = pars$k, a = pars$a, G0 = pars$G0,
    obj = -lLik / nrow(y), nlopt = best$nlopt,
    y = y, xs = xs, xe = NULL, pred = pred,
    rresids_I = rresids_I,
    rresids_G0 = rresids_G0,
    rresids_std = rresids_std,
    dists = acos(dots),
    DoF = DoF, AIC = 2 * DoF - 2 * lLik, lLik = lLik,
    initial = initial,
    component_fits = lapply(fits, function(z) {
      list(det_sign = z$det_sign, status = z$nlopt$status,
           message = z$nlopt$message, lLik_approx = z$lLik_approx)
    }),
    linktype = list(
      type = "Prop4i", det_constraint = det_constraint,
      G01behaviour = G01behaviour,
      intercept = FALSE, euclidean = FALSE
    )
  )
  class(out) <- c("mobius_SvMF_prop4i", "mobius_SvMF", "list")
  out
}

.prop4i_euc_joint_spec <- function(
    mean_ref, prelim, G0reference, G01behaviour) {

  mean_ref <- as_mobius_link_prop4i(mean_ref, check = FALSE)
  p <- nrow(mean_ref$Rtilde0)

  mspec <- .prop4i_euc_mean_spec(
    Rref = mean_ref$Rtilde0,
    Pref = mean_ref$P,
    Vref = mean_ref$Qe_star,
    bstart = mean_ref$Be,
    beta_start = mean_ref$beta_s,
    xe_center = mean_ref$xe_center,
    xe_scale = mean_ref$xe_scale,
    intercept = mean_ref$intercept
  )

  dA <- max(p - 2L, 0L)
  dG <- if (G01behaviour == "free") {
    .prop4_skew_dim(p)
  } else {
    .prop4_skew_dim(p - 1L)
  }

  pos <- length(mspec$par0)
  take <- function(n) {
    if (n == 0L) return(integer(0))
    z <- pos + seq_len(n)
    pos <<- pos + n
    z
  }

  idxK <- take(1L)
  idxA <- take(dA)
  idxG <- take(dG)

  Gstart <- .prop4_nearest_orthogonal(prelim$G0, det_sign = 1)
  Gcoord <- if (is.null(G0reference)) {
    Gstart
  } else {
    .prop4_nearest_orthogonal(G0reference, det_sign = 1)
  }

  if (G01behaviour == "free") {
    Gref <- Gcoord
    thetaG0 <- .prop4_inverse_cayley(Gstart %*% t(Gref))
    Ganchor <- NULL
  } else {
    target <- if (G01behaviour == "p1") {
      mean_ref$P[, 1L]
    } else {
      Gstart[, 1L]
    }

    Ganchor <- .prop4_rebase_G0(Gcoord, target)
    if (is.null(Ganchor)) {
      Ganchor <- .prop4_rebase_G0(Gstart, target)
    }
    if (is.null(Ganchor)) {
      stop("Could not construct the G0 reference frame.")
    }

    Gstart2 <- .prop4_rebase_G0(Gstart, target)
    H0 <- crossprod(
      Ganchor[, -1L, drop = FALSE],
      Gstart2[, -1L, drop = FALSE]
    )
    H0 <- .prop4_nearest_orthogonal(H0, det_sign = 1)
    thetaG0 <- .prop4_inverse_cayley(H0)
    Gref <- NULL
  }

  par0 <- c(mspec$par0, rep(0, pos - length(mspec$par0)))
  par0[idxK] <- log(max(prelim$k, 1e-8))
  if (dA > 0L) {
    par0[idxA] <- .prop4_scales_to_eta(prelim$a)
  }
  if (dG > 0L) {
    par0[idxG] <- thetaG0
  }

  lower <- c(mspec$lower, rep(-Inf, pos - length(mspec$lower)))
  upper <- c(mspec$upper, rep( Inf, pos - length(mspec$upper)))
  lower[idxK] <- log(1e-8)
  upper[idxK] <- log(1e5)
  if (dA > 0L) {
    lower[idxA] <- -10
    upper[idxA] <- 10
  }

  list(
    mean_spec = mspec,
    mean_length = length(mspec$par0),
    p = p,
    a1 = prelim$a[1L],
    G01behaviour = G01behaviour,
    Gref = Gref,
    Ganchor = Ganchor,
    idxK = idxK,
    idxA = idxA,
    idxG = idxG,
    par0 = par0,
    lower = lower,
    upper = upper
  )
}

.prop4i_euc_unpack_joint <- function(theta, spec) {
  mp <- .prop4i_euc_unpack_mean(
    theta[seq_len(spec$mean_length)], spec$mean_spec
  )

  kappa <- exp(theta[spec$idxK])
  eta <- if (length(spec$idxA)) theta[spec$idxA] else numeric(0)
  a <- .prop4_eta_to_scales(eta, spec$a1, spec$p)

  if (spec$G01behaviour == "free") {
    G0 <- .prop4_cayley(theta[spec$idxG], spec$p) %*% spec$Gref
  } else {
    Gbase <- spec$Ganchor
    if (spec$G01behaviour == "p1") {
      Gbase <- .prop4_rebase_G0(
        spec$Ganchor, mp$P[, 1L]
      )
      if (is.null(Gbase)) return(NULL)
    }

    H <- if (spec$p == 2L) {
      matrix(1, 1, 1)
    } else {
      .prop4_cayley(theta[spec$idxG], spec$p - 1L)
    }

    B <- diag(spec$p)
    B[-1L, -1L] <- H
    G0 <- Gbase %*% B
  }

  c(mp, list(k = kappa, a = a, G0 = G0))
}

.prop4i_euc_joint_loglik <- function(
    y, xs, xe, mean, k, a, G0, approximate = TRUE) {

  if (!is.finite(k) || k <= 0 ||
      any(!is.finite(a)) || any(a <= 0)) {
    return(-Inf)
  }

  pred <- try(
    mobius_link(xs = xs, xe = xe, param = mean, check = FALSE),
    silent = TRUE
  )
  if (inherits(pred, "try-error") || any(!is.finite(pred))) {
    return(-Inf)
  }

  yrot <- try(
    undo_partransport(y, pred, G01 = G0[, 1L]),
    silent = TRUE
  )
  if (inherits(yrot, "try-error") || any(!is.finite(yrot))) {
    return(-Inf)
  }

  ld <- NULL
  if (approximate) {
    z <- try(
      ldSvMF_cann(yrot, k = k, a = a, G = G0),
      silent = TRUE
    )
    if (!inherits(z, "try-error") && all(is.finite(z))) {
      ld <- z
    }
  }

  if (is.null(ld)) {
    z <- try(
      SvMF_log_lik_cann(
        yrot, SvMF_cann(k = k, a = a, G = G0)
      ),
      silent = TRUE
    )
    if (!inherits(z, "try-error") && all(is.finite(z))) {
      ld <- z
    }
  }

  if (is.null(ld)) return(-Inf)
  sum(ld)
}

.prop4i_euc_fit_joint_component <- function(
    y, xs, xe, mean_ref, prelim, G0reference,
    G01behaviour, algorithm, opts) {

  spec <- .prop4i_euc_joint_spec(
    mean_ref = mean_ref,
    prelim = prelim,
    G0reference = G0reference,
    G01behaviour = G01behaviour
  )

  eval_f <- function(theta) {
    pars <- .prop4i_euc_unpack_joint(theta, spec)
    if (is.null(pars)) return(1e100)

    ll <- .prop4i_euc_joint_loglik(
      y, xs, xe,
      pars$mean, pars$k, pars$a, pars$G0,
      approximate = TRUE
    )
    if (!is.finite(ll)) return(1e100)
    -ll / nrow(y)
  }

  default_opts <- list(
    algorithm = algorithm,
    xtol_rel = 1e-9,
    maxeval = 10000,
    print_level = 0
  )
  combined_opts <- utils::modifyList(default_opts, opts)
  uses_gradient <- grepl("^NLOPT_LD_", combined_opts$algorithm)

  eval_nloptr <- if (uses_gradient) {
    function(theta) {
      list(
        objective = eval_f(theta),
        gradient = .prop4_fd_gradient(
          eval_f, theta, spec$lower, spec$upper
        )
      )
    }
  } else {
    eval_f
  }

  fit <- try(
    nloptr::nloptr(
      x0 = spec$par0,
      eval_f = eval_nloptr,
      lb = spec$lower,
      ub = spec$upper,
      opts = combined_opts
    ),
    silent = TRUE
  )

  # Preserve the valid preliminary point if the optimiser fails.
  theta_best <- spec$par0
  f_best <- eval_f(theta_best)
  if (!inherits(fit, "try-error") &&
      is.finite(fit$objective) && fit$objective < f_best) {
    theta_best <- fit$solution
    f_best <- fit$objective
  }

  pars <- .prop4i_euc_unpack_joint(theta_best, spec)
  ll <- if (is.null(pars)) {
    -Inf
  } else {
    .prop4i_euc_joint_loglik(
      y, xs, xe,
      pars$mean, pars$k, pars$a, pars$G0,
      approximate = TRUE
    )
  }

  if (inherits(fit, "try-error")) {
    fit <- list(
      status = -999L,
      message = as.character(fit),
      objective = f_best,
      solution = theta_best
    )
  } else {
    fit$objective <- f_best
    fit$solution <- theta_best
  }

  list(
    nlopt = fit,
    pars = pars,
    lLik_approx = ll,
    spec = spec
  )
}

mobius_SvMF_joint_fit_prop4i_euc <- function(
    y, xs, xe, mean, k, a, G0,
    G0reference = NULL, G01behaviour = "p1",
    det_constraint = c("orthogonal", "rotation"),
    intercept = TRUE,
    algorithm = "NLOPT_LD_SLSQP", opts = list(),
    check = TRUE) {

  det_constraint <- match.arg(det_constraint)
  G01behaviour <- match.arg(
    G01behaviour, c("p1", "free", "fixed")
  )

  y <- as.matrix(y)
  xs <- as.matrix(xs)
  xe <- as.matrix(xe)
  p <- ncol(y)

  if (nrow(y) != nrow(xs) || nrow(y) != nrow(xe) ||
      ncol(xs) != p) {
    stop("For Proposition 4(i)+Euclidean, y and xs must be n by p and xe must have n rows.")
  }
  if (check) {
    .prop4_check_unit_rows(y, "y")
    .prop4_check_unit_rows(xs, "xs")
  }

  start_mean <- as_mobius_link_prop4i(mean, check = FALSE)
  start_a <- .prop4_normalise_scales(a, p, a1 = a[1L])
  start_G0 <- .prop4_nearest_orthogonal(G0, det_sign = 1)

  initial <- list(
    mean = start_mean,
    Rtilde0 = start_mean$Rtilde0,
    P = start_mean$P,
    beta_s = start_mean$beta_s,
    Qe_star = start_mean$Qe_star,
    Be = start_mean$Be,
    k = k,
    a = start_a,
    G0 = start_G0
  )

  target_signs <- if (det_constraint == "rotation") 1 else c(1, -1)
  fits <- vector("list", length(target_signs))

  for (j in seq_along(target_signs)) {
    sgn <- target_signs[j]

    Rref <- if ((if (det(start_mean$Rtilde0) >= 0) 1 else -1) == sgn) {
      start_mean$Rtilde0
    } else {
      .prop4_procrustes(y, xs, det_sign = sgn)$R
    }

    mean_ref <- mobius_link_prop4i(
      Rtilde0 = Rref,
      P = start_mean$P,
      beta_s = start_mean$beta_s,
      Qe_star = start_mean$Qe_star,
      Be = start_mean$Be,
      xe_center = start_mean$xe_center,
      xe_scale = start_mean$xe_scale,
      intercept = start_mean$intercept,
      check = FALSE
    )

    pre_j <- initial
    pre_j$mean <- mean_ref
    pre_j$Rtilde0 <- Rref

    if (G01behaviour == "p1") {
      rebased <- .prop4_rebase_G0(
        pre_j$G0, mean_ref$P[, 1L]
      )
      if (is.null(rebased)) {
        stop("Could not rebase G0 to P[,1].")
      }
      pre_j$G0 <- rebased
    }

    fits[[j]] <- .prop4i_euc_fit_joint_component(
      y = y, xs = xs, xe = xe,
      mean_ref = mean_ref,
      prelim = pre_j,
      G0reference = if (is.null(G0reference)) pre_j$G0 else G0reference,
      G01behaviour = G01behaviour,
      algorithm = algorithm,
      opts = opts
    )
    fits[[j]]$det_sign <- sgn
  }

  lla <- vapply(fits, function(z) z$lLik_approx, numeric(1))
  if (all(!is.finite(lla))) {
    stop("All Proposition 4(i)+Euclidean optimisations failed.")
  }

  best <- fits[[which.max(lla)]]
  pars <- best$pars

  # Same residual-axis identifiability convention as the existing models.
  aord <- order(pars$a[-1L], decreasing = TRUE)
  pars$a[-1L] <- pars$a[-1L][aord]
  pars$G0[, -1L] <- pars$G0[, -1L, drop = FALSE][, aord, drop = FALSE]
  pars$G0[, -1L] <- standardise_col_signs(
    pars$G0[, -1L, drop = FALSE]
  )
  if (det(pars$G0) < 0) {
    pars$G0[, p] <- -pars$G0[, p]
  }

  pred <- mobius_link(
    xs = xs, xe = xe, param = pars$mean, check = FALSE
  )

  if (p != 3L) {
    kres <- mobius_SvMF_konly(
      y = y, ymean = pred,
      a = pars$a, G0 = pars$G0
    )
    pars$k <- kres$k
    lLik <- kres$lLik
  } else {
    yrot <- undo_partransport(
      y, pred, G01 = pars$G0[, 1L]
    )
    lLik <- sum(
      SvMF_log_lik_cann(
        yrot,
        SvMF_cann(k = pars$k, a = pars$a, G = pars$G0)
      )
    )
  }

  SvMF_cann_check(
    SvMF_cann(k = pars$k, a = pars$a, G = pars$G0)
  )

  dots <- pmax(-1, pmin(1, rowSums(y * pred)))
  rrI <- rotated_resid(y, pred, north_pole(p))
  rresids_I <- rrI[, -1L, drop = FALSE]
  attr(rresids_I, "samehemisphere") <- attr(rrI, "samehemisphere")
  colnames(rresids_I) <- paste0("r", seq_len(ncol(rresids_I)))

  rresids_G0 <- resid_SvMF_partransport(
    y, pred, G0 = pars$G0, scale = FALSE
  )
  rresids_std <- resid_SvMF_partransport(
    y, pred, pars$k, pars$a, pars$G0, scale = TRUE
  )

  ep_m <- nrow(pars$mean$Qe_star)
  kmean <- p - 1L

  dmean <- .prop4_skew_dim(p) +       # Rtilde0
    .prop4_skew_dim(p) +              # P
    DoF_Stiefel(ep_m, kmean) + # Qe_star
    kmean +                             # Be
    1L                                  # beta_s

  dscale <- max(p - 2L, 0L)
  dG0 <- if (G01behaviour == "free") {
    .prop4_skew_dim(p)
  } else {
    .prop4_skew_dim(p - 1L)
  }
  DoF <- dmean + 1L + dscale + dG0

  mean_obj <- pars$mean
  out <- list(
    mean = mean_obj, est = mean_obj,
    Rtilde0 = mean_obj$Rtilde0,
    P = mean_obj$P,
    beta_s = mean_obj$beta_s,
    phi = mean_obj$phi,
    psi = mean_obj$psi,
    rs1 = mean_obj$rs1,
    Qe_star = mean_obj$Qe_star,
    Be = mean_obj$Be,
    k = pars$k,
    a = pars$a,
    G0 = pars$G0,
    obj = -lLik / nrow(y),
    nlopt = best$nlopt,
    y = y, xs = xs, xe = xe,
    pred = pred,
    rresids_I = rresids_I,
    rresids_G0 = rresids_G0,
    rresids_std = rresids_std,
    dists = acos(dots),
    DoF = DoF,
    AIC = 2 * DoF - 2 * lLik,
    lLik = lLik,
    initial = initial,
    component_fits = lapply(fits, function(z) {
      list(
        det_sign = z$det_sign,
        status = z$nlopt$status,
        message = z$nlopt$message,
        lLik_approx = z$lLik_approx
      )
    }),
    linktype = list(
      type = "Prop4i",
      det_constraint = det_constraint,
      G01behaviour = G01behaviour,
      intercept = intercept,
      euclidean = TRUE,
      spherical_constraint = "Bs = beta_s I"
    )
  )

  class(out) <- c(
    "mobius_SvMF_prop4i_euc",
    "mobius_SvMF_prop4i",
    "mobius_SvMF",
    "list"
  )
  out
}

mobius_SvMF_prop4i <- function(
    y, xs, xe = NULL, mean = NULL, k = NULL, a = NULL,
    G0 = NULL, G0reference = NULL, G01behaviour = "p1",
    doprelim = TRUE,
    det_constraint = c("orthogonal", "rotation"),
    algorithm = "NLOPT_LD_SLSQP", opts = list(), check = TRUE, ...) {

  if (!is.null(xe) && ncol(as.matrix(xe)) > 0L) {
    stop("Strict Proposition 4(i) has qe = 0; xe must be NULL.")
  }
  det_constraint <- match.arg(det_constraint)
  G01behaviour <- match.arg(G01behaviour, c("p1", "free", "fixed"))
  dots <- list(...)
  if (length(dots)) opts <- utils::modifyList(opts, dots)

  if (doprelim) {
    preest <- mobius_SvMF_partransport_prelim_prop4i(
      y = y, xs = xs, xe = NULL, mean = mean, G0 = G0,
      G01behaviour = G01behaviour,
      det_constraint = det_constraint,
      algorithm = algorithm, opts = opts, check = check
    )
  } else {
    if (is.null(mean) || is.null(k) || is.null(a) || is.null(G0)) {
      stop("When doprelim = FALSE, mean, k, a, and G0 must all be supplied.")
    }
    mean <- as_mobius_link_prop4i(mean)
    preest <- list(
      mean = mean, Rtilde0 = mean$Rtilde0, psi = mean$psi,
      beta_s = mean$beta_s, k = k, a = a, G0 = G0
    )
  }

  finalest <- mobius_SvMF_joint_fit_prop4i(
    y = y, xs = xs, xe = NULL,
    mean = preest$mean,
    k = if (is.null(k)) preest$k else k,
    a = if (is.null(a)) preest$a else a,
    G0 = preest$G0,
    G0reference = if (is.null(G0reference)) preest$G0 else G0reference,
    G01behaviour = G01behaviour,
    det_constraint = det_constraint,
    algorithm = algorithm, opts = opts, check = check
  )
  finalest$preest <- preest
  finalest
}

mobius_SvMF_prop4i_euc <- function(
    y, xs, xe, mean = NULL, k = NULL, a = NULL,
    G0 = NULL, G0reference = NULL, G01behaviour = "p1",
    doprelim = TRUE,
    det_constraint = c("orthogonal", "rotation"),
    intercept = TRUE,
    algorithm = "NLOPT_LD_SLSQP", opts = list(),
    check = TRUE, ...) {

  det_constraint <- match.arg(det_constraint)
  G01behaviour <- match.arg(
    G01behaviour, c("p1", "free", "fixed")
  )

  dots <- list(...)
  if (length(dots)) {
    opts <- utils::modifyList(opts, dots)
  }

  if (doprelim) {
    preest <- mobius_SvMF_partransport_prelim_prop4i_euc(
      y = y, xs = xs, xe = xe,
      mean = mean, G0 = G0,
      G01behaviour = G01behaviour,
      det_constraint = det_constraint,
      intercept = intercept,
      algorithm = algorithm,
      opts = opts,
      check = check
    )
  } else {
    if (is.null(mean) || is.null(k) ||
        is.null(a) || is.null(G0)) {
      stop("When doprelim=FALSE, mean, k, a, and G0 must all be supplied.")
    }

    mean <- as_mobius_link_prop4i(mean, check = FALSE)
    preest <- list(
      mean = mean,
      Rtilde0 = mean$Rtilde0,
      P = mean$P,
      beta_s = mean$beta_s,
      Qe_star = mean$Qe_star,
      Be = mean$Be,
      k = k,
      a = a,
      G0 = G0
    )
  }

  finalest <- mobius_SvMF_joint_fit_prop4i_euc(
    y = y, xs = xs, xe = xe,
    mean = preest$mean,
    k = if (is.null(k)) preest$k else k,
    a = if (is.null(a)) preest$a else a,
    G0 = preest$G0,
    G0reference = if (is.null(G0reference)) preest$G0 else G0reference,
    G01behaviour = G01behaviour,
    det_constraint = det_constraint,
    intercept = intercept,
    algorithm = algorithm,
    opts = opts,
    check = check
  )

  finalest$preest <- preest
  finalest
}

#' @export
print.mobius_vMF_prop4i <- function(x, ...) {
  cat("vMF regression with Proposition 4(i) spherical M\u00f6bius link\n")
  cat("  beta_s =", format(x$beta_s, digits = 6), "\n")
  cat("  phi =", format(x$phi, digits = 6), "\n")
  cat("  det(Rtilde0) =", format(det(x$Rtilde0), digits = 6), "\n")
  cat("  k =", format(x$k, digits = 6), "\n")
  cat("  log-likelihood =", format(x$lLik, digits = 7), "\n")
  cat("  AIC =", format(x$AIC, digits = 7), "\n")
  invisible(x)
}

#' @export
print.mobius_SvMF_prop4i <- function(x, ...) {
  cat("SvMF regression with Proposition 4(i) spherical M\u00f6bius link\n")
  cat("  G0 estimation: parallel-transport/Cayley procedure\n")
  cat("  G01behaviour =", x$linktype$G01behaviour, "\n")
  cat("  beta_s =", format(x$beta_s, digits = 6), "\n")
  cat("  phi =", format(x$phi, digits = 6), "\n")
  cat("  det(Rtilde0) =", format(det(x$Rtilde0), digits = 6), "\n")
  cat("  k =", format(x$k, digits = 6), "\n")
  cat("  a =", paste(format(x$a, digits = 5), collapse = ", "), "\n")
  cat("  log-likelihood =", format(x$lLik, digits = 7), "\n")
  cat("  AIC =", format(x$AIC, digits = 7), "\n")
  invisible(x)
}

#' @export
predict.mobius_vMF_prop4i <- function(object, newdata = NULL,
                                       xs = newdata, ...) {
  if (is.null(xs)) xs <- object$xs
  mobius_link(xs = as.matrix(xs), xe = NULL, param = object$mean)
}

#' @export
predict.mobius_SvMF_prop4i <- function(object, newdata = NULL,
                                        xs = newdata, ...) {
  predict.mobius_vMF_prop4i(object, newdata = newdata, xs = xs, ...)
}

#' @export
print.mobius_vMF_prop4i_euc <- function(x, ...) {
  cat("vMF regression with Proposition 4(i)-constrained spherical link and LinEuc covariates\n")
  cat("  beta_s =", format(x$beta_s, digits = 6), "\n")
  cat("  phi =", format(x$phi, digits = 6), "\n")
  cat("  det(Rtilde0) =", format(det(x$Rtilde0), digits = 6), "\n")
  cat("  Be =", paste(format(x$Be, digits = 5), collapse = ", "), "\n")
  cat("  k =", format(x$k, digits = 6), "\n")
  cat("  log-likelihood =", format(x$lLik, digits = 7), "\n")
  cat("  AIC =", format(x$AIC, digits = 7), "\n")
  invisible(x)
}

#' @export
print.mobius_SvMF_prop4i_euc <- function(x, ...) {
  cat("SvMF regression with Proposition 4(i)-constrained spherical link and LinEuc covariates\n")
  cat("  G0 estimation: parallel-transport/Cayley procedure\n")
  cat("  G01behaviour =", x$linktype$G01behaviour, "\n")
  cat("  beta_s =", format(x$beta_s, digits = 6), "\n")
  cat("  phi =", format(x$phi, digits = 6), "\n")
  cat("  det(Rtilde0) =", format(det(x$Rtilde0), digits = 6), "\n")
  cat("  Be =", paste(format(x$Be, digits = 5), collapse = ", "), "\n")
  cat("  k =", format(x$k, digits = 6), "\n")
  cat("  a =", paste(format(x$a, digits = 5), collapse = ", "), "\n")
  cat("  log-likelihood =", format(x$lLik, digits = 7), "\n")
  cat("  AIC =", format(x$AIC, digits = 7), "\n")
  invisible(x)
}

#' @export
predict.mobius_vMF_prop4i_euc <- function(
    object, newdata = NULL, xs = NULL, xe = NULL, ...) {

  if (is.null(xs)) xs <- if (is.null(newdata)) object$xs else newdata
  if (is.null(xe)) xe <- object$xe

  mobius_link(
    xs = as.matrix(xs),
    xe = as.matrix(xe),
    param = object$mean
  )
}

#' @export
predict.mobius_SvMF_prop4i_euc <- function(
    object, newdata = NULL, xs = NULL, xe = NULL, ...) {

  predict.mobius_vMF_prop4i_euc(
    object, newdata = newdata, xs = xs, xe = xe, ...
  )
}
