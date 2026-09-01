# Proposition 4(ii): the scaled Moebius link with Bs = I_(p-1).
#
# Spherical covariate only (p = qs, qe = 0):
#     mu(xs) = Rtilde0 xs,   Rtilde0 in O(p)
# With Euclidean covariates the LinEuc construction is retained while Bs stays fixed at
# the identity:
#     mu(xs, xe) = P iSp{ Sp(Qs^T xs) + Be Qe*^T xe_tilde },   Rtilde0 = P Qs^T
# where xe_tilde holds the standardised Euclidean covariates and (when intercept = TRUE)
# an intercept, Qe* has orthonormal columns and Be is diagonal. Setting the Euclidean
# contribution to zero recovers the spherical-only link exactly.
#
# G0 is estimated exactly as in the unconstrained code: preliminary moment estimation from
# parallel-transported residuals, then joint likelihood optimisation in Cayley coordinates
# relative to G0reference, then permutation/sign standardisation of (a_j, G0[,j]), j >= 2.
#
# Reached through mobius_vMF(type = "Prop4ii") / mobius_SvMF(type = "Prop4ii"); the
# dispatchers in Mobius_vMF.R and Mobius_SvMF.R pick the _euc variant when xe is supplied.
# Shared internals live in Mobius_prop4_helpers.R.
#
# Future unification candidate: the spherical-only and _euc joint-fit engines below are
# structurally parallel but parameterise different manifolds. Merging them would change
# optimiser coordinates, so it is deliberately not done here.

#' @name mobius_link_prop4ii
#' @title Proposition 4(ii) Mean Link Parameters
#' @description Parameters of the scaled Mobius mean link constrained by
#' \eqn{p = q_s} and \eqn{B_s = I_{p-1}} --- Proposition 4(ii) of the manuscript.
#' Without Euclidean covariates the link is simply \eqn{\mu(x_s) = \tilde R_0 x_s}.
#' @param Rtilde0 A `p` x `p` orthonormal matrix.
#' @param P Final rotation on the response sphere: a `p` x `p` rotation matrix. Supplied
#'   only when the fit has Euclidean covariates, in which case \eqn{\tilde R_0 = P Q_s^\top}.
#' @param Qe_star An `m` x `(p-1)` matrix with orthonormal columns acting on the
#'   standardised Euclidean covariates, where `m` counts them plus any intercept.
#' @param Be The `p-1` diagonal entries of the Euclidean scaling matrix, each in `(0, 1)`.
#' @param xe_center,xe_scale The centre and scale used to standardise the Euclidean
#'   covariates during fitting, so that predictions reproduce them.
#' @param intercept `TRUE` if an intercept column was appended to the Euclidean covariates.
#' @param check Whether to validate the result.
#' @param obj An object to coerce or validate.
#' @return `mobius_link_prop4ii()` and `as_mobius_link_prop4ii()` return an object of class
#'   "mobius_link_prop4ii". `dim()` returns the named vector `c(p, qs, qe)`.
#' @details
#' `P`, `Qe_star` and `Be` must be supplied together or not at all. When they are supplied
#' the link is \eqn{\mu(x_s, x_e) = P\,\mathcal{S}^{[-1]}\{\mathcal{S}(Q_s^\top x_s) +
#' B_e Q_e^{*\top} \tilde x_e\}}, matching the `LinEuc` construction with \eqn{B_s} held
#' at the identity.
#'
#' Use [`as_mobius_link_cann()`] to convert to the canonical parameterisation for
#' interpretation; see [`mobius_link_params`].
#' @seealso [`mobius_link_prop4i`] for the isotropic-scale case, [`mobius_link()`],
#'   [`mobius_SvMF()`].
#' @family link-function
#' @export
mobius_link_prop4ii <- function(Rtilde0, P = NULL, Qe_star = NULL,
                                Be = NULL, xe_center = NULL,
                                xe_scale = NULL, intercept = TRUE,
                                check = TRUE) {
  Rtilde0 <- as.matrix(Rtilde0)
  storage.mode(Rtilde0) <- "double"
  has_euc <- !is.null(P) || !is.null(Qe_star) || !is.null(Be)
  # mobius_link_prop4ii_check() enforces this too, but as.matrix(NULL) below would fail
  # first with a much less helpful message. Matches mobius_link_prop4i().
  if (has_euc && (is.null(P) || is.null(Qe_star) || is.null(Be))) {
    stop("P, Qe_star, and Be must either all be supplied or all be NULL.")
  }

  obj <- list(
    Rtilde0 = Rtilde0,
    P = if (has_euc) as.matrix(P) else NULL,
    Qe_star = if (has_euc) as.matrix(Qe_star) else NULL,
    Be = if (has_euc) as.numeric(Be) else NULL,
    xe_center = if (has_euc) as.numeric(xe_center) else NULL,
    xe_scale = if (has_euc) as.numeric(xe_scale) else NULL,
    intercept = if (has_euc) isTRUE(intercept) else FALSE
  )
  class(obj) <- c("mobius_link_prop4ii", "list")
  if (check) mobius_link_prop4ii_check(obj)
  obj
}

mobius_link_prop4ii_check <- function(obj,
                                      tol = 100 * sqrt(.Machine$double.eps)) {
  if (!inherits(obj, "mobius_link_prop4ii")) {
    stop("obj must inherit from 'mobius_link_prop4ii'.")
  }
  R <- obj$Rtilde0
  if (!is.matrix(R) || nrow(R) != ncol(R) || nrow(R) < 2L) {
    stop("Rtilde0 must be a square p by p matrix with p >= 2.")
  }
  if (any(!is.finite(R))) stop("Rtilde0 contains non-finite values.")
  if (max(abs(crossprod(R) - diag(ncol(R)))) > tol) {
    stop("Rtilde0 is not orthogonal.")
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
  V <- obj$Qe_star
  b <- obj$Be
  if (!all(dim(P) == c(p, p))) stop("P must be p by p.")
  if (max(abs(crossprod(P) - diag(p))) > tol || det(P) < 0) {
    stop("P must be a proper orthogonal p by p matrix.")
  }
  if (!is.matrix(V) || ncol(V) != k || nrow(V) < k) {
    stop("Qe_star must be m by (p-1), with m >= p-1.")
  }
  if (max(abs(crossprod(V) - diag(k))) > tol) {
    stop("Qe_star must have orthonormal columns.")
  }
  if (length(b) != k || any(!is.finite(b)) || any(b <= 0) || any(b >= 1)) {
    stop("Be must contain p-1 values strictly between zero and one.")
  }
  q <- nrow(V) - as.integer(isTRUE(obj$intercept))
  if (q < 0 || length(obj$xe_center) != q || length(obj$xe_scale) != q) {
    stop("Euclidean centering/scaling information has the wrong dimension.")
  }
  invisible(NULL)
}

#' @describeIn mobius_link_prop4ii Coerce a compatible object to the Proposition 4(ii)
#'   parameterisation. Accepts a `p` x `p` matrix, a list with an `Rtilde0` element, or a
#'   canonical object with `Bs = I` and `p = qs`.
#' @export
as_mobius_link_prop4ii <- function(obj, check = TRUE) {
  if (inherits(obj, "mobius_link_prop4ii")) {
    if (check) mobius_link_prop4ii_check(obj)
    return(obj)
  }
  if (is.matrix(obj)) return(mobius_link_prop4ii(obj, check = check))
  if (is.list(obj) && !is.null(obj$Rtilde0)) {
    args <- obj[intersect(
      c("Rtilde0", "P", "Qe_star", "Be", "xe_center", "xe_scale", "intercept"),
      names(obj)
    )]
    args$check <- check
    return(do.call(mobius_link_prop4ii, args))
  }
  if (inherits(obj, "mobius_link_cann")) {
    p <- nrow(obj$P)
    if (is.null(obj$Qs) || !all(dim(obj$Qs) == c(p, p))) {
      stop("The canonical object does not have p = qs.")
    }
    if (max(abs(obj$Bs - diag(p - 1L))) >
        100 * sqrt(.Machine$double.eps)) {
      stop("The spherical part requires Bs = I_(p-1).")
    }
    R <- obj$P %*% t(obj$Qs)
    if (is.null(obj$Qe)) return(mobius_link_prop4ii(R, check = check))

    if (!is_LinEuc(obj)) {
      stop("Only the LinEuc Euclidean form is supported with type='Prop4ii'.")
    }
    V <- obj$Qe[-1, -1, drop = FALSE]
    return(mobius_link_prop4ii(
      Rtilde0 = R,
      P = obj$P,
      Qe_star = V,
      Be = diag(obj$Be),
      xe_center = rep(0, nrow(V) - 1L),
      xe_scale = rep(1, nrow(V) - 1L),
      intercept = TRUE,
      check = check
    ))
  }
  stop("obj is not a compatible Proposition 4(ii) parameter object.")
}

#' @export
print.mobius_link_prop4ii <- function(x, ...) {
  cat("Proposition 4(ii) spherical link")
  if (!is.null(x$P)) cat(" with optional LinEuc covariates")
  cat("\n")
  cat("  p = qs =", nrow(x$Rtilde0), "\n")
  cat("  qe =", if (is.null(x$P)) 0L else length(x$xe_center), "\n")
  cat("  det(Rtilde0) =", format(det(x$Rtilde0), digits = 6), "\n")
  if (!is.null(x$Be)) {
    cat("  Be =", paste(format(x$Be, digits = 5), collapse = ", "), "\n")
  }
  invisible(x)
}

#' @export
dim.mobius_link_prop4ii <- function(x) {
  p <- nrow(x$Rtilde0)
  c(p = p, qs = p, qe = if (is.null(x$P)) 0L else length(x$xe_center))
}

#' @noRd
#' @description Evaluate the Proposition 4(ii) mean link. Called by [`mobius_link()`].
mobius_link_pred_prop4ii <- function(xs = NULL, xe = NULL, param = NULL, check = TRUE) {
  if (check) mobius_link_prop4ii_check(param)
  if (is.null(xs)) stop("The Proposition 4(ii) link requires xs.")
  xs <- as.matrix(xs)
  if (ncol(xs) != nrow(param$Rtilde0)) {
    stop("ncol(xs) must equal nrow(Rtilde0).")
  }

  if (is.null(param$P)) {
    if (.is_xe_nonempty(xe)) {
      stop("This fitted Proposition 4(ii) link has no Euclidean component.")
    }
    return(xs %*% t(param$Rtilde0))
  }

  if (is.null(xe)) stop("xe must be supplied for this fitted model.")
  ep <- .prop4_prepare_euc(
    xe, intercept = param$intercept,
    center = param$xe_center, scale = param$xe_scale,
    fitting = FALSE
  )
  if (nrow(xs) != nrow(ep$model)) stop("xs and xe must have the same number of rows.")

  # Rtilde0 = P Qs^T, hence Qs = Rtilde0^T P.
  Qs <- t(param$Rtilde0) %*% param$P
  spherical_part <- Sp(xs %*% Qs)
  euclidean_part <- sweep(ep$model %*% param$Qe_star, 2, param$Be, "*")
  iSp(spherical_part + euclidean_part) %*% t(param$P)
}

mobius_vMF_prop4ii <- function(y, xs, xe = NULL,
                               det_constraint = c("orthogonal", "rotation"),
                               start = NULL, check = TRUE) {
  det_constraint <- match.arg(det_constraint)
  y <- as.matrix(y)
  xs <- as.matrix(xs)
  storage.mode(y) <- storage.mode(xs) <- "double"

  if (!is.null(xe) && ncol(as.matrix(xe)) > 0L) {
    stop("Proposition 4(ii) requires qe = 0; xe must be NULL.")
  }
  if (nrow(y) != nrow(xs) || ncol(y) != ncol(xs)) {
    stop("For Proposition 4(ii), y and xs must have the same n by p dimensions.")
  }
  if (ncol(y) < 2L) stop("p must be at least 2.")
  if (check) {
    .prop4_check_unit_rows(y, "y")
    .prop4_check_unit_rows(xs, "xs")
  }
  if (!is.null(start)) {
    warning("start is ignored because the vMF estimate of Rtilde0 is available in closed form.")
  }

  det_sign <- if (det_constraint == "rotation") 1 else NULL
  pro <- .prop4_procrustes(y, xs, det_sign = det_sign)
  Rhat <- pro$R
  est <- mobius_link_prop4ii(Rhat)
  pred <- xs %*% t(Rhat)
  dots <- pmax(-1, pmin(1, rowSums(y * pred)))
  rbar <- mean(dots)
  kfit <- .prop4_estimate_k(rbar, ncol(y))

  p <- ncol(y)
  n <- nrow(y)
  DoF <- p * (p - 1L) / 2L + 1L
  lLik <- n * kfit$loglik_per_obs
  dists <- acos(dots)
  rresids <- NULL
  rr <- try(rotated_resid(y, pred, north_pole(p)), silent = TRUE)
  if (!inherits(rr, "try-error")) {
    rresids <- rr[, -1, drop = FALSE]
    attr(rresids, "samehemisphere") <- attr(rr, "samehemisphere")
    colnames(rresids) <- paste0("r", seq_len(ncol(rresids)))
  }

  out <- list(
    est = est,
    mean = est,
    Rtilde0 = Rhat,
    k = kfit$k,
    obj = -rbar,
    solution = est,
    nlopt = list(
      status = 1L,
      message = "Closed-form orthogonal Procrustes solution.",
      objective = -rbar,
      solution = as.vector(Rhat)
    ),
    y = y,
    xs = xs,
    xe = NULL,
    pred = pred,
    rresids = rresids,
    dists = dists,
    DoF = DoF,
    AIC = 2 * DoF - 2 * lLik,
    lLik = lLik,
    start = start,
    singular_values = pro$singular_values,
    linktype = list(type = "Prop4ii", det_constraint = det_constraint,
                    intercept = FALSE)
  )
  class(out) <- c("mobius_vMF_prop4ii", "mobius_vMF", "list")
  out
}

mobius_vMF_prop4ii_euc <- function(
    y, xs, xe = NULL,
    det_constraint = c("orthogonal", "rotation"),
    start = NULL, check = TRUE, intercept = TRUE,
    algorithm = "NLOPT_LD_SLSQP", opts = list(), ...) {

  det_constraint <- match.arg(det_constraint)
  y <- as.matrix(y)
  xs <- as.matrix(xs)
  xe <- as.matrix(xe)
  storage.mode(y) <- storage.mode(xs) <- storage.mode(xe) <- "double"
  if (nrow(y) != nrow(xs) || nrow(y) != nrow(xe) || ncol(y) != ncol(xs)) {
    stop("For type='Prop4ii', y and xs must be n by p and xe must have n rows.")
  }
  p <- ncol(y)
  k <- p - 1L
  if (check) {
    .prop4_check_unit_rows(y, "y")
    .prop4_check_unit_rows(xs, "xs")
  }
  ep <- .prop4_prepare_euc(xe, intercept = intercept, fitting = TRUE)
  if (ep$m < k) {
    stop("With Euclidean covariates, ncol(xe) + intercept must be at least p-1.")
  }

  start_obj <- if (is.null(start)) NULL else as_mobius_link_prop4ii(start)
  target_signs <- if (det_constraint == "rotation") 1 else c(1, -1)
  fits <- vector("list", length(target_signs))

  for (j in seq_along(target_signs)) {
    sgn <- target_signs[j]
    Rref <- if (!is.null(start_obj) &&
                (if (det(start_obj$Rtilde0) >= 0) 1 else -1) == sgn) {
      start_obj$Rtilde0
    } else {
      .prop4_procrustes(y, xs, det_sign = sgn)$R
    }
    Pref <- if (!is.null(start_obj) && !is.null(start_obj$P)) {
      start_obj$P
    } else {
      .prop4_nearest_orthogonal(Rref, det_sign = 1)
    }
    Vref <- if (!is.null(start_obj) && !is.null(start_obj$Qe_star) &&
                all(dim(start_obj$Qe_star) == c(ep$m, k))) {
      start_obj$Qe_star
    } else {
      diag(ep$m)[, seq_len(k), drop = FALSE]
    }
    bstart <- if (!is.null(start_obj) && !is.null(start_obj$Be)) {
      start_obj$Be
    } else {
      rep(0.05, k)
    }

    spec <- .prop4ii_mean_spec_euc(
      Rref = Rref, Pref = Pref, Vref = Vref, bstart = bstart,
      xe_center = ep$center, xe_scale = ep$scale,
      intercept = intercept
    )
    fits[[j]] <- .prop4ii_fit_vmf_euc_component(
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
  dV <- DoF_Stiefel(ep$m, k)
  dB <- k
  DoF <- dR + dP + dV + dB + 1L
  lLik <- n * kfit$loglik_per_obs

  rr <- try(rotated_resid(y, pred, north_pole(p)), silent = TRUE)
  rresids <- if (inherits(rr, "try-error")) NULL else rr[, -1, drop = FALSE]
  if (!is.null(rresids)) {
    attr(rresids, "samehemisphere") <- attr(rr, "samehemisphere")
    colnames(rresids) <- paste0("r", seq_len(ncol(rresids)))
  }

  out <- list(
    est = pars$mean,
    mean = pars$mean,
    Rtilde0 = pars$Rtilde0,
    P = pars$P,
    Qe_star = pars$Qe_star,
    Be = pars$Be,
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
      list(det_sign = z$det_sign, status = z$nlopt$status,
           message = z$nlopt$message, objective = z$objective)
    }),
    linktype = list(
      type = "Prop4ii", det_constraint = det_constraint,
      intercept = intercept, euclidean = TRUE
    )
  )
  class(out) <- c("mobius_vMF_prop4ii", "mobius_vMF", "list")
  out
}

.prop4ii_mean_spec_euc <- function(Rref, Pref, Vref, bstart,
                                   xe_center, xe_scale, intercept) {
  p <- nrow(Rref)
  k <- p - 1L
  m <- nrow(Vref)
  dR <- .prop4_skew_dim(p)
  dP <- .prop4_skew_dim(p)
  dV <- DoF_Stiefel(m, k)
  dB <- k

  pos <- 0L
  take <- function(n) {
    if (n == 0L) return(integer(0))
    z <- pos + seq_len(n)
    pos <<- pos + n
    z
  }
  idxR <- take(dR)
  idxP <- take(dP)
  idxV <- take(dV)
  idxB <- take(dB)

  bstart <- pmin(pmax(as.numeric(bstart), 1e-5), 1 - 1e-5)
  par0 <- numeric(pos)
  par0[idxB] <- qlogis(bstart)
  lower <- rep(-Inf, pos)
  upper <- rep( Inf, pos)
  lower[idxB] <- -12
  upper[idxB] <- 12

  list(
    p = p, k = k, m = m,
    Rref = Rref,
    Pref = .prop4_nearest_orthogonal(Pref, det_sign = 1),
    Vref = qr.Q(qr(Vref), complete = FALSE)[, seq_len(k), drop = FALSE],
    xe_center = xe_center,
    xe_scale = xe_scale,
    intercept = intercept,
    idxR = idxR, idxP = idxP, idxV = idxV, idxB = idxB,
    par0 = par0, lower = lower, upper = upper
  )
}

.prop4ii_unpack_mean_euc <- function(theta, spec) {
  R <- .prop4_cayley(theta[spec$idxR], spec$p) %*% spec$Rref
  P <- .prop4_cayley(theta[spec$idxP], spec$p) %*% spec$Pref
  V <- .prop4_stiefel_cayley(theta[spec$idxV], spec$Vref)
  b <- plogis(theta[spec$idxB])
  mean <- mobius_link_prop4ii(
    Rtilde0 = R, P = P, Qe_star = V, Be = b,
    xe_center = spec$xe_center, xe_scale = spec$xe_scale,
    intercept = spec$intercept, check = FALSE
  )
  list(mean = mean, Rtilde0 = R, P = P, Qe_star = V, Be = b)
}

.prop4ii_fit_vmf_euc_component <- function(y, xs, xe, spec,
                                           algorithm, opts) {
  eval_f <- function(theta) {
    pars <- .prop4ii_unpack_mean_euc(theta, spec)
    pred <- try(mobius_link(xs = xs, xe = xe, param = pars$mean,
                            check = FALSE), silent = TRUE)
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

  fit <- nloptr::nloptr(
    x0 = spec$par0,
    eval_f = eval_nloptr,
    lb = spec$lower,
    ub = spec$upper,
    opts = combined_opts
  )
  pars <- .prop4ii_unpack_mean_euc(fit$solution, spec)
  pred <- mobius_link(xs = xs, xe = xe, param = pars$mean, check = FALSE)
  list(nlopt = fit, pars = pars, pred = pred,
       objective = -mean(rowSums(y * pred)))
}

mobius_SvMF_partransport_prelim_prop4ii <- function(
    y, xs, xe = NULL, mean = NULL, G0 = NULL,
    G01behaviour = "p1",
    det_constraint = c("orthogonal", "rotation"), check = TRUE, ...) {

  det_constraint <- match.arg(det_constraint)
  G01behaviour <- match.arg(G01behaviour, c("p1", "free", "fixed"))
  y <- as.matrix(y)
  xs <- as.matrix(xs)

  if (!is.null(xe) && ncol(as.matrix(xe)) > 0L) {
    stop("Proposition 4(ii) requires qe = 0; xe must be NULL.")
  }
  if (nrow(y) != nrow(xs) || ncol(y) != ncol(xs)) {
    stop("For Proposition 4(ii), y and xs must have identical dimensions.")
  }
  if (check) {
    .prop4_check_unit_rows(y, "y")
    .prop4_check_unit_rows(xs, "xs")
  }

  prelim <- mobius_vMF_prop4ii(
    y = y, xs = xs, xe = NULL, start = mean,
    det_constraint = det_constraint, check = FALSE
  )
  mean <- prelim$est
  k <- prelim$k
  Rtilde0 <- prelim$Rtilde0

  # This is exactly the original switch, with prelim$est$p1 replaced by the
  # identifiable representative p1 = Rtilde0[,1].
  if (G01behaviour == "fixed" && is.null(G0)) {
    stop("At least the first column of G0 must be supplied when G01behaviour = 'fixed'.")
  }
  G01 <- switch(
    G01behaviour,
    p1 = Rtilde0[, 1],
    free = if (is.null(G0)) Rtilde0[, 1] else G0[, 1],
    fixed = G0[, 1]
  )

  rresid <- rotated_resid(y, prelim$pred, base = G01)

  if (!is.null(G0) && all(!is.na(G0))) {
    if (G01behaviour == "p1") {
      # Original re-expression of supplied axes at the new p1 base.
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
    Rtilde0 = Rtilde0,
    k = k,
    a = c(1, aremaining),
    G0 = G0,
    nlopt = prelim$nlopt,
    pred = prelim$pred
  )
}

mobius_SvMF_partransport_prelim_prop4ii_euc <- function(
    y, xs, xe = NULL, mean = NULL, G0 = NULL,
    G01behaviour = "p1",
    det_constraint = c("orthogonal", "rotation"), check = TRUE,
    intercept = TRUE, algorithm = "NLOPT_LD_SLSQP",
    opts = list(), ...) {

  det_constraint <- match.arg(det_constraint)
  G01behaviour <- match.arg(G01behaviour, c("p1", "free", "fixed"))
  prelim <- mobius_vMF_prop4ii_euc(
    y = y, xs = xs, xe = xe, start = mean,
    det_constraint = det_constraint, check = check,
    intercept = intercept, algorithm = algorithm, opts = opts
  )

  # With Euclidean covariates, p1 is identified as the first column of P.
  p1 <- prelim$mean$P[, 1]
  if (G01behaviour == "fixed" && is.null(G0)) {
    stop("At least the first column of G0 must be supplied when G01behaviour='fixed'.")
  }
  G01 <- switch(
    G01behaviour,
    p1 = p1,
    free = if (is.null(G0)) p1 else G0[, 1],
    fixed = G0[, 1]
  )

  rresid <- rotated_resid(as.matrix(y), prelim$pred, base = G01)
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
    mean = prelim$mean,
    Rtilde0 = prelim$Rtilde0,
    P = prelim$P,
    Qe_star = prelim$Qe_star,
    Be = prelim$Be,
    k = prelim$k,
    a = c(1, aremaining),
    G0 = G0,
    nlopt = prelim$nlopt,
    pred = prelim$pred
  )
}

.prop4ii_joint_loglik <- function(y, xs, Rtilde0, k, a, G0,
                                  approximate = TRUE) {
  if (!is.finite(k) || k <= 0 || any(!is.finite(a)) || any(a <= 0)) {
    return(-Inf)
  }
  pred <- xs %*% t(Rtilde0)
  yrot <- try(undo_partransport(y, pred, G01 = G0[, 1]), silent = TRUE)
  if (inherits(yrot, "try-error") || any(!is.finite(yrot))) return(-Inf)

  # Try the same fast approximate density used by the original optimiser first.
  # For some otherwise valid parameter values ldSvMF_cann() can return a
  # non-finite value.  In that case fall back to the exact R implementation
  # rather than treating the whole determinant component as having -Inf
  # likelihood.  This does not change the model; it is only a numerical
  # safeguard for Proposition 4(ii).
  ld <- NULL
  if (approximate) {
    ld_try <- try(ldSvMF_cann(yrot, k = k, a = a, G = G0), silent = TRUE)
    if (!inherits(ld_try, "try-error") && all(is.finite(ld_try))) {
      ld <- ld_try
    }
  }
  if (is.null(ld)) {
    ld_try <- try(
      SvMF_log_lik_cann(yrot, SvMF_cann(k = k, a = a, G = G0)),
      silent = TRUE
    )
    if (!inherits(ld_try, "try-error") && all(is.finite(ld_try))) {
      ld <- ld_try
    }
  }
  if (is.null(ld)) return(-Inf)
  sum(ld)
}

.prop4ii_build_joint_spec <- function(Rref, prelim, G0reference,
                                      G01behaviour) {
  p <- nrow(Rref)
  G01behaviour <- match.arg(G01behaviour, c("p1", "free", "fixed"))

  dR <- .prop4_skew_dim(p)
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
    Q0 <- Gstart %*% t(Gref)
    thetaG0 <- .prop4_inverse_cayley(Q0)
    Ganchor <- NULL
    fixed_base <- NULL
  } else {
    target_base <- if (G01behaviour == "p1") Rref[, 1] else Gstart[, 1]

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
  par0[idxK] <- log(max(prelim$k, 1e-8))
  if (dA > 0L) par0[idxA] <- .prop4_scales_to_eta(prelim$a)
  if (dG > 0L) par0[idxG] <- thetaG0

  # As in the original Cayley parameterisation, rotation coordinates are
  # unconstrained. Bounds are used only for k and the log-scale coordinates.
  lower <- rep(-Inf, pos)
  upper <- rep( Inf, pos)
  lower[idxK] <- log(1e-8)
  upper[idxK] <- log(1e5)
  if (dA > 0L) {
    lower[idxA] <- -10
    upper[idxA] <- 10
  }

  list(
    p = p,
    Rref = Rref,
    a1 = prelim$a[1],
    G01behaviour = G01behaviour,
    Gref = Gref,
    Ganchor = Ganchor,
    fixed_base = fixed_base,
    idxR = idxR,
    idxK = idxK,
    idxA = idxA,
    idxG = idxG,
    par0 = par0,
    lower = lower,
    upper = upper
  )
}

.prop4ii_unpack_joint <- function(theta, spec) {
  p <- spec$p
  Rtilde0 <- .prop4_cayley(theta[spec$idxR], p) %*% spec$Rref
  k <- exp(theta[spec$idxK])
  eta <- if (length(spec$idxA)) theta[spec$idxA] else numeric(0)
  a <- .prop4_eta_to_scales(eta, spec$a1, p)

  if (spec$G01behaviour == "free") {
    G0 <- .prop4_cayley(theta[spec$idxG], p) %*% spec$Gref
  } else {
    Gbase <- spec$Ganchor
    if (spec$G01behaviour == "p1") {
      # Same identification as the original code: G0[,1] is p1. Here the
      # identifiable representative is p1 = Rtilde0[,1].
      Gbase <- .prop4_rebase_G0(spec$Ganchor, Rtilde0[, 1])
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

  list(Rtilde0 = Rtilde0, k = k, a = a, G0 = G0)
}

.prop4ii_fit_component <- function(y, xs, Rref, prelim, G0reference,
                                   G01behaviour, algorithm, opts,
                                   approximate = TRUE) {
  spec <- .prop4ii_build_joint_spec(
    Rref = Rref,
    prelim = prelim,
    G0reference = G0reference,
    G01behaviour = G01behaviour
  )

  eval_f <- function(theta) {
    pars <- .prop4ii_unpack_joint(theta, spec)
    if (is.null(pars)) return(1e100)
    ll <- .prop4ii_joint_loglik(
      y, xs, pars$Rtilde0, pars$k, pars$a, pars$G0,
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
      value <- eval_f(theta)
      gradient <- .prop4_fd_gradient(
        eval_f, theta, spec$lower, spec$upper
      )
      list(objective = value, gradient = gradient)
    }
  } else {
    eval_f
  }

  fit <- nloptr::nloptr(
    x0 = spec$par0,
    eval_f = eval_nloptr,
    lb = spec$lower,
    ub = spec$upper,
    opts = combined_opts
  )

  pars <- .prop4ii_unpack_joint(fit$solution, spec)
  ll <- if (is.null(pars)) -Inf else
    .prop4ii_joint_loglik(
      y, xs, pars$Rtilde0, pars$k, pars$a, pars$G0,
      approximate = approximate
    )

  list(nlopt = fit, pars = pars, lLik_approx = ll, spec = spec)
}

mobius_SvMF_joint_fit_prop4ii <- function(
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
    stop("Proposition 4(ii) requires qe = 0; xe must be NULL.")
  }
  if (nrow(y) != nrow(xs) || ncol(y) != ncol(xs)) {
    stop("For Proposition 4(ii), y and xs must have identical dimensions.")
  }
  p <- ncol(y)
  if (check) {
    .prop4_check_unit_rows(y, "y")
    .prop4_check_unit_rows(xs, "xs")
  }

  start_mean <- as_mobius_link_prop4ii(mean)
  start_R <- start_mean$Rtilde0
  start_a <- .prop4_normalise_scales(a, p, a1 = a[1])
  start_G0 <- .prop4_nearest_orthogonal(G0, det_sign = 1)

  initial <- list(mean = start_mean, Rtilde0 = start_R,
                  k = k, a = start_a, G0 = start_G0)

  target_signs <- if (det_constraint == "rotation") 1 else c(1, -1)
  fits <- vector("list", length(target_signs))

  for (j in seq_along(target_signs)) {
    sgn <- target_signs[j]
    if ((if (det(start_R) >= 0) 1 else -1) == sgn) {
      Rref <- start_R
    } else {
      Rref <- .prop4_procrustes(y, xs, det_sign = sgn)$R
    }

    pre_j <- initial
    pre_j$Rtilde0 <- Rref
    pre_j$mean <- mobius_link_prop4ii(Rref)

    # Under G01behaviour="p1", the original G0 base must follow p1.
    if (G01behaviour == "p1") {
      rebased <- .prop4_rebase_G0(pre_j$G0, Rref[, 1])
      if (is.null(rebased)) stop("Could not rebase G0 to Rtilde0[,1].")
      pre_j$G0 <- rebased
    }

    fits[[j]] <- .prop4ii_fit_component(
      y = y,
      xs = xs,
      Rref = Rref,
      prelim = pre_j,
      G0reference = if (is.null(G0reference)) pre_j$G0 else G0reference,
      G01behaviour = G01behaviour,
      algorithm = algorithm,
      opts = opts,
      approximate = TRUE
    )
    fits[[j]]$det_sign <- sgn
  }

  ll_approx <- vapply(fits, function(z) z$lLik_approx, numeric(1))
  if (all(!is.finite(ll_approx))) stop("All Proposition 4(ii) optimisations failed.")
  best <- fits[[which.max(ll_approx)]]
  pars <- best$pars

  # Original identifiability convention for (a_j, G0[,j]), j >= 2.
  aord <- order(pars$a[-1], decreasing = TRUE)
  pars$a[-1] <- pars$a[-1][aord]
  pars$G0[, -1] <- pars$G0[, -1, drop = FALSE][, aord, drop = FALSE]
  pars$G0[, -1] <- standardise_col_signs(pars$G0[, -1, drop = FALSE])
  if (det(pars$G0) < 0) pars$G0[, p] <- -pars$G0[, p]

  est <- mobius_link_prop4ii(pars$Rtilde0)
  pred <- xs %*% t(pars$Rtilde0)

  # The original code refines k with the exact likelihood when p != 3.
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

  # Same final validity/unimodality check as in the original code.
  SvMF_cann_check(SvMF_cann(k = pars$k, a = pars$a, G = pars$G0))

  dots <- pmax(-1, pmin(1, rowSums(y * pred)))
  dists <- acos(dots)

  rrI <- rotated_resid(y, pred, north_pole(p))
  rresids_I <- rrI[, -1, drop = FALSE]
  attr(rresids_I, "samehemisphere") <- attr(rrI, "samehemisphere")
  colnames(rresids_I) <- paste0("r", seq_len(ncol(rresids_I)))

  rresids_G0 <- resid_SvMF_partransport(
    y, pred, G0 = pars$G0, scale = FALSE
  )
  rresids_std <- resid_SvMF_partransport(
    y, pred, pars$k, pars$a, pars$G0, scale = TRUE
  )

  dmean <- .prop4_skew_dim(p)
  dscale <- max(p - 2L, 0L)
  dG0 <- if (G01behaviour == "free") {
    .prop4_skew_dim(p)
  } else {
    .prop4_skew_dim(p - 1L)
  }
  DoF <- dmean + 1L + dscale + dG0

  out <- list(
    mean = est,
    est = est,
    Rtilde0 = pars$Rtilde0,
    k = pars$k,
    a = pars$a,
    G0 = pars$G0,
    obj = -lLik / nrow(y),
    nlopt = best$nlopt,
    y = y,
    xs = xs,
    xe = NULL,
    pred = pred,
    rresids_I = rresids_I,
    rresids_G0 = rresids_G0,
    rresids_std = rresids_std,
    dists = dists,
    DoF = DoF,
    AIC = 2 * DoF - 2 * lLik,
    lLik = lLik,
    initial = initial,
    component_fits = lapply(fits, function(z) {
      list(det_sign = z$det_sign,
           status = z$nlopt$status,
           message = z$nlopt$message,
           lLik_approx = z$lLik_approx)
    }),
    linktype = list(
      type = "Prop4ii",
      det_constraint = det_constraint,
      G01behaviour = G01behaviour,
      intercept = FALSE
    )
  )
  class(out) <- c("mobius_SvMF_prop4ii", "mobius_SvMF", "list")
  out
}

.prop4ii_joint_spec_euc <- function(mean_ref, prelim, G0reference,
                                    G01behaviour) {
  mean_ref <- as_mobius_link_prop4ii(mean_ref)
  p <- nrow(mean_ref$Rtilde0)
  k <- p - 1L
  mspec <- .prop4ii_mean_spec_euc(
    Rref = mean_ref$Rtilde0,
    Pref = mean_ref$P,
    Vref = mean_ref$Qe_star,
    bstart = mean_ref$Be,
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
  Gcoord <- if (is.null(G0reference)) Gstart else
    .prop4_nearest_orthogonal(G0reference, det_sign = 1)

  if (G01behaviour == "free") {
    Gref <- Gcoord
    thetaG0 <- .prop4_inverse_cayley(Gstart %*% t(Gref))
    Ganchor <- NULL
  } else {
    target <- if (G01behaviour == "p1") mean_ref$P[, 1] else Gstart[, 1]
    Ganchor <- .prop4_rebase_G0(Gcoord, target)
    if (is.null(Ganchor)) Ganchor <- .prop4_rebase_G0(Gstart, target)
    if (is.null(Ganchor)) stop("Could not construct the G0 reference frame.")
    Gstart2 <- .prop4_rebase_G0(Gstart, target)
    H0 <- crossprod(Ganchor[, -1, drop = FALSE],
                    Gstart2[, -1, drop = FALSE])
    H0 <- .prop4_nearest_orthogonal(H0, det_sign = 1)
    thetaG0 <- .prop4_inverse_cayley(H0)
    Gref <- NULL
  }

  par0 <- c(mspec$par0, rep(0, pos - length(mspec$par0)))
  par0[idxK] <- log(max(prelim$k, 1e-8))
  if (dA > 0L) par0[idxA] <- .prop4_scales_to_eta(prelim$a)
  if (dG > 0L) par0[idxG] <- thetaG0
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
    a1 = prelim$a[1],
    G01behaviour = G01behaviour,
    Gref = Gref,
    Ganchor = Ganchor,
    idxK = idxK, idxA = idxA, idxG = idxG,
    par0 = par0, lower = lower, upper = upper
  )
}

.prop4ii_unpack_joint_euc <- function(theta, spec) {
  mp <- .prop4ii_unpack_mean_euc(
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
      Gbase <- .prop4_rebase_G0(spec$Ganchor, mp$P[, 1])
      if (is.null(Gbase)) return(NULL)
    }
    H <- if (spec$p == 2L) matrix(1, 1, 1) else
      .prop4_cayley(theta[spec$idxG], spec$p - 1L)
    B <- diag(spec$p)
    B[-1, -1] <- H
    G0 <- Gbase %*% B
  }

  c(mp, list(k = kappa, a = a, G0 = G0))
}

.prop4ii_joint_loglik_euc <- function(y, xs, xe, mean, k, a, G0,
                                      approximate = TRUE) {
  if (!is.finite(k) || k <= 0 || any(!is.finite(a)) || any(a <= 0)) {
    return(-Inf)
  }
  pred <- try(mobius_link(xs = xs, xe = xe, param = mean,
                          check = FALSE), silent = TRUE)
  if (inherits(pred, "try-error") || any(!is.finite(pred))) return(-Inf)
  yrot <- try(undo_partransport(y, pred, G01 = G0[, 1]), silent = TRUE)
  if (inherits(yrot, "try-error") || any(!is.finite(yrot))) return(-Inf)
  # Same numerical safeguard as in the spherical-only Proposition 4(ii)
  # likelihood: use the fast approximation when finite and otherwise fall back
  # to the exact R density.
  ld <- NULL
  if (approximate) {
    ld_try <- try(ldSvMF_cann(yrot, k = k, a = a, G = G0), silent = TRUE)
    if (!inherits(ld_try, "try-error") && all(is.finite(ld_try))) ld <- ld_try
  }
  if (is.null(ld)) {
    ld_try <- try(
      SvMF_log_lik_cann(yrot, SvMF_cann(k = k, a = a, G = G0)),
      silent = TRUE
    )
    if (!inherits(ld_try, "try-error") && all(is.finite(ld_try))) ld <- ld_try
  }
  if (is.null(ld)) return(-Inf)
  sum(ld)
}

.prop4ii_fit_joint_euc_component <- function(
    y, xs, xe, mean_ref, prelim, G0reference, G01behaviour,
    algorithm, opts) {

  spec <- .prop4ii_joint_spec_euc(
    mean_ref = mean_ref, prelim = prelim,
    G0reference = G0reference,
    G01behaviour = G01behaviour
  )
  eval_f <- function(theta) {
    pars <- .prop4ii_unpack_joint_euc(theta, spec)
    if (is.null(pars)) return(1e100)
    ll <- .prop4ii_joint_loglik_euc(
      y, xs, xe, pars$mean, pars$k, pars$a, pars$G0,
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

  fit <- nloptr::nloptr(
    x0 = spec$par0,
    eval_f = eval_nloptr,
    lb = spec$lower,
    ub = spec$upper,
    opts = combined_opts
  )
  pars <- .prop4ii_unpack_joint_euc(fit$solution, spec)
  ll <- if (is.null(pars)) -Inf else .prop4ii_joint_loglik_euc(
    y, xs, xe, pars$mean, pars$k, pars$a, pars$G0,
    approximate = TRUE
  )
  list(nlopt = fit, pars = pars, lLik_approx = ll, spec = spec)
}

mobius_SvMF_joint_fit_prop4ii_euc <- function(
    y, xs, xe = NULL, mean, k, a, G0,
    G0reference = NULL, G01behaviour = "p1",
    det_constraint = c("orthogonal", "rotation"),
    algorithm = "NLOPT_LD_SLSQP", opts = list(), check = TRUE,
    intercept = TRUE) {

  det_constraint <- match.arg(det_constraint)
  G01behaviour <- match.arg(G01behaviour, c("p1", "free", "fixed"))
  y <- as.matrix(y)
  xs <- as.matrix(xs)
  xe <- as.matrix(xe)
  p <- ncol(y)
  if (nrow(y) != nrow(xs) || nrow(y) != nrow(xe) || ncol(xs) != p) {
    stop("For type='Prop4ii', y and xs must be n by p and xe must have n rows.")
  }
  if (check) {
    .prop4_check_unit_rows(y, "y")
    .prop4_check_unit_rows(xs, "xs")
  }

  start_mean <- as_mobius_link_prop4ii(mean)
  if (is.null(start_mean$P)) stop("The starting mean has no Euclidean component.")
  start_a <- .prop4_normalise_scales(a, p, a1 = a[1])
  start_G0 <- .prop4_nearest_orthogonal(G0, det_sign = 1)
  initial <- list(mean = start_mean, Rtilde0 = start_mean$Rtilde0,
                  P = start_mean$P, Qe_star = start_mean$Qe_star,
                  Be = start_mean$Be, k = k, a = start_a, G0 = start_G0)

  target_signs <- if (det_constraint == "rotation") 1 else c(1, -1)
  fits <- vector("list", length(target_signs))
  for (j in seq_along(target_signs)) {
    sgn <- target_signs[j]
    Rref <- if ((if (det(start_mean$Rtilde0) >= 0) 1 else -1) == sgn) {
      start_mean$Rtilde0
    } else {
      .prop4_procrustes(y, xs, det_sign = sgn)$R
    }
    mean_ref <- mobius_link_prop4ii(
      Rtilde0 = Rref,
      P = start_mean$P,
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
      rebased <- .prop4_rebase_G0(pre_j$G0, mean_ref$P[, 1])
      if (is.null(rebased)) stop("Could not rebase G0 to P[,1].")
      pre_j$G0 <- rebased
    }
    fits[[j]] <- .prop4ii_fit_joint_euc_component(
      y = y, xs = xs, xe = xe,
      mean_ref = mean_ref, prelim = pre_j,
      G0reference = if (is.null(G0reference)) pre_j$G0 else G0reference,
      G01behaviour = G01behaviour,
      algorithm = algorithm, opts = opts
    )
    fits[[j]]$det_sign <- sgn
  }

  lla <- vapply(fits, function(z) z$lLik_approx, numeric(1))
  if (all(!is.finite(lla))) stop("All Proposition 4(ii)+Euclidean optimisations failed.")
  best <- fits[[which.max(lla)]]
  pars <- best$pars

  # Original residual-axis identifiability convention.
  aord <- order(pars$a[-1], decreasing = TRUE)
  pars$a[-1] <- pars$a[-1][aord]
  pars$G0[, -1] <- pars$G0[, -1, drop = FALSE][, aord, drop = FALSE]
  pars$G0[, -1] <- standardise_col_signs(pars$G0[, -1, drop = FALSE])
  if (det(pars$G0) < 0) pars$G0[, p] <- -pars$G0[, p]

  pred <- mobius_link(xs = xs, xe = xe, param = pars$mean, check = FALSE)
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

  ep_m <- nrow(pars$mean$Qe_star)
  kmean <- p - 1L
  dmean <- .prop4_skew_dim(p) +       # Rtilde0
    .prop4_skew_dim(p) +              # P
    DoF_Stiefel(ep_m, kmean) + # Qe_star
    kmean                               # Be
  dscale <- max(p - 2L, 0L)
  dG0 <- if (G01behaviour == "free") {
    .prop4_skew_dim(p)
  } else {
    .prop4_skew_dim(p - 1L)
  }
  DoF <- dmean + 1L + dscale + dG0

  out <- list(
    mean = pars$mean,
    est = pars$mean,
    Rtilde0 = pars$Rtilde0,
    P = pars$P,
    Qe_star = pars$Qe_star,
    Be = pars$Be,
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
      list(det_sign = z$det_sign, status = z$nlopt$status,
           message = z$nlopt$message, lLik_approx = z$lLik_approx)
    }),
    linktype = list(
      type = "Prop4ii", det_constraint = det_constraint,
      G01behaviour = G01behaviour,
      intercept = intercept, euclidean = TRUE
    )
  )
  class(out) <- c("mobius_SvMF_prop4ii", "mobius_SvMF", "list")
  out
}

mobius_SvMF_prop4ii <- function(
    y, xs, xe = NULL, mean = NULL, k = NULL, a = NULL,
    G0 = NULL, G0reference = NULL, G01behaviour = "p1",
    doprelim = TRUE,
    det_constraint = c("orthogonal", "rotation"),
    algorithm = "NLOPT_LD_SLSQP", opts = list(), check = TRUE, ...) {

  det_constraint <- match.arg(det_constraint)
  G01behaviour <- match.arg(G01behaviour, c("p1", "free", "fixed"))
  dots <- list(...)
  if (length(dots)) opts <- utils::modifyList(opts, dots)

  if (doprelim) {
    preest <- mobius_SvMF_partransport_prelim_prop4ii(
      y = y, xs = xs, xe = xe, mean = mean, G0 = G0,
      G01behaviour = G01behaviour,
      det_constraint = det_constraint, check = check
    )
  } else {
    if (is.null(mean) || is.null(k) || is.null(a) || is.null(G0)) {
      stop("When doprelim = FALSE, mean, k, a, and G0 must all be supplied.")
    }
    mean <- as_mobius_link_prop4ii(mean)
    preest <- list(
      mean = mean,
      Rtilde0 = mean$Rtilde0,
      k = k,
      a = a,
      G0 = G0
    )
  }

  finalest <- mobius_SvMF_joint_fit_prop4ii(
    y = y,
    xs = xs,
    xe = xe,
    mean = preest$mean,
    k = if (is.null(k)) preest$k else k,
    a = if (is.null(a)) preest$a else a,
    G0 = preest$G0,
    G0reference = if (is.null(G0reference)) preest$G0 else G0reference,
    G01behaviour = G01behaviour,
    det_constraint = det_constraint,
    algorithm = algorithm,
    opts = opts,
    check = check
  )
  finalest$preest <- preest
  finalest
}

mobius_SvMF_prop4ii_euc <- function(
    y, xs, xe = NULL, mean = NULL, k = NULL, a = NULL,
    G0 = NULL, G0reference = NULL, G01behaviour = "p1",
    doprelim = TRUE,
    det_constraint = c("orthogonal", "rotation"),
    algorithm = "NLOPT_LD_SLSQP", opts = list(), check = TRUE,
    intercept = TRUE, ...) {

  det_constraint <- match.arg(det_constraint)
  G01behaviour <- match.arg(G01behaviour, c("p1", "free", "fixed"))
  dots <- list(...)
  if (length(dots)) opts <- utils::modifyList(opts, dots)

  if (doprelim) {
    preest <- mobius_SvMF_partransport_prelim_prop4ii_euc(
      y = y, xs = xs, xe = xe, mean = mean, G0 = G0,
      G01behaviour = G01behaviour,
      det_constraint = det_constraint, check = check,
      intercept = intercept, algorithm = algorithm, opts = opts
    )
  } else {
    if (is.null(mean) || is.null(k) || is.null(a) || is.null(G0)) {
      stop("When doprelim=FALSE, mean, k, a, and G0 must all be supplied.")
    }
    mean <- as_mobius_link_prop4ii(mean)
    if (is.null(mean$P)) stop("mean must contain the Euclidean component.")
    preest <- list(
      mean = mean, Rtilde0 = mean$Rtilde0,
      P = mean$P, Qe_star = mean$Qe_star, Be = mean$Be,
      k = k, a = a, G0 = G0
    )
  }

  finalest <- mobius_SvMF_joint_fit_prop4ii_euc(
    y = y, xs = xs, xe = xe,
    mean = preest$mean,
    k = if (is.null(k)) preest$k else k,
    a = if (is.null(a)) preest$a else a,
    G0 = preest$G0,
    G0reference = if (is.null(G0reference)) preest$G0 else G0reference,
    G01behaviour = G01behaviour,
    det_constraint = det_constraint,
    algorithm = algorithm,
    opts = opts,
    check = check,
    intercept = intercept
  )
  finalest$preest <- preest
  finalest
}

#' @export
print.mobius_vMF_prop4ii <- function(x, ...) {
  cat("vMF regression with Proposition 4(ii) spherical link")
  if (isTRUE(x$linktype$euclidean)) cat(" and LinEuc covariates")
  cat("\n")
  cat("  det(Rtilde0) =", format(det(x$Rtilde0), digits = 6), "\n")
  if (!is.null(x$Be)) cat("  Be =", paste(format(x$Be, digits = 5), collapse = ", "), "\n")
  cat("  k =", format(x$k, digits = 6), "\n")
  cat("  log-likelihood =", format(x$lLik, digits = 7), "\n")
  cat("  AIC =", format(x$AIC, digits = 7), "\n")
  invisible(x)
}

#' @export
print.mobius_SvMF_prop4ii <- function(x, ...) {
  cat("SvMF regression with Proposition 4(ii) spherical link")
  if (isTRUE(x$linktype$euclidean)) cat(" and LinEuc covariates")
  cat("\n")
  cat("  G0 estimation: original parallel-transport/Cayley procedure\n")
  cat("  G01behaviour =", x$linktype$G01behaviour, "\n")
  cat("  det(Rtilde0) =", format(det(x$Rtilde0), digits = 6), "\n")
  if (!is.null(x$Be)) cat("  Be =", paste(format(x$Be, digits = 5), collapse = ", "), "\n")
  cat("  k =", format(x$k, digits = 6), "\n")
  cat("  a =", paste(format(x$a, digits = 5), collapse = ", "), "\n")
  cat("  log-likelihood =", format(x$lLik, digits = 7), "\n")
  cat("  AIC =", format(x$AIC, digits = 7), "\n")
  invisible(x)
}

#' @export
predict.mobius_vMF_prop4ii <- function(object, newdata = NULL,
                                        xs = NULL, xe = NULL, ...) {
  if (is.list(newdata) && is.null(xs)) {
    xs <- newdata$xs
    if (is.null(xe)) xe <- newdata$xe
  } else if (!is.null(newdata) && is.null(xs)) {
    xs <- newdata
  }
  if (is.null(xs)) xs <- object$xs
  if (is.null(xe) && !is.null(object$mean$P)) xe <- object$xe
  mobius_link(xs = as.matrix(xs), xe = xe, param = object$mean)
}

#' @export
predict.mobius_SvMF_prop4ii <- function(object, newdata = NULL,
                                         xs = NULL, xe = NULL, ...) {
  predict.mobius_vMF_prop4ii(object, newdata = newdata, xs = xs, xe = xe, ...)
}
