# Shared internals for the Proposition 4(i) and 4(ii) constrained links.
#
# Proposition 4 of "Regression for spherical responses with linear and spherical
# covariates using a scaled link function" gives two special cases of the scaled
# Moebius link:
#   4(ii)  Bs = I_(p-1)          -- mu(xs) = Rtilde0 xs
#   4(i)   Bs = beta_s I_(p-1)   -- the isotropic Moebius map, 0 < beta_s <= 1
# Each is fitted with and without Euclidean covariates; the four resulting submodels
# are the ones reported in Section 5.1 of the manuscript. Implementations live in
# Mobius_prop4i.R and Mobius_prop4ii.R; this file holds only what they share.
#
# These fitters optimise with central finite-difference gradients (.prop4_fd_gradient)
# rather than the C++ AD tapes used by the SpEuc/LinEuc paths, so they are slower.
#
# Naming: internals carry a leading dot and a .prop4_ / .prop4i_ / .prop4ii_ prefix.
#
# Future unification candidates (deliberately NOT merged, because doing so would change
# optimiser coordinates and shift fitted values): the four parallel
# _mean_spec / _unpack_* / _joint_loglik / _fit_component / _joint_fit families in
# Mobius_prop4i.R and Mobius_prop4ii.R, and .prop4_inverse_cayley versus the C++
# inverse_cayley_transform() (which uses the lower-triangle convention).

#' @noRd
#' @title Recognise the Proposition 4 link types
#' @description Accepts "Prop4i"/"Proposition4i" and "Prop4ii"/"Proposition4ii" in any
#'   capitalisation or punctuation, matching how `type` is written in the vignettes.

.is_prop4ii_type <- function(type) {
  if (length(type) != 1L || is.na(type)) return(FALSE)
  key <- tolower(gsub("[^[:alnum:]]", "", as.character(type)))
  key %in% c("prop4ii", "proposition4ii")
}

#' @noRd
.is_prop4i_type <- function(type) {
  if (length(type) != 1L || is.na(type)) return(FALSE)
  key <- tolower(gsub("[^[:alnum:]]", "", as.character(type)))
  key %in% c("prop4i", "proposition4i")
}

#' @noRd
#' @description TRUE when `xe` supplies at least one Euclidean covariate.
.prop4_has_euc <- function(xe) {
  !is.null(xe) && ncol(as.matrix(xe)) > 0L
}

#' @noRd
.prop4_check_unit_rows <- function(x, name,
                                     tol = 1000 * sqrt(.Machine$double.eps)) {
  if (!is.matrix(x) || !is.numeric(x)) {
    stop(name, " must be a numeric matrix.")
  }
  if (any(!is.finite(x))) stop(name, " contains non-finite values.")
  nr <- sqrt(rowSums(x^2))
  err <- max(abs(nr - 1))
  if (err > tol) {
    stop(name, " must contain unit vectors as rows; maximum norm error is ",
         format(err, digits = 4), ".")
  }
  invisible(NULL)
}

#' @noRd
#' @description Nearest orthogonal matrix to `M` in Frobenius norm (polar/SVD
#' projection), optionally forced to a given determinant sign. Unlike `as_rotation_mat()`
#' this projects a non-orthogonal matrix rather than flipping a column sign.
.prop4_nearest_orthogonal <- function(M, det_sign = NULL) {
  M <- as.matrix(M)
  if (nrow(M) != ncol(M)) stop("M must be square.")
  p <- nrow(M)
  z <- svd(M, nu = p, nv = p)
  D <- diag(p)
  R0 <- z$u %*% t(z$v)
  if (!is.null(det_sign)) {
    requested <- if (det_sign >= 0) 1 else -1
    current <- if (det(R0) >= 0) 1 else -1
    D[p, p] <- requested / current
  }
  z$u %*% D %*% t(z$v)
}

#' @noRd
.prop4_unit <- function(x, name = "vector") {
  x <- as.numeric(x)
  nrm <- sqrt(sum(x^2))
  if (!is.finite(nrm) || nrm <= sqrt(.Machine$double.eps)) {
    stop(name, " must be a finite nonzero vector.")
  }
  x / nrm
}

#' @noRd
#' @description Centre and scale the Euclidean covariates and append an intercept.
#' Deliberately separate from `standardise_Euc()`/`addEuccovars()`: those PCA-rotate
#' without scaling (the SpEuc link is not equivariant to scaling of xe), whereas the
#' Proposition 4 Euclidean term is linear and is fitted in centred, sd-scaled coordinates.
#' The centre and scale are stored on the fitted link so predictions can reproduce them.
.prop4_prepare_euc <- function(xe, intercept = TRUE,
                                 center = NULL, scale = NULL,
                                 fitting = FALSE) {
  xe <- as.matrix(xe)
  storage.mode(xe) <- "double"
  if (any(!is.finite(xe))) stop("xe contains non-finite values.")
  q <- ncol(xe)

  if (fitting) {
    if (intercept) {
      center <- if (q) colMeans(xe) else numeric(0)
      scale <- if (q) apply(xe, 2, stats::sd) else numeric(0)
      scale[!is.finite(scale) | scale <= sqrt(.Machine$double.eps)] <- 1
    } else {
      center <- rep(0, q)
      scale <- rep(1, q)
    }
  } else {
    if (is.null(center) || length(center) != q) {
      stop("The number of columns of xe does not match the fitted model.")
    }
    if (is.null(scale) || length(scale) != q) {
      stop("The Euclidean scaling information is invalid.")
    }
  }

  xstd <- if (q) {
    sweep(sweep(xe, 2, center, "-"), 2, scale, "/")
  } else {
    matrix(numeric(0), nrow(xe), 0L)
  }
  xmodel <- if (intercept) cbind(xstd, `(Intercept)` = 1) else xstd

  list(
    raw = xe,
    standardised = xstd,
    model = xmodel,
    center = as.numeric(center),
    scale = as.numeric(scale),
    q = q,
    m = ncol(xmodel),
    intercept = isTRUE(intercept)
  )
}

#' @noRd
.prop4_skew_dim <- function(p) as.integer(p * (p - 1L) / 2L)

#' @noRd
#' @title Cayley transform, Proposition 4 sign convention
#' @description Maps the upper-triangle entries of a skew-symmetric `A` to the orthogonal
#' matrix \eqn{(I - A)^{-1}(I + A)}. The package's `cayley()` computes
#' \eqn{(I - A)(I + A)^{-1}}; for skew-symmetric `A` the two are transposes of each other,
#' so this is `t(cayley(theta))` (equivalently `cayley(-theta)`). The Proposition 4 code
#' uses this orientation throughout, paired with `.prop4_inverse_cayley()`.
#' @param theta numeric vector of length `p*(p-1)/2`, the upper triangle of `A` filled
#'   column by column -- the same ordering `cayley()` uses.
#' @param p dimension of the resulting orthogonal matrix.
.prop4_cayley <- function(theta, p) {
  stopifnot(length(theta) == p * (p - 1L) / 2L)
  t(cayley(theta))
}

#' @noRd
.prop4_inverse_cayley <- function(Q, tol = 1e-10) {
  Q <- as.matrix(Q)
  p <- nrow(Q)
  I <- diag(p)
  if (ncol(Q) != p || rcond(Q + I) < tol) {
    return(rep(0, .prop4_skew_dim(p)))
  }
  A <- (Q - I) %*% solve(Q + I)
  A <- (A - t(A)) / 2
  A[upper.tri(A)]
}

#' @noRd
.prop4_complete_orthogonal <- function(V) {
  V <- as.matrix(V)
  m <- nrow(V)
  k <- ncol(V)
  if (k > m) stop("A Stiefel matrix cannot have more columns than rows.")
  V <- qr.Q(qr(V), complete = FALSE)[, seq_len(k), drop = FALSE]
  if (k == m) return(V)
  z <- svd(t(V), nu = 0, nv = m)$v
  W <- z[, (k + 1L):m, drop = FALSE]
  Q <- cbind(V, W)
  if (det(Q) < 0) Q[, m] <- -Q[, m]
  Q
}

#' @noRd
.prop4_stiefel_cayley <- function(theta, Vref) {
  Vref <- as.matrix(Vref)
  m <- nrow(Vref)
  k <- ncol(Vref)
  dA <- k * (k - 1L) / 2L
  dB <- (m - k) * k
  if (length(theta) != dA + dB) stop("Incorrect Stiefel-coordinate length.")

  A <- matrix(0, k, k)
  if (dA > 0L) {
    ij <- which(upper.tri(A), arr.ind = TRUE)
    A[ij] <- theta[seq_len(dA)]
    A[cbind(ij[, 2], ij[, 1])] <- -theta[seq_len(dA)]
  }
  B <- if (dB > 0L) {
    matrix(theta[dA + seq_len(dB)], nrow = m - k, ncol = k)
  } else {
    matrix(numeric(0), 0L, k)
  }

  K <- matrix(0, m, m)
  K[seq_len(k), seq_len(k)] <- A
  if (m > k) {
    K[k + seq_len(m - k), seq_len(k)] <- B
    K[seq_len(k), k + seq_len(m - k)] <- -t(B)
  }

  Qref <- .prop4_complete_orthogonal(Vref)
  Qref %*% .prop4_cayley(K[upper.tri(K)], m)[, seq_len(k), drop = FALSE]
}

#' @noRd
.prop4_frame_with_base <- function(base, tangent, tol = 1e-10) {
  base <- .prop4_unit(base, "base")
  p <- length(base)
  tangent <- as.matrix(tangent)
  if (!all(dim(tangent) == c(p, p - 1L))) {
    stop("tangent must be p by (p-1).")
  }
  base_col <- matrix(base, ncol = 1L)
  tangent <- tangent - base_col %*% crossprod(base, tangent)
  q <- qr(tangent, tol = tol)
  if (q$rank < p - 1L) {
    seed <- diag(p)[, -which.max(abs(base)), drop = FALSE]
    seed <- seed - base_col %*% crossprod(base, seed)
    tangent <- qr.Q(qr(seed), complete = FALSE)[, seq_len(p - 1L), drop = FALSE]
  } else {
    tangent <- qr.Q(q, complete = FALSE)[, seq_len(p - 1L), drop = FALSE]
  }
  G <- cbind(base, tangent)
  if (det(G) < 0) G[, p] <- -G[, p]
  G
}

#' @noRd
#' @description Re-express G0 about a new base direction by parallel transport.
.prop4_rebase_G0 <- function(G0, newbase, tol = 1e-8) {
  G0 <- .prop4_nearest_orthogonal(G0, det_sign = 1)
  oldbase <- G0[, 1]
  newbase <- .prop4_unit(newbase, "newbase")
  if (sum(oldbase * newbase) <= -1 + tol) return(NULL)

  # jupp_Rmat() is the preferred transport; parallel_transport_mat() is the fallback
  # when it is singular (antipodal base directions).
  tangent <- try(-jupp_Rmat(oldbase, newbase) %*%
                   G0[, -1, drop = FALSE], silent = TRUE)
  if (inherits(tangent, "try-error")) {
    tangent <- try(parallel_transport_mat(oldbase, newbase) %*%
                     G0[, -1, drop = FALSE], silent = TRUE)
  }
  if (inherits(tangent, "try-error") || any(!is.finite(tangent))) return(NULL)
  .prop4_frame_with_base(newbase, tangent)
}

#' @noRd
.prop4_normalise_scales <- function(a, p, a1 = NULL) {
  a <- as.numeric(a)
  if (length(a) != p || any(!is.finite(a)) || any(a <= 0)) {
    stop("a must be a positive vector of length p.")
  }
  if (!is.null(a1)) a[1] <- a1
  if (p > 1L) a[-1] <- a[-1] / exp(mean(log(a[-1])))
  a
}

#' @noRd
.prop4_scales_to_eta <- function(a) {
  p <- length(a)
  if (p <= 2L) return(numeric(0))
  log(a[2:(p - 1L)])
}

#' @noRd
.prop4_eta_to_scales <- function(eta, a1, p) {
  if (p == 2L) return(c(a1, 1))
  if (length(eta) != p - 2L) stop("Incorrect number of scale coordinates.")
  c(a1, exp(eta), exp(-sum(eta)))
}

#' @noRd
.prop4_procrustes <- function(y, xs, det_sign = NULL) {
  p <- ncol(y)
  C <- crossprod(y, xs)             # Y^T X
  z <- svd(C, nu = p, nv = p)
  D <- diag(p)
  R0 <- z$u %*% t(z$v)
  if (!is.null(det_sign)) {
    requested <- if (det_sign >= 0) 1 else -1
    current <- if (det(R0) >= 0) 1 else -1
    D[p, p] <- requested / current
  }
  list(R = z$u %*% D %*% t(z$v), singular_values = z$d)
}

#' @noRd
.prop4_estimate_k <- function(rbar, p) {
  rbar <- min(max(as.numeric(rbar), -1 + 1e-12), 1 - 1e-12)
  objective <- function(k) {
    -vMF_log_norm_const_exact(k, p) + k * rbar
  }
  fit <- optimise(objective, lower = 1e-8, upper = 1e5, maximum = TRUE)
  list(k = fit$maximum, loglik_per_obs = fit$objective)
}

#' @noRd
.prop4_fd_gradient <- function(fn, x, lower, upper,
                                 rel_step = .Machine$double.eps^(1/3)) {
  g <- numeric(length(x))
  f0 <- fn(x)
  for (j in seq_along(x)) {
    h <- rel_step * (1 + abs(x[j]))
    xp <- xm <- x
    xp[j] <- min(x[j] + h, upper[j])
    xm[j] <- max(x[j] - h, lower[j])
    if (xp[j] > x[j] && xm[j] < x[j]) {
      g[j] <- (fn(xp) - fn(xm)) / (xp[j] - xm[j])
    } else if (xp[j] > x[j]) {
      g[j] <- (fn(xp) - f0) / (xp[j] - x[j])
    } else if (xm[j] < x[j]) {
      g[j] <- (f0 - fn(xm)) / (x[j] - xm[j])
    } else {
      g[j] <- 0
    }
  }
  g
}

#' @noRd
#' @title Canonical form of a Proposition 4 mean link
#' @description Re-express a `mobius_link_prop4i` or `mobius_link_prop4ii` parameter object
#' as a [`mobius_link_cann()`] object, so that fitted Proposition 4 models can be read and
#' plotted with the same code as unconstrained fits. Reached through
#' [`as_mobius_link_cann()`].
#' @details
#' The canonical link is
#' \deqn{\mu(x) = P\,\mathcal{S}^{[-1]}\!\left(B_s \mathcal{S}(Q_s^\top x_s) +
#'   \frac{B_e Q_e[,-1]^\top x_e}{Q_e[,1]^\top x_e + c_e}\right).}
#'
#' * Proposition 4(ii) fixes \eqn{B_s = I_{p-1}}; Proposition 4(i) fixes
#'   \eqn{B_s = \beta_s I_{p-1}}.
#' * Without Euclidean covariates, `Qs` is any rotation whose first column is the
#'   reference direction \eqn{r_{s1}} and `P = Rtilde0 Qs`. Scaling the stereographic
#'   coordinates by \eqn{\beta_s} about \eqn{r_{s1}} is exactly the isotropic Moebius map
#'   that `.prop4i_map()` applies -- this is Proposition 4(i) read backwards. At
#'   \eqn{\beta_s = 1} the direction \eqn{r_{s1}} is not identified and `Qs = I` is used.
#' * With Euclidean covariates, `Rtilde0 = P Qs^T` gives `Qs`, and the leading zero row of
#'   `Qe` is the LinEuc dummy-zero covariate that `addEuccovars()` also uses, so the
#'   denominator is 1 and the Euclidean term reduces to the fitted `xe_model Qe_star Be`.
#'
#' NOTE: with Euclidean covariates the result is expressed in the fit's INTERNAL
#' standardised covariates -- `.prop4_prepare_euc()` centres by `xe_center`, scales by
#' `xe_scale` and appends an intercept -- because that standardisation cannot be absorbed
#' into an orthonormal `Qe`. Those three are attached as attributes. Evaluating
#' [`mobius_link()`] on the result therefore needs `cbind(dummyzero = 0, <standardised xe>)`,
#' not the raw `xe`. The conversion is intended for interpretation and display.
.prop4_as_cann <- function(obj) {
  p <- nrow(obj$Rtilde0)
  beta_s <- if (is.null(obj$beta_s)) 1 else obj$beta_s   # Prop 4(ii) has Bs = I
  Bs <- diag(beta_s, p - 1L)

  if (is.null(obj$P)) {
    unidentified <- beta_s >= 1 - 1e-12 || is.null(obj$rs1) || any(!is.finite(obj$rs1))
    Qs <- if (unidentified) {
      diag(p)
    } else {
      .prop4_frame_with_base(obj$rs1, diag(p)[, -1L, drop = FALSE])
    }
    return(mobius_link_cann(P = obj$Rtilde0 %*% Qs, Bs = Bs, Qs = Qs, check = FALSE))
  }

  # Rtilde0 = P Qs^T, hence Qs = Rtilde0^T P.
  Qs <- t(obj$Rtilde0) %*% obj$P
  m <- nrow(obj$Qe_star)
  Qe <- matrix(0, m + 1L, p)
  Qe[1L, 1L] <- 1                      # dummy-zero covariate: keeps the denominator at ce
  Qe[-1L, -1L] <- obj$Qe_star
  if (!is.null(rownames(obj$Qe_star))) {
    rownames(Qe) <- c("dummyzero", rownames(obj$Qe_star))
  }
  out <- mobius_link_cann(P = obj$P, Bs = Bs, Qs = Qs,
                          Be = diag(obj$Be, p - 1L), Qe = Qe, ce = 1,
                          check = FALSE)
  attr(out, "xe_center") <- obj$xe_center
  attr(out, "xe_scale") <- obj$xe_scale
  attr(out, "intercept") <- obj$intercept
  out
}

#' @noRd
#' @title Reject Proposition 4 fits in the sign-flip / multistart helpers
#' @description These helpers explore sign flips of the pole directions p1, qs1 and qe1 in
#' the canonical parameterisation (see `signflip_starts()`), which the constrained links
#' of Proposition 4 do not have -- their mean is an orthogonal `Rtilde0` (plus, optionally,
#' `beta_s` and a Euclidean block). Fail with a clear message rather than deep inside.
#' @param x A fitted model, or a `type` string.
#' @param what Name of the calling function, used in the error message.
.prop4_reject <- function(x, what){
  is_prop4 <- if (is.character(x)){
    .is_prop4i_type(x) || .is_prop4ii_type(x)
  } else {
    inherits(x, c("mobius_link_prop4i", "mobius_link_prop4ii",
                  "mobius_vMF_prop4i", "mobius_vMF_prop4ii",
                  "mobius_SvMF_prop4i", "mobius_SvMF_prop4ii")) ||
      (is.list(x) && inherits(x$mean, c("mobius_link_prop4i", "mobius_link_prop4ii"))) ||
      (is.list(x) && inherits(x$est, c("mobius_link_prop4i", "mobius_link_prop4ii")))
  }
  if (isTRUE(is_prop4)){
    stop(what, "() does not support the constrained links of Proposition 4: they have no ",
         "qs1/qe1 pole directions to sign-flip. Refit with mobius_vMF()/mobius_SvMF() and ",
         "type = \"Prop4i\" or \"Prop4ii\" instead.")
  }
  invisible(NULL)
}
