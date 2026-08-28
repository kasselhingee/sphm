# =============================================================================
# REFERENCE COPY -- DO NOT TIDY, DO NOT CALL FROM PACKAGE CODE
# =============================================================================
#
# Shogo Kato's Proposition 4(i)/4(ii) implementation exactly as delivered, kept in
# the package under `obsolete_`-prefixed names so that the refactored implementation
# in Mobius_prop4_helpers.R / Mobius_prop4i.R / Mobius_prop4ii.R can be checked
# against it by tests/testthat/test-prop4_vs_obsolete.R on every test run.
#
# It is DELIBERATELY not tidied. In particular it still contains:
#   * obsolete_mobius_link / obsolete_mobius_vMF / obsolete_mobius_SvMF defined four
#     times each, chained through the .obsolete_mobius_*_before_* shims, so that
#     last-definition-wins reproduces what source() produced; and
#   * fourteen Proposition 4(ii) functions defined twice (spherical-only, then a
#     superset handling optional Euclidean covariates).
# Removing either would defeat the point of the file.
#
# Differences from the delivered file, all mechanical:
#   1. Every identifier it defines is prefixed `obsolete_` (`.foo` -> `.obsolete_foo`),
#      as is every S3 class it creates, so nothing collides with the live package.
#      Names it only *consumes* -- undo_partransport(), mobius_SvMF_konly(),
#      is_LinEuc(), the "mobius_link_cann" class, and the package internals it calls --
#      are left bare and resolve to the live package.
#   2. The first 1466 lines of the delivered file, a verbatim copy of an older snapshot
#      of the package's own R files, are dropped. The Proposition 4 code reached into
#      them for only undo_partransport(), mobius_SvMF_konly() and is_LinEuc(), all
#      three byte-identical to the live definitions.
#   3. .prop4ii_import_sphm_internal() and its lapply() call are dropped: they assigned
#      package internals into .GlobalEnv, which R CMD check rejects and which is
#      unnecessary inside a namespace.
#   4. The three `_definition1` shims are call-time wrappers rather than load-time
#      copies, so this file does not depend on R file collation order.
#   5. .obsolete_prop4ii_vmf_log_norm_const_exact() has its exists() guard removed --
#      see the comment at that function.
#
# Delete this file together with test-prop4_vs_obsolete.R once the Proposition 4
# results have been through review and the equivalence is no longer in question.
#
# @noRd
# =============================================================================

# =============================================================================
# Proposition 4(ii) extension for the spherical regression code
#
# Mean link:
#   mu(xs) = Rtilde0 %*% xs,  Rtilde0 in O(p)
#
# The estimation of G0 is kept in the same form as in the original code:
#   * preliminary moment estimation from parallel-transported residuals;
#   * G01behaviour = "p1", "free", or "fixed";
#   * final joint likelihood estimation using Cayley coordinates relative to
#     G0reference;
#   * permutation/sign standardisation of (a_j, G0[,j]), j >= 2.
#
# Source this file AFTER the original package R files have been loaded.
# =============================================================================


# The three entry points as the live package defines them. Call-time wrappers rather
# than load-time copies so this file does not depend on R file collation order.
.obsolete_mobius_link_definition1 <- function(...) obsolete_mobius_link(...)
.obsolete_mobius_vMF_definition1 <- function(...) obsolete_mobius_vMF(...)
.obsolete_mobius_SvMF_definition1 <- function(...) obsolete_mobius_SvMF(...)

.obsolete_is_prop4ii_type <- function(type) {
  if (length(type) != 1L || is.na(type)) return(FALSE)
  key <- tolower(gsub("[^[:alnum:]]", "", as.character(type)))
  key %in% c("prop4ii", "proposition4ii")
}

.obsolete_prop4ii_check_unit_rows <- function(x, name,
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

.obsolete_prop4ii_nearest_orthogonal <- function(M, det_sign = NULL) {
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

# -----------------------------------------------------------------------------
# Proposition 4(ii) parameter class and mean calculation
# -----------------------------------------------------------------------------

obsolete_mobius_link_prop4ii <- function(Rtilde0, check = TRUE) {
  Rtilde0 <- as.matrix(Rtilde0)
  storage.mode(Rtilde0) <- "double"
  obj <- list(Rtilde0 = Rtilde0)
  class(obj) <- c("obsolete_mobius_link_prop4ii", "list")
  if (check) obsolete_mobius_link_prop4ii_check(obj)
  obj
}

obsolete_mobius_link_prop4ii_check <- function(obj,
                                      tol = 100 * sqrt(.Machine$double.eps)) {
  if (!inherits(obj, "obsolete_mobius_link_prop4ii")) {
    stop("obj must inherit from 'obsolete_mobius_link_prop4ii'.")
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
  invisible(NULL)
}

obsolete_as_mobius_link_prop4ii <- function(obj, check = TRUE) {
  if (inherits(obj, "obsolete_mobius_link_prop4ii")) {
    if (check) obsolete_mobius_link_prop4ii_check(obj)
    return(obj)
  }
  if (is.matrix(obj)) return(obsolete_mobius_link_prop4ii(obj, check = check))
  if (is.list(obj) && !is.null(obj$Rtilde0)) {
    return(obsolete_mobius_link_prop4ii(obj$Rtilde0, check = check))
  }
  if (inherits(obj, "mobius_link_cann")) {
    p <- nrow(obj$P)
    if (!is.null(obj$Qe) || is.null(obj$Qs) ||
        !all(dim(obj$Qs) == c(p, p))) {
      stop("The canonical object is not a Proposition 4(ii) link.")
    }
    if (max(abs(obj$Bs - diag(p - 1L))) >
        100 * sqrt(.Machine$double.eps)) {
      stop("Proposition 4(ii) requires Bs = I_(p-1).")
    }
    return(obsolete_mobius_link_prop4ii(obj$P %*% t(obj$Qs), check = check))
  }
  stop("obj must be a obsolete_mobius_link_prop4ii object, a square matrix, or a compatible canonical object.")
}

#' @noRd
#' @exportS3Method
print.obsolete_mobius_link_prop4ii <- function(x, ...) {
  cat("Proposition 4(ii) spherical link\n")
  cat("  p = qs =", nrow(x$Rtilde0), ", qe = 0\n")
  cat("  det(Rtilde0) =", format(det(x$Rtilde0), digits = 6), "\n")
  print(x$Rtilde0, ...)
  invisible(x)
}

#' @noRd
#' @exportS3Method
dim.obsolete_mobius_link_prop4ii <- function(x) {
  p <- nrow(x$Rtilde0)
  c(p = p, qs = p, qe = 0L)
}

# Wrapper: retain the original Definition 1 link for its original classes.
obsolete_mobius_link <- function(xs = NULL, xe = NULL, param = NULL, check = TRUE) {
  if (inherits(param, "obsolete_mobius_link_prop4ii")) {
    if (check) obsolete_mobius_link_prop4ii_check(param)
    if (is.null(xs)) stop("The Proposition 4(ii) link requires xs.")
    xs <- as.matrix(xs)
    if (!is.null(xe) && ncol(as.matrix(xe)) > 0L) {
      stop("Proposition 4(ii) requires qe = 0; xe must be NULL.")
    }
    if (ncol(xs) != nrow(param$Rtilde0)) {
      stop("ncol(xs) must equal nrow(Rtilde0).")
    }
    # Rows contain transposed column vectors: mu_i^T = xs_i^T Rtilde0^T.
    return(xs %*% t(param$Rtilde0))
  }
  .obsolete_mobius_link_definition1(xs = xs, xe = xe, param = param, check = check)
}

# -----------------------------------------------------------------------------
# Direct vMF estimation of Rtilde0 by orthogonal Procrustes
# -----------------------------------------------------------------------------

.obsolete_prop4ii_procrustes <- function(y, xs, det_sign = NULL) {
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

.obsolete_prop4ii_estimate_k <- function(rbar, p) {
  rbar <- min(max(as.numeric(rbar), -1 + 1e-12), 1 - 1e-12)
  objective <- function(k) {
    -vMF_log_norm_const_exact(k, p) + k * rbar
  }
  fit <- optimise(objective, lower = 1e-8, upper = 1e5, maximum = TRUE)
  list(k = fit$maximum, loglik_per_obs = fit$objective)
}

obsolete_mobius_vMF_prop4ii <- function(y, xs, xe = NULL,
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
    .obsolete_prop4ii_check_unit_rows(y, "y")
    .obsolete_prop4ii_check_unit_rows(xs, "xs")
  }
  if (!is.null(start)) {
    warning("start is ignored because the vMF estimate of Rtilde0 is available in closed form.")
  }

  det_sign <- if (det_constraint == "rotation") 1 else NULL
  pro <- .obsolete_prop4ii_procrustes(y, xs, det_sign = det_sign)
  Rhat <- pro$R
  est <- obsolete_mobius_link_prop4ii(Rhat)
  pred <- xs %*% t(Rhat)
  dots <- pmax(-1, pmin(1, rowSums(y * pred)))
  rbar <- mean(dots)
  kfit <- .obsolete_prop4ii_estimate_k(rbar, ncol(y))

  p <- ncol(y)
  n <- nrow(y)
  DoF <- p * (p - 1L) / 2L + 1L
  lLik <- n * kfit$loglik_per_obs
  dists <- acos(dots)
  rresids <- NULL
  if (exists("rotated_resid", mode = "function")) {
    rr <- try(rotated_resid(y, pred, north_pole(p)), silent = TRUE)
    if (!inherits(rr, "try-error")) {
      rresids <- rr[, -1, drop = FALSE]
      attr(rresids, "samehemisphere") <- attr(rr, "samehemisphere")
      colnames(rresids) <- paste0("r", seq_len(ncol(rresids)))
    }
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
  class(out) <- c("obsolete_mobius_vMF_prop4ii", "obsolete_mobius_vMF", "list")
  out
}

#' @noRd
#' @exportS3Method
print.obsolete_mobius_vMF_prop4ii <- function(x, ...) {
  cat("vMF regression with Proposition 4(ii) link\n")
  cat("  det(Rtilde0) =", format(det(x$Rtilde0), digits = 6), "\n")
  cat("  k =", format(x$k, digits = 6), "\n")
  cat("  log-likelihood =", format(x$lLik, digits = 7), "\n")
  cat("  AIC =", format(x$AIC, digits = 7), "\n")
  invisible(x)
}

#' @noRd
#' @exportS3Method
predict.obsolete_mobius_vMF_prop4ii <- function(object, newdata = NULL,
                                        xs = newdata, ...) {
  if (is.null(xs)) xs <- object$xs
  obsolete_mobius_link(xs = as.matrix(xs), param = object$mean)
}

# Wrapper retaining the original vMF implementation for the original link types.
obsolete_mobius_vMF <- function(y, xs = NULL, xe = NULL, start = NULL,
                       type = "SpEuc", fix_qs1 = FALSE,
                       fix_qe1 = (type == "LinEuc"), intercept = TRUE,
                       lb = NULL, ub = NULL,
                       det_constraint = c("orthogonal", "rotation"), ...) {
  if (.obsolete_is_prop4ii_type(type)) {
    return(obsolete_mobius_vMF_prop4ii(
      y = y, xs = xs, xe = xe, start = start,
      det_constraint = det_constraint
    ))
  }
  .obsolete_mobius_vMF_definition1(
    y = y, xs = xs, xe = xe, start = start, type = type,
    fix_qs1 = fix_qs1, fix_qe1 = fix_qe1, intercept = intercept,
    lb = lb, ub = ub, ...
  )
}

# -----------------------------------------------------------------------------
# Cayley-coordinate helpers
# -----------------------------------------------------------------------------

.obsolete_prop4ii_skew_dim <- function(p) as.integer(p * (p - 1L) / 2L)

.obsolete_prop4ii_vec_to_skew <- function(theta, p) {
  d <- .obsolete_prop4ii_skew_dim(p)
  if (length(theta) != d) stop("Incorrect skew-vector length.")
  A <- matrix(0, p, p)
  if (d > 0L) {
    ij <- which(upper.tri(A), arr.ind = TRUE)
    A[ij] <- theta
    A[cbind(ij[, 2], ij[, 1])] <- -theta
  }
  A
}

.obsolete_prop4ii_cayley <- function(theta, p) {
  A <- .obsolete_prop4ii_vec_to_skew(theta, p)
  I <- diag(p)
  solve(I - A, I + A)
}

.obsolete_prop4ii_inverse_cayley <- function(Q, tol = 1e-10) {
  Q <- as.matrix(Q)
  p <- nrow(Q)
  I <- diag(p)
  if (ncol(Q) != p || rcond(Q + I) < tol) {
    return(rep(0, .obsolete_prop4ii_skew_dim(p)))
  }
  A <- (Q - I) %*% solve(Q + I)
  A <- (A - t(A)) / 2
  A[upper.tri(A)]
}

.obsolete_prop4ii_unit <- function(x, name = "vector") {
  x <- as.numeric(x)
  nrm <- sqrt(sum(x^2))
  if (!is.finite(nrm) || nrm <= sqrt(.Machine$double.eps)) {
    stop(name, " must be a finite nonzero vector.")
  }
  x / nrm
}

.obsolete_prop4ii_frame_with_base <- function(base, tangent, tol = 1e-10) {
  base <- .obsolete_prop4ii_unit(base, "base")
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

# Same parallel-transport re-expression used in the original preliminary code.
.obsolete_prop4ii_rebase_G0 <- function(G0, newbase, tol = 1e-8) {
  G0 <- .obsolete_prop4ii_nearest_orthogonal(G0, det_sign = 1)
  oldbase <- G0[, 1]
  newbase <- .obsolete_prop4ii_unit(newbase, "newbase")
  if (sum(oldbase * newbase) <= -1 + tol) return(NULL)

  tangent <- NULL
  if (exists("jupp_Rmat", mode = "function")) {
    tangent <- try(-jupp_Rmat(oldbase, newbase) %*%
                     G0[, -1, drop = FALSE], silent = TRUE)
  }
  if (is.null(tangent) || inherits(tangent, "try-error")) {
    tangent <- try(parallel_transport_mat(oldbase, newbase) %*%
                     G0[, -1, drop = FALSE], silent = TRUE)
  }
  if (inherits(tangent, "try-error") || any(!is.finite(tangent))) return(NULL)
  .obsolete_prop4ii_frame_with_base(newbase, tangent)
}

.obsolete_prop4ii_normalise_scales <- function(a, p, a1 = NULL) {
  a <- as.numeric(a)
  if (length(a) != p || any(!is.finite(a)) || any(a <= 0)) {
    stop("a must be a positive vector of length p.")
  }
  if (!is.null(a1)) a[1] <- a1
  if (p > 1L) a[-1] <- a[-1] / exp(mean(log(a[-1])))
  a
}

.obsolete_prop4ii_scales_to_eta <- function(a) {
  p <- length(a)
  if (p <= 2L) return(numeric(0))
  log(a[2:(p - 1L)])
}

.obsolete_prop4ii_eta_to_scales <- function(eta, a1, p) {
  if (p == 2L) return(c(a1, 1))
  if (length(eta) != p - 2L) stop("Incorrect number of scale coordinates.")
  c(a1, exp(eta), exp(-sum(eta)))
}

# -----------------------------------------------------------------------------
# Preliminary SvMF step: G0 estimation is the original procedure
# -----------------------------------------------------------------------------

obsolete_mobius_SvMF_partransport_prelim_prop4ii <- function(
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
    .obsolete_prop4ii_check_unit_rows(y, "y")
    .obsolete_prop4ii_check_unit_rows(xs, "xs")
  }

  prelim <- obsolete_mobius_vMF_prop4ii(
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
  } else {
    # Original moment estimator of the axes and scales.
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

# -----------------------------------------------------------------------------
# Joint SvMF likelihood. G0 uses the original Cayley-coordinate structure.
# -----------------------------------------------------------------------------

.obsolete_prop4ii_joint_loglik <- function(y, xs, Rtilde0, k, a, G0,
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
  if (approximate && exists("ldSvMF_cann", mode = "function")) {
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

.obsolete_prop4ii_build_joint_spec <- function(Rref, prelim, G0reference,
                                      G01behaviour) {
  p <- nrow(Rref)
  G01behaviour <- match.arg(G01behaviour, c("p1", "free", "fixed"))

  dR <- .obsolete_prop4ii_skew_dim(p)
  dA <- max(p - 2L, 0L)
  dG <- if (G01behaviour == "free") {
    .obsolete_prop4ii_skew_dim(p)
  } else {
    .obsolete_prop4ii_skew_dim(p - 1L)
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

  Gstart <- .obsolete_prop4ii_nearest_orthogonal(prelim$G0, det_sign = 1)
  Gcoord <- if (is.null(G0reference)) {
    Gstart
  } else {
    .obsolete_prop4ii_nearest_orthogonal(G0reference, det_sign = 1)
  }

  if (G01behaviour == "free") {
    Gref <- Gcoord
    Q0 <- Gstart %*% t(Gref)
    thetaG0 <- .obsolete_prop4ii_inverse_cayley(Q0)
    Ganchor <- NULL
    fixed_base <- NULL
  } else {
    target_base <- if (G01behaviour == "p1") Rref[, 1] else Gstart[, 1]

    Ganchor <- .obsolete_prop4ii_rebase_G0(Gcoord, target_base)
    if (is.null(Ganchor)) Ganchor <- .obsolete_prop4ii_rebase_G0(Gstart, target_base)
    if (is.null(Ganchor)) stop("Could not construct the G0 reference frame.")

    Gstart2 <- .obsolete_prop4ii_rebase_G0(Gstart, target_base)
    if (is.null(Gstart2)) stop("Could not rebase the starting G0 frame.")

    H0 <- crossprod(Ganchor[, -1, drop = FALSE],
                    Gstart2[, -1, drop = FALSE])
    H0 <- .obsolete_prop4ii_nearest_orthogonal(H0, det_sign = 1)
    thetaG0 <- .obsolete_prop4ii_inverse_cayley(H0)
    Gref <- NULL
    fixed_base <- if (G01behaviour == "fixed") target_base else NULL
  }

  par0 <- numeric(pos)
  par0[idxR] <- 0
  par0[idxK] <- log(max(prelim$k, 1e-8))
  if (dA > 0L) par0[idxA] <- .obsolete_prop4ii_scales_to_eta(prelim$a)
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

.obsolete_prop4ii_unpack_joint <- function(theta, spec) {
  p <- spec$p
  Rtilde0 <- .obsolete_prop4ii_cayley(theta[spec$idxR], p) %*% spec$Rref
  k <- exp(theta[spec$idxK])
  eta <- if (length(spec$idxA)) theta[spec$idxA] else numeric(0)
  a <- .obsolete_prop4ii_eta_to_scales(eta, spec$a1, p)

  if (spec$G01behaviour == "free") {
    G0 <- .obsolete_prop4ii_cayley(theta[spec$idxG], p) %*% spec$Gref
  } else {
    Gbase <- spec$Ganchor
    if (spec$G01behaviour == "p1") {
      # Same identification as the original code: G0[,1] is p1. Here the
      # identifiable representative is p1 = Rtilde0[,1].
      Gbase <- .obsolete_prop4ii_rebase_G0(spec$Ganchor, Rtilde0[, 1])
      if (is.null(Gbase)) return(NULL)
    }
    H <- if (p == 2L) {
      matrix(1, 1, 1)
    } else {
      .obsolete_prop4ii_cayley(theta[spec$idxG], p - 1L)
    }
    B <- diag(p)
    B[-1, -1] <- H
    G0 <- Gbase %*% B
  }

  list(Rtilde0 = Rtilde0, k = k, a = a, G0 = G0)
}


.obsolete_prop4ii_fd_gradient <- function(fn, x, lower, upper,
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

.obsolete_prop4ii_fit_component <- function(y, xs, Rref, prelim, G0reference,
                                   G01behaviour, algorithm, opts,
                                   approximate = TRUE) {
  spec <- .obsolete_prop4ii_build_joint_spec(
    Rref = Rref,
    prelim = prelim,
    G0reference = G0reference,
    G01behaviour = G01behaviour
  )

  eval_f <- function(theta) {
    pars <- .obsolete_prop4ii_unpack_joint(theta, spec)
    if (is.null(pars)) return(1e100)
    ll <- .obsolete_prop4ii_joint_loglik(
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
      gradient <- .obsolete_prop4ii_fd_gradient(
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

  pars <- .obsolete_prop4ii_unpack_joint(fit$solution, spec)
  ll <- if (is.null(pars)) -Inf else
    .obsolete_prop4ii_joint_loglik(
      y, xs, pars$Rtilde0, pars$k, pars$a, pars$G0,
      approximate = approximate
    )

  list(nlopt = fit, pars = pars, lLik_approx = ll, spec = spec)
}

obsolete_mobius_SvMF_joint_fit_prop4ii <- function(
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
    .obsolete_prop4ii_check_unit_rows(y, "y")
    .obsolete_prop4ii_check_unit_rows(xs, "xs")
  }

  start_mean <- obsolete_as_mobius_link_prop4ii(mean)
  start_R <- start_mean$Rtilde0
  start_a <- .obsolete_prop4ii_normalise_scales(a, p, a1 = a[1])
  start_G0 <- .obsolete_prop4ii_nearest_orthogonal(G0, det_sign = 1)

  initial <- list(mean = start_mean, Rtilde0 = start_R,
                  k = k, a = start_a, G0 = start_G0)

  target_signs <- if (det_constraint == "rotation") 1 else c(1, -1)
  fits <- vector("list", length(target_signs))

  for (j in seq_along(target_signs)) {
    sgn <- target_signs[j]
    if ((if (det(start_R) >= 0) 1 else -1) == sgn) {
      Rref <- start_R
    } else {
      Rref <- .obsolete_prop4ii_procrustes(y, xs, det_sign = sgn)$R
    }

    pre_j <- initial
    pre_j$Rtilde0 <- Rref
    pre_j$mean <- obsolete_mobius_link_prop4ii(Rref)

    # Under G01behaviour="p1", the original G0 base must follow p1.
    if (G01behaviour == "p1") {
      rebased <- .obsolete_prop4ii_rebase_G0(pre_j$G0, Rref[, 1])
      if (is.null(rebased)) stop("Could not rebase G0 to Rtilde0[,1].")
      pre_j$G0 <- rebased
    }

    fits[[j]] <- .obsolete_prop4ii_fit_component(
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
  if (exists("standardise_col_signs", mode = "function")) {
    pars$G0[, -1] <- standardise_col_signs(pars$G0[, -1, drop = FALSE])
  }
  if (det(pars$G0) < 0) pars$G0[, p] <- -pars$G0[, p]

  est <- obsolete_mobius_link_prop4ii(pars$Rtilde0)
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

  dmean <- .obsolete_prop4ii_skew_dim(p)
  dscale <- max(p - 2L, 0L)
  dG0 <- if (G01behaviour == "free") {
    .obsolete_prop4ii_skew_dim(p)
  } else {
    .obsolete_prop4ii_skew_dim(p - 1L)
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
  class(out) <- c("obsolete_mobius_SvMF_prop4ii", "obsolete_mobius_SvMF", "list")
  out
}

obsolete_mobius_SvMF_prop4ii <- function(
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
    preest <- obsolete_mobius_SvMF_partransport_prelim_prop4ii(
      y = y, xs = xs, xe = xe, mean = mean, G0 = G0,
      G01behaviour = G01behaviour,
      det_constraint = det_constraint, check = check
    )
  } else {
    if (is.null(mean) || is.null(k) || is.null(a) || is.null(G0)) {
      stop("When doprelim = FALSE, mean, k, a, and G0 must all be supplied.")
    }
    mean <- obsolete_as_mobius_link_prop4ii(mean)
    preest <- list(
      mean = mean,
      Rtilde0 = mean$Rtilde0,
      k = k,
      a = a,
      G0 = G0
    )
  }

  finalest <- obsolete_mobius_SvMF_joint_fit_prop4ii(
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

#' @noRd
#' @exportS3Method
print.obsolete_mobius_SvMF_prop4ii <- function(x, ...) {
  cat("SvMF regression with Proposition 4(ii) link\n")
  cat("  G0 estimation: original parallel-transport/Cayley procedure\n")
  cat("  G01behaviour =", x$linktype$G01behaviour, "\n")
  cat("  det(Rtilde0) =", format(det(x$Rtilde0), digits = 6), "\n")
  cat("  k =", format(x$k, digits = 6), "\n")
  cat("  a =", paste(format(x$a, digits = 5), collapse = ", "), "\n")
  cat("  log-likelihood =", format(x$lLik, digits = 7), "\n")
  cat("  AIC =", format(x$AIC, digits = 7), "\n")
  invisible(x)
}

#' @noRd
#' @exportS3Method
predict.obsolete_mobius_SvMF_prop4ii <- function(object, newdata = NULL,
                                         xs = newdata, ...) {
  if (is.null(xs)) xs <- object$xs
  obsolete_mobius_link(xs = as.matrix(xs), param = object$mean)
}

# Wrapper retaining the original SvMF implementation for the original link types.
obsolete_mobius_SvMF <- function(
    y, xs = NULL, xe = NULL, mean = NULL, k = NULL, a = NULL,
    G0 = NULL, G0reference = NULL, G01behaviour = "p1",
    type = "LinEuc", fix_qs1 = FALSE,
    fix_qe1 = (type == "LinEuc"), intercept = TRUE,
    doprelim = TRUE,
    det_constraint = c("orthogonal", "rotation"),
    prop4ii_algorithm = "NLOPT_LD_SLSQP",
    prop4ii_opts = list(), ...) {

  if (.obsolete_is_prop4ii_type(type)) {
    return(obsolete_mobius_SvMF_prop4ii(
      y = y, xs = xs, xe = xe, mean = mean, k = k, a = a,
      G0 = G0, G0reference = G0reference,
      G01behaviour = G01behaviour, doprelim = doprelim,
      det_constraint = det_constraint,
      algorithm = prop4ii_algorithm,
      opts = prop4ii_opts,
      ...
    ))
  }

  .obsolete_mobius_SvMF_definition1(
    y = y, xs = xs, xe = xe, mean = mean, k = k, a = a,
    G0 = G0, G0reference = G0reference,
    G01behaviour = G01behaviour,
    type = type, fix_qs1 = fix_qs1, fix_qe1 = fix_qe1,
    intercept = intercept, doprelim = doprelim, ...
  )
}
# =============================================================================
# MINIMAL PATCH: Proposition 4(ii) spherical part + optional LinEuc covariates
#
# The spherical part is constrained by p = qs and Bs = I_(p-1).  With no
# Euclidean covariates, the link is exactly
#
#     mu(xs) = Rtilde0 xs.
#
# If Euclidean covariates are supplied, the original LinEuc construction is
# retained while Bs is fixed at the identity:
#
#     mu(xs, xe) = P iSp{ Sp(Qs^T xs) + Be Qe*^T xe_tilde },
#
# where Rtilde0 = P Qs^T, xe_tilde contains standardised Euclidean covariates
# and (when intercept = TRUE) an intercept, Qe* has orthonormal columns, and Be
# is diagonal.  Thus setting the Euclidean contribution to zero gives the
# Proposition 4(ii) link exactly.
#
# This section is intentionally appended at the end of the combined source file
# so that the original code is not rewritten.  Only the Proposition 4(ii)
# wrappers are overridden.
# =============================================================================

if (!exists(".obsolete_mobius_vMF_prop4ii_spherical_only", inherits = FALSE)) {
  .obsolete_mobius_vMF_prop4ii_spherical_only <- obsolete_mobius_vMF_prop4ii
}
if (!exists(".obsolete_mobius_SvMF_prelim_prop4ii_spherical_only", inherits = FALSE)) {
  .obsolete_mobius_SvMF_prelim_prop4ii_spherical_only <-
    obsolete_mobius_SvMF_partransport_prelim_prop4ii
}
if (!exists(".obsolete_mobius_SvMF_joint_prop4ii_spherical_only", inherits = FALSE)) {
  .obsolete_mobius_SvMF_joint_prop4ii_spherical_only <- obsolete_mobius_SvMF_joint_fit_prop4ii
}
if (!exists(".obsolete_mobius_SvMF_prop4ii_spherical_only", inherits = FALSE)) {
  .obsolete_mobius_SvMF_prop4ii_spherical_only <- obsolete_mobius_SvMF_prop4ii
}

.obsolete_prop4ii_has_euc <- function(xe) {
  !is.null(xe) && ncol(as.matrix(xe)) > 0L
}

.obsolete_prop4ii_vmf_log_norm_const_exact <- function(k, p) {
  # FIDELITY EDIT (see file header): the delivered code guarded this with
  #   if (exists("vMF_log_norm_const_exact", mode = "function")) return(...)
  # Under source() into the global environment that guard was FALSE, because sphm does
  # not export vMF_log_norm_const_exact -- so the fallback below is what actually ran.
  # Inside the namespace the guard would flip to TRUE and silently switch implementation,
  # so it is removed to keep this file a faithful reference. The two agree analytically
  # for p = 3: both reduce to log(2*pi) + log(1 - exp(-2k)) + k - log(k).
  k <- as.numeric(k)
  nu <- p / 2 - 1
  out <- numeric(length(k))
  small <- k < 1e-7
  out[small] <- log(2) + (p / 2) * log(pi) - lgamma(p / 2)
  if (any(!small)) {
    kk <- k[!small]
    logI <- log(besselI(kk, nu, expon.scaled = TRUE)) + kk
    out[!small] <- (p / 2) * log(2 * pi) + logI - nu * log(kk)
  }
  out
}

# Replace the earlier helper so the pure-spherical and optional-Euclidean paths
# both work even when the non-exported package helper is unavailable.
.obsolete_prop4ii_estimate_k <- function(rbar, p) {
  rbar <- min(max(as.numeric(rbar), -1 + 1e-12), 1 - 1e-12)
  objective <- function(k) {
    -.obsolete_prop4ii_vmf_log_norm_const_exact(k, p) + k * rbar
  }
  fit <- optimise(objective, lower = 1e-8, upper = 1e5, maximum = TRUE)
  list(k = fit$maximum, loglik_per_obs = fit$objective)
}

.obsolete_prop4ii_prepare_euc <- function(xe, intercept = TRUE,
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

.obsolete_prop4ii_stiefel_dim <- function(m, k) {
  as.integer(m * k - k * (k + 1L) / 2L)
}

.obsolete_prop4ii_complete_orthogonal <- function(V) {
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

.obsolete_prop4ii_stiefel_cayley <- function(theta, Vref) {
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

  Qref <- .obsolete_prop4ii_complete_orthogonal(Vref)
  Qref %*% .obsolete_prop4ii_cayley(K[upper.tri(K)], m)[, seq_len(k), drop = FALSE]
}

# Extended parameter class.  The original call obsolete_mobius_link_prop4ii(Rtilde0)
# remains valid; Euclidean fields are added only when xe is used.
obsolete_mobius_link_prop4ii <- function(Rtilde0, P = NULL, Qe_star = NULL,
                                Be = NULL, xe_center = NULL,
                                xe_scale = NULL, intercept = TRUE,
                                check = TRUE) {
  Rtilde0 <- as.matrix(Rtilde0)
  storage.mode(Rtilde0) <- "double"
  has_euc <- !is.null(P) || !is.null(Qe_star) || !is.null(Be)

  obj <- list(
    Rtilde0 = Rtilde0,
    P = if (has_euc) as.matrix(P) else NULL,
    Qe_star = if (has_euc) as.matrix(Qe_star) else NULL,
    Be = if (has_euc) as.numeric(Be) else NULL,
    xe_center = if (has_euc) as.numeric(xe_center) else NULL,
    xe_scale = if (has_euc) as.numeric(xe_scale) else NULL,
    intercept = if (has_euc) isTRUE(intercept) else FALSE
  )
  class(obj) <- c("obsolete_mobius_link_prop4ii", "list")
  if (check) obsolete_mobius_link_prop4ii_check(obj)
  obj
}

obsolete_mobius_link_prop4ii_check <- function(obj,
                                      tol = 100 * sqrt(.Machine$double.eps)) {
  if (!inherits(obj, "obsolete_mobius_link_prop4ii")) {
    stop("obj must inherit from 'obsolete_mobius_link_prop4ii'.")
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

obsolete_as_mobius_link_prop4ii <- function(obj, check = TRUE) {
  if (inherits(obj, "obsolete_mobius_link_prop4ii")) {
    if (check) obsolete_mobius_link_prop4ii_check(obj)
    return(obj)
  }
  if (is.matrix(obj)) return(obsolete_mobius_link_prop4ii(obj, check = check))
  if (is.list(obj) && !is.null(obj$Rtilde0)) {
    args <- obj[intersect(
      c("Rtilde0", "P", "Qe_star", "Be", "xe_center", "xe_scale", "intercept"),
      names(obj)
    )]
    args$check <- check
    return(do.call(obsolete_mobius_link_prop4ii, args))
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
    if (is.null(obj$Qe)) return(obsolete_mobius_link_prop4ii(R, check = check))

    if (!is_LinEuc(obj)) {
      stop("Only the LinEuc Euclidean form is supported with type='Prop4ii'.")
    }
    V <- obj$Qe[-1, -1, drop = FALSE]
    return(obsolete_mobius_link_prop4ii(
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

#' @noRd
#' @exportS3Method
print.obsolete_mobius_link_prop4ii <- function(x, ...) {
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

#' @noRd
#' @exportS3Method
dim.obsolete_mobius_link_prop4ii <- function(x) {
  p <- nrow(x$Rtilde0)
  c(p = p, qs = p, qe = if (is.null(x$P)) 0L else length(x$xe_center))
}

# Wrapper: Definition 1 classes continue to use the original implementation.
obsolete_mobius_link <- function(xs = NULL, xe = NULL, param = NULL, check = TRUE) {
  if (!inherits(param, "obsolete_mobius_link_prop4ii")) {
    return(.obsolete_mobius_link_definition1(xs = xs, xe = xe,
                                    param = param, check = check))
  }
  if (check) obsolete_mobius_link_prop4ii_check(param)
  if (is.null(xs)) stop("The Proposition 4(ii) link requires xs.")
  xs <- as.matrix(xs)
  if (ncol(xs) != nrow(param$Rtilde0)) {
    stop("ncol(xs) must equal nrow(Rtilde0).")
  }

  if (is.null(param$P)) {
    if (.obsolete_prop4ii_has_euc(xe)) {
      stop("This fitted Proposition 4(ii) link has no Euclidean component.")
    }
    return(xs %*% t(param$Rtilde0))
  }

  if (is.null(xe)) stop("xe must be supplied for this fitted model.")
  ep <- .obsolete_prop4ii_prepare_euc(
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

.obsolete_prop4ii_mean_spec_euc <- function(Rref, Pref, Vref, bstart,
                                   xe_center, xe_scale, intercept) {
  p <- nrow(Rref)
  k <- p - 1L
  m <- nrow(Vref)
  dR <- .obsolete_prop4ii_skew_dim(p)
  dP <- .obsolete_prop4ii_skew_dim(p)
  dV <- .obsolete_prop4ii_stiefel_dim(m, k)
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
    Pref = .obsolete_prop4ii_nearest_orthogonal(Pref, det_sign = 1),
    Vref = qr.Q(qr(Vref), complete = FALSE)[, seq_len(k), drop = FALSE],
    xe_center = xe_center,
    xe_scale = xe_scale,
    intercept = intercept,
    idxR = idxR, idxP = idxP, idxV = idxV, idxB = idxB,
    par0 = par0, lower = lower, upper = upper
  )
}

.obsolete_prop4ii_unpack_mean_euc <- function(theta, spec) {
  R <- .obsolete_prop4ii_cayley(theta[spec$idxR], spec$p) %*% spec$Rref
  P <- .obsolete_prop4ii_cayley(theta[spec$idxP], spec$p) %*% spec$Pref
  V <- .obsolete_prop4ii_stiefel_cayley(theta[spec$idxV], spec$Vref)
  b <- plogis(theta[spec$idxB])
  mean <- obsolete_mobius_link_prop4ii(
    Rtilde0 = R, P = P, Qe_star = V, Be = b,
    xe_center = spec$xe_center, xe_scale = spec$xe_scale,
    intercept = spec$intercept, check = FALSE
  )
  list(mean = mean, Rtilde0 = R, P = P, Qe_star = V, Be = b)
}

.obsolete_prop4ii_fit_vmf_euc_component <- function(y, xs, xe, spec,
                                           algorithm, opts) {
  eval_f <- function(theta) {
    pars <- .obsolete_prop4ii_unpack_mean_euc(theta, spec)
    pred <- try(obsolete_mobius_link(xs = xs, xe = xe, param = pars$mean,
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
        gradient = .obsolete_prop4ii_fd_gradient(
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
  pars <- .obsolete_prop4ii_unpack_mean_euc(fit$solution, spec)
  pred <- obsolete_mobius_link(xs = xs, xe = xe, param = pars$mean, check = FALSE)
  list(nlopt = fit, pars = pars, pred = pred,
       objective = -mean(rowSums(y * pred)))
}

obsolete_mobius_vMF_prop4ii <- function(
    y, xs, xe = NULL,
    det_constraint = c("orthogonal", "rotation"),
    start = NULL, check = TRUE, intercept = TRUE,
    algorithm = "NLOPT_LD_SLSQP", opts = list(), ...) {

  if (!.obsolete_prop4ii_has_euc(xe)) {
    return(.obsolete_mobius_vMF_prop4ii_spherical_only(
      y = y, xs = xs, xe = NULL,
      det_constraint = det_constraint, start = start, check = check
    ))
  }

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
    .obsolete_prop4ii_check_unit_rows(y, "y")
    .obsolete_prop4ii_check_unit_rows(xs, "xs")
  }
  ep <- .obsolete_prop4ii_prepare_euc(xe, intercept = intercept, fitting = TRUE)
  if (ep$m < k) {
    stop("With Euclidean covariates, ncol(xe) + intercept must be at least p-1.")
  }

  start_obj <- if (is.null(start)) NULL else obsolete_as_mobius_link_prop4ii(start)
  target_signs <- if (det_constraint == "rotation") 1 else c(1, -1)
  fits <- vector("list", length(target_signs))

  for (j in seq_along(target_signs)) {
    sgn <- target_signs[j]
    Rref <- if (!is.null(start_obj) &&
                (if (det(start_obj$Rtilde0) >= 0) 1 else -1) == sgn) {
      start_obj$Rtilde0
    } else {
      .obsolete_prop4ii_procrustes(y, xs, det_sign = sgn)$R
    }
    Pref <- if (!is.null(start_obj) && !is.null(start_obj$P)) {
      start_obj$P
    } else {
      .obsolete_prop4ii_nearest_orthogonal(Rref, det_sign = 1)
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

    spec <- .obsolete_prop4ii_mean_spec_euc(
      Rref = Rref, Pref = Pref, Vref = Vref, bstart = bstart,
      xe_center = ep$center, xe_scale = ep$scale,
      intercept = intercept
    )
    fits[[j]] <- .obsolete_prop4ii_fit_vmf_euc_component(
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
  kfit <- .obsolete_prop4ii_estimate_k(mean(dots), p)
  n <- nrow(y)

  dR <- .obsolete_prop4ii_skew_dim(p)
  dP <- .obsolete_prop4ii_skew_dim(p)
  dV <- .obsolete_prop4ii_stiefel_dim(ep$m, k)
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
  class(out) <- c("obsolete_mobius_vMF_prop4ii", "obsolete_mobius_vMF", "list")
  out
}

# Wrapper retaining the original vMF implementation for its original link types.
obsolete_mobius_vMF <- function(
    y, xs = NULL, xe = NULL, start = NULL,
    type = "SpEuc", fix_qs1 = FALSE,
    fix_qe1 = (type == "LinEuc"), intercept = TRUE,
    lb = NULL, ub = NULL,
    det_constraint = c("orthogonal", "rotation"),
    prop4ii_algorithm = "NLOPT_LD_SLSQP",
    prop4ii_opts = list(), ...) {

  if (.obsolete_is_prop4ii_type(type)) {
    dots <- list(...)
    if (length(dots)) prop4ii_opts <- utils::modifyList(prop4ii_opts, dots)
    return(obsolete_mobius_vMF_prop4ii(
      y = y, xs = xs, xe = xe, start = start,
      det_constraint = det_constraint, intercept = intercept,
      algorithm = prop4ii_algorithm, opts = prop4ii_opts
    ))
  }
  .obsolete_mobius_vMF_definition1(
    y = y, xs = xs, xe = xe, start = start, type = type,
    fix_qs1 = fix_qs1, fix_qe1 = fix_qe1, intercept = intercept,
    lb = lb, ub = ub, ...
  )
}

obsolete_mobius_SvMF_partransport_prelim_prop4ii <- function(
    y, xs, xe = NULL, mean = NULL, G0 = NULL,
    G01behaviour = "p1",
    det_constraint = c("orthogonal", "rotation"), check = TRUE,
    intercept = TRUE, algorithm = "NLOPT_LD_SLSQP",
    opts = list(), ...) {

  if (!.obsolete_prop4ii_has_euc(xe)) {
    return(.obsolete_mobius_SvMF_prelim_prop4ii_spherical_only(
      y = y, xs = xs, xe = NULL, mean = mean, G0 = G0,
      G01behaviour = G01behaviour,
      det_constraint = det_constraint, check = check
    ))
  }

  det_constraint <- match.arg(det_constraint)
  G01behaviour <- match.arg(G01behaviour, c("p1", "free", "fixed"))
  prelim <- obsolete_mobius_vMF_prop4ii(
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

.obsolete_prop4ii_joint_spec_euc <- function(mean_ref, prelim, G0reference,
                                    G01behaviour) {
  mean_ref <- obsolete_as_mobius_link_prop4ii(mean_ref)
  p <- nrow(mean_ref$Rtilde0)
  k <- p - 1L
  mspec <- .obsolete_prop4ii_mean_spec_euc(
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
    .obsolete_prop4ii_skew_dim(p)
  } else {
    .obsolete_prop4ii_skew_dim(p - 1L)
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

  Gstart <- .obsolete_prop4ii_nearest_orthogonal(prelim$G0, det_sign = 1)
  Gcoord <- if (is.null(G0reference)) Gstart else
    .obsolete_prop4ii_nearest_orthogonal(G0reference, det_sign = 1)

  if (G01behaviour == "free") {
    Gref <- Gcoord
    thetaG0 <- .obsolete_prop4ii_inverse_cayley(Gstart %*% t(Gref))
    Ganchor <- NULL
  } else {
    target <- if (G01behaviour == "p1") mean_ref$P[, 1] else Gstart[, 1]
    Ganchor <- .obsolete_prop4ii_rebase_G0(Gcoord, target)
    if (is.null(Ganchor)) Ganchor <- .obsolete_prop4ii_rebase_G0(Gstart, target)
    if (is.null(Ganchor)) stop("Could not construct the G0 reference frame.")
    Gstart2 <- .obsolete_prop4ii_rebase_G0(Gstart, target)
    H0 <- crossprod(Ganchor[, -1, drop = FALSE],
                    Gstart2[, -1, drop = FALSE])
    H0 <- .obsolete_prop4ii_nearest_orthogonal(H0, det_sign = 1)
    thetaG0 <- .obsolete_prop4ii_inverse_cayley(H0)
    Gref <- NULL
  }

  par0 <- c(mspec$par0, rep(0, pos - length(mspec$par0)))
  par0[idxK] <- log(max(prelim$k, 1e-8))
  if (dA > 0L) par0[idxA] <- .obsolete_prop4ii_scales_to_eta(prelim$a)
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

.obsolete_prop4ii_unpack_joint_euc <- function(theta, spec) {
  mp <- .obsolete_prop4ii_unpack_mean_euc(
    theta[seq_len(spec$mean_length)], spec$mean_spec
  )
  kappa <- exp(theta[spec$idxK])
  eta <- if (length(spec$idxA)) theta[spec$idxA] else numeric(0)
  a <- .obsolete_prop4ii_eta_to_scales(eta, spec$a1, spec$p)

  if (spec$G01behaviour == "free") {
    G0 <- .obsolete_prop4ii_cayley(theta[spec$idxG], spec$p) %*% spec$Gref
  } else {
    Gbase <- spec$Ganchor
    if (spec$G01behaviour == "p1") {
      Gbase <- .obsolete_prop4ii_rebase_G0(spec$Ganchor, mp$P[, 1])
      if (is.null(Gbase)) return(NULL)
    }
    H <- if (spec$p == 2L) matrix(1, 1, 1) else
      .obsolete_prop4ii_cayley(theta[spec$idxG], spec$p - 1L)
    B <- diag(spec$p)
    B[-1, -1] <- H
    G0 <- Gbase %*% B
  }

  c(mp, list(k = kappa, a = a, G0 = G0))
}

.obsolete_prop4ii_joint_loglik_euc <- function(y, xs, xe, mean, k, a, G0,
                                      approximate = TRUE) {
  if (!is.finite(k) || k <= 0 || any(!is.finite(a)) || any(a <= 0)) {
    return(-Inf)
  }
  pred <- try(obsolete_mobius_link(xs = xs, xe = xe, param = mean,
                          check = FALSE), silent = TRUE)
  if (inherits(pred, "try-error") || any(!is.finite(pred))) return(-Inf)
  yrot <- try(undo_partransport(y, pred, G01 = G0[, 1]), silent = TRUE)
  if (inherits(yrot, "try-error") || any(!is.finite(yrot))) return(-Inf)
  # Same numerical safeguard as in the spherical-only Proposition 4(ii)
  # likelihood: use the fast approximation when finite and otherwise fall back
  # to the exact R density.
  ld <- NULL
  if (approximate && exists("ldSvMF_cann", mode = "function")) {
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

.obsolete_prop4ii_fit_joint_euc_component <- function(
    y, xs, xe, mean_ref, prelim, G0reference, G01behaviour,
    algorithm, opts) {

  spec <- .obsolete_prop4ii_joint_spec_euc(
    mean_ref = mean_ref, prelim = prelim,
    G0reference = G0reference,
    G01behaviour = G01behaviour
  )
  eval_f <- function(theta) {
    pars <- .obsolete_prop4ii_unpack_joint_euc(theta, spec)
    if (is.null(pars)) return(1e100)
    ll <- .obsolete_prop4ii_joint_loglik_euc(
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
        gradient = .obsolete_prop4ii_fd_gradient(
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
  pars <- .obsolete_prop4ii_unpack_joint_euc(fit$solution, spec)
  ll <- if (is.null(pars)) -Inf else .obsolete_prop4ii_joint_loglik_euc(
    y, xs, xe, pars$mean, pars$k, pars$a, pars$G0,
    approximate = TRUE
  )
  list(nlopt = fit, pars = pars, lLik_approx = ll, spec = spec)
}

obsolete_mobius_SvMF_joint_fit_prop4ii <- function(
    y, xs, xe = NULL, mean, k, a, G0,
    G0reference = NULL, G01behaviour = "p1",
    det_constraint = c("orthogonal", "rotation"),
    algorithm = "NLOPT_LD_SLSQP", opts = list(), check = TRUE,
    intercept = TRUE) {

  if (!.obsolete_prop4ii_has_euc(xe)) {
    return(.obsolete_mobius_SvMF_joint_prop4ii_spherical_only(
      y = y, xs = xs, xe = NULL, mean = mean, k = k, a = a, G0 = G0,
      G0reference = G0reference, G01behaviour = G01behaviour,
      det_constraint = det_constraint, algorithm = algorithm,
      opts = opts, check = check
    ))
  }

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
    .obsolete_prop4ii_check_unit_rows(y, "y")
    .obsolete_prop4ii_check_unit_rows(xs, "xs")
  }

  start_mean <- obsolete_as_mobius_link_prop4ii(mean)
  if (is.null(start_mean$P)) stop("The starting mean has no Euclidean component.")
  start_a <- .obsolete_prop4ii_normalise_scales(a, p, a1 = a[1])
  start_G0 <- .obsolete_prop4ii_nearest_orthogonal(G0, det_sign = 1)
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
      .obsolete_prop4ii_procrustes(y, xs, det_sign = sgn)$R
    }
    mean_ref <- obsolete_mobius_link_prop4ii(
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
      rebased <- .obsolete_prop4ii_rebase_G0(pre_j$G0, mean_ref$P[, 1])
      if (is.null(rebased)) stop("Could not rebase G0 to P[,1].")
      pre_j$G0 <- rebased
    }
    fits[[j]] <- .obsolete_prop4ii_fit_joint_euc_component(
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
  if (exists("standardise_col_signs", mode = "function")) {
    pars$G0[, -1] <- standardise_col_signs(pars$G0[, -1, drop = FALSE])
  }
  if (det(pars$G0) < 0) pars$G0[, p] <- -pars$G0[, p]

  pred <- obsolete_mobius_link(xs = xs, xe = xe, param = pars$mean, check = FALSE)
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
  dmean <- .obsolete_prop4ii_skew_dim(p) +       # Rtilde0
    .obsolete_prop4ii_skew_dim(p) +              # P
    .obsolete_prop4ii_stiefel_dim(ep_m, kmean) + # Qe_star
    kmean                               # Be
  dscale <- max(p - 2L, 0L)
  dG0 <- if (G01behaviour == "free") {
    .obsolete_prop4ii_skew_dim(p)
  } else {
    .obsolete_prop4ii_skew_dim(p - 1L)
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
  class(out) <- c("obsolete_mobius_SvMF_prop4ii", "obsolete_mobius_SvMF", "list")
  out
}

obsolete_mobius_SvMF_prop4ii <- function(
    y, xs, xe = NULL, mean = NULL, k = NULL, a = NULL,
    G0 = NULL, G0reference = NULL, G01behaviour = "p1",
    doprelim = TRUE,
    det_constraint = c("orthogonal", "rotation"),
    algorithm = "NLOPT_LD_SLSQP", opts = list(), check = TRUE,
    intercept = TRUE, ...) {

  if (!.obsolete_prop4ii_has_euc(xe)) {
    return(.obsolete_mobius_SvMF_prop4ii_spherical_only(
      y = y, xs = xs, xe = NULL, mean = mean, k = k, a = a,
      G0 = G0, G0reference = G0reference,
      G01behaviour = G01behaviour, doprelim = doprelim,
      det_constraint = det_constraint, algorithm = algorithm,
      opts = utils::modifyList(opts, list(...)), check = check
    ))
  }

  det_constraint <- match.arg(det_constraint)
  G01behaviour <- match.arg(G01behaviour, c("p1", "free", "fixed"))
  dots <- list(...)
  if (length(dots)) opts <- utils::modifyList(opts, dots)

  if (doprelim) {
    preest <- obsolete_mobius_SvMF_partransport_prelim_prop4ii(
      y = y, xs = xs, xe = xe, mean = mean, G0 = G0,
      G01behaviour = G01behaviour,
      det_constraint = det_constraint, check = check,
      intercept = intercept, algorithm = algorithm, opts = opts
    )
  } else {
    if (is.null(mean) || is.null(k) || is.null(a) || is.null(G0)) {
      stop("When doprelim=FALSE, mean, k, a, and G0 must all be supplied.")
    }
    mean <- obsolete_as_mobius_link_prop4ii(mean)
    if (is.null(mean$P)) stop("mean must contain the Euclidean component.")
    preest <- list(
      mean = mean, Rtilde0 = mean$Rtilde0,
      P = mean$P, Qe_star = mean$Qe_star, Be = mean$Be,
      k = k, a = a, G0 = G0
    )
  }

  finalest <- obsolete_mobius_SvMF_joint_fit_prop4ii(
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

#' @noRd
#' @exportS3Method
print.obsolete_mobius_vMF_prop4ii <- function(x, ...) {
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

#' @noRd
#' @exportS3Method
print.obsolete_mobius_SvMF_prop4ii <- function(x, ...) {
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

#' @noRd
#' @exportS3Method
predict.obsolete_mobius_vMF_prop4ii <- function(object, newdata = NULL,
                                        xs = NULL, xe = NULL, ...) {
  if (is.list(newdata) && is.null(xs)) {
    xs <- newdata$xs
    if (is.null(xe)) xe <- newdata$xe
  } else if (!is.null(newdata) && is.null(xs)) {
    xs <- newdata
  }
  if (is.null(xs)) xs <- object$xs
  if (is.null(xe) && !is.null(object$mean$P)) xe <- object$xe
  obsolete_mobius_link(xs = as.matrix(xs), xe = xe, param = object$mean)
}

#' @noRd
#' @exportS3Method
predict.obsolete_mobius_SvMF_prop4ii <- function(object, newdata = NULL,
                                         xs = NULL, xe = NULL, ...) {
  predict.obsolete_mobius_vMF_prop4ii(object, newdata = newdata, xs = xs, xe = xe, ...)
}

# Final wrapper.  The user-facing command remains type = "Prop4ii".  Passing
# xe includes the optional LinEuc component; xe = NULL fits the original
# spherical-only Proposition 4(ii) model.
obsolete_mobius_SvMF <- function(
    y, xs = NULL, xe = NULL, mean = NULL, k = NULL, a = NULL,
    G0 = NULL, G0reference = NULL, G01behaviour = "p1",
    type = "LinEuc", fix_qs1 = FALSE,
    fix_qe1 = (type == "LinEuc"), intercept = TRUE,
    doprelim = TRUE,
    det_constraint = c("orthogonal", "rotation"),
    prop4ii_algorithm = "NLOPT_LD_SLSQP",
    prop4ii_opts = list(), ...) {

  if (.obsolete_is_prop4ii_type(type)) {
    return(obsolete_mobius_SvMF_prop4ii(
      y = y, xs = xs, xe = xe, mean = mean, k = k, a = a,
      G0 = G0, G0reference = G0reference,
      G01behaviour = G01behaviour, doprelim = doprelim,
      det_constraint = det_constraint,
      algorithm = prop4ii_algorithm,
      opts = prop4ii_opts,
      intercept = intercept,
      ...
    ))
  }

  .obsolete_mobius_SvMF_definition1(
    y = y, xs = xs, xe = xe, mean = mean, k = k, a = a,
    G0 = G0, G0reference = G0reference,
    G01behaviour = G01behaviour,
    type = type, fix_qs1 = fix_qs1, fix_qe1 = fix_qe1,
    intercept = intercept, doprelim = doprelim, ...
  )
}

# =============================================================================
# Proposition 4(i) extension: spherical-only isotropic scaled Moebius link
#
# Paper parameterisation (strict Proposition 4(i)):
#   qe = 0, p = qs, Bs = beta_s I_(p-1), 0 < beta_s <= 1,
#
#   mu_s(xs) = M_S(xs : Rtilde0, psi),
#   psi = phi r_s1,   phi = (1 - beta_s)/(1 + beta_s),
#
# where Rtilde0 is orthogonal and ||psi|| < 1.  beta_s = 1 gives
# Proposition 4(ii).  We optimise psi through an unconstrained coordinate u,
#   psi = u / sqrt(1 + ||u||^2),
# so that ||psi|| < 1 automatically.
#
# This block is appended after the Proposition 4(ii) implementation so the
# existing type="Prop4ii" behaviour is unchanged.  The new user-facing type is
# type="Prop4i" (or "Proposition4i").
# =============================================================================

.obsolete_is_prop4i_type <- function(type) {
  if (length(type) != 1L || is.na(type)) return(FALSE)
  key <- tolower(gsub("[^[:alnum:]]", "", as.character(type)))
  key %in% c("prop4i", "proposition4i")
}

.obsolete_prop4i_u_to_psi <- function(u) {
  u <- as.numeric(u)
  den <- sqrt(1 + sum(u^2))
  u / den
}

.obsolete_prop4i_psi_to_u <- function(psi, tol = 1e-12) {
  psi <- as.numeric(psi)
  r2 <- sum(psi^2)
  if (!is.finite(r2) || r2 >= 1) {
    r2 <- min(r2, 1 - tol)
    psi <- psi * sqrt(r2 / max(sum(psi^2), tol))
  }
  psi / sqrt(max(1 - sum(psi^2), tol))
}

.obsolete_prop4i_beta_from_psi <- function(psi) {
  phi <- sqrt(sum(as.numeric(psi)^2))
  (1 - phi) / (1 + phi)
}

.obsolete_prop4i_reference_direction <- function(mean, tol = 1e-10) {
  mean <- obsolete_as_mobius_link_prop4i(mean, check = FALSE)
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

obsolete_mobius_link_prop4i <- function(Rtilde0, psi = NULL,
                               beta_s = NULL, rs1 = NULL,
                               check = TRUE) {
  Rtilde0 <- as.matrix(Rtilde0)
  storage.mode(Rtilde0) <- "double"
  p <- nrow(Rtilde0)

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
  rs1_out <- if (phi > 1e-10) psi / phi else rep(NA_real_, p)
  obj <- list(
    Rtilde0 = Rtilde0,
    psi = psi,
    phi = phi,
    beta_s = (1 - phi) / (1 + phi),
    rs1 = rs1_out
  )
  class(obj) <- c("obsolete_mobius_link_prop4i", "list")
  if (check) obsolete_mobius_link_prop4i_check(obj)
  obj
}

obsolete_mobius_link_prop4i_check <- function(obj,
                                     tol = 100 * sqrt(.Machine$double.eps)) {
  if (!inherits(obj, "obsolete_mobius_link_prop4i")) {
    stop("obj must inherit from 'obsolete_mobius_link_prop4i'.")
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
  invisible(NULL)
}

obsolete_as_mobius_link_prop4i <- function(obj, check = TRUE) {
  if (inherits(obj, "obsolete_mobius_link_prop4i")) {
    if (check) obsolete_mobius_link_prop4i_check(obj)
    return(obj)
  }
  if (is.list(obj) && !is.null(obj$Rtilde0)) {
    if (!is.null(obj$psi)) {
      return(obsolete_mobius_link_prop4i(obj$Rtilde0, psi = obj$psi, check = check))
    }
    if (!is.null(obj$beta_s)) {
      return(obsolete_mobius_link_prop4i(obj$Rtilde0, beta_s = obj$beta_s,
                                rs1 = obj$rs1, check = check))
    }
  }
  stop("obj must be a obsolete_mobius_link_prop4i object or a compatible list.")
}

#' @noRd
#' @exportS3Method
print.obsolete_mobius_link_prop4i <- function(x, ...) {
  cat("Proposition 4(i) spherical M\u00f6bius link\n")
  cat("  p = qs =", nrow(x$Rtilde0), ", qe = 0\n")
  cat("  beta_s =", format(x$beta_s, digits = 6), "\n")
  cat("  phi =", format(x$phi, digits = 6), "\n")
  cat("  det(Rtilde0) =", format(det(x$Rtilde0), digits = 6), "\n")
  invisible(x)
}

#' @noRd
#' @exportS3Method
dim.obsolete_mobius_link_prop4i <- function(x) {
  p <- nrow(x$Rtilde0)
  c(p = p, qs = p, qe = 0L)
}

.obsolete_prop4i_map <- function(xs, Rtilde0, psi) {
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

# Save the fully working Proposition 4(ii)+optional-Euclidean dispatcher before
# adding the new Proposition 4(i) class.
.obsolete_mobius_link_before_prop4i <- obsolete_mobius_link
.obsolete_mobius_vMF_before_prop4i <- obsolete_mobius_vMF
.obsolete_mobius_SvMF_before_prop4i <- obsolete_mobius_SvMF

obsolete_mobius_link <- function(xs = NULL, xe = NULL, param = NULL, check = TRUE) {
  if (inherits(param, "obsolete_mobius_link_prop4i")) {
    if (check) obsolete_mobius_link_prop4i_check(param)
    if (is.null(xs)) stop("The Proposition 4(i) link requires xs.")
    if (!is.null(xe) && ncol(as.matrix(xe)) > 0L) {
      stop("Strict Proposition 4(i) has qe = 0; xe must be NULL.")
    }
    xs <- as.matrix(xs)
    if (ncol(xs) != nrow(param$Rtilde0)) {
      stop("ncol(xs) must equal nrow(Rtilde0).")
    }
    return(.obsolete_prop4i_map(xs, param$Rtilde0, param$psi))
  }
  .obsolete_mobius_link_before_prop4i(xs = xs, xe = xe, param = param, check = check)
}

# -----------------------------------------------------------------------------
# vMF preliminary fit for Proposition 4(i)
# -----------------------------------------------------------------------------

.obsolete_prop4i_mean_spec <- function(Rref, psi_start = NULL) {
  p <- nrow(Rref)
  dR <- .obsolete_prop4ii_skew_dim(p)
  idxR <- seq_len(dR)
  idxPsi <- dR + seq_len(p)
  par0 <- numeric(dR + p)
  if (!is.null(psi_start)) par0[idxPsi] <- .obsolete_prop4i_psi_to_u(psi_start)
  list(
    p = p, Rref = Rref,
    idxR = idxR, idxPsi = idxPsi,
    par0 = par0,
    lower = rep(-Inf, length(par0)),
    upper = rep( Inf, length(par0))
  )
}

.obsolete_prop4i_unpack_mean <- function(theta, spec) {
  R <- .obsolete_prop4ii_cayley(theta[spec$idxR], spec$p) %*% spec$Rref
  psi <- .obsolete_prop4i_u_to_psi(theta[spec$idxPsi])
  mean <- obsolete_mobius_link_prop4i(R, psi = psi, check = FALSE)
  list(Rtilde0 = R, psi = psi, mean = mean)
}

.obsolete_prop4i_fit_vmf_component <- function(y, xs, spec,
                                      algorithm = "NLOPT_LD_SLSQP",
                                      opts = list()) {
  eval_f <- function(theta) {
    pars <- .obsolete_prop4i_unpack_mean(theta, spec)
    pred <- try(.obsolete_prop4i_map(xs, pars$Rtilde0, pars$psi), silent = TRUE)
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
        gradient = .obsolete_prop4ii_fd_gradient(
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
  pars <- .obsolete_prop4i_unpack_mean(theta_best, spec)
  pred <- .obsolete_prop4i_map(xs, pars$Rtilde0, pars$psi)

  if (inherits(fit, "try-error")) {
    fit <- list(status = -999L, message = as.character(fit),
                objective = f_best, solution = theta_best)
  } else {
    fit$objective <- f_best
    fit$solution <- theta_best
  }

  list(nlopt = fit, pars = pars, pred = pred, objective = f_best)
}

obsolete_mobius_vMF_prop4i <- function(
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
    .obsolete_prop4ii_check_unit_rows(y, "y")
    .obsolete_prop4ii_check_unit_rows(xs, "xs")
  }

  start_obj <- if (is.null(start)) NULL else obsolete_as_mobius_link_prop4i(start, check = FALSE)
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
      .obsolete_prop4ii_procrustes(y, xs, det_sign = sgn)$R
    }
    spec <- .obsolete_prop4i_mean_spec(Rref, psi_start = start_psi)
    fits[[j]] <- .obsolete_prop4i_fit_vmf_component(
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
  kfit <- .obsolete_prop4ii_estimate_k(mean(dots), ncol(y))
  n <- nrow(y)
  p <- ncol(y)
  dmean <- .obsolete_prop4ii_skew_dim(p) + p
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
  class(out) <- c("obsolete_mobius_vMF_prop4i", "obsolete_mobius_vMF", "list")
  out
}

# -----------------------------------------------------------------------------
# SvMF fit for Proposition 4(i)
# -----------------------------------------------------------------------------

obsolete_mobius_SvMF_partransport_prelim_prop4i <- function(
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
    .obsolete_prop4ii_check_unit_rows(y, "y")
    .obsolete_prop4ii_check_unit_rows(xs, "xs")
  }

  prelim <- obsolete_mobius_vMF_prop4i(
    y = y, xs = xs, xe = NULL, start = mean,
    det_constraint = det_constraint,
    algorithm = algorithm, opts = opts, check = FALSE
  )
  mean <- prelim$mean
  k <- prelim$k
  p1 <- .obsolete_prop4i_reference_direction(mean)

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

.obsolete_prop4i_joint_loglik <- function(y, xs, mean, k, a, G0,
                                 approximate = TRUE) {
  if (!is.finite(k) || k <= 0 || any(!is.finite(a)) || any(a <= 0)) {
    return(-Inf)
  }
  pred <- try(obsolete_mobius_link(xs = xs, xe = NULL, param = mean, check = FALSE),
              silent = TRUE)
  if (inherits(pred, "try-error") || any(!is.finite(pred))) return(-Inf)
  yrot <- try(undo_partransport(y, pred, G01 = G0[, 1]), silent = TRUE)
  if (inherits(yrot, "try-error") || any(!is.finite(yrot))) return(-Inf)

  ld <- NULL
  if (approximate && exists("ldSvMF_cann", mode = "function")) {
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

.obsolete_prop4i_build_joint_spec <- function(Rref, prelim, G0reference,
                                     G01behaviour) {
  p <- nrow(Rref)
  G01behaviour <- match.arg(G01behaviour, c("p1", "free", "fixed"))
  dR <- .obsolete_prop4ii_skew_dim(p)
  dPsi <- p
  dA <- max(p - 2L, 0L)
  dG <- if (G01behaviour == "free") {
    .obsolete_prop4ii_skew_dim(p)
  } else {
    .obsolete_prop4ii_skew_dim(p - 1L)
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

  start_mean <- obsolete_as_mobius_link_prop4i(prelim$mean, check = FALSE)
  Gstart <- .obsolete_prop4ii_nearest_orthogonal(prelim$G0, det_sign = 1)
  Gcoord <- if (is.null(G0reference)) {
    Gstart
  } else {
    .obsolete_prop4ii_nearest_orthogonal(G0reference, det_sign = 1)
  }

  mean_ref <- obsolete_mobius_link_prop4i(Rref, psi = start_mean$psi, check = FALSE)
  target_base <- if (G01behaviour == "p1") {
    .obsolete_prop4i_reference_direction(mean_ref)
  } else {
    Gstart[, 1]
  }

  if (G01behaviour == "free") {
    Gref <- Gcoord
    Q0 <- Gstart %*% t(Gref)
    thetaG0 <- .obsolete_prop4ii_inverse_cayley(Q0)
    Ganchor <- NULL
    fixed_base <- NULL
  } else {
    Ganchor <- .obsolete_prop4ii_rebase_G0(Gcoord, target_base)
    if (is.null(Ganchor)) Ganchor <- .obsolete_prop4ii_rebase_G0(Gstart, target_base)
    if (is.null(Ganchor)) stop("Could not construct the G0 reference frame.")

    Gstart2 <- .obsolete_prop4ii_rebase_G0(Gstart, target_base)
    if (is.null(Gstart2)) stop("Could not rebase the starting G0 frame.")
    H0 <- crossprod(Ganchor[, -1, drop = FALSE],
                    Gstart2[, -1, drop = FALSE])
    H0 <- .obsolete_prop4ii_nearest_orthogonal(H0, det_sign = 1)
    thetaG0 <- .obsolete_prop4ii_inverse_cayley(H0)
    Gref <- NULL
    fixed_base <- if (G01behaviour == "fixed") target_base else NULL
  }

  par0 <- numeric(pos)
  par0[idxR] <- 0
  par0[idxPsi] <- .obsolete_prop4i_psi_to_u(start_mean$psi)
  par0[idxK] <- log(max(prelim$k, 1e-8))
  if (dA > 0L) par0[idxA] <- .obsolete_prop4ii_scales_to_eta(prelim$a)
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

.obsolete_prop4i_unpack_joint <- function(theta, spec) {
  p <- spec$p
  R <- .obsolete_prop4ii_cayley(theta[spec$idxR], p) %*% spec$Rref
  psi <- .obsolete_prop4i_u_to_psi(theta[spec$idxPsi])
  mean <- obsolete_mobius_link_prop4i(R, psi = psi, check = FALSE)
  k <- exp(theta[spec$idxK])
  eta <- if (length(spec$idxA)) theta[spec$idxA] else numeric(0)
  a <- .obsolete_prop4ii_eta_to_scales(eta, spec$a1, p)

  if (spec$G01behaviour == "free") {
    G0 <- .obsolete_prop4ii_cayley(theta[spec$idxG], p) %*% spec$Gref
  } else {
    Gbase <- spec$Ganchor
    if (spec$G01behaviour == "p1") {
      target <- .obsolete_prop4i_reference_direction(mean)
      Gbase <- .obsolete_prop4ii_rebase_G0(spec$Ganchor, target)
      if (is.null(Gbase)) return(NULL)
    }
    H <- if (p == 2L) {
      matrix(1, 1, 1)
    } else {
      .obsolete_prop4ii_cayley(theta[spec$idxG], p - 1L)
    }
    B <- diag(p)
    B[-1, -1] <- H
    G0 <- Gbase %*% B
  }

  list(mean = mean, Rtilde0 = R, psi = psi, k = k, a = a, G0 = G0)
}

.obsolete_prop4i_fit_component <- function(y, xs, Rref, prelim, G0reference,
                                  G01behaviour, algorithm, opts,
                                  approximate = TRUE) {
  spec <- .obsolete_prop4i_build_joint_spec(
    Rref = Rref, prelim = prelim,
    G0reference = G0reference,
    G01behaviour = G01behaviour
  )

  eval_f <- function(theta) {
    pars <- .obsolete_prop4i_unpack_joint(theta, spec)
    if (is.null(pars)) return(1e100)
    ll <- .obsolete_prop4i_joint_loglik(
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
        gradient = .obsolete_prop4ii_fd_gradient(
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
  pars <- .obsolete_prop4i_unpack_joint(theta_best, spec)
  ll <- if (is.null(pars)) -Inf else .obsolete_prop4i_joint_loglik(
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

obsolete_mobius_SvMF_joint_fit_prop4i <- function(
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
    .obsolete_prop4ii_check_unit_rows(y, "y")
    .obsolete_prop4ii_check_unit_rows(xs, "xs")
  }

  start_mean <- obsolete_as_mobius_link_prop4i(mean, check = FALSE)
  start_R <- start_mean$Rtilde0
  start_a <- .obsolete_prop4ii_normalise_scales(a, p, a1 = a[1])
  start_G0 <- .obsolete_prop4ii_nearest_orthogonal(G0, det_sign = 1)
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
      .obsolete_prop4ii_procrustes(y, xs, det_sign = sgn)$R
    }
    pre_j <- initial
    pre_j$mean <- obsolete_mobius_link_prop4i(
      Rref, psi = start_mean$psi, check = FALSE
    )
    pre_j$Rtilde0 <- Rref
    pre_j$psi <- start_mean$psi
    if (G01behaviour == "p1") {
      newbase <- .obsolete_prop4i_reference_direction(pre_j$mean)
      rebased <- .obsolete_prop4ii_rebase_G0(pre_j$G0, newbase)
      if (is.null(rebased)) stop("Could not rebase G0 for Proposition 4(i).")
      pre_j$G0 <- rebased
    }
    fits[[j]] <- .obsolete_prop4i_fit_component(
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
  if (exists("standardise_col_signs", mode = "function")) {
    pars$G0[, -1] <- standardise_col_signs(pars$G0[, -1, drop = FALSE])
  }
  if (det(pars$G0) < 0) pars$G0[, p] <- -pars$G0[, p]

  pred <- obsolete_mobius_link(xs = xs, xe = NULL, param = pars$mean, check = FALSE)
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

  dmean <- .obsolete_prop4ii_skew_dim(p) + p
  dscale <- max(p - 2L, 0L)
  dG0 <- if (G01behaviour == "free") {
    .obsolete_prop4ii_skew_dim(p)
  } else {
    .obsolete_prop4ii_skew_dim(p - 1L)
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
  class(out) <- c("obsolete_mobius_SvMF_prop4i", "obsolete_mobius_SvMF", "list")
  out
}

obsolete_mobius_SvMF_prop4i <- function(
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
    preest <- obsolete_mobius_SvMF_partransport_prelim_prop4i(
      y = y, xs = xs, xe = NULL, mean = mean, G0 = G0,
      G01behaviour = G01behaviour,
      det_constraint = det_constraint,
      algorithm = algorithm, opts = opts, check = check
    )
  } else {
    if (is.null(mean) || is.null(k) || is.null(a) || is.null(G0)) {
      stop("When doprelim = FALSE, mean, k, a, and G0 must all be supplied.")
    }
    mean <- obsolete_as_mobius_link_prop4i(mean)
    preest <- list(
      mean = mean, Rtilde0 = mean$Rtilde0, psi = mean$psi,
      beta_s = mean$beta_s, k = k, a = a, G0 = G0
    )
  }

  finalest <- obsolete_mobius_SvMF_joint_fit_prop4i(
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

#' @noRd
#' @exportS3Method
print.obsolete_mobius_vMF_prop4i <- function(x, ...) {
  cat("vMF regression with Proposition 4(i) spherical M\u00f6bius link\n")
  cat("  beta_s =", format(x$beta_s, digits = 6), "\n")
  cat("  phi =", format(x$phi, digits = 6), "\n")
  cat("  det(Rtilde0) =", format(det(x$Rtilde0), digits = 6), "\n")
  cat("  k =", format(x$k, digits = 6), "\n")
  cat("  log-likelihood =", format(x$lLik, digits = 7), "\n")
  cat("  AIC =", format(x$AIC, digits = 7), "\n")
  invisible(x)
}

#' @noRd
#' @exportS3Method
print.obsolete_mobius_SvMF_prop4i <- function(x, ...) {
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

#' @noRd
#' @exportS3Method
predict.obsolete_mobius_vMF_prop4i <- function(object, newdata = NULL,
                                       xs = newdata, ...) {
  if (is.null(xs)) xs <- object$xs
  obsolete_mobius_link(xs = as.matrix(xs), xe = NULL, param = object$mean)
}

#' @noRd
#' @exportS3Method
predict.obsolete_mobius_SvMF_prop4i <- function(object, newdata = NULL,
                                        xs = newdata, ...) {
  predict.obsolete_mobius_vMF_prop4i(object, newdata = newdata, xs = xs, ...)
}

# Final user-facing dispatchers.  Existing Proposition 4(ii) and Definition 1
# behaviour is delegated unchanged to the saved working wrappers.
obsolete_mobius_vMF <- function(
    y, xs = NULL, xe = NULL, start = NULL,
    type = "SpEuc", fix_qs1 = FALSE,
    fix_qe1 = (type == "LinEuc"), intercept = TRUE,
    lb = NULL, ub = NULL,
    det_constraint = c("orthogonal", "rotation"),
    prop4i_algorithm = "NLOPT_LD_SLSQP",
    prop4i_opts = list(),
    prop4ii_algorithm = "NLOPT_LD_SLSQP",
    prop4ii_opts = list(), ...) {

  if (.obsolete_is_prop4i_type(type)) {
    dots <- list(...)
    if (length(dots)) prop4i_opts <- utils::modifyList(prop4i_opts, dots)
    return(obsolete_mobius_vMF_prop4i(
      y = y, xs = xs, xe = xe, start = start,
      det_constraint = det_constraint,
      algorithm = prop4i_algorithm, opts = prop4i_opts
    ))
  }

  .obsolete_mobius_vMF_before_prop4i(
    y = y, xs = xs, xe = xe, start = start,
    type = type, fix_qs1 = fix_qs1, fix_qe1 = fix_qe1,
    intercept = intercept, lb = lb, ub = ub,
    det_constraint = det_constraint,
    prop4ii_algorithm = prop4ii_algorithm,
    prop4ii_opts = prop4ii_opts, ...
  )
}

obsolete_mobius_SvMF <- function(
    y, xs = NULL, xe = NULL, mean = NULL, k = NULL, a = NULL,
    G0 = NULL, G0reference = NULL, G01behaviour = "p1",
    type = "LinEuc", fix_qs1 = FALSE,
    fix_qe1 = (type == "LinEuc"), intercept = TRUE,
    doprelim = TRUE,
    det_constraint = c("orthogonal", "rotation"),
    prop4i_algorithm = "NLOPT_LD_SLSQP",
    prop4i_opts = list(),
    prop4ii_algorithm = "NLOPT_LD_SLSQP",
    prop4ii_opts = list(), ...) {

  if (.obsolete_is_prop4i_type(type)) {
    return(obsolete_mobius_SvMF_prop4i(
      y = y, xs = xs, xe = xe, mean = mean, k = k, a = a,
      G0 = G0, G0reference = G0reference,
      G01behaviour = G01behaviour, doprelim = doprelim,
      det_constraint = det_constraint,
      algorithm = prop4i_algorithm,
      opts = prop4i_opts,
      ...
    ))
  }

  .obsolete_mobius_SvMF_before_prop4i(
    y = y, xs = xs, xe = xe, mean = mean, k = k, a = a,
    G0 = G0, G0reference = G0reference,
    G01behaviour = G01behaviour,
    type = type, fix_qs1 = fix_qs1, fix_qe1 = fix_qe1,
    intercept = intercept, doprelim = doprelim,
    det_constraint = det_constraint,
    prop4ii_algorithm = prop4ii_algorithm,
    prop4ii_opts = prop4ii_opts,
    ...
  )
}

# =============================================================================
# Proposition 4(i)-constrained spherical component + optional Euclidean terms
#
# IMPORTANT:
# Proposition 4(i) in the paper itself is stated for qe = 0.  The extension
# below keeps the Proposition 4(i) restriction on the SPHERICAL component,
#
#     p = qs,   Bs = beta_s I_(p-1),   0 < beta_s <= 1,
#
# while adding the same linear Euclidean term used by the existing LinEuc
# implementation:
#
#   mu(xs,xe) = P S^{-1}{
#                 beta_s S(Qs^T xs) + Be Qe_star^T x_e,model
#               },
#
# with Rtilde0 = P Qs^T.  Hence beta_s = 1 reduces exactly to the existing
# Proposition 4(ii) + Euclidean model.
#
# User-facing usage:
#   type = "Prop4i", xe = NULL     : strict Proposition 4(i)
#   type = "Prop4i", xe != NULL    : Proposition 4(i)-constrained spherical
#                                    component + LinEuc covariates
# =============================================================================

# Save the working dispatchers before adding the optional-Euclidean Prop. 4(i)
# extension.  The spherical-only Prop. 4(i), Prop. 4(ii), and Definition 1
# behaviours remain unchanged.
.obsolete_mobius_link_before_prop4i_euc <- obsolete_mobius_link
.obsolete_mobius_vMF_before_prop4i_euc <- obsolete_mobius_vMF
.obsolete_mobius_SvMF_before_prop4i_euc <- obsolete_mobius_SvMF


# -----------------------------------------------------------------------------
# Parameter object and mean link
# -----------------------------------------------------------------------------

obsolete_mobius_link_prop4i_euc <- function(
    Rtilde0, P, beta_s, Qe_star, Be,
    xe_center, xe_scale, intercept = TRUE, check = TRUE) {

  Rtilde0 <- as.matrix(Rtilde0)
  P <- as.matrix(P)
  Qe_star <- as.matrix(Qe_star)
  storage.mode(Rtilde0) <- "double"
  storage.mode(P) <- "double"
  storage.mode(Qe_star) <- "double"

  beta_s <- as.numeric(beta_s)[1L]
  Be <- as.numeric(Be)
  xe_center <- as.numeric(xe_center)
  xe_scale <- as.numeric(xe_scale)

  # Rtilde0 = P Qs^T, hence Qs = Rtilde0^T P.
  Qs <- t(Rtilde0) %*% P
  rs1 <- Qs[, 1L]
  phi <- (1 - beta_s) / (1 + beta_s)
  psi <- phi * rs1

  obj <- list(
    Rtilde0 = Rtilde0,
    P = P,
    beta_s = beta_s,
    phi = phi,
    psi = psi,
    rs1 = rs1,
    Qe_star = Qe_star,
    Be = Be,
    xe_center = xe_center,
    xe_scale = xe_scale,
    intercept = isTRUE(intercept)
  )
  class(obj) <- c("obsolete_mobius_link_prop4i_euc", "list")
  if (check) obsolete_mobius_link_prop4i_euc_check(obj)
  obj
}


obsolete_mobius_link_prop4i_euc_check <- function(
    obj, tol = 100 * sqrt(.Machine$double.eps)) {

  if (!inherits(obj, "obsolete_mobius_link_prop4i_euc")) {
    stop("obj must inherit from 'obsolete_mobius_link_prop4i_euc'.")
  }

  R <- obj$Rtilde0
  P <- obj$P
  if (!is.matrix(R) || nrow(R) != ncol(R) || nrow(R) < 2L) {
    stop("Rtilde0 must be a square p by p matrix with p >= 2.")
  }
  p <- nrow(R)
  k <- p - 1L
  if (any(!is.finite(R)) ||
      max(abs(crossprod(R) - diag(p))) > tol) {
    stop("Rtilde0 must be orthogonal.")
  }

  if (!all(dim(P) == c(p, p)) || any(!is.finite(P)) ||
      max(abs(crossprod(P) - diag(p))) > tol || det(P) < 0) {
    stop("P must be a proper orthogonal p by p matrix.")
  }

  if (!is.finite(obj$beta_s) || obj$beta_s <= 0 || obj$beta_s > 1) {
    stop("beta_s must satisfy 0 < beta_s <= 1.")
  }

  V <- obj$Qe_star
  if (!is.matrix(V) || ncol(V) != k || nrow(V) < k ||
      any(!is.finite(V))) {
    stop("Qe_star must be m by (p-1), with m >= p-1.")
  }
  if (max(abs(crossprod(V) - diag(k))) > tol) {
    stop("Qe_star must have orthonormal columns.")
  }

  b <- obj$Be
  if (length(b) != k || any(!is.finite(b)) ||
      any(b <= 0) || any(b >= 1)) {
    stop("Be must contain p-1 values strictly between zero and one.")
  }

  q <- nrow(V) - as.integer(isTRUE(obj$intercept))
  if (q < 0L || length(obj$xe_center) != q ||
      length(obj$xe_scale) != q) {
    stop("Euclidean centering/scaling information has the wrong dimension.")
  }
  if (any(!is.finite(obj$xe_center)) ||
      any(!is.finite(obj$xe_scale)) ||
      any(obj$xe_scale <= 0)) {
    stop("Euclidean centering/scaling information is invalid.")
  }

  invisible(NULL)
}


obsolete_as_mobius_link_prop4i_euc <- function(obj, check = TRUE) {
  if (inherits(obj, "obsolete_mobius_link_prop4i_euc")) {
    if (check) obsolete_mobius_link_prop4i_euc_check(obj)
    return(obj)
  }
  if (is.list(obj) &&
      all(c("Rtilde0", "P", "beta_s", "Qe_star", "Be",
            "xe_center", "xe_scale") %in% names(obj))) {
    return(obsolete_mobius_link_prop4i_euc(
      Rtilde0 = obj$Rtilde0,
      P = obj$P,
      beta_s = obj$beta_s,
      Qe_star = obj$Qe_star,
      Be = obj$Be,
      xe_center = obj$xe_center,
      xe_scale = obj$xe_scale,
      intercept = if (is.null(obj$intercept)) TRUE else obj$intercept,
      check = check
    ))
  }
  stop("obj is not a compatible Proposition 4(i)+Euclidean parameter object.")
}


#' @noRd
#' @exportS3Method
print.obsolete_mobius_link_prop4i_euc <- function(x, ...) {
  cat("Proposition 4(i)-constrained spherical link with LinEuc covariates\n")
  cat("  p = qs =", nrow(x$Rtilde0), "\n")
  cat("  qe =", length(x$xe_center), "\n")
  cat("  beta_s =", format(x$beta_s, digits = 6), "\n")
  cat("  phi =", format(x$phi, digits = 6), "\n")
  cat("  det(Rtilde0) =", format(det(x$Rtilde0), digits = 6), "\n")
  cat("  Be =", paste(format(x$Be, digits = 5), collapse = ", "), "\n")
  invisible(x)
}


#' @noRd
#' @exportS3Method
dim.obsolete_mobius_link_prop4i_euc <- function(x) {
  p <- nrow(x$Rtilde0)
  c(p = p, qs = p, qe = length(x$xe_center))
}


# Mean-link dispatcher.  beta_s = 1 gives exactly the existing Prop. 4(ii)
# optional-Euclidean link.
obsolete_mobius_link <- function(xs = NULL, xe = NULL, param = NULL, check = TRUE) {
  if (!inherits(param, "obsolete_mobius_link_prop4i_euc")) {
    return(.obsolete_mobius_link_before_prop4i_euc(
      xs = xs, xe = xe, param = param, check = check
    ))
  }

  if (check) obsolete_mobius_link_prop4i_euc_check(param)
  if (is.null(xs)) stop("The Proposition 4(i)+Euclidean link requires xs.")
  if (is.null(xe)) stop("xe must be supplied for this fitted model.")

  xs <- as.matrix(xs)
  if (ncol(xs) != nrow(param$Rtilde0)) {
    stop("ncol(xs) must equal nrow(Rtilde0).")
  }

  ep <- .obsolete_prop4ii_prepare_euc(
    xe,
    intercept = param$intercept,
    center = param$xe_center,
    scale = param$xe_scale,
    fitting = FALSE
  )
  if (nrow(xs) != nrow(ep$model)) {
    stop("xs and xe must have the same number of rows.")
  }

  Qs <- t(param$Rtilde0) %*% param$P
  spherical_part <- param$beta_s * Sp(xs %*% Qs)
  euclidean_part <- sweep(
    ep$model %*% param$Qe_star, 2L, param$Be, "*"
  )

  iSp(spherical_part + euclidean_part) %*% t(param$P)
}


# -----------------------------------------------------------------------------
# Mean-parameter optimisation
# -----------------------------------------------------------------------------

.obsolete_prop4i_euc_mean_spec <- function(
    Rref, Pref, Vref, bstart, beta_start,
    xe_center, xe_scale, intercept) {

  base <- .obsolete_prop4ii_mean_spec_euc(
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


.obsolete_prop4i_euc_unpack_mean <- function(theta, spec) {
  mp <- .obsolete_prop4ii_unpack_mean_euc(
    theta[seq_len(spec$base_length)], spec$base
  )
  beta_s <- plogis(theta[spec$idxBeta])

  mean <- obsolete_mobius_link_prop4i_euc(
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


.obsolete_prop4i_euc_fit_vmf_component <- function(
    y, xs, xe, spec,
    algorithm = "NLOPT_LD_SLSQP", opts = list()) {

  eval_f <- function(theta) {
    pars <- .obsolete_prop4i_euc_unpack_mean(theta, spec)
    pred <- try(
      obsolete_mobius_link(xs = xs, xe = xe, param = pars$mean, check = FALSE),
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
        gradient = .obsolete_prop4ii_fd_gradient(
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

  pars <- .obsolete_prop4i_euc_unpack_mean(theta_best, spec)
  pred <- obsolete_mobius_link(
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


# Construct a proper P starting value from a strict Proposition 4(i) mean.
# This is used only for numerical initialisation; the fitted P is free.
.obsolete_prop4i_euc_P_from_strict <- function(mean) {
  mean <- obsolete_as_mobius_link_prop4i(mean, check = FALSE)
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
  .obsolete_prop4ii_nearest_orthogonal(P, det_sign = 1)
}


obsolete_mobius_vMF_prop4i_euc <- function(
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
    .obsolete_prop4ii_check_unit_rows(y, "y")
    .obsolete_prop4ii_check_unit_rows(xs, "xs")
  }

  ep <- .obsolete_prop4ii_prepare_euc(xe, intercept = intercept, fitting = TRUE)
  if (ep$m < kmean) {
    stop("With Euclidean covariates, ncol(xe) + intercept must be at least p-1.")
  }

  if (!is.null(start) && is.list(start) && !is.null(start$mean)) {
    start <- start$mean
  }

  start_euc <- NULL
  start_strict <- NULL
  if (!is.null(start)) {
    if (inherits(start, "obsolete_mobius_link_prop4i_euc")) {
      start_euc <- obsolete_as_mobius_link_prop4i_euc(start, check = FALSE)
    } else if (inherits(start, "obsolete_mobius_link_prop4i")) {
      start_strict <- obsolete_as_mobius_link_prop4i(start, check = FALSE)
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
    .obsolete_prop4i_euc_P_from_strict(start_strict)
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
      .obsolete_prop4ii_procrustes(y, xs, det_sign = sgn)$R
    }

    Pref <- if (!is.null(start_P)) {
      .obsolete_prop4ii_nearest_orthogonal(start_P, det_sign = 1)
    } else {
      .obsolete_prop4ii_nearest_orthogonal(Rref, det_sign = 1)
    }

    spec <- .obsolete_prop4i_euc_mean_spec(
      Rref = Rref,
      Pref = Pref,
      Vref = start_V,
      bstart = start_Be,
      beta_start = start_beta,
      xe_center = ep$center,
      xe_scale = ep$scale,
      intercept = intercept
    )

    fits[[j]] <- .obsolete_prop4i_euc_fit_vmf_component(
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
  kfit <- .obsolete_prop4ii_estimate_k(mean(dots), p)
  n <- nrow(y)

  dR <- .obsolete_prop4ii_skew_dim(p)
  dP <- .obsolete_prop4ii_skew_dim(p)
  dV <- .obsolete_prop4ii_stiefel_dim(ep$m, kmean)
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
    "obsolete_mobius_vMF_prop4i_euc", "obsolete_mobius_vMF_prop4i", "obsolete_mobius_vMF", "list"
  )
  out
}


# -----------------------------------------------------------------------------
# SvMF fit
# -----------------------------------------------------------------------------

obsolete_mobius_SvMF_partransport_prelim_prop4i_euc <- function(
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

  prelim <- obsolete_mobius_vMF_prop4i_euc(
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


.obsolete_prop4i_euc_joint_spec <- function(
    mean_ref, prelim, G0reference, G01behaviour) {

  mean_ref <- obsolete_as_mobius_link_prop4i_euc(mean_ref, check = FALSE)
  p <- nrow(mean_ref$Rtilde0)

  mspec <- .obsolete_prop4i_euc_mean_spec(
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
    .obsolete_prop4ii_skew_dim(p)
  } else {
    .obsolete_prop4ii_skew_dim(p - 1L)
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

  Gstart <- .obsolete_prop4ii_nearest_orthogonal(prelim$G0, det_sign = 1)
  Gcoord <- if (is.null(G0reference)) {
    Gstart
  } else {
    .obsolete_prop4ii_nearest_orthogonal(G0reference, det_sign = 1)
  }

  if (G01behaviour == "free") {
    Gref <- Gcoord
    thetaG0 <- .obsolete_prop4ii_inverse_cayley(Gstart %*% t(Gref))
    Ganchor <- NULL
  } else {
    target <- if (G01behaviour == "p1") {
      mean_ref$P[, 1L]
    } else {
      Gstart[, 1L]
    }

    Ganchor <- .obsolete_prop4ii_rebase_G0(Gcoord, target)
    if (is.null(Ganchor)) {
      Ganchor <- .obsolete_prop4ii_rebase_G0(Gstart, target)
    }
    if (is.null(Ganchor)) {
      stop("Could not construct the G0 reference frame.")
    }

    Gstart2 <- .obsolete_prop4ii_rebase_G0(Gstart, target)
    H0 <- crossprod(
      Ganchor[, -1L, drop = FALSE],
      Gstart2[, -1L, drop = FALSE]
    )
    H0 <- .obsolete_prop4ii_nearest_orthogonal(H0, det_sign = 1)
    thetaG0 <- .obsolete_prop4ii_inverse_cayley(H0)
    Gref <- NULL
  }

  par0 <- c(mspec$par0, rep(0, pos - length(mspec$par0)))
  par0[idxK] <- log(max(prelim$k, 1e-8))
  if (dA > 0L) {
    par0[idxA] <- .obsolete_prop4ii_scales_to_eta(prelim$a)
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


.obsolete_prop4i_euc_unpack_joint <- function(theta, spec) {
  mp <- .obsolete_prop4i_euc_unpack_mean(
    theta[seq_len(spec$mean_length)], spec$mean_spec
  )

  kappa <- exp(theta[spec$idxK])
  eta <- if (length(spec$idxA)) theta[spec$idxA] else numeric(0)
  a <- .obsolete_prop4ii_eta_to_scales(eta, spec$a1, spec$p)

  if (spec$G01behaviour == "free") {
    G0 <- .obsolete_prop4ii_cayley(theta[spec$idxG], spec$p) %*% spec$Gref
  } else {
    Gbase <- spec$Ganchor
    if (spec$G01behaviour == "p1") {
      Gbase <- .obsolete_prop4ii_rebase_G0(
        spec$Ganchor, mp$P[, 1L]
      )
      if (is.null(Gbase)) return(NULL)
    }

    H <- if (spec$p == 2L) {
      matrix(1, 1, 1)
    } else {
      .obsolete_prop4ii_cayley(theta[spec$idxG], spec$p - 1L)
    }

    B <- diag(spec$p)
    B[-1L, -1L] <- H
    G0 <- Gbase %*% B
  }

  c(mp, list(k = kappa, a = a, G0 = G0))
}


.obsolete_prop4i_euc_joint_loglik <- function(
    y, xs, xe, mean, k, a, G0, approximate = TRUE) {

  if (!is.finite(k) || k <= 0 ||
      any(!is.finite(a)) || any(a <= 0)) {
    return(-Inf)
  }

  pred <- try(
    obsolete_mobius_link(xs = xs, xe = xe, param = mean, check = FALSE),
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
  if (approximate && exists("ldSvMF_cann", mode = "function")) {
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


.obsolete_prop4i_euc_fit_joint_component <- function(
    y, xs, xe, mean_ref, prelim, G0reference,
    G01behaviour, algorithm, opts) {

  spec <- .obsolete_prop4i_euc_joint_spec(
    mean_ref = mean_ref,
    prelim = prelim,
    G0reference = G0reference,
    G01behaviour = G01behaviour
  )

  eval_f <- function(theta) {
    pars <- .obsolete_prop4i_euc_unpack_joint(theta, spec)
    if (is.null(pars)) return(1e100)

    ll <- .obsolete_prop4i_euc_joint_loglik(
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
        gradient = .obsolete_prop4ii_fd_gradient(
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

  pars <- .obsolete_prop4i_euc_unpack_joint(theta_best, spec)
  ll <- if (is.null(pars)) {
    -Inf
  } else {
    .obsolete_prop4i_euc_joint_loglik(
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


obsolete_mobius_SvMF_joint_fit_prop4i_euc <- function(
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
    .obsolete_prop4ii_check_unit_rows(y, "y")
    .obsolete_prop4ii_check_unit_rows(xs, "xs")
  }

  start_mean <- obsolete_as_mobius_link_prop4i_euc(mean, check = FALSE)
  start_a <- .obsolete_prop4ii_normalise_scales(a, p, a1 = a[1L])
  start_G0 <- .obsolete_prop4ii_nearest_orthogonal(G0, det_sign = 1)

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
      .obsolete_prop4ii_procrustes(y, xs, det_sign = sgn)$R
    }

    mean_ref <- obsolete_mobius_link_prop4i_euc(
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
      rebased <- .obsolete_prop4ii_rebase_G0(
        pre_j$G0, mean_ref$P[, 1L]
      )
      if (is.null(rebased)) {
        stop("Could not rebase G0 to P[,1].")
      }
      pre_j$G0 <- rebased
    }

    fits[[j]] <- .obsolete_prop4i_euc_fit_joint_component(
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
  if (exists("standardise_col_signs", mode = "function")) {
    pars$G0[, -1L] <- standardise_col_signs(
      pars$G0[, -1L, drop = FALSE]
    )
  }
  if (det(pars$G0) < 0) {
    pars$G0[, p] <- -pars$G0[, p]
  }

  pred <- obsolete_mobius_link(
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

  dmean <- .obsolete_prop4ii_skew_dim(p) +       # Rtilde0
    .obsolete_prop4ii_skew_dim(p) +              # P
    .obsolete_prop4ii_stiefel_dim(ep_m, kmean) + # Qe_star
    kmean +                             # Be
    1L                                  # beta_s

  dscale <- max(p - 2L, 0L)
  dG0 <- if (G01behaviour == "free") {
    .obsolete_prop4ii_skew_dim(p)
  } else {
    .obsolete_prop4ii_skew_dim(p - 1L)
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
    "obsolete_mobius_SvMF_prop4i_euc",
    "obsolete_mobius_SvMF_prop4i",
    "obsolete_mobius_SvMF",
    "list"
  )
  out
}


obsolete_mobius_SvMF_prop4i_euc <- function(
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
    preest <- obsolete_mobius_SvMF_partransport_prelim_prop4i_euc(
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

    mean <- obsolete_as_mobius_link_prop4i_euc(mean, check = FALSE)
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

  finalest <- obsolete_mobius_SvMF_joint_fit_prop4i_euc(
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


#' @noRd
#' @exportS3Method
print.obsolete_mobius_vMF_prop4i_euc <- function(x, ...) {
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


#' @noRd
#' @exportS3Method
print.obsolete_mobius_SvMF_prop4i_euc <- function(x, ...) {
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


#' @noRd
#' @exportS3Method
predict.obsolete_mobius_vMF_prop4i_euc <- function(
    object, newdata = NULL, xs = NULL, xe = NULL, ...) {

  if (is.null(xs)) xs <- if (is.null(newdata)) object$xs else newdata
  if (is.null(xe)) xe <- object$xe

  obsolete_mobius_link(
    xs = as.matrix(xs),
    xe = as.matrix(xe),
    param = object$mean
  )
}


#' @noRd
#' @exportS3Method
predict.obsolete_mobius_SvMF_prop4i_euc <- function(
    object, newdata = NULL, xs = NULL, xe = NULL, ...) {

  predict.obsolete_mobius_vMF_prop4i_euc(
    object, newdata = newdata, xs = xs, xe = xe, ...
  )
}


# -----------------------------------------------------------------------------
# Final user-facing dispatchers
# -----------------------------------------------------------------------------

obsolete_mobius_vMF <- function(
    y, xs = NULL, xe = NULL, start = NULL,
    type = "SpEuc", fix_qs1 = FALSE,
    fix_qe1 = (type == "LinEuc"), intercept = TRUE,
    lb = NULL, ub = NULL,
    det_constraint = c("orthogonal", "rotation"),
    prop4i_algorithm = "NLOPT_LD_SLSQP",
    prop4i_opts = list(),
    prop4ii_algorithm = "NLOPT_LD_SLSQP",
    prop4ii_opts = list(), ...) {

  if (.obsolete_is_prop4i_type(type) && .obsolete_prop4ii_has_euc(xe)) {
    dots <- list(...)
    if (length(dots)) {
      prop4i_opts <- utils::modifyList(prop4i_opts, dots)
    }
    return(obsolete_mobius_vMF_prop4i_euc(
      y = y, xs = xs, xe = xe, start = start,
      det_constraint = det_constraint,
      intercept = intercept,
      algorithm = prop4i_algorithm,
      opts = prop4i_opts
    ))
  }

  .obsolete_mobius_vMF_before_prop4i_euc(
    y = y, xs = xs, xe = xe, start = start,
    type = type,
    fix_qs1 = fix_qs1,
    fix_qe1 = fix_qe1,
    intercept = intercept,
    lb = lb, ub = ub,
    det_constraint = det_constraint,
    prop4i_algorithm = prop4i_algorithm,
    prop4i_opts = prop4i_opts,
    prop4ii_algorithm = prop4ii_algorithm,
    prop4ii_opts = prop4ii_opts,
    ...
  )
}


obsolete_mobius_SvMF <- function(
    y, xs = NULL, xe = NULL, mean = NULL, k = NULL, a = NULL,
    G0 = NULL, G0reference = NULL, G01behaviour = "p1",
    type = "LinEuc", fix_qs1 = FALSE,
    fix_qe1 = (type == "LinEuc"), intercept = TRUE,
    doprelim = TRUE,
    det_constraint = c("orthogonal", "rotation"),
    prop4i_algorithm = "NLOPT_LD_SLSQP",
    prop4i_opts = list(),
    prop4ii_algorithm = "NLOPT_LD_SLSQP",
    prop4ii_opts = list(), ...) {

  if (.obsolete_is_prop4i_type(type) && .obsolete_prop4ii_has_euc(xe)) {
    return(obsolete_mobius_SvMF_prop4i_euc(
      y = y, xs = xs, xe = xe,
      mean = mean, k = k, a = a,
      G0 = G0, G0reference = G0reference,
      G01behaviour = G01behaviour,
      doprelim = doprelim,
      det_constraint = det_constraint,
      intercept = intercept,
      algorithm = prop4i_algorithm,
      opts = prop4i_opts,
      ...
    ))
  }

  .obsolete_mobius_SvMF_before_prop4i_euc(
    y = y, xs = xs, xe = xe,
    mean = mean, k = k, a = a,
    G0 = G0, G0reference = G0reference,
    G01behaviour = G01behaviour,
    type = type,
    fix_qs1 = fix_qs1,
    fix_qe1 = fix_qe1,
    intercept = intercept,
    doprelim = doprelim,
    det_constraint = det_constraint,
    prop4i_algorithm = prop4i_algorithm,
    prop4i_opts = prop4i_opts,
    prop4ii_algorithm = prop4ii_algorithm,
    prop4ii_opts = prop4ii_opts,
    ...
  )
}

# =============================================================================
# End Proposition 4(i) + Euclidean extension
# =============================================================================
