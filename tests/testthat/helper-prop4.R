# Shared fixtures for the Proposition 4 tests.

# A random orthogonal p x p matrix with determinant +1.
prop4_rot <- function(p) {
  R <- qr.Q(qr(matrix(stats::rnorm(p * p), p, p)))
  if (det(R) < 0) R[, p] <- -R[, p]
  R
}

# An m x k matrix with orthonormal columns. qr.Q() of an m x k matrix already returns
# m x k, so no subsetting is needed.
prop4_stiefel <- function(m, k) qr.Q(qr(matrix(stats::rnorm(m * k), m, k)))

# vMF noise about each row of `ymean`, the idiom used throughout test-mobius_vMF.R.
# movMF is an Imports dependency, so it is always available.
prop4_rvMF <- function(ymean, k) {
  t(apply(ymean, 1, function(mn) movMF::rmovMF(1, k * mn)))
}

# The centre and scale a Proposition 4 fit will compute for `xe`.
#
# This is load-bearing for the recovery tests. .prop4_prepare_euc() derives the
# standardisation from the data at fit time (colMeans / sd) and the fitted link stores it,
# and the Euclidean term cannot absorb an affine change of coordinates because Qe_star is
# constrained to have orthonormal columns. A truth built with any other centre and scale
# therefore lies outside the model the fitter searches, and nothing will match.
prop4_euc_std <- function(xe, intercept = TRUE) {
  ep <- .prop4_prepare_euc(xe, intercept = intercept, fitting = TRUE)
  list(center = ep$center, scale = ep$scale, m = ep$m)
}

# Root-mean-square geodesic distance between two sets of unit vectors. Clamped before
# acos(), as the package does throughout (e.g. R/Mobius_prop4i.R:475).
prop4_rms_angle <- function(a, b) {
  sqrt(mean(acos(pmax(-1, pmin(1, rowSums(a * b))))^2))
}

# Check that a fit recovered the truth it was simulated from.
#
# Which quantities may be compared elementwise is dictated by an exact symmetry of the
# Euclidean form of the link -- see the header of test-prop4_recovery.R. In short:
# predictions, Rtilde0, beta_s, psi and rs1 are invariant; P, Qe_star and Be are not, so
# the Euclidean block is checked through the invariants M %*% t(M) and sort(Be).
#
# Every tolerance is supplied by the caller and derived from a standard error there; this
# function chooses none of them.
expect_prop4_recovers <- function(fit, truth, xs, xe = NULL,
                                  tol_pred, tol_R,
                                  tol_beta = NULL, tol_psi = NULL,
                                  tol_gram = NULL, tol_Be = NULL,
                                  k_true = NULL, tol_k = NULL,
                                  a_true = NULL, tol_a = NULL,
                                  G0_true = NULL, tol_G0 = NULL,
                                  info = NULL) {
  pred_fit <- mobius_link(xs = xs, xe = xe, param = fit$mean)
  pred_true <- mobius_link(xs = xs, xe = xe, param = truth)
  testthat::expect_lt(prop4_rms_angle(pred_fit, pred_true), tol_pred)

  testthat::expect_lt(sqrt(sum((fit$mean$Rtilde0 - truth$Rtilde0)^2)), tol_R)

  if (!is.null(tol_beta)) {
    testthat::expect_lt(abs(fit$mean$beta_s - truth$beta_s), tol_beta)
  }
  if (!is.null(tol_psi)) {
    testthat::expect_lt(sqrt(sum((fit$mean$psi - truth$psi)^2)), tol_psi)
    # rs1 is NA whenever a fit lands at beta_s = 1, where it is not identified
    if (all(is.finite(fit$mean$rs1)) && all(is.finite(truth$rs1))) {
      testthat::expect_lt(sqrt(sum((fit$mean$rs1 - truth$rs1)^2)), tol_psi / truth$phi)
    }
  }

  if (!is.null(tol_gram)) {
    M_fit <- sweep(fit$mean$Qe_star, 2, fit$mean$Be, "*")
    M_true <- sweep(truth$Qe_star, 2, truth$Be, "*")
    testthat::expect_lt(max(abs(tcrossprod(M_fit) - tcrossprod(M_true))), tol_gram)
    testthat::expect_lt(max(abs(sort(fit$mean$Be) - sort(truth$Be))), tol_Be)
  }

  # k and a are compared on the absolute scale, which is what the Monte Carlo study
  # measured; expect_equal()'s relative tolerance would not correspond to those SDs.
  if (!is.null(k_true)) {
    testthat::expect_lt(abs(fit$k - k_true), tol_k)
  }
  if (!is.null(a_true)) {
    testthat::expect_lt(max(abs(fit$a - a_true)), tol_a)
  }
  if (!is.null(G0_true)) {
    testthat::expect_lt(
      max(axis_distance(acos(pmax(-1, pmin(1, colSums(fit$G0 * G0_true)))))), tol_G0)
  }
  invisible(fit)
}

# Tolerances for the two non-negative distance statistics.
#
# `pred` (RMS angular prediction error) and `dR` (Frobenius norm of the Rtilde0 error) are
# norms of an error, so their expectation is strictly positive even for a perfect
# estimator: 3 sd on its own would sit at or below the mean and fail most of the time.
# Both are therefore centred on a THEORETICAL expectation -- never on the simulation mean,
# which would absorb any bias the estimator has -- and 3 sd is added to that.

# E[mean d_i^2] = d / (n k) for a correctly specified fit with d free mean-link parameters.
# Validated against the Monte Carlo: 0.0100 predicted vs 0.0103 measured at d = 3.
# `d` is the mean-link parameter count only, which is fit$DoF - 1 for the vMF fits but not
# for the SvMF ones (whose DoF also counts k, a and G0), so it is passed explicitly.
prop4_tol_pred <- function(d, n, k, sd_pred) sqrt(d / (n * k)) + 3 * sd_pred

# For a rotation error exp(A) about R0 with A skew built from delta ~ N_3(0, sigma^2 I),
# ||R0 exp(A) - R0||_F = ||A||_F = sqrt(2) ||delta||, and ||delta|| ~ sigma * chi_3. So
# E[dR] and sd(dR) are both fixed multiples of sigma, and their ratio E/sd = 2.370 is a
# property of the chi_3 distribution alone. That converts the measured sd into the
# theoretical mean without ever touching the simulation mean. Checked against the Monte
# Carlo: predicts 0.0175 / 0.0206 / 0.0239 against measured 0.0179 / 0.0200 / 0.0276.
prop4_tol_R <- function(sd_R) {
  m <- 2 * sqrt(2 / pi)      # E[chi_3]
  s <- sqrt(3 - m^2)         # sd[chi_3]
  (m / s + 3) * sd_R
}

# `a` and the G0 axis errors are compared as max_j |estimate_j - truth_j|, which is again a
# non-negative statistic with a strictly positive mean, so it needs the same treatment as
# pred and dR rather than a bare 3 sd.
#
# For errors that are marginally normal, max_j |error_j| is distributed as the maximum of m
# absolute standard normals scaled by sigma, so E and sd are both fixed multiples of sigma
# and E/sd is a distributional constant: 1.323, 1.872, 2.264 for m = 1, 2, 3 (computed from
# 4e6 standard normal draws -- a property of the normal, not of any fitter). As with
# prop4_tol_R(), that converts the measured sd into the theoretical mean without using the
# simulation mean.
#
# The m components are not exactly independent here -- a[-1] is constrained to have product
# one, and the G0 axes are orthonormal -- so the constant is an approximation. It errs on
# the conservative side for positively dependent components, whose maximum is less spread
# out than the independent case.
prop4_tol_maxabs <- function(m, sd_stat) {
  ratio <- c(1.323, 1.872, 2.264)[m]
  (ratio + 3) * sd_stat
}
