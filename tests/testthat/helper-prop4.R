# Shared fixtures for the Proposition 4 tests.

# The midatlantic data, built exactly as vignettes/reproduce_midatlantic.qmd builds it.
# Returns NULL when the vignette data is not available (e.g. an installed-package check).
prop4_midatlantic <- function() {
  f <- testthat::test_path("..", "..", "vignettes", "midatlantic.mat")
  if (!file.exists(f) || !requireNamespace("R.matlab", quietly = TRUE)) return(NULL)
  raw <- R.matlab::readMat(f)
  Y <- t(raw$X); colnames(Y) <- paste0("Y", 1:3)
  Xs <- t(raw$Y); colnames(Xs) <- paste0("X", 1:3)
  list(Y = Y, Xs = Xs,
       westedge = matrix(as.numeric(Y[, 2] < -0.3), ncol = 1,
                         dimnames = list(NULL, "westedge")))
}

# A random orthogonal p x p matrix with determinant +1.
prop4_rot <- function(p) {
  R <- qr.Q(qr(matrix(stats::rnorm(p * p), p, p)))
  if (det(R) < 0) R[, p] <- -R[, p]
  R
}

# Compare two fits field by field.
#
# Two tolerances, because the fits differ in two quite different ways.
#
# The reported quantities -- lLik, AIC, DoF, obj -- agree to ~1e-14 relative, so they are
# checked tightly. A logic change would move them far more than that.
#
# The fitted parameters (k, a, G0, and the mean link) agree to ~1e-7 relative rather than
# exactly, because the merge made two numerically-equivalent-but-not-bit-identical
# substitutions: vMF_log_norm_const_exact() in place of the delivered fallback (they agree
# to 1e-16 relative), and .prop4_cayley() built on cayley() (1e-15). Those perturbations
# move where the optimisers land, not what they land on: k is found by optimise(), whose
# default convergence tolerance is .Machine$double.eps^0.25 ~ 1.2e-4 relative, and the
# likelihood is very flat near the optimum -- for the midatlantic fits k ~ 490 shifts by
# 2e-5 while lLik moves by 2e-12. Testing those fields tighter than the optimiser's own
# tolerance would be testing noise, so they get 1e-4.
expect_fit_equal <- function(new, old, tol_reported = 1e-10, tol_estimated = 1e-4) {
  for (f in c("lLik", "AIC", "DoF", "obj")) {
    testthat::expect_equal(new[[f]], old[[f]], tolerance = tol_reported,
                           ignore_attr = TRUE, info = f)
  }
  for (f in c("k", "a", "G0", "pred", "dists")) {
    testthat::expect_equal(new[[f]], old[[f]], tolerance = tol_estimated,
                           ignore_attr = TRUE, info = f)
  }
  for (f in c("Rtilde0", "P", "beta_s", "phi", "psi", "rs1",
              "Qe_star", "Be", "xe_center", "xe_scale")) {
    testthat::expect_equal(new$mean[[f]], old$mean[[f]], tolerance = tol_estimated,
                           ignore_attr = TRUE, info = paste0("mean$", f))
  }
}
