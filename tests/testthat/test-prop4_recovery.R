# Simulate-and-recover checks for the Proposition 4 fitters: draw data from a known
# Proposition 4 link and confirm the fitter gets that link back. This is what
# test-prop4_vs_obsolete.R cannot do -- that file pins the implementation against the code
# as delivered, which shows the merge changed nothing, not that the estimator is right.
#
# All eight fitters are covered: {vMF, SvMF} x {Prop4i, Prop4ii} x {xs only, xs and xe}.
#
#
# WHAT MAY BE COMPARED TO TRUTH
#
# The Euclidean form of the link carries an exact finite symmetry, so P, Qe_star and Be
# are not individually identified -- see the header of R/Mobius_prop4_helpers.R for the
# derivation. Predictions, Rtilde0, beta_s, psi and rs1 are invariant and are compared
# directly; the Euclidean block is compared only through the invariants M %*% t(M) and
# sort(Be), where M = Qe_star diag(Be). Asserting on P or on Be in the fitted order would
# fail on roughly half of all seeds.
#
#
# WHERE THE TOLERANCES COME FROM
#
# Every threshold is 3 sd of the statistic's sampling distribution, estimated by Monte
# Carlo over replicate noise draws from a fixed truth and fixed covariates (30 replicates
# for the vMF fits, 8 for the SvMF fits, whose cost is prohibitive). The simulation MEAN is
# deliberately never used: a tolerance centred on the observed mean would absorb whatever
# bias the estimator has and stop the test from ever seeing it.
#
# For the parameters whose deviation from truth is centred at zero under a correct fitter
# -- k, beta_s, and the norm-like Euclidean invariants -- the threshold is simply 3 sd.
#
# `pred`, `dR`, `a` and the G0 axis errors are norms or maxima of an error, so they have a
# strictly positive expectation even for a perfect estimator and 3 sd alone would sit at or
# below their mean -- for the larger models that fails ~90% of the time on a correct fit.
# These are centred on a THEORETICAL expectation instead, via prop4_tol_pred(),
# prop4_tol_R() and prop4_tol_maxabs() in helper-prop4.R, each of which derives that
# expectation from a closed form or a distributional constant rather than from the
# simulation. The measured SDs are recorded in the comment above each call.
#
# Two closed forms corroborate the Monte Carlo and can be used to rescale if n or k
# changes:
#   RMS angular prediction error   ~ sqrt(d / (n k)),  relative sd sqrt(1 / (2 d))
#   sd(k_hat)                      ~ k / sqrt(n)
# where d = fit$DoF - 1 is the number of free mean-link parameters. At n = 500, k = 60
# these give 0.0100 and 2.68 for the spherical-only Prop 4(ii) fit, against measured
# 0.0103 and 2.37.
#
#
# RUNTIME
#
# The Proposition 4 fitters use central finite-difference gradients rather than the C++ AD
# tapes (R/Mobius_prop4_helpers.R), so the SvMF fits are slow: the Euclidean ones dominate
# this file's runtime by a wide margin. n is deliberately kept high while the fitters are
# being validated; once they are trusted, drop n or gate the SvMF+Euclidean test behind an
# environment variable.

# The common truth: covariates, rotations and Euclidean blocks, built under a fixed seed so
# the tolerances above apply. Only the noise draw varies between this and the study.
prop4_recovery_fixture <- function() {
  set.seed(101)
  p <- 3L; n <- 500L; k <- 60; q <- 2L
  xs <- matrix(stats::rnorm(n * p), n, p); xs <- xs / sqrt(rowSums(xs^2))
  xe <- matrix(stats::rnorm(n * q), n, q, dimnames = list(NULL, c("e1", "e2")))
  std <- prop4_euc_std(xe, intercept = TRUE)
  R0 <- prop4_rot(p); P0 <- prop4_rot(p); V0 <- prop4_stiefel(std$m, p - 1L)
  G0 <- prop4_rot(p)
  Be <- c(0.3, 0.15)
  rs1 <- c(0.3, -0.5, 0.81); rs1 <- rs1 / sqrt(sum(rs1^2))
  euc <- list(P = P0, Qe_star = V0, Be = Be,
              xe_center = std$center, xe_scale = std$scale, intercept = TRUE)
  list(
    p = p, n = n, k = k, xs = xs, xe = xe, R0 = R0, G0 = G0,
    a = c(1, 1.4, 1 / 1.4),
    ii_sph = mobius_link_prop4ii(R0),
    i_sph  = mobius_link_prop4i(R0, beta_s = 0.6, rs1 = rs1),
    ii_euc = do.call(mobius_link_prop4ii, c(list(R0), euc)),
    i_euc  = do.call(mobius_link_prop4i, c(list(R0, beta_s = 0.6), euc))
  )
}


test_that("vMF Prop4ii recovers a spherical Prop 4(ii) truth", {
  skip_on_cran()
  f <- prop4_recovery_fixture()
  set.seed(1001)
  y <- prop4_rvMF(mobius_link(xs = f$xs, param = f$ii_sph), f$k)

  fit <- mobius_vMF(y, xs = f$xs, type = "Prop4ii", det_constraint = "orthogonal")
  expect_s3_class(fit$mean, "mobius_link_prop4ii")
  expect_true(all(is.finite(c(fit$lLik, fit$AIC, fit$DoF, fit$k))))

  # MC over 30 noise draws (sd): pred 0.0042, dR 0.0074, k 2.37
  expect_prop4_recovers(fit, f$ii_sph, f$xs,
                        tol_pred = prop4_tol_pred(3, f$n, f$k, 0.0042),
                        tol_R = prop4_tol_R(0.0074),
                        k_true = f$k, tol_k = 7.1)
})


test_that("vMF Prop4i recovers a spherical Prop 4(i) truth with beta_s < 1", {
  skip_on_cran()
  f <- prop4_recovery_fixture()
  set.seed(1001)
  y <- prop4_rvMF(mobius_link(xs = f$xs, param = f$i_sph), f$k)

  fit <- mobius_vMF(y, xs = f$xs, type = "Prop4i", det_constraint = "orthogonal")
  expect_s3_class(fit$mean, "mobius_link_prop4i")

  # MC over 30 noise draws (sd): pred 0.0045, dR 0.0087, beta_s 0.0040, dpsi 0.0029, k 2.34
  expect_prop4_recovers(fit, f$i_sph, f$xs,
                        tol_pred = prop4_tol_pred(6, f$n, f$k, 0.0045),
                        tol_R = prop4_tol_R(0.0087),
                        tol_beta = 0.0120, tol_psi = 0.0086,
                        k_true = f$k, tol_k = 7.0)

  # Proposition 4(i) contains Proposition 4(ii) at beta_s = 1, and its start set includes
  # the Procrustes solution at psi = 0, so this inequality is structural rather than luck.
  # Both fitters minimise -mean(rowSums(y * pred)), so the objectives are comparable.
  fit_ii <- mobius_vMF(y, xs = f$xs, type = "Prop4ii", det_constraint = "orthogonal")
  expect_lte(fit$obj, fit_ii$obj + 1e-8)
  # and beta_s is doing real work on data generated with beta_s = 0.6
  expect_lt(fit$obj, fit_ii$obj)
})


test_that("vMF Prop4ii recovers a Prop 4(ii) truth with Euclidean covariates", {
  skip_on_cran()
  f <- prop4_recovery_fixture()
  set.seed(1001)
  y <- prop4_rvMF(mobius_link(xs = f$xs, xe = f$xe, param = f$ii_euc), f$k)

  fit <- mobius_vMF(y, xs = f$xs, xe = f$xe, type = "Prop4ii",
                    det_constraint = "orthogonal", intercept = TRUE)
  expect_s3_class(fit$mean, "mobius_link_prop4ii")
  expect_false(is.null(fit$mean$P))

  # MC over 30 noise draws (sd): pred 0.0046, dR 0.0101, gram 0.0026, sort(Be) 0.0042, k 2.31
  expect_prop4_recovers(fit, f$ii_euc, f$xs, f$xe,
                        tol_pred = prop4_tol_pred(11, f$n, f$k, 0.0046),
                        tol_R = prop4_tol_R(0.0101),
                        tol_gram = 0.0077, tol_Be = 0.0127,
                        k_true = f$k, tol_k = 6.9)
})


test_that("vMF Prop4i recovers a Prop 4(i) truth with Euclidean covariates", {
  skip_on_cran()
  f <- prop4_recovery_fixture()
  set.seed(1001)
  y <- prop4_rvMF(mobius_link(xs = f$xs, xe = f$xe, param = f$i_euc), f$k)

  fit <- mobius_vMF(y, xs = f$xs, xe = f$xe, type = "Prop4i",
                    det_constraint = "orthogonal", intercept = TRUE)
  expect_s3_class(fit$mean, "mobius_link_prop4i")
  expect_false(is.null(fit$mean$P))

  # MC over 30 noise draws (sd): pred 0.0051, dR 0.0225, beta_s 0.0080, dpsi 0.0048,
  # gram 0.0064, sort(Be) 0.0093, k 2.26. Rtilde0 is the least well determined quantity in
  # this model, hence the comparatively wide tol_R.
  expect_prop4_recovers(fit, f$i_euc, f$xs, f$xe,
                        tol_pred = prop4_tol_pred(12, f$n, f$k, 0.0051),
                        tol_R = prop4_tol_R(0.0225),
                        tol_beta = 0.0240, tol_psi = 0.0144,
                        tol_gram = 0.0191, tol_Be = 0.0280,
                        k_true = f$k, tol_k = 6.8)
})


# ---------------------------------------------------------------------------------------
# SvMF. These fits are slow (see RUNTIME above): roughly 80 s and 130 s for the two
# spherical-only fits, and 19 min and 31 min for the two Euclidean ones, at n = 500.
#
# Cold starts throughout. The Monte Carlo study measured cold starts, so the tolerances
# below apply to them; warm-starting from a spherical-only fit (as the midatlantic vignette
# does for Prop4i+Euclidean) would be a different estimator with unmeasured sampling error.
# ---------------------------------------------------------------------------------------

test_that("SvMF Prop4ii recovers the mean link from a spherical Prop 4(ii) truth", {
  skip_on_cran()
  f <- prop4_recovery_fixture()
  set.seed(2001)
  y <- rmobius_SvMF(xs = f$xs, xe = NULL, mean = f$ii_sph,
                    k = f$k, a = f$a, G0 = f$G0)[, seq_len(f$p)]

  fit <- mobius_SvMF(y, xs = f$xs, type = "Prop4ii",
                     det_constraint = "orthogonal", G01behaviour = "free")
  expect_s3_class(fit$mean, "mobius_link_prop4ii")
  expect_true(all(is.finite(c(fit$lLik, fit$AIC, fit$DoF, fit$k))))

  # MC over 8 noise draws (sd): pred 0.0028, dR 0.0050.
  # k, a and G0 are deliberately NOT checked here -- see the next test.
  expect_prop4_recovers(fit, f$ii_sph, f$xs,
                        tol_pred = prop4_tol_pred(3, f$n, f$k, 0.0028),
                        tol_R = prop4_tol_R(0.0050))
})


test_that("SvMF Prop4ii recovers the error structure from a spherical Prop 4(ii) truth", {
  skip_on_cran()
  # THIS TEST IS EXPECTED TO FAIL: it documents an open defect in the cold-start path of
  # mobius_SvMF_prop4ii() (spherical-only, i.e. xe = NULL). It is deliberately left failing
  # rather than skipped, so the suite stays red until the prelim is fixed.
  #
  # From the default cold start the fit lands about 74 log-likelihood units BELOW the
  # likelihood at the true parameters (532.3 vs 606.0), with a collapsed to isotropic
  # (1.03, 0.97) against a true (1.4, 0.714) and the G0 axes about 84 degrees out.
  # Starting the same fit at the truth, or supplying the true G0, reaches 612.0 -- above
  # the truth, as an MLE should -- and recovers a to 0.05 and the axes to 0.02 rad. The
  # objective and the joint optimiser are therefore sound; the starting values from
  # mobius_SvMF_partransport_prelim_prop4ii() are at fault.
  #
  # The failure is truth-dependent rather than universal, which is what makes it a
  # starting-value problem rather than a broken objective: it occurs for this fixture's G0
  # in all 8 Monte Carlo replicates, and for a true G0 of diag(3) (the configuration the
  # previous Rtilde0-only test used, lLik 72.5 below truth), but at least one other random
  # rotation recovers correctly (+5.7 above truth, axes to 0.016 rad). It occurs for
  # G01behaviour "p1" and "free" alike.
  #
  # Rtilde0 is unaffected either way -- it comes from the closed-form Procrustes prelim --
  # which is why the previous Rtilde0-only test never caught this. Prop4ii+Euclidean and
  # both Prop4i fitters recover correctly.

  f <- prop4_recovery_fixture()
  set.seed(2001)
  y <- rmobius_SvMF(xs = f$xs, xe = NULL, mean = f$ii_sph,
                    k = f$k, a = f$a, G0 = f$G0)[, seq_len(f$p)]
  fit <- mobius_SvMF(y, xs = f$xs, type = "Prop4ii",
                     det_constraint = "orthogonal", G01behaviour = "free")
  # 3 sd; a, G0 and k SDs carried over from the Prop4i spherical MC, the closest measured
  # case (this fit's own SDs describe the broken optimum, not a working estimator)
  expect_prop4_recovers(fit, f$ii_sph, f$xs,
                        tol_pred = prop4_tol_pred(3, f$n, f$k, 0.0028),
                        tol_R = prop4_tol_R(0.0050),
                        k_true = f$k, tol_k = 7.1,
                        a_true = f$a, tol_a = prop4_tol_maxabs(2, 0.0149),
                        G0_true = f$G0, tol_G0 = prop4_tol_maxabs(3, 0.0347))
})


test_that("SvMF Prop4i recovers a spherical Prop 4(i) truth with beta_s < 1", {
  skip_on_cran()
  f <- prop4_recovery_fixture()
  set.seed(2001)
  y <- rmobius_SvMF(xs = f$xs, xe = NULL, mean = f$i_sph,
                    k = f$k, a = f$a, G0 = f$G0)[, seq_len(f$p)]

  fit <- mobius_SvMF(y, xs = f$xs, type = "Prop4i",
                     det_constraint = "orthogonal", G01behaviour = "free")
  expect_s3_class(fit$mean, "mobius_link_prop4i")

  # MC over 8 noise draws (sd): pred 0.0043, dR 0.0092, beta_s 0.0046, k 2.36,
  # max|a - a_true| 0.0149, G0 axis error 0.0347 rad. dpsi was not recorded for the SvMF
  # arm, so its tolerance is the vMF Prop4i spherical measurement.
  expect_prop4_recovers(fit, f$i_sph, f$xs,
                        tol_pred = prop4_tol_pred(6, f$n, f$k, 0.0043),
                        tol_R = prop4_tol_R(0.0092),
                        tol_beta = 0.0138, tol_psi = 0.0086,
                        k_true = f$k, tol_k = 7.1,
                        a_true = f$a, tol_a = prop4_tol_maxabs(2, 0.0149),
                        G0_true = f$G0, tol_G0 = prop4_tol_maxabs(3, 0.0347))

  # Proposition 4(i) nests Proposition 4(ii), so it cannot fit worse
  fit_ii <- mobius_SvMF(y, xs = f$xs, type = "Prop4ii",
                        det_constraint = "orthogonal", G01behaviour = "free")
  expect_gte(fit$lLik, fit_ii$lLik - 1e-6)
})


test_that("SvMF Prop4ii recovers a Prop 4(ii) truth with Euclidean covariates", {
  skip_on_cran()
  f <- prop4_recovery_fixture()
  set.seed(2001)
  y <- rmobius_SvMF(xs = f$xs, xe = f$xe, mean = f$ii_euc,
                    k = f$k, a = f$a, G0 = f$G0)[, seq_len(f$p)]

  fit <- mobius_SvMF(y, xs = f$xs, xe = f$xe, type = "Prop4ii", intercept = TRUE,
                     det_constraint = "orthogonal", G01behaviour = "free")
  expect_s3_class(fit$mean, "mobius_link_prop4ii")
  expect_false(is.null(fit$mean$P))

  # Only one replicate is affordable here (~19 min), so the thresholds are the measured
  # vMF+Euclidean ones: on the sphere the SvMF sampling SDs are 0.67-1.06 times the vMF
  # ones, so reusing them is conservative. a and G0 use the Prop4i spherical SvMF MC.
  # The single SvMF replicate gave pred 0.0240, dR 0.0406, gram 0.0030, sort(Be) 0.0062,
  # k 58.45, max|a - a_true| 0.0421, G0 0.0217 rad -- all comfortably inside.
  expect_prop4_recovers(fit, f$ii_euc, f$xs, f$xe,
                        tol_pred = prop4_tol_pred(11, f$n, f$k, 0.0046),
                        tol_R = prop4_tol_R(0.0101),
                        tol_gram = 0.0077, tol_Be = 0.0127,
                        k_true = f$k, tol_k = 6.9,
                        a_true = f$a, tol_a = prop4_tol_maxabs(2, 0.0149),
                        G0_true = f$G0, tol_G0 = prop4_tol_maxabs(3, 0.0347))
})


test_that("SvMF Prop4i recovers a Prop 4(i) truth with Euclidean covariates", {
  skip_on_cran()
  f <- prop4_recovery_fixture()
  set.seed(2001)
  y <- rmobius_SvMF(xs = f$xs, xe = f$xe, mean = f$i_euc,
                    k = f$k, a = f$a, G0 = f$G0)[, seq_len(f$p)]

  fit <- mobius_SvMF(y, xs = f$xs, xe = f$xe, type = "Prop4i", intercept = TRUE,
                     det_constraint = "orthogonal", G01behaviour = "free")
  expect_s3_class(fit$mean, "mobius_link_prop4i")
  expect_false(is.null(fit$mean$P))

  # As above, one replicate (~31 min); thresholds from the vMF+Euclidean MC. The single
  # SvMF replicate gave pred 0.0255, dR 0.0385, beta_s 0.5980, gram 0.0045,
  # sort(Be) 0.0062, k 58.64, max|a - a_true| 0.0398, G0 0.0450 rad.
  expect_prop4_recovers(fit, f$i_euc, f$xs, f$xe,
                        tol_pred = prop4_tol_pred(12, f$n, f$k, 0.0051),
                        tol_R = prop4_tol_R(0.0225),
                        tol_beta = 0.0240, tol_psi = 0.0144,
                        tol_gram = 0.0191, tol_Be = 0.0280,
                        k_true = f$k, tol_k = 6.8,
                        a_true = f$a, tol_a = prop4_tol_maxabs(2, 0.0149),
                        G0_true = f$G0, tol_G0 = prop4_tol_maxabs(3, 0.0347))
})
