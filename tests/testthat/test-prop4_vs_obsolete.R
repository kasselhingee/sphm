# Equivalence of the refactored Proposition 4 implementation against the code as Shogo
# Kato delivered it, preserved verbatim in R/obsolete_Mobius_prop4.R under obsolete_-
# prefixed names. See that file's header.
#
# This is the check that the merge preserved behaviour; test-mobius_prop4.R pins the
# individual substitutions so that a failure here can be localised quickly.

test_that("the helper substitutions match what they replaced", {
  # The Cayley pair now goes through the package's C++ cayley_transform() /
  # inverse_cayley_transform() and their lower-triangle vectorisers, where the delivered
  # code filled the UPPER triangle in R. The C++ row-major-lower traversal transposes onto
  # R's column-major-upper order, so the coordinate ORDER is unchanged and only the sign
  # differs: theta_new = -theta_old. Forward and inverse flip together, so every
  # inverse -> optimise -> forward composition is preserved.
  set.seed(7)
  for (p in 2:6) {
    theta <- stats::rnorm(p * (p - 1) / 2, sd = 2)
    expect_equal(.prop4_cayley(theta, p),
                 .obsolete_prop4ii_cayley(-theta, p), tolerance = 1e-12)

    Q <- prop4_rot(p)
    expect_equal(.prop4_inverse_cayley(Q),
                 -.obsolete_prop4ii_inverse_cayley(Q), tolerance = 1e-10)

    # the invariant that holds whatever the convention
    expect_equal(.prop4_cayley(.prop4_inverse_cayley(Q), p), Q, tolerance = 1e-10)
    expect_equal(crossprod(.prop4_cayley(theta, p)), diag(p), tolerance = 1e-12)
  }

  # inverse_cayley_transform() returns NaN when Q + I is singular; the guard in
  # .prop4_inverse_cayley() is the only thing between the callers and that.
  Qsing <- diag(c(-1, -1, 1))
  expect_equal(rcond(Qsing + diag(3)), 0)
  expect_equal(.prop4_inverse_cayley(Qsing), rep(0, 3))

  # .prop4_stiefel_cayley consumes the same convention, so its skew sub-block flips too
  # (its coordinates always start at zero, so nothing downstream can observe the sign)
  set.seed(9)
  for (m in 3:5) for (k in 1:(m - 1)) {
    V <- qr.Q(qr(matrix(stats::rnorm(m * k), m, k)))[, seq_len(k), drop = FALSE]
    dA <- k * (k - 1) / 2
    dB <- (m - k) * k
    theta <- stats::rnorm(dA + dB, sd = 0.5)
    theta_old <- c(if (dA > 0) -theta[seq_len(dA)] else numeric(0), theta[dA + seq_len(dB)])
    expect_equal(.prop4_stiefel_cayley(theta, V),
                 .obsolete_prop4ii_stiefel_cayley(theta_old, V), tolerance = 1e-11)
  }

  # DoF_Stiefel replaces .prop4ii_stiefel_dim
  for (m in 2:8) for (k in 1:m) {
    expect_identical(as.integer(DoF_Stiefel(m, k)),
                     as.integer(.obsolete_prop4ii_stiefel_dim(m, k)))
  }

  # vMF_log_norm_const_exact replaces the delivered fallback
  kk <- c(0.5, 1, 5, 20, 100, 500)
  for (p in c(3, 5)) {
    expect_equal(vMF_log_norm_const_exact(kk, p),
                 .obsolete_prop4ii_vmf_log_norm_const_exact(kk, p), tolerance = 1e-12)
  }
})

test_that("the mean link matches the delivered implementation", {
  set.seed(3)
  p <- 3; n <- 30; q <- 2
  xs <- matrix(stats::rnorm(n * p), n, p); xs <- xs / sqrt(rowSums(xs^2))
  xe <- matrix(stats::rnorm(n * q), n, q, dimnames = list(NULL, c("a", "b")))
  R0 <- prop4_rot(p); P <- prop4_rot(p)
  ep <- .prop4_prepare_euc(xe, intercept = TRUE, fitting = TRUE)
  V <- qr.Q(qr(matrix(stats::rnorm(ep$m * (p - 1)), ep$m, p - 1)))[, seq_len(p - 1), drop = FALSE]
  Be <- c(0.4, 0.25)

  expect_equal(mobius_link(xs = xs, param = mobius_link_prop4ii(R0)),
               obsolete_mobius_link(xs = xs, param = obsolete_mobius_link_prop4ii(R0)))

  expect_equal(
    mobius_link(xs = xs, xe = xe,
                param = mobius_link_prop4ii(R0, P = P, Qe_star = V, Be = Be,
                                            xe_center = ep$center, xe_scale = ep$scale)),
    obsolete_mobius_link(xs = xs, xe = xe,
                param = obsolete_mobius_link_prop4ii(R0, P = P, Qe_star = V, Be = Be,
                                            xe_center = ep$center, xe_scale = ep$scale)))

  for (beta_s in c(1, 0.6)) {
    rs1 <- c(0.3, -0.5, 0.81); rs1 <- rs1 / sqrt(sum(rs1^2))
    expect_equal(
      mobius_link(xs = xs, param = mobius_link_prop4i(R0, beta_s = beta_s, rs1 = rs1)),
      obsolete_mobius_link(xs = xs,
        param = obsolete_mobius_link_prop4i(R0, beta_s = beta_s, rs1 = rs1)))

    expect_equal(
      mobius_link(xs = xs, xe = xe,
        param = mobius_link_prop4i(R0, beta_s = beta_s, P = P, Qe_star = V, Be = Be,
                                   xe_center = ep$center, xe_scale = ep$scale)),
      obsolete_mobius_link(xs = xs, xe = xe,
        param = obsolete_mobius_link_prop4i_euc(R0, P = P, beta_s = beta_s, Qe_star = V,
                                   Be = Be, xe_center = ep$center, xe_scale = ep$scale)))
  }
})

test_that("vMF fits match the delivered implementation", {
  skip_on_cran()
  d <- prop4_midatlantic(); skip_if(is.null(d), "midatlantic.mat not available")
  for (ty in c("Prop4i", "Prop4ii")) {
    new <- mobius_vMF(d$Y, xs = d$Xs, type = ty, det_constraint = "orthogonal")
    old <- obsolete_mobius_vMF(d$Y, xs = d$Xs, type = ty, det_constraint = "orthogonal")
    expect_equal(new$mean$Rtilde0, old$mean$Rtilde0, tolerance = 1e-4, ignore_attr = TRUE)
    expect_equal(new$k, old$k, tolerance = 1e-4)
    expect_equal(new$obj, old$obj, tolerance = 1e-10)
    expect_equal(new$lLik, old$lLik, tolerance = 1e-10)
  }
})

test_that("SvMF fits match the delivered implementation", {
  skip_on_cran()
  d <- prop4_midatlantic(); skip_if(is.null(d), "midatlantic.mat not available")
  args <- list(det_constraint = "orthogonal", G01behaviour = "free")

  new_4i  <- do.call(mobius_SvMF, c(list(d$Y, xs = d$Xs, xe = NULL, type = "Prop4i"), args))
  old_4i  <- do.call(obsolete_mobius_SvMF, c(list(d$Y, xs = d$Xs, xe = NULL, type = "Prop4i"), args))
  expect_fit_equal(new_4i, old_4i)

  new_4ii <- do.call(mobius_SvMF, c(list(d$Y, xs = d$Xs, xe = NULL, type = "Prop4ii"), args))
  old_4ii <- do.call(obsolete_mobius_SvMF, c(list(d$Y, xs = d$Xs, xe = NULL, type = "Prop4ii"), args))
  expect_fit_equal(new_4ii, old_4ii)

  expect_fit_equal(
    do.call(mobius_SvMF, c(list(d$Y, xs = d$Xs, xe = d$westedge, type = "Prop4ii",
                                intercept = TRUE), args)),
    do.call(obsolete_mobius_SvMF, c(list(d$Y, xs = d$Xs, xe = d$westedge, type = "Prop4ii",
                                intercept = TRUE), args)))

  # started from the spherical-only fit, as the vignette does
  expect_fit_equal(
    do.call(mobius_SvMF, c(list(d$Y, xs = d$Xs, xe = d$westedge, type = "Prop4i",
                                mean = new_4i$mean, intercept = TRUE), args)),
    do.call(obsolete_mobius_SvMF, c(list(d$Y, xs = d$Xs, xe = d$westedge, type = "Prop4i",
                                mean = old_4i$mean, intercept = TRUE), args)))
})
