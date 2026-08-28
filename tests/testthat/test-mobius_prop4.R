# Unit tests for the Proposition 4 API. The behaviour-preservation check against the
# code as delivered lives in test-prop4_vs_obsolete.R.

test_that("the Proposition 4(ii) link is a plain rotation of xs", {
  set.seed(11)
  p <- 4; n <- 25
  xs <- matrix(stats::rnorm(n * p), n, p); xs <- xs / sqrt(rowSums(xs^2))
  R0 <- prop4_rot(p)
  expect_equal(mobius_link(xs = xs, param = mobius_link_prop4ii(R0)), xs %*% t(R0))
})

test_that("Proposition 4(i) at beta_s = 1 reduces to Proposition 4(ii)", {
  set.seed(12)
  p <- 3; n <- 25; q <- 2
  xs <- matrix(stats::rnorm(n * p), n, p); xs <- xs / sqrt(rowSums(xs^2))
  xe <- matrix(stats::rnorm(n * q), n, q)
  R0 <- prop4_rot(p); P <- prop4_rot(p)

  expect_equal(mobius_link(xs = xs, param = mobius_link_prop4i(R0, beta_s = 1)),
               mobius_link(xs = xs, param = mobius_link_prop4ii(R0)))

  ep <- .prop4_prepare_euc(xe, intercept = TRUE, fitting = TRUE)
  V <- qr.Q(qr(matrix(stats::rnorm(ep$m * (p - 1)), ep$m, p - 1)))[, seq_len(p - 1), drop = FALSE]
  Be <- c(0.4, 0.25)
  euc <- list(P = P, Qe_star = V, Be = Be, xe_center = ep$center, xe_scale = ep$scale)
  expect_equal(
    mobius_link(xs = xs, xe = xe, param = do.call(mobius_link_prop4i, c(list(R0, beta_s = 1), euc))),
    mobius_link(xs = xs, xe = xe, param = do.call(mobius_link_prop4ii, c(list(R0), euc))))
})

test_that("as_mobius_link_cann() reproduces the Proposition 4 links", {
  set.seed(13)
  p <- 3; n <- 30; q <- 2
  xs <- matrix(stats::rnorm(n * p), n, p); xs <- xs / sqrt(rowSums(xs^2))
  xe <- matrix(stats::rnorm(n * q), n, q, dimnames = list(NULL, c("a", "b")))
  R0 <- prop4_rot(p); P <- prop4_rot(p)
  ep <- .prop4_prepare_euc(xe, intercept = TRUE, fitting = TRUE)
  V <- qr.Q(qr(matrix(stats::rnorm(ep$m * (p - 1)), ep$m, p - 1)))[, seq_len(p - 1), drop = FALSE]
  euc <- list(P = P, Qe_star = V, Be = c(0.4, 0.25),
              xe_center = ep$center, xe_scale = ep$scale, intercept = TRUE)
  rs1 <- c(0.3, -0.5, 0.81); rs1 <- rs1 / sqrt(sum(rs1^2))

  params <- list(
    "prop4ii"            = mobius_link_prop4ii(R0),
    "prop4ii + Euc"      = do.call(mobius_link_prop4ii, c(list(R0), euc)),
    "prop4i, beta_s = 1" = mobius_link_prop4i(R0, beta_s = 1),
    "prop4i, beta_s < 1" = mobius_link_prop4i(R0, beta_s = 0.6, rs1 = rs1),
    "prop4i + Euc"       = do.call(mobius_link_prop4i, c(list(R0, beta_s = 0.6), euc))
  )
  for (nm in names(params)) {
    param <- params[[nm]]
    cann <- as_mobius_link_cann(param)
    has_euc <- !is.null(param$P)
    direct <- mobius_link(xs = xs, xe = if (has_euc) xe else NULL, param = param)
    viacann <- if (!has_euc) {
      mobius_link(xs = xs, param = cann)
    } else {
      # the canonical form is in the fit's internal standardised coordinates, preceded by
      # the LinEuc dummy-zero covariate -- see .prop4_as_cann()
      m <- .prop4_prepare_euc(xe, intercept = attr(cann, "intercept"),
                              center = attr(cann, "xe_center"),
                              scale = attr(cann, "xe_scale"), fitting = FALSE)$model
      mobius_link(xs = xs, xe = cbind(dummyzero = 0, m), param = cann)
    }
    expect_equal(direct, viacann, tolerance = 1e-10, info = nm)
  }
})

test_that("the parameter class round-trips and reports its dimensions", {
  set.seed(14)
  p <- 3
  R0 <- prop4_rot(p); P <- prop4_rot(p)
  V <- qr.Q(qr(matrix(stats::rnorm(3 * 2), 3, 2)))[, 1:2, drop = FALSE]
  euc <- list(P = P, Qe_star = V, Be = c(0.4, 0.25),
              xe_center = c(0, 0), xe_scale = c(1, 1), intercept = TRUE)

  for (obj in list(mobius_link_prop4ii(R0),
                   do.call(mobius_link_prop4ii, c(list(R0), euc)),
                   mobius_link_prop4i(R0, beta_s = 0.6),
                   do.call(mobius_link_prop4i, c(list(R0, beta_s = 0.6), euc)))) {
    coerce <- if (inherits(obj, "mobius_link_prop4i")) as_mobius_link_prop4i else as_mobius_link_prop4ii
    expect_equal(coerce(obj), obj)
    expect_equal(unname(dim(obj)[c("p", "qs")]), c(p, p))
    expect_equal(unname(dim(obj)["qe"]), if (is.null(obj$P)) 0L else 2L)
  }

  # partially supplied Euclidean blocks are rejected
  expect_error(mobius_link_prop4ii(R0, P = P), "all be supplied or all be NULL")
  expect_error(mobius_link_prop4i(R0, beta_s = 0.6, P = P), "all be supplied or all be NULL")
  # Proposition 4(i) needs 0 < beta_s <= 1
  expect_error(mobius_link_prop4i(R0, beta_s = 1.5), "beta_s")
})

test_that("the sign-flip and multistart helpers reject Proposition 4", {
  expect_error(mobius_vMF_multistart(matrix(0), type = "Prop4i"), "Proposition 4")
  expect_error(mobius_SvMF_multistart(matrix(0), type = "Prop4ii"), "Proposition 4")
  fake <- structure(list(mean = mobius_link_prop4ii(diag(3))), class = "list")
  expect_error(mobius_SvMF_signflip_refit(fake), "Proposition 4")
  expect_error(mobius_vMF_signflip_refit(fake), "Proposition 4")
})

test_that("the four types fit simulated data and recover Rtilde0", {
  skip_on_cran()
  set.seed(15)
  p <- 3; n <- 300
  R0 <- prop4_rot(p)
  xs <- matrix(stats::rnorm(n * p), n, p); xs <- xs / sqrt(rowSums(xs^2))
  sim <- rmobius_SvMF(xs = xs, xe = NULL, mean = mobius_link_prop4ii(R0),
                      k = 60, a = c(1, 1.4, 1 / 1.4), G0 = diag(p))
  y <- sim[, seq_len(p)]

  fit <- mobius_SvMF(y, xs = xs, type = "Prop4ii",
                     det_constraint = "orthogonal", G01behaviour = "free")
  expect_s3_class(fit$mean, "mobius_link_prop4ii")
  expect_true(all(is.finite(c(fit$lLik, fit$AIC, fit$DoF, fit$k))))
  expect_lt(max(abs(fit$mean$Rtilde0 - R0)), 0.15)

  # Proposition 4(i) nests Proposition 4(ii), so it cannot fit worse
  fit4i <- mobius_SvMF(y, xs = xs, type = "Prop4i",
                       det_constraint = "orthogonal", G01behaviour = "free")
  expect_s3_class(fit4i$mean, "mobius_link_prop4i")
  expect_gte(fit4i$lLik, fit$lLik - 1e-6)
  expect_lte(fit4i$mean$beta_s, 1)
})
