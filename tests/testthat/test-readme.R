# Runs the Quick Start block of README.md verbatim, so that a renamed or re-argued
# exported function breaks a test instead of silently rotting the README.
# If you change this test, change the README block to match (and vice versa).

test_that("the README Quick Start runs", {
  # The spherical-only fits warn that `intercept = TRUE` is ignored because there are no
  # Euclidean covariates. That is expected here; it is not what this test is checking.
  suppressWarnings({

    # Spherical covariates on S^2
    set.seed(1)
    n <- 50;  p <- 3
    xs <- matrix(rnorm(n * p), n, p)
    xs <- xs / sqrt(rowSums(xs^2))

    # Simulate responses from the model
    mean_params <- rand_mobius_link_cann(p = p, qs = p, qe = 0)
    sim <- rmobius_SvMF(xs = xs, xe = NULL, mean = mean_params,
                        k = 5, a = c(1, 2, 0.5), G0 = diag(p))
    y <- sim[, 1:p]

    # Step 1: fast preliminary vMF regression (provides good starting values)
    fit_vMF <- mobius_vMF(y = y, xs = xs)

    # Step 2: full SvMF regression with anisotropic error
    fit_SvMF <- mobius_SvMF(y = y, xs = xs, mean = fit_vMF$mean)

    # Inspect the fitted mean link (canonical parameterisation: P, Bs, Qs)
    cann <- as_mobius_link_cann(fit_SvMF$mean)

  })

  # rmobius_SvMF() returns p+1 columns: the simulated response then its log-density,
  # so sim[, 1:p] in the README is the response.
  expect_equal(dim(sim), c(n, p + 1))
  expect_equal(rowSums(y^2), rep(1, n))

  expect_s3_class(fit_vMF$mean, "mobius_link_Omega")
  expect_s3_class(cann, "mobius_link_cann")
  expect_equal(dim(cann$P), c(p, p))
  expect_equal(dim(cann$Bs), c(p - 1, p - 1))
  expect_equal(dim(cann$Qs), c(p, p))
})

test_that("every function named in the README is exported", {
  # The README's "Key functions" table and Quick Start block only help if the names in
  # them are the names the package actually exports.
  readme_functions <- c(
    "rand_mobius_link_cann", "rmobius_SvMF", "mobius_vMF", "mobius_SvMF",
    "as_mobius_link_cann", "mobius_link", "rotated_resid", "mobius_dof",
    "mobius_link_cann", "mobius_link_Omega", "as_mobius_link_Omega")
  expect_true(all(readme_functions %in% getNamespaceExports("sphm")))
})
