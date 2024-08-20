test_that("rSvMF runs", {
  p <- 5
  set.seed(1)
  a <- runif(p)
  a[p] <- a[p]/prod(a[-1])
  set.seed(1); Gamma <- mclust::randomOrthogonalMatrix(p, p)
  obj <- SvMFcann(0.5, a, Gamma)
  set.seed(2)
  y <- rSvMF(10, obj)
  expect_equal(dim(y), c(10, p))
})
