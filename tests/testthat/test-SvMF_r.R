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

test_that("simple estimators correct at high concentration and data", {
  p <- 5
  set.seed(1)
  a <- c(1, sort(runif(p-1), decreasing = TRUE))
  a[-1] <- a[-1] / (prod(a[-1])^(1/(p-1)))
  set.seed(2)
  Gamma <- mclust::randomOrthogonalMatrix(p, p)
  Gamma <- toBigPosElRot_keepfirst(Gamma) # make Gamma consistentish signs, and a rotation matrix
  k <- 200
  
  obj <- SvMFcann(k, a, Gamma)
  set.seed(3)
  y <- rSvMF(1E4, obj)
  
  Gest <- SvMF_mom_axes(y = y, mu = Gamma[,1])
  diff <- acos(t(Gest[,-1]) %*% Gamma[,-1]) #pairwise geodesic differences, ignoring mu, which isnt estimated
  expect_equal(diag(diff), rep(0, p-1), tolerance = 1E-2)
  expect_equal(diff[upper.tri(diff)], rep(pi/2, (p-1)*(p-2)/2), tolerance = 1E-2)
  
  scales1 <- SvMF_prelim_scales(y, Gamma)
  expect_equal(scales1, a[-1], tolerance = 1E-2)
  
  scales2 <- SvMF_prelim_scales(y, Gest)
  expect_equal(scales2, a[-1], tolerance = 1E-2, ignore_attr = "names")
})
