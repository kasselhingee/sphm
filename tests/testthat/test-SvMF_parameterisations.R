test_that("SvMF_cann() creates objects that pass checks", {
  p <- 5
  set.seed(1)
  a <- runif(p)
  a[p] <- a[p]/prod(a[-1])
  set.seed(1); Gamma <- mclust::randomOrthogonalMatrix(p, p)
  obj <- SvMF_cann(0.5, a, Gamma)
  expect_silent(SvMF_cann_check(obj))
})

test_that("SvMF_cann2muV and reverse passes checks", {
  p <- 5
  set.seed(1)
  a <- runif(p)
  a[p] <- a[p]/prod(a[-1])
  a <- c(a[1], sort(a[-1], decreasing = TRUE))
  set.seed(1); Gamma <- mclust::randomOrthogonalMatrix(p, p)
  obj <- SvMF_cann(0.5, a, Gamma)
  
  obj2 <- SvMF_cann2muV(obj)
  expect_silent(SvMF_muV_check(obj2))
  
  expect_equal(SvMF_muV2cann(obj2)[1:2], obj[1:2]) #up to directions of G[, -1]!
  expect_equal(abs(SvMF_muV2cann(obj2)[[3]]), abs(obj[[3]]), ignore_attr = TRUE) #up to directions of G[, -1]!
})



