test_that("SvMF_ll is the same using either parameterisatiion", {
  p <- 3
  set.seed(1)
  a <- runif(p)
  a[p] <- a[p]/prod(a[-1])
  set.seed(1); Gamma <- mclust::randomOrthogonalMatrix(p, p)
  obj <- SvMFcann(0.5, a, Gamma)
  set.seed(2)
  y <-  rSvMF(10, obj)
  ll_cann <- SvMF_ll_cann(y, obj)
  
  obj2 <- SvMF_cann2muV(obj)
  ll_muV <- SvMF_ll_muV(y, obj2)
  expect_equal(ll_muV, ll_cann)
})