test_that("SvMFcann() creates objects that pass checks", {
  p <- 5
  set.seed(1)
  a <- runif(p)
  a[p] <- a[p]/prod(a[-1])
  set.seed(1); Gamma <- mclust::randomOrthogonalMatrix(p, p)
  obj <- SvMFcann(0.5, a, Gamma)
  expect_silent(SvMFcann_check(obj))
})

test_that("SvMF_cann2muV and reversible passes checks", {
  p <- 5
  set.seed(1)
  a <- runif(p)
  a[p] <- a[p]/prod(a[-1])
  a <- c(a[1], sort(a[-1], decreasing = TRUE))
  set.seed(1); Gamma <- mclust::randomOrthogonalMatrix(p, p)
  obj <- SvMFcann(0.5, a, Gamma)
  
  obj2 <- SvMF_cann2muV(obj)
  expect_silent(SvMFmuV_check(obj2))
  
  expect_equal(SvMF_muV2cann(obj2)[1:2], obj[1:2]) #up to directions of G[, -1]!
  expect_equal(abs(SvMF_muV2cann(obj2)[3]), abs(obj[3])) #up to directions of G[, -1]!
})


test_that("alignedG() gives an orthogonal matrix", {
  set.seed(1)
  p <- 5
  P <- mclust::randomOrthogonalMatrix(p, p)
  set.seed(3)
  m <- runif(p, -1, 1)
  m <- m/vnorm(m)
  G <- alignedG(m, P)
  expect_equal(G[,1], m)
  expect_equal(G %*% t(G), diag(1, p))
})

