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
  ll_cpp <- do.call(ldSvMF_cann, c(list(y = y), obj))
  expect_equal(ll_cpp, ll_cann)
  
  obj2 <- SvMF_cann2muV(obj)
  ll_muV <- SvMF_ll_muV(y, obj2)
  expect_equal(ll_muV, ll_cann)
  
  expect_equal(dSvMF(y, obj), dSvMF(y, obj2))
  
  ll_cpp <- do.call(ldSvMF_muV, c(list(y = y), obj2))
  expect_equal(ll_cpp, ll_muV)
})

test_that("ll using mean link in C++ matches R", {
  p <- 3
  q <- 5
  # data generating parameters:
  set.seed(1)
  P <- mclust::randomOrthogonalMatrix(p, p)
  set.seed(2)
  Q <- mclust::randomOrthogonalMatrix(q, p)
  set.seed(3)
  B <- diag(sort(runif(p-1), decreasing = TRUE))
  omegapar <- as_OmegaS2S(cannS2S(P,Q,B))
  
  #generate covariates uniformly on the sphere
  set.seed(4)
  x <- matrix(rnorm(1000*q), nrow = 1000)
  x <- sweep(x, 1, apply(x, 1, vnorm), FUN = "/")
  
  ymean <- meanlinkS2S(x = x, paramobj = omegapar)
  
  # generate noise
  if (!requireNamespace("movMF", quietly = TRUE)){skip("Need movMF package")}
  set.seed(5)
  y <- t(apply(ymean, 1, function(mn){movMF::rmovMF(1, 10*mn)}))
  
  # log-likelihood via R
  k <- 10
  a <- c(1, seq(5, 0.2, length.out = p-1))
  # fix a P matrix
  P <- diag(1, p)
  ld <- rep(NA, nrow(y))
  for (i in 1:nrow(y)){
    G <- alignedG(ymean[i, ], P)
    ld[i] <- ldSvMF_cann(y[i, , drop = FALSE], k, a, G)
  }
  
  
})
