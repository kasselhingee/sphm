test_that("ll using alignedG_mean link in C++ matches R", {
  p <- 3
  q <- 5
  # data generating parameters:
  set.seed(1)
  P <- mclust::randomOrthogonalMatrix(p, p)
  set.seed(2)
  Q <- mclust::randomOrthogonalMatrix(q, p)
  set.seed(3)
  B <- diag(sort(runif(p-1), decreasing = TRUE))
  omegapar <- as_mnlink_Omega(cannS2S(P,Q,B))
  
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
  # fix a P matrix, may as well be the true one
  ld <- rep(NA, nrow(y))
  for (i in 1:nrow(y)){
    G <- alignedG(ymean[i, ], P)
    ld[i] <- uldSvMF_cann(y[i, , drop = FALSE], k, a, G)
  }
  
  ldcpp <- ull_S2S_alignedG_mean(mnlink_Omega_vec(omegapar), dyn = c(k, a, as.vector(P)), p, cbind(y, x))
  expect_equal(ld, ldcpp)
  
  # compute likelihood when a2, ... is the independent vector and P, k is fixed
  ldcpp <- ull_S2S_alignedG_a(log(a[-c(1,2)]), c(k, a[1]), c(p, mnlink_Omega_vec(omegapar), as.vector(P)), cbind(y, x))
  expect_equal(ld, ldcpp)

  # compute likelihood when k is the only independent vector
  ldcpp <- ull_S2S_alignedG_k(k, c(mnlink_Omega_vec(omegapar), a, as.vector(P)), p, cbind(y, x))
  expect_equal(ld, ldcpp)
})


test_that("maximum likelihood for alignedG link", {
  p <- 3
  q <- 5
  # data generating parameters:
  set.seed(1)
  P <- mclust::randomOrthogonalMatrix(p, p)
  set.seed(2)
  Q <- mclust::randomOrthogonalMatrix(q, p)
  set.seed(3)
  B <- diag(sort(runif(p-1), decreasing = TRUE))
  omegapar <- as_mnlink_Omega(cannS2S(P,Q,B))
  
  #generate covariates uniformly on the sphere
  set.seed(4)
  x <- matrix(rnorm(1000*q), nrow = 1000)
  x <- sweep(x, 1, apply(x, 1, vnorm), FUN = "/")
  
  ymean <- meanlinkS2S(x = x, paramobj = omegapar)
  
  # generate noise according to SvMF
  k <- 20
  a <- c(1, seq(5, 0.2, length.out = p-1))
  set.seed(5)
  y <- t(apply(ymean, 1, function(mn){
    G <- alignedG(mn, P)
    rSvMF(1, SvMFcann(k, a, G))
  }))

  # when at the starting guess, expect it to work well
  out <- optim_alignedG(y, x, a[1], omegapar, k, a[-1], xtol_rel = 1E-4)
  
  expect_equal(out$solution, list(mean = mnlink_Omega_vec(omegapar), k = k, aremaining = a[-1]),
               tolerance = 2E-2)

})

test_that("maximum likelihood for alignedG link p = 4", {
  p <- 4
  q <- 5
  # data generating parameters:
  set.seed(1)
  P <- mclust::randomOrthogonalMatrix(p, p)
  set.seed(2)
  Q <- mclust::randomOrthogonalMatrix(q, p)
  set.seed(3)
  B <- diag(sort(runif(p-1), decreasing = TRUE))
  omegapar <- as_mnlink_Omega(cannS2S(P,Q,B))
  
  #generate covariates uniformly on the sphere
  set.seed(4)
  x <- matrix(rnorm(1000*q), nrow = 1000)
  x <- sweep(x, 1, apply(x, 1, vnorm), FUN = "/")
  
  ymean <- meanlinkS2S(x = x, paramobj = omegapar)
  
  # generate noise according to SvMF
  k <- 30
  a <- c(1, seq(5, 0.2, length.out = p-1))
  a[-1] <- a[-1]/prod(a[-1])^(1/(p-1))
  set.seed(5)
  y <- t(apply(ymean, 1, function(mn){
    G <- alignedG(mn, P)
    rSvMF(1, SvMFcann(k, a, G))
  }))

  # when at the starting guess, expect it to work well
  out <- optim_alignedG(y, x, a[1], omegapar, k, a[-1], xtol_rel = 1E-4)
  
  expect_equal(out$solution, list(mean = mnlink_Omega_vec(omegapar), k = k, aremaining = a[-1]),
               tolerance = 2E-2)

})
