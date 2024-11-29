test_that("Concatenating Omegas behaves as my maths suggests", {
  p <- 4
  qs <- 5
  qe <- 7
  # data generating parameters:
  set.seed(2)
  P <- mclust::randomOrthogonalMatrix(p, p)
  set.seed(3)
  Qs <- mclust::randomOrthogonalMatrix(qs, p)
  set.seed(4)
  Bs <- diag(sort(runif(p-1), decreasing = TRUE))
  set.seed(5)
  Be <- diag(sort(runif(p-1), decreasing = TRUE))
  set.seed(6)
  Qe <- mclust::randomOrthogonalMatrix(qe, p -1)
  
  Pstar <- P[, -1]
  Qsstar <- Qs[, -1]
  
  Omegas <- Pstar %*% Bs %*% t(Qsstar)
  Omegae <- Pstar %*% Be %*% t(Qe)
 
  Omega <- cbind(Omegas, Omegae)
  
  #Omega can act on spherical covariates
  set.seed(6)
  xs = runif(qs, -1, 1)
  xs <- xs / vnorm(xs)
  enlarges <- rbind(diag(1, qs), matrix(0, qe, qs)) 
  expect_equal(Omega %*% enlarges %*% xs,
               Omegas %*% xs)
  
  # Omega can act on Euc covariates
  set.seed(7)
  xe <- runif(qe, -1, 1)
  enlargee <- rbind(matrix(0, qs, qe), diag(1, qe))
  expect_equal(Omega %*% enlargee %*% xe,
               Omegae %*% xe)
  
  # check manual SVD returns Omega
  manSV <- sqrt(diag(Bs)^2 + diag(Be)^2)
  manrvecs <- rbind(Qsstar %*% Bs, Qe %*% Be)
  manrvecs <- t(t(manrvecs)/manSV)
  expect_equal(Pstar %*% diag(manSV) %*% t(manrvecs), Omega)
  expect_equal(apply(manrvecs, 2, vnorm), rep(1, p-1))
  
  # check SVD matches manual SVD
  Omega_svd <- svd(Omega, nu = nrow(Omega) - 1, nv = nrow(Omega) - 1)
  expect_equal(topos1strow(Omega_svd$u), topos1strow(Pstar))
  expect_equal(Omega_svd$d[1:(p-1)], manSV)
  expect_equal(topos1strow(Omega_svd$v), topos1strow(manrvecs)) 
  
  # recovering Qs and Qe
  expect_equal(topos1strow(apply(Omega_svd$v[1:qs, ], 2, function(v){v/vnorm(v)})),
               topos1strow(Qsstar))
  expect_equal(topos1strow(apply(Omega_svd$v[qs + (1:qe), ], 2, function(v){v/vnorm(v)})),
               topos1strow(Qe))
  
  #recovering Bs and Be
  expect_equal(apply(Omega_svd$v[1:qs, ], 2, vnorm) * Omega_svd$d[1:(p-1)], diag(Bs))
  expect_equal(apply(Omega_svd$v[qs + (1:qe), ], 2, vnorm) * Omega_svd$d[1:(p-1)], diag(Be))
})