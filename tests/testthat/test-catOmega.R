test_that("Concatenating Omegas behaves as my maths suggests", {
  p <- 4
  qs <- 5
  qe <- 7
  # data generating parameters:
  set.seed(1)
  P <- mclust::randomOrthogonalMatrix(p, p)
  set.seed(2)
  Qs <- mclust::randomOrthogonalMatrix(qs, p)
  set.seed(3)
  Bs <- diag(sort(runif(p-1), decreasing = TRUE))
  set.seed(4)
  Be <- diag(sort(runif(p-1), decreasing = TRUE))
  set.seed(5)
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
  
  # SVD recovers Pstar
  Omega_svd <- svd(Omega)#, nu = nrow(Omega) - 1, nv = nrow(Omega) - 1)
  expect_equal(topos1strow(Omega_svd$u[, -p]), topos1strow(Pstar))
  
  # partial SVD gets back action on covariates
  Omega_part <- diag(Omega_svd$d) %*% t(Omega_svd$v)
  expect_equal(topos1strow(Omega_part %*% enlarges %*% xs)[-p, , drop = FALSE],
               topos1strow(Bs %*% t(Qsstar) %*% xs))
  expect_equal(topos1strow(Omega_part %*% enlargee %*% xe)[-p, , drop = FALSE],
               topos1strow(Be %*% t(Qe) %*% xe))
  
  # trying to get back Be Bs, Qe and Qs
  BQs <- (diag(Omega_svd$d) %*% t(Omega_svd$v))[-p, ] %*% enlarges
  expect_equal(topos1strow(t(BQs)),
               topos1strow(t(Bs %*% t(Qsstar))))
  BQe <- (diag(Omega_svd$d) %*% t(Omega_svd$v))[-p, ] %*% enlargee
  expect_equal(topos1strow(t(BQe)),
               topos1strow(t(Be %*% t(Qe))))
  
  normss <- apply(BQs, 1, vnorm)
  normse <- apply(BQe, 1, vnorm)
  expect_equal(diag(normss), Bs)
  expect_equal(diag(normse), Be)
  expect_equal(topos1strow(t(diag(1/normss) %*% BQs)), topos1strow(Qsstar))
  expect_equal(topos1strow(t(diag(1/normse) %*% BQe)), topos1strow(Qe))
})