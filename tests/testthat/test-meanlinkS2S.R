test_that("iSp inverses Sp", {
  # for a vector on the sphere
  set.seed(64)
  x <- rnorm(5)
  x <- x / vnorm(x)
  expect_equal(iSp(Sp(x)), x)
  
  #-e1
  expect_equal(Sp(c(1, rep(0, 3))), rep(1E9, 3))
  
  # for a vector inside the disk, iSp is not an inverse
  x2 <- x/1.5
  expect_equal(Sp(x2), x2[-1]/(1 + x2[1]))
  expect_false(any(iSp(Sp(x2)) == x2))
  
  #a matrix
  expect_equal(Sp(rbind(x, x2)), rbind(Sp(x), Sp(x2)), ignore_attr = "dimnames")
  expect_equal(iSp(Sp(rbind(x, x2))), rbind(x, iSp(Sp(x2))), ignore_attr = "dimnames")
})

test_that("Omega and cannonical versions match manual version, and give same result to each other", {
  set.seed(1)
  p <- 3
  P <- mclust::randomOrthogonalMatrix(p, p)
  qs <- 5
  set.seed(2)
  Qs <- mclust::randomOrthogonalMatrix(qs, p)
  set.seed(3)
  Bs <- diag(sort(runif(p-1), decreasing = TRUE))
  qe <- 4
  set.seed(12)
  Qe <- mclust::randomOrthogonalMatrix(qe, p)
  set.seed(13)
  Be <- diag(sort(runif(p-1), decreasing = TRUE))
  set.seed(14)
  ce <- runif(p)
  paramobj <- mnlink_cann(P, Bs = Bs, Qs = Qs, Be = Be, Qe = Qe, ce = ce, check = TRUE)
  
  set.seed(4)
  xs <- matrix(rnorm(2*qs), nrow = 2)
  xs <- sweep(xs, 1, apply(xs, 1, vnorm), FUN = "/")
  xs <- rbind(xs, xs[1, ])
  
  set.seed(5)
  xe <- matrix(rnorm(2*qe), nrow = 2)
  xe <- rbind(xe, xe[1, ])
  
  # direct canonical 
  mnA <- mnlink_pred_cann(xs, xe, paramobj)
  # via reparameterisation
  Om <- as_mnlink_Omega(paramobj)
  mnB <- meanlinkS2Scpp(xs, xe, mnlink_Omega_vec(Om), p)
  expect_equal(mnA, mnB)
  expect_equal(mnA[1, ], mnA[3, ]) #a good way to check that row vectors are being treated as such.
  
  # Sph only
  pSph <- do.call(mnlink_cann, paramobj[c("P", "Bs", "Qs")])
  Sph_Om <- as_mnlink_Omega(pSph)
  mnA <- mnlink_pred_cann(xs, xe, pSph)
  mnB <- meanlinkS2Scpp(xs, xe = matrix(ncol = 0, nrow = nrow(xs)), mnlink_Omega_vec(Sph_Om), p = p)
  expect_equal(mnA, mnB)
  expect_equal(mnB[1, ], drop(P %*% iSp(drop(Bs %*% Sp(drop(t(Qs) %*% xs[1, ]))))))

  # Euc only
  pEuc <- do.call(mnlink_cann, paramobj[c("P", "Be", "Qe", "ce")])
  mnA <- mnlink_pred_cann(xs, xe, pEuc)
  mnB <- meanlinkS2Scpp(xs = matrix(ncol = 0, nrow = nrow(xe)), xe, mnlink_Omega_vec(as_mnlink_Omega(pEuc)), p = p)
  expect_equal(mnA, mnB)
  expect_equal(mnA[1, ], drop(P %*% iSp(drop( Be %*% (t(Qe[,-1]) %*% xe[1, ] + ce[-1]) / drop(Qe[, 1] %*% xe[1, ]  + ce[1]) ) )))

})

test_that("meanlinkS2S() works with a variety of inputs", {
  set.seed(1)
  p <- 3
  q <- 5
  P <- mclust::randomOrthogonalMatrix(p, p)
  set.seed(2)
  Q <- mclust::randomOrthogonalMatrix(q, p)
  set.seed(3)
  B <- diag(sort(runif(p-1), decreasing = TRUE))
  paramobj <- cannS2S(P, Q, B)
  set.seed(4)
  x <- matrix(rnorm(4*q), nrow = 4)
  x <- sweep(x, 1, apply(x, 1, vnorm), FUN = "/")
  expect_equal(meanlinkS2S(x[1, ], P = P, Q = Q, B = B), meanlinkS2S(x[1, ], paramobj = paramobj))
  expect_equal(meanlinkS2S(x[1, ], paramobj = as_mnlink_Omega(paramobj)), meanlinkS2S(x[1, ], paramobj = paramobj))
  expect_equal(meanlinkS2S(x, P = P, Q = Q, B = B), meanlinkS2S(x, paramobj = as_mnlink_Omega(paramobj)))
  expect_equal(meanlinkS2S(x, paramobj = paramobj), meanlinkS2S(x, paramobj = as_mnlink_Omega(paramobj)))
})

test_that("meanlinkS2Scpp() works and matches R version", {
  set.seed(1)
  p <- 3
  q <- 5
  P <- mclust::randomOrthogonalMatrix(p, p)
  set.seed(2)
  Q <- mclust::randomOrthogonalMatrix(q, p)
  set.seed(3)
  B <- diag(sort(runif(p-1), decreasing = TRUE))
  paramobj <- cannS2S(P, Q, B)
  set.seed(4)
  x <- matrix(rnorm(4*q), nrow = 4)
  x <- sweep(x, 1, apply(x, 1, vnorm), FUN = "/")
  ymeanR <- meanlinkS2S(x, paramobj = as_mnlink_Omega(paramobj))
  ymeanCpp <- meanlinkS2Scpp(x, vec = mnlink_Omega_vec(as_mnlink_Omega(paramobj)), p)
  expect_equal(ymeanCpp, ymeanR)
})
  
