test_that("iSp inverses Sp", {
  # for a vector on the sphere
  set.seed(64)
  x <- rnorm(5)
  x <- x / vnorm(x)
  expect_equal(iSp(Sp(x)), x)
  
  # for a vector inside the disk, iSp is not an inverse
  x2 <- x/1.5
  expect_false(any(iSp(Sp(x2)) == x2))
})


test_that("Omega and cannonical versions give same result", {
  set.seed(1)
  p <- 3
  q <- 5
  P <- mclust::randomOrthogonalMatrix(p, p)
  set.seed(2)
  Q <- mclust::randomOrthogonalMatrix(q, p)
  set.seed(3)
  B <- diag(sort(runif(p-1), decreasing = TRUE))

  
  set.seed(4)
  x <- matrix(rnorm(2*q), nrow = 2)
  x <- sweep(x, 1, apply(x, 1, vnorm), FUN = "/")
  x <- rbind(x, x[1, ])
  
  # direct canonical 
  mnA <- meanlinkS2S(x[1, ], P, Q, B)
  
  # via reparameterisation
  mnB <- do.call(meanlinkS2S_Omega, c(list(x = x), param_cann2omega(P, Q, B)))
  # back - remember many signs get ignored, only the result of the mean link matters
  newcann <- do.call(param_omega2cann, param_cann2omega(P, Q, B))
  mnC <- do.call(meanlinkS2S, c(list(x = x[1, ]), newcann))
  
  expect_equal(mnA, mnC)
  expect_equal(drop(mnA), mnB[1, ])
  expect_equal(mnB[1,], mnB[3, ])
  expect_equal(do.call(meanlinkS2S_Omega, c(list(x = x[2, ]), param_cann2omega(P, Q, B))), mnB[2, ])
})
  