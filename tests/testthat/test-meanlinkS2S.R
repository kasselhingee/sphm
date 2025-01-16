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


test_that("Omega and cannonical versions give same result", {
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
  x <- matrix(rnorm(2*q), nrow = 2)
  x <- sweep(x, 1, apply(x, 1, vnorm), FUN = "/")
  x <- rbind(x, x[1, ])
  
  # direct canonical 
  mnA <- meanlinkS2S_cann(x[1, ], paramobj)
  
  # via reparameterisation
  mnB <- meanlinkS2S_Omega(x, as_mnlink_Omega(paramobj))
  # back - remember many signs get ignored, only the result of the mean link matters
  newcann <- as_mnlink_cann(as_mnlink_Omega(paramobj))
  mnC <- meanlinkS2S_cann(x[1, ], newcann)
  
  expect_equal(mnA, mnC)
  expect_equal(drop(mnA), mnB[1, ])
  expect_equal(mnB[1,], mnB[3, ])
  expect_equal(meanlinkS2S_Omega(x[2, ], as_mnlink_Omega(paramobj)), mnB[2, ])
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
  
