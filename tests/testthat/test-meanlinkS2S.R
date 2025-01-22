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
  rmnlink_cann__place_in_env()
  
  set.seed(4)
  xs <- matrix(rnorm(2*qs), nrow = 2)
  xs <- sweep(xs, 1, apply(xs, 1, vnorm), FUN = "/")
  xs <- rbind(xs, xs[1, ])
  
  set.seed(5)
  xe <- matrix(rnorm(2*qe), nrow = 2)
  xe <- rbind(xe, xe[1, ])
  
  mnA <- mnlink(xs, xe, paramobj)
  mnB <- mnlink(xs, xe, as_mnlink_Omega(paramobj))
  expect_equal(mnA[1, ], mnA[3, ]) #a good way to check that row vectors are being treated as such.
  
  # Sph only
  pSph <- do.call(mnlink_cann, paramobj[c("P", "Bs", "Qs")])
  mnA <- mnlink(xs, param = pSph)
  mnB <- mnlink(xs, param = as_mnlink_Omega(pSph))
  expect_equal(mnA, mnB)
  expect_equal(mnB[1, ], drop(P %*% iSp(drop(Bs %*% Sp(drop(t(Qs) %*% xs[1, ]))))))

  # Euc only
  pEuc <- do.call(mnlink_cann, paramobj[c("P", "Be", "Qe", "ce")])
  mnA <- mnlink(xe = xe, param = pEuc)
  mnB <- mnlink(xe = xe, param = as_mnlink_Omega(pEuc))
  expect_equal(mnA, mnB)
  expect_equal(mnA[1, ], drop(P %*% iSp(drop( Be %*% (t(Qe[,-1]) %*% xe[1, ] + ce[-1]) / drop(Qe[, 1] %*% xe[1, ]  + ce[1]) ) )))

  # errors
  expect_error(mnlink(xs, as_mnlink_Omega(pSph)), "parameter object")
  expect_error(mnlink(xe = xe, pEuc), "parameter object")
  expect_error(mnlink(xe = xe, param = list(P, Be, Qs)), "class")
})

  
