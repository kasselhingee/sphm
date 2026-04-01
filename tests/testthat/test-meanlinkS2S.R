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
  rand_mobius_link_cann__place_in_env()
  
  set.seed(4)
  xs <- matrix(rnorm(2*qs), nrow = 2)
  xs <- sweep(xs, 1, apply(xs, 1, vnorm), FUN = "/")
  xs <- rbind(xs, xs[1, ])
  
  set.seed(5)
  xe <- matrix(rnorm(2*qe), nrow = 2)
  xe <- rbind(xe, xe[1, ])
  
  mnA <- mobius_link(xs, xe, paramobj)
  mnB <- mobius_link(xs, xe, as_mobius_link_Omega(paramobj))
  expect_equal(mnA[1, ], mnA[3, ]) #a good way to check that row vectors are being treated as such.
  
  # Sph only
  pSph <- do.call(mobius_link_cann, paramobj[c("P", "Bs", "Qs")])
  mnA <- mobius_link(xs, param = pSph)
  mnB <- mobius_link(xs, param = as_mobius_link_Omega(pSph))
  expect_equal(mnA, mnB)
  expect_equal(mnB[1, ], drop(P %*% iSp(drop(Bs %*% Sp(drop(t(Qs) %*% xs[1, ]))))))

  # Euc only
  pEuc <- do.call(mobius_link_cann, paramobj[c("P", "Be", "Qe", "ce")])
  mnA <- mobius_link(xe = xe, param = pEuc)
  mnB <- mobius_link(xe = xe, param = as_mobius_link_Omega(pEuc))
  expect_equal(mnA, mnB)
  expect_equal(mnA[1, ], drop(P %*% iSp(drop( Be %*% (t(Qe[,-1]) %*% xe[1, ]) / drop(Qe[, 1] %*% xe[1, ]  + ce) ) )))

  # errors
  expect_error(mobius_link(xs, as_mobius_link_Omega(pSph)), "parameter object")
  expect_error(mobius_link(xe = xe, pEuc), "parameter object")
  expect_error(mobius_link(xe = xe, param = list(P, Be, Qs)), "class")
})

  
