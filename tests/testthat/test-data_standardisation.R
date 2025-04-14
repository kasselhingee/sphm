test_that("standardise spherical data gives data with correct moments", {
  rmnlink_cann__place_in_env(3, 5, 4)
  
  #generate covariates Gaussianly
  set.seed(4)
  xe <- matrix(rnorm(100*qe), nrow = 100)
  #generate covariates on the sphere
  set.seed(4)
  xs <- matrix(rnorm(100*qs), nrow = 100)
  xs <- sweep(xs, 1, apply(xs, 1, vnorm), FUN = "/")
  
  ymean <- mnlink(xs = xs, xe = xe, param = paramobj)
  
  # generate noise
  if (!requireNamespace("movMF", quietly = TRUE)){skip("Need movMF package")}
  set.seed(5)
  y <- t(apply(ymean, 1, function(mn){movMF::rmovMF(1, 2*mn)}))
  
  mat1 <- standardise_mat(y)
  y2 <- standardise(y, tG = t(mat1))
  expect_equal(colMeans(y2)[-1], rep(0, p-1))
  # standardisation should do nothing up to sign:
  expect_equal(standardise_mat(y2), diag(diag(standardise_mat(y2))), ignore_attr = TRUE)
  
  # undo standardisation
  expect_equal(destandardise(y2, t(mat1)), y)
})

test_that("standadise_Euc works seemlessly when covariates are all 1", {
  rmnlink_cann__place_in_env(3, 5, 4, preseed = 3)
  
  #generate covariates Gaussianly
  set.seed(4)
  xe <- matrix(rnorm(20*qe), nrow = 20)
  colnames(xe) <- c("x1", "x2", "", "x4")
  
  # check standardisation and destandardisation for xe without constants
  xestd <- standardise_Euc(xe)
  expect_equal(cor(xestd), diag(qe), ignore_attr = TRUE) #check correlation only because not standardising scale
  xe2 <- destandardise_Euc(xestd, center = attr(xestd, "std_center"), rotation = attr(xestd, "std_rotation"))
  expect_equal(xe2, xe)
  
  # check with constants
  xe <- cbind("const1" = 1, "const2" = 1.5, xe)
  xestd <- standardise_Euc(xe)
  expect_equal(xestd[, c(1, 2)], xe[, c(1,2)], ignore_attr = FALSE)
  expect_equal(xestd[, -c(1, 2)], standardise_Euc(xe[, -c(1, 2)]), ignore_attr = TRUE)
  xe2 <- destandardise_Euc(xestd, center = attr(xestd, "std_center"), rotation = attr(xestd, "std_rotation"))
  expect_equal(xe2, xe)
})

test_that("standardisations applied to parameters work", {
  rmnlink_cann__place_in_env(3, 5, 4)
  
  #generate covariates Gaussianly
  set.seed(4)
  xe <- cbind(1, matrix(rnorm(20*(qe-1)), nrow = 20))
  #generate covariates on the sphere
  set.seed(4)
  xs <- matrix(rnorm(20*qs), nrow = 20)
  xs <- sweep(xs, 1, apply(xs, 1, vnorm), FUN = "/")
  
  ymean <- mnlink(xs = xs, xe = xe, param = paramobj)
  
  #standardise xs, xe and response and compute paramobj that includes the reverse
  xsstd <- standardise(xs)
  xestd <- standardise_Euc(xe)
  ystd <- standardise(ymean)
  
  # standardise cann version:
  paramstd <- paramobj
  paramstd$Qs <- attr(xsstd, "std_rotation") %*% paramobj$Qs
  paramstd$Qe <- attr(xestd, "std_rotation") %*% paramobj$Qe
  paramstd$ce <- drop(t(paramobj$Qe) %*% attr(xestd, "std_center")) + paramobj$ce
  expect_equal(mnlink(xsstd, xestd, param = paramstd), ymean)
  # extra for ystd:
  paramstd$P <-  attr(ystd, "std_rotation") %*% paramobj$P
  expect_equal(mnlink(xsstd, xestd, param = paramstd), ystd, ignore_attr = TRUE)
  
  # standardise omega version
  om <- as_mnlink_Omega(paramobj)
  omstd <- om
  omstd$qs1 <- drop(attr(xsstd, "std_rotation") %*% om$qs1)
  omstd$qe1 <- drop(attr(xestd, "std_rotation") %*% om$qe1)
  omstd$ce1 <- drop(t(om$qe1) %*% attr(xestd, "std_center")) + om$ce1
  omstd$PBce <- drop(om$Omega[, seq.int(qs + 1, length.out = qe)] %*% attr(xestd, "std_center") + om$PBce)
  omstd$Omega[, seq.int(1, length.out = qs)] <- om$Omega[, seq.int(1, length.out = qs)] %*% t(attr(xsstd, "std_rotation"))
  omstd$Omega[, seq.int(qs + 1, length.out = qe)] <- om$Omega[, seq.int(qs + 1, length.out = qe)] %*% t(attr(xestd, "std_rotation"))
  # add std of P
  omstd$PBce <- drop(attr(ystd, "std_rotation") %*% omstd$PBce)
  omstd$Omega <- attr(ystd, "std_rotation") %*% omstd$Omega
  omstd$p1 <- drop(attr(ystd, "std_rotation") %*% omstd$p1)
  expect_equal(omstd, as_mnlink_Omega(paramstd), ignore_attr = TRUE)
  
  # given parameters for standardised data, solve for the corresponding parameters of the non-standard data
})
