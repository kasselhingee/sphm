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
  paramstd <- recoordinate_cann(paramobj, #yrot = diag(nrow(param$P)), 
                              xsrot = attr(xsstd, "std_rotation"),
                              xerot = attr(xestd, "std_rotation"), 
                              xecenter = attr(xestd, "std_center"))
  expect_equal(mnlink(xsstd, xestd, param = paramstd), ymean)
  # extra for ystd:
  paramstd <- recoordinate_cann(paramstd, yrot = attr(ystd, "std_rotation"))
  expect_equal(mnlink(xsstd, xestd, param = paramstd), ystd, ignore_attr = TRUE)
  
  # standardise omega version
  omstd <- recoordinate_Omega(as_mnlink_Omega(paramobj), 
                                yrot = attr(ystd, "std_rotation"), 
                                xsrot = attr(xsstd, "std_rotation"),
                                xerot = attr(xestd, "std_rotation"), 
                                xecenter = attr(xestd, "std_center"))
  expect_equal(omstd, as_mnlink_Omega(paramstd), ignore_attr = TRUE)
  
  # given parameters for standardised data, solve for the corresponding parameters of the non-standard data
})
