test_that("standardise spherical data gives data with correct moments", {
  rand_mobius_link_cann__place_in_env(3, 5, 4)
  
  #generate covariates Gaussianly
  set.seed(4)
  xe <- matrix(rnorm(100*qe), nrow = 100)
  #generate covariates on the sphere
  set.seed(4)
  xs <- matrix(rnorm(100*qs), nrow = 100)
  xs <- sweep(xs, 1, apply(xs, 1, vnorm), FUN = "/")
  
  ymean <- mobius_link(xs = xs, xe = xe, param = paramobj)
  
  # generate noise
  if (!requireNamespace("movMF", quietly = TRUE)){skip("Need movMF package")}
  set.seed(5)
  y <- t(apply(ymean, 1, function(mn){movMF::rmovMF(1, 2*mn)}))
  
  mat1 <- second_moment_mat(y)
  y2 <- standardise_sph(y, rotation = t(mat1))
  expect_equal(colMeans(y2)[-1], rep(0, p-1))
  # standardisation should do nothing up to sign:
  expect_equal(second_moment_mat(y2), diag(diag(second_moment_mat(y2))), ignore_attr = TRUE)
  
  # undo standardisation
  expect_equal(destandardise_sph(y2, t(mat1)), y)
})

test_that("standadise_Euc works seemlessly when covariates are all 0 or 1", {
  rand_mobius_link_cann__place_in_env(3, 5, 4, preseed = 3)
  
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
  xe <- cbind("const1" = 0, "const2" = 1, xe)
  xestd <- standardise_Euc(xe)
  expect_equal(xestd[, c(1, 2)], xe[, c(1,2)], ignore_attr = FALSE)
  expect_equal(xestd[, -c(1, 2)], standardise_Euc(xe[, -c(1, 2)]), ignore_attr = TRUE)
  xe2 <- destandardise_Euc(xestd, center = attr(xestd, "std_center"), rotation = attr(xestd, "std_rotation"))
  expect_equal(xe2, xe)
})

test_that("defaultstart() sets correct default starting parameters", {
  rand_mobius_link_cann__place_in_env(3, 5, 4)

  set.seed(4)
  xe <- matrix(rnorm(100 * qe), nrow = 100)
  set.seed(4)
  xs <- matrix(rnorm(100 * qs), nrow = 100)
  xs <- sweep(xs, 1, apply(xs, 1, vnorm), FUN = "/")
  set.seed(5)
  y <- matrix(rnorm(100 * p), nrow = 100)
  y <- sweep(y, 1, apply(y, 1, vnorm), FUN = "/")

  preplist <- list(y = y, xs = xs, xe = xe, start = NULL)
  preplist <- addEuccovars(preplist, type = "SpEuc", intercept = FALSE)
  preplist <- standardise_data(preplist, intercept = FALSE)

  # xe is not standardised when intercept = FALSE, so preplist$xe == xe
  expected_ce <- max(-preplist$xe) + 0.1 * IQR(preplist$xe)

  preplist <- defaultstart(preplist, type = "SpEuc", estimatescales = FALSE)

  start <- preplist$start
  mobius_link_cann_check(start)
  expect_equal(start$P,  diag(p),        ignore_attr = TRUE)
  expect_equal(start$Bs, diag(0.9, p-1), ignore_attr = TRUE)
  expect_equal(start$Qs, diag(1, qs, p), ignore_attr = TRUE)
  expect_equal(start$Be, diag(0.9, p-1), ignore_attr = TRUE)
  expect_equal(start$Qe, diag(1, qe, p), ignore_attr = TRUE)
  expect_equal(start$ce, expected_ce)
})

test_that("defaultstart() is a no-op when start is already set", {
  rand_mobius_link_cann__place_in_env(3, 5, 4)

  set.seed(4)
  xe <- matrix(rnorm(100 * qe), nrow = 100)
  set.seed(4)
  xs <- matrix(rnorm(100 * qs), nrow = 100)
  xs <- sweep(xs, 1, apply(xs, 1, vnorm), FUN = "/")
  set.seed(5)
  y <- matrix(rnorm(100 * p), nrow = 100)
  y <- sweep(y, 1, apply(y, 1, vnorm), FUN = "/")

  preplist <- list(y = y, xs = xs, xe = xe, start = paramobj)
  preplist <- addEuccovars(preplist, type = "SpEuc", intercept = FALSE)
  preplist <- standardise_data(preplist, intercept = FALSE)

  start_before <- preplist$start  # already recoordinated by standardise_data
  preplist <- defaultstart(preplist, type = "SpEuc")

  expect_identical(preplist$start, start_before)
})

test_that("defaultstart() with estimatescales=TRUE gives non-negative Bs/Be and consistent Qs/Qe", {
  rand_mobius_link_cann__place_in_env(3, 5, 4)

  set.seed(4)
  xe <- matrix(rnorm(100 * qe), nrow = 100)
  set.seed(4)
  xs <- matrix(rnorm(100 * qs), nrow = 100)
  xs <- sweep(xs, 1, apply(xs, 1, vnorm), FUN = "/")
  set.seed(5)
  y <- matrix(rnorm(100 * p), nrow = 100)
  y <- sweep(y, 1, apply(y, 1, vnorm), FUN = "/")

  preplist <- list(y = y, xs = xs, xe = xe, start = NULL)
  preplist <- addEuccovars(preplist, type = "SpEuc", intercept = FALSE)
  preplist <- standardise_data(preplist, intercept = FALSE)
  startce <- max(-preplist$xe) + 0.1 * IQR(preplist$xe)

  preplist <- defaultstart(preplist, type = "SpEuc", estimatescales = TRUE)

  start <- preplist$start
  mobius_link_cann_check(start)
  # Scale diagonals must be non-negative after sign absorption
  expect_true(all(diag(start$Bs) >= 0))
  expect_true(all(diag(start$Be) >= 0))
  # Re-run OLS with the sign-corrected Qs/Qe: all coefficients should now be non-negative
  olsY  <- Sp(preplist$y %*% t(solve(diag(p))))
  olsX1 <- Sp(preplist$xs %*% start$Qs)
  olsX2 <- (preplist$xe %*% start$Qe[, -1]) /
            drop(preplist$xe %*% start$Qe[, 1, drop = FALSE] + startce)
  check_fits <- lapply(seq_len(p - 1), function(i) lm(olsY[, i] ~ olsX1[, i] + olsX2[, i] - 1))
  expect_true(all(sapply(check_fits, function(f) coef(f)[1]) >= 0))
  expect_true(all(sapply(check_fits, function(f) coef(f)[2]) >= 0))
})

test_that("defaultstart() works with xs only (no xe)", {
  rand_mobius_link_cann__place_in_env(3, 5, 4)

  set.seed(4)
  xs <- matrix(rnorm(100 * qs), nrow = 100)
  xs <- sweep(xs, 1, apply(xs, 1, vnorm), FUN = "/")
  set.seed(5)
  y <- matrix(rnorm(100 * p), nrow = 100)
  y <- sweep(y, 1, apply(y, 1, vnorm), FUN = "/")

  preplist <- list(y = y, xs = xs, xe = NULL, start = NULL)
  preplist <- addEuccovars(preplist, type = "SpEuc", intercept = FALSE)
  preplist <- standardise_data(preplist, intercept = FALSE)
  preplist <- defaultstart(preplist, type = "SpEuc")

  start <- preplist$start
  mobius_link_cann_check(start)
  expect_true(all(diag(start$Bs) >= 0))
  expect_null(start$Be)
  expect_null(start$Qe)
})

test_that("defaultstart() works with xe only (no xs)", {
  rand_mobius_link_cann__place_in_env(3, 5, 4)

  set.seed(4)
  xe <- matrix(rnorm(100 * qe), nrow = 100)
  set.seed(5)
  y <- matrix(rnorm(100 * p), nrow = 100)
  y <- sweep(y, 1, apply(y, 1, vnorm), FUN = "/")

  preplist <- list(y = y, xs = NULL, xe = xe, start = NULL)
  preplist <- addEuccovars(preplist, type = "SpEuc", intercept = FALSE)
  preplist <- standardise_data(preplist, intercept = FALSE)
  preplist <- defaultstart(preplist, type = "SpEuc")

  start <- preplist$start
  mobius_link_cann_check(start)
  expect_true(all(diag(start$Be) >= 0))
  expect_null(start$Bs)
  expect_null(start$Qs)
})

test_that("recoordination of parameters work", {
  rand_mobius_link_cann__place_in_env(3, 5, 4)
  
  #generate covariates Gaussianly
  set.seed(4)
  xe <- cbind(matrix(rnorm(20*qe), nrow = 20))
  #generate covariates on the sphere
  set.seed(4)
  xs <- matrix(rnorm(20*qs), nrow = 20)
  xs <- sweep(xs, 1, apply(xs, 1, vnorm), FUN = "/")
  
  # add dummy covariates
  xe <- cbind(zeros = 0, ones = 1, xe)
  paramobj$Qe <- rbind(0, 0, paramobj$Qe)
  mobius_link_cann_check(paramobj)
  
  ymean <- mobius_link(xs = xs, xe = xe, param = paramobj)
  
  
  #standardise xs, xe and response and compute paramobj that includes the reverse
  xsstd <- standardise_sph(xs)
  xestd <- standardise_Euc(xe)
  ystd <- standardise_sph(ymean)
  
  # standardise omega version
  omstd <- recoordinate_Omega(as_mobius_link_Omega(paramobj), 
                                yrot = attr(ystd, "std_rotation"),
                                xsrot = attr(xsstd, "std_rotation"),
                                xerot = attr(xestd, "std_rotation"),
                                xecenter = attr(xestd, "std_center"),
                                onescovaridx = 2)
  expect_equal(mobius_link(xsstd,  xestd, param = omstd), ystd, ignore_attr = "std_rotation")
  
  # given parameters for standardised data, solve for the corresponding parameters of the non-standard data
  om2 <- undo_recoordinate_Omega(omstd,  
                          yrot = attr(ystd, "std_rotation"), 
                          xsrot = attr(xsstd, "std_rotation"),
                          xerot = attr(xestd, "std_rotation"), 
                          xecenter = attr(xestd, "std_center"),
                          onescovaridx = 2)
  expect_equal(om2, as_mobius_link_Omega(paramobj), ignore_attr = TRUE)
})

