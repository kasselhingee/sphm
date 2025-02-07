test_that("mnlink_cann() objects pass check: Sph only", {
  set.seed(1)
  p <- 3
  qs <- 5
  P <- mclust::randomOrthogonalMatrix(p, p)
  set.seed(2)
  Qs <- mclust::randomOrthogonalMatrix(qs, p)
  set.seed(3)
  Bs <- diag(sort(runif(p-1), decreasing = TRUE))
  obj <- mnlink_cann(P, Bs = Bs, Qs = Qs, check = FALSE)
  expect_silent(mnlink_cann_check(obj))
})

test_that("mnlink_cann() objects pass check: Euc only", {
  set.seed(1)
  p <- 3
  qe <- 5
  P <- mclust::randomOrthogonalMatrix(p, p)
  set.seed(2)
  Qe <- mclust::randomOrthogonalMatrix(qe, p)
  set.seed(3)
  Be <- diag(sort(runif(p-1), decreasing = TRUE))
  set.seed(4)
  ce <- runif(p)
  obj <- mnlink_cann(P, Be = Be, Qe = Qe, ce = ce, check = FALSE)
  expect_silent(mnlink_cann_check(obj))
})

test_that("mnlink_cann() objects pass check: Sph + Euc only", {
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
  obj <- mnlink_cann(P, Bs = Bs, Qs = Qs, Be = Be, Qe = Qe, ce = ce, check = FALSE)
  expect_silent(mnlink_cann_check(obj))
})

test_that("mnlink_cann(): common mistakes", {
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
  
  expect_error(mnlink_cann(P, Bs = Bs, Qs = Qs, Be = Be, Qe = Qe, ce = NULL, check = TRUE))
  expect_error(mnlink_cann(P, Bs = Bs, Qs = Qs, Be = Be, Qe = Qe, ce = ce[-1], check = TRUE))
  expect_error(mnlink_cann(P, Bs = Bs, Qe = Qs, check = TRUE))
  expect_error(mnlink_cann(P, Be = Be, Qs = Qe, check = TRUE))
})

test_that("Conversions work: Sph + Euc", {
  rmnlink_cann__place_in_env(3, 5, 4)
  cann <- paramobj
  
  Om <- cann2Omega(cann, check = FALSE)
  

  
  # check that convert back matches
  cann2 <- Omega2cann(Om, check = FALSE)
  expect_silent(mnlink_cann_check(cann2))
  cann2 <- P_signswitch(cann2, sign(cann2$P[1, ]) != sign(cann$P[1, ]))
  expect_equal(cann2, cann)
  
  # check properties of Om
  expect_equal(sum(diag(t(Om$Omega) %*% Om$Omega)), sum(Be^2 + Bs^2))
  expect_silent(mnlink_Omega_check(Om))
  # further constraints
  Is_tilde <- diag(1, qs + qe, qs)
  Ie_tilde <- matrix(0, qs + qe, qe)
  Ie_tilde[qs + (1:qe), ] <- diag(1, qe)
  expect_equal(Om$Omega %*% Is_tilde, Om$Omega[, 1:qs])
  eigfull <- eigen(Om$Omega %*% t(Om$Omega))
  eigs <- eigen(Om$Omega %*% (Is_tilde %*% t(Is_tilde)) %*% t(Om$Omega))
  eige <- eigen(Om$Omega %*% (Ie_tilde %*% t(Ie_tilde)) %*% t(Om$Omega))
  expect_equal(eigfull$vectors, eigs$vectors)
  expect_equal(eigfull$vectors, eige$vectors)

  # commutivity doesn't apply to general matrices:
  dummy <- matrix(rnorm(p * (qs + qe)), p, qs + qe)
  Om$Omega <- dummy
  expect_true(all(mnlink_Omega_check_numerical(Om)[c("Omega_comm_s", "Omega_comm_e")] > rep(sqrt(.Machine$double.eps), 2)))
})

test_that("Conversions work: Sph only", {
  rmnlink_cann__place_in_env(3, 5, 0)
  cann <- paramobj
  
  Om <- cann2Omega(cann, check = FALSE)
  
  # check properties of Om
  expect_equal(sum(diag(t(Om$Omega) %*% Om$Omega)), sum(Bs^2))
  expect_silent(mnlink_Omega_check(Om))
  
  # check that convert back matches
  cann2 <- Omega2cann(Om, check = FALSE)
  expect_silent(mnlink_cann_check(cann2))
  cann2 <- P_signswitch(cann2, sign(cann2$P[1, ]) != sign(cann$P[1, ]))
  expect_equal(cann2, cann)
})

test_that("Conversions work: Euc only", {
  rmnlink_cann__place_in_env(3, 0, 4)
  cann <- paramobj
  
  Om <- cann2Omega(cann, check = FALSE)
  
  # check properties of Om
  expect_equal(sum(diag(t(Om$Omega) %*% Om$Omega)), sum(Be^2))
  expect_silent(mnlink_Omega_check(Om))
  
  # check that convert back matches
  cann2 <- Omega2cann(Om, check = FALSE)
  expect_silent(mnlink_cann_check(cann2))
  cann2 <- P_signswitch(cann2, sign(cann2$P[1, ]) != sign(cann$P[1, ]))
  expect_equal(cann2, cann)
})


test_that("mnlink_Omega works directly", {
  rmnlink_cann__place_in_env(3, 5, 4)
  cann <- paramobj
  Om <- cann2Omega(cann)
  
  Om2 <- mnlink_Omega(p1 = Om$p1, qs1 = Om$qs1, qe1 = Om$qe1, Omega = Om$Omega, ce1 = Om$ce1, PBce = Om$PBce, check = TRUE)
  expect_equal(Om2, Om)
})
  
test_that("projections orthogonal to p1 and qe1, qs1", {
  rmnlink_cann__place_in_env(3, 5, 4)
  cann <- paramobj
  Om <- cann2Omega(cann)
  
  # project Omega perpendicular to p1 and q1
  expect_equal(Omega_proj(Om), Om)
  Ommod <- Om
  Ommod$Omega <- Om$Omega * (1 + 1E-2 * matrix(runif(length(Om$Omega), -1, 1), p, qs + qe))
  expect_error(mnlink_Omega_check(Ommod), "checks failed")
  expect_equal(Omega_proj(Ommod)$Omega, Om$Omega, tolerance = 1E-2)
  expect_equal(Omega_proj(Ommod)[c("qs1", "qe1", "p1")], Om[c("qs1", "qe1", "p1")])
  
  # just Sph
  cann <- mnlink_cann(P, Bs = Bs, Qs = Qs)
  Om <- cann2Omega(cann)
  expect_equal(Omega_proj(Om), Om)
  Ommod <- Om
  Ommod$Omega <- Om$Omega * (1 + 1E-2 * matrix(runif(length(Om$Omega), -1, 1), p, qs))
  expect_error(mnlink_Omega_check(Ommod), "checks failed")
  expect_equal(Omega_proj(Ommod)$Omega, Om$Omega, tolerance = 1E-2)
  expect_equal(Omega_proj(Ommod)[c("qs1", "qe1", "p1")], Om[c("qs1", "qe1", "p1")])
  
  # just Euc
  cann <- mnlink_cann(P, Be = Be, Qe = Qe, ce = ce)
  Om <- cann2Omega(cann)
  expect_equal(Omega_proj(Om), Om)
  Ommod <- Om
  Ommod$Omega <- Om$Omega * (1 + 1E-2 * matrix(runif(length(Om$Omega), -1, 1), p, qe))
  expect_error(mnlink_Omega_check(Ommod), "checks failed")
  expect_equal(Omega_proj(Ommod)$Omega, Om$Omega, tolerance = 1E-2)
  expect_equal(Omega_proj(Ommod)[c("qs1", "qe1", "p1")], Om[c("qs1", "qe1", "p1")])
}) 
  
test_that("vec and unvec", {
  rmnlink_cann__place_in_env(3, 5, 4)
  cann <- paramobj
  Om <- cann2Omega(cann)
  
  #vec and unvec
  expect_equal(mnlink_Omega_unvec(mnlink_Omega_vec(Om), p, qe = qe), Om)
  
  #sph only
  Om <- cann2Omega(mnlink_cann(P, Bs = Bs, Qs = Qs))
  expect_equal(mnlink_Omega_unvec(mnlink_Omega_vec(Om), p, qe = 0), Om)
  
  #Euc only
  Om <- cann2Omega(mnlink_cann(P, Be = Be, Qe = Qe, ce = ce))
  expect_equal(mnlink_Omega_unvec(mnlink_Omega_vec(Om), p, qe = qe), Om)
})

test_that("check Shogo conversion", {
  rmnlink_cann__place_in_env(3, 0, 4)
  # convert to Shogo form:
  bigQe <- rbind(0, Qe)
  bigQe[, 1] <- 0
  bigQe[1,1] <- 1
  bigce <- ce
  bigce[1] <- 1
  ce <- ce[-1]
  paramobj <- mnlink_cann(P, Be = Be, Qe = bigQe, ce = bigce)
  expect_true(is_Shogo(paramobj))
  
  # check manual
  x <- runif(qe)
  expect_equal(drop(mnlink(xe = matrix(c(0, x), nrow = 1), param = paramobj)),
               drop(P %*% iSp(drop(Be %*% (t(bigQe[-1,-1]) %*% x + ce)))))
  
  # check Shogo:
  # need to check linearity of P^{-1} S(y) wrt x_e and ce or even Qe
  # expect differences between the intermediate output to be purely due to the difference between inputs
  out <- Sp(mnlink(xe = rbind(0, c(0,x), c(0, 2*x)), param = paramobj) %*% P)
  expect_equal(out[2, ] - out[1, ], out[3, ] - out[2, ])
})

