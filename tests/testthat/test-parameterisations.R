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
  
  cann <- mnlink_cann(P, Bs = Bs, Qs = Qs, Be = Be, Qe = Qe, ce = ce)
  
  Om <- cann2Omega(cann, check = FALSE)
  
  # check properties of Om
  expect_equal(sum(diag(t(Om$Omega) %*% Om$Omega)), sum(Be^2 + Bs^2))
  expect_silent(mnlink_Omega_check(Om))
  
  # check that convert back matches
  cann2 <- Omega2cann(Om, check = FALSE)
  expect_silent(mnlink_cann_check(cann2))
  expect_equal(topos1strow(cann2$P), topos1strow(cann$P))
  # Qs and Qe
  expect_equal(topos1strow(cann2$Qs),topos1strow(cann$Qs))
  expect_equal(topos1strow(cann2$Qe),topos1strow(cann$Qe))
  
  #recovering Bs and Be
  expect_equal(cann2$Bs,cann$Bs)
  expect_equal(cann2$Be,cann$Be)
  expect_equal(cann2$ce,cann$ce)
})

test_that("Conversions work: Sph only", {
  set.seed(1)
  p <- 3
  P <- mclust::randomOrthogonalMatrix(p, p)
  qs <- 5
  set.seed(2)
  Qs <- mclust::randomOrthogonalMatrix(qs, p)
  set.seed(3)
  Bs <- diag(sort(runif(p-1), decreasing = TRUE))
  cann <- mnlink_cann(P, Bs = Bs, Qs = Qs)
  
  Om <- cann2Omega(cann, check = FALSE)
  
  # check properties of Om
  expect_equal(sum(diag(t(Om$Omega) %*% Om$Omega)), sum(Bs^2))
  expect_silent(mnlink_Omega_check(Om))
  
  # check that convert back matches
  cann2 <- Omega2cann(Om, check = FALSE)
  expect_silent(mnlink_cann_check(cann2))
  expect_equal(topos1strow(cann2$P), topos1strow(cann$P))
  # Qs
  expect_equal(topos1strow(cann2$Qs),topos1strow(cann$Qs))
  #recovering Bs and Be
  expect_equal(cann2$Bs,cann$Bs)
})

test_that("Conversions work: Euc only", {
  set.seed(1)
  p <- 3
  P <- mclust::randomOrthogonalMatrix(p, p)
  qe <- 5
  set.seed(2)
  Qe <- mclust::randomOrthogonalMatrix(qe, p)
  set.seed(3)
  Be <- diag(sort(runif(p-1), decreasing = TRUE))
  set.seed(4)
  ce <- runif(p)
  cann <- mnlink_cann(P, Be = Be, Qe = Qe, ce = ce)
  
  Om <- cann2Omega(cann, check = FALSE)
  
  # check properties of Om
  expect_equal(sum(diag(t(Om$Omega) %*% Om$Omega)), sum(Be^2))
  expect_silent(mnlink_Omega_check(Om))
  
  # check that convert back matches
  cann2 <- Omega2cann(Om, check = FALSE)
  expect_silent(mnlink_cann_check(cann2))
  expect_equal(topos1strow(cann2$P), topos1strow(cann$P))
  # Qe
  expect_equal(topos1strow(cann2$Qe),topos1strow(cann$Qe))
  #recovering Bs and Be
  expect_equal(cann2$Be,cann$Be)
  expect_equal(cann2$ce,cann$ce)
})


test_that("mnlink_Omega works directly", {
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
  cann <- mnlink_cann(P, Bs = Bs, Qs = Qs, Be = Be, Qe = Qe, ce = ce)
  Om <- cann2Omega(cann)
  
  Om2 <- mnlink_Omega(p1 = Om$p1, qs1 = Om$qs1, qe1 = Om$qe1, Omega = Om$Omega, ce = Om$ce, check = TRUE)
  expect_equal(Om2, Om)
})
  
test_that("projections orthogonal to p1 and qe1, qs1", {
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
  cann <- mnlink_cann(P, Bs = Bs, Qs = Qs, Be = Be, Qe = Qe, ce = ce)
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
  
test_that("vec and unvec - to finish", {
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
  cann <- mnlink_cann(P, Bs = Bs, Qs = Qs, Be = Be, Qe = Qe, ce = ce)
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


