test_that("cannS2S() creates objects that pass check", {
  set.seed(1)
  p <- 3
  q <- 5
  P <- mclust::randomOrthogonalMatrix(p, p)
  set.seed(2)
  Q <- mclust::randomOrthogonalMatrix(q, p)
  set.seed(3)
  B <- diag(sort(runif(p-1), decreasing = TRUE))
  obj <- cannS2S(P, Q, B)
  expect_silent(mnlink_cann_check(obj))
})

test_that("mnlink_cann() creates objects that pass check", {
  set.seed(1)
  p <- 3
  q <- 5
  P <- mclust::randomOrthogonalMatrix(p, p)
  set.seed(2)
  Qs <- mclust::randomOrthogonalMatrix(q, p)
  Qe <- mclust::randomOrthogonalMatrix(q+2, p)[,-1]
  set.seed(3)
  Bs <- diag(sort(runif(p-1), decreasing = TRUE))
  Be <- diag(sort(runif(p-1), decreasing = TRUE))
  obj <- mnlink_cann(P = P, Qs = Qs, Bs = Bs, Be = Be, Qe = Qe)
  expect_silent(mnlink_cann_check(obj))
})

test_that("OmegaS2S works and conversions", {
  set.seed(1)
  p <- 3
  q <- 5
  P <- mclust::randomOrthogonalMatrix(p, p)
  set.seed(2)
  Q <- mclust::randomOrthogonalMatrix(q, p)
  set.seed(3)
  B <- diag(sort(runif(p-1), decreasing = TRUE))
  cann <- cannS2S(P, Q, B)
  om <- as_mnlink_Omega(cann)  
  expect_silent(mnlink_Omega_check(om))
  expect_equal(as_mnlink_Omega(as_mnlink_cann(om)), om)
  expect_equal(OmegaS2S(om$p1, om$q1, om$Omega), om)
  expect_equal(sum(diag(t(obj$Omega) %*% obj$Omega)), sum(Be^2 + Bs^2))
  
  #vec and unvec
  expect_equal(OmegaS2S_unvec(OmegaS2S_vec(om), p), om)
  
  # project p1 and q1 perpendicular to Omega
  expect_equal(OmegaS2S_proj(om), om)
  ommod <- om
  ommod$p1 <- om$p1 * (1 + 1E-2 * runif(length(om$p1), -1, 1))
  ommod$q1 <- om$q1 * (1 + 1E-2 * runif(length(om$q1), -1, 1))
  expect_error(mnlink_Omega_check(ommod), "checks failed")
  expect_equal(OmegaS2S_proj(ommod, method = "p1q1"), om, tolerance = 1E-2)
  expect_equal(OmegaS2S_proj(ommod, method = "p1q1")$Omega, om$Omega)
  
  # project Omega perpendicular to p1 and q1
  expect_equal(OmegaS2S_proj(om, method = "Omega"), om)
  ommod <- om
  ommod$Omega <- om$Omega * (1 + 1E-2 * matrix(runif(length(om$Omega), -1, 1), p, q))
  expect_error(mnlink_Omega_check(ommod), "checks failed")
  expect_equal(OmegaS2S_proj(ommod, method = "Omega"), om, tolerance = 1E-2)
  expect_equal(OmegaS2S_proj(ommod, method = "Omega")[c("q1", "p1")], om[c("q1", "p1")])
})


