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
  expect_silent(cannS2S_check(obj))
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
  om <- as_OmegaS2S(cann)  
  expect_silent(OmegaS2S_check(om))
  expect_equal(as_OmegaS2S(as_cannS2S(om)), om)
  expect_equal(OmegaS2S(om$p1, om$q1, om$Omega), om)
  
  #vec and unvec
  expect_equal(OmegaS2S_unvec(OmegaS2S_vec(om), p), om)
  
  # project
  expect_equal(OmegaS2S_proj(om), om)
  ommod <- om
  ommod$p1 <- om$p1 * (1 + 1E-2 * runif(length(om$p1), -1, 1))
  ommod$q1 <- om$q1 * (1 + 1E-2 * runif(length(om$q1), -1, 1))
  expect_error(OmegaS2S_check(ommod), "checks failed")
  expect_equal(OmegaS2S_proj(ommod), om, tolerance = 1E-2)
  expect_equal(OmegaS2S_proj(ommod)$Omega, om$Omega)
})


