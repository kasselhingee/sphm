test_that("maximum likelihood for parallel axes per Jupp's path", {
  p <- 3
  q <- 5
  # data generating parameters:
  set.seed(1)
  P <- mclust::randomOrthogonalMatrix(p, p)
  set.seed(2)
  Q <- mclust::randomOrthogonalMatrix(q, p)
  set.seed(3)
  B <- diag(sort(runif(p-1), decreasing = TRUE))
  omegapar <- as_OmegaS2S(cannS2S(P,Q,B))
  
  #generate covariates uniformly on the sphere
  set.seed(4)
  x <- matrix(rnorm(1000*q), nrow = 1000)
  x <- sweep(x, 1, apply(x, 1, vnorm), FUN = "/")
  
  ymean <- meanlinkS2S(x = x, paramobj = omegapar)
 
  # generate noise
  # step 1: axes defined at P[, 1], orthogonal to P[, 1]
  set.seed(5)
  Kstar <- mclust::randomOrthogonalMatrix(p-1, p-1)
  Gstar <- getHstar(P[ ,1]) %*% Kstar
  # step 2: assume concentration and scales are constant
  k <- 30
  a <- c(1, seq(5, 0.2, length.out = p-1))
  a[-1] <- a[-1]/prod(a[-1])^(1/(p-1))
  #step 3: for each location define G as mean, Jupp transport of Gstar columns
  # then simulate from a SvMF, and evaluate density at the noise
  set.seed(6)
  y_ld <- t(apply(ymean, 1, function(mn){
    G <- cbind(mn, JuppRmat(P[,1], mn) %*% Gstar)
    obs <- rSvMF(1, SvMFcann(k, a, G))
    ld <- uldSvMF_cann(obs, k = k, a = a, G = G)
    return(c(obs, ld))
  }))
  
  # check ull_S2S_constV in C++
  ldCpp <- ull_S2S_constV_forR(y = y_ld[, 1:3], x = x, omvec = OmegaS2S_vec(omegapar), k = k,
                      a1 = a[1], aremaining = a[-1], Gstar = Gstar)
  expect_equal(ldCpp, y_ld[, 4])
                      
})


test_that("Cayley transform and vectorisation works", {
  p <- 3
  set.seed(100)
  M <- mclust::randomOrthogonalMatrix(p, p)
  A <- inverseCayleyTransform(M)
  expect_equal(A[upper.tri(A)], -A[lower.tri(A)])
  expect_equal(diag(A), rep(0, p))
  expect_equal(cayleyTransform(A), M)
  
  expect_equal(vectorizeLowerTriangle(A), A[lower.tri(A)])
  expect_equal(inverseVectorizeLowerTriangle(vectorizeLowerTriangle(A)), A)
})