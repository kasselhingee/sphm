test_that("rotatedresiduals() rotates residuals to the northpole along a geodesic correctly", {
  eps <- 0.1
  ypred <- rbind(c(0, 1, 0),
             c(0, 1, 0),
             c(0, 1, 0),
             c(0, 0, 1),
             c(0, 1/sqrt(2), 1/sqrt(2)))
  y <- rbind(c(0, 1, 0) + eps * c(1, 0, 0), #tangent to geodesic
                 c(0, 1, 0) + eps * c(0, 0, 1), #orthogonal to geodesic
                 c(0, 1, 0) + eps * c(1/sqrt(2), 0, 1/sqrt(2)), #pi/4 off from geodesic
                 c(0, 0, 1) + eps * c(1, 0, 0), #tangent to geodesic
                 c(0, 1/sqrt(2), 1/sqrt(2)) + eps * c(1, 0, 0)) #tangent to geodesic
  rresid <- rotatedresid(y = y, ypred = ypred, base = c(1, 0, 0))
  expect_equal(rresid[1, ], eps * c(0, -1, 0))
  expect_equal(rresid[2, ], eps * c(0, 0, 1))
  c(1/sqrt(2), 0, 1/sqrt(2)) %*% c(0, 0, 1)
  expect_equal(rresid[2, ] %*% rresid[3, ], eps^2 * c(1/sqrt(2), 0, 1/sqrt(2)) %*% c(0, 0, 1))
  expect_equal(rresid[4, ], eps * c(0, 0, -1))
  expect_equal(rresid[5, ], eps * c(0, -1, -1)/sqrt(2))
  
  #and actually Jupp's transport seems to give the reflected residual of geodesic transport
  rresid2 <- rotatedresid(y = y, ypred = ypred, base = c(1, 0, 0), path = "Jupp")
  expect_equal(rresid2, -rresid)
})

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
  Kstar[, 1] <- det(Kstar) * Kstar[,1]
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
                      a1 = a[1], aremaining = a[-1], Kstar = Kstar)
  expect_equal(ldCpp, y_ld[, 4])
  
  #check tape:
  S2S_constV_nota1_tovecparams(omvec = OmegaS2S_vec(omegapar), k = k,
                               aremaining = a[-1], Kstar = Kstar)
  ulltape <- tape_ull_S2S_constV_nota1(omvec = OmegaS2S_vec(omegapar), k = k,
                            a1 = a[1], aremaining = a[-1], Kstar = Kstar,
                            p, cbind(y_ld[, 1:3], x))
  expect_equal(ulltape$forward(0, ulltape$xtape), y_ld[, 4])
  
  exactll <- sum(ulltape$forward(0, ulltape$xtape))
  
  set.seed(7)
  Kstardifferent <- mclust::randomOrthogonalMatrix(p-1, p-1)
  Kstardifferent[, 1] <- det(Kstardifferent) * Kstardifferent[,1]
  badll <- sum(ulltape$forward(0, S2S_constV_nota1_tovecparams(omvec = OmegaS2S_vec(omegapar), k = k,
                               aremaining = a[-1], Kstar = Kstardifferent)))
  expect_lt(badll, exactll)
  
  # now try optimisation!
  omvec <- OmegaS2S_vec(omegapar)
  ll_mean_constraint <- tape_namedfun("wrap_OmegaS2S_constraints", omvec, vector(mode = "numeric"), p, matrix(nrow = 0, ncol = 0), check_for_nan = FALSE)
  est <- nloptr::nloptr(
    x0 = S2S_constV_nota1_tovecparams(omvec = omvec, k = k, aremaining = a[-1], Kstar = Kstardifferent),
    eval_f = function(theta){-sum(ulltape$eval(theta, a[1]))},
    eval_grad_f = function(theta){-colSums(matrix(ulltape$Jac(theta, a[1]), byrow = TRUE, ncol = length(theta)))},
    eval_g_eq =  function(theta){ll_mean_constraint$eval(theta[1:length(omvec)], vector(mode = "numeric"))},
    eval_jac_g_eq =  function(theta){
      cbind(matrix(ll_mean_constraint$Jac(theta[1:length(omvec)], vector(mode = "numeric")), byrow = TRUE, ncol = length(omvec)),
      matrix(0, nrow = 2, ncol = length(ulltape$xtape) - length(ll_mean_constraint$xtape)))
      },
    opts = list(algorithm = "NLOPT_LD_SLSQP", tol_constraints_eq = rep(1E-1, 2), xtol_rel = 1E-4, maxeval = 200, check_derivatives = TRUE, check_derivatives_print = "errors", print_level = 3)
  )
  cbind(est$solution, S2S_constV_nota1_tovecparams(omvec = omvec, k = k, aremaining = a[-1], Kstar = Kstar))
  expect_equal(est$solution, S2S_constV_nota1_tovecparams(omvec = omvec, k = k, aremaining = a[-1], Kstar = Kstar),
               tolerance = 1E-1)
  
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
