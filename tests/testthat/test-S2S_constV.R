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

test_that("maximum likelihood for parallel axes per geodesic path", {
  rmnlink_cann__place_in_env(3, 5, 4)
  omegapar <- as_mnlink_Omega(paramobj)
  
  #generate covariates Gaussianly
  set.seed(4)
  xe <- matrix(rnorm(1000*qe), nrow = 1000)
  #generate covariates on the sphere
  set.seed(4)
  xs <- matrix(rnorm(1000*qs), nrow = 1000)
  xs <- sweep(xs, 1, apply(xs, 1, vnorm), FUN = "/")
  
  ymean <- mnlink(xs = xs, xe = xe, param = paramobj)
 
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
  #step 3: for each location define G as mean and geodesic transport of Gstar columns
  # then simulate from a SvMF, and evaluate density at the noise
  set.seed(6)
  y_ld <- t(apply(ymean, 1, function(mn){
    G <- cbind(mn, rotationmat_amaral(P[,1], mn) %*% Gstar)
    obs <- rSvMF(1, SvMFcann(k, a, G))
    ld <- uldSvMF_cann(obs, k = k, a = a, G = G)
    return(c(obs, ld))
  }))
  
  # check ull_S2S_constV in C++
  ldCpp <- ull_S2S_constV_forR(y = y_ld[, 1:p], xs = xs, xe = xe, omvec = mnlink_Omega_vec(omegapar), k = k,
                      a1 = a[1], aremaining = a[-1], Kstar = Kstar)
  expect_equal(ldCpp, y_ld[, p+1])
  
  # check vectorisation and reverse
  vecparams <- S2S_constV_nota1_tovecparams(omvec = mnlink_Omega_vec(omegapar), k = k,
                               aremaining = a[-1], Kstar = Kstar)
  expect_equal(S2S_constV_nota1_fromvecparamsR(vecparams, p, qs, qe),
               list(omvec = mnlink_Omega_vec(omegapar),
                    k = k,
                    aremaining = a[-1],
                    Kstar = Kstar), ignore_attr = TRUE)
  
  #check tape:
  ulltape <- tape_ull_S2S_constV_nota1(omvec = mnlink_Omega_vec(omegapar), k = k,
                            a1 = a[1], aremaining = a[-1], Kstar = Kstar,
                            yx = cbind(y_ld[, 1:p], xs, xe)[1, ],
                            p, qe)
  expect_equal(ulltape$forward(0, ulltape$xtape), y_ld[1, 4])
  
  exactll <- sum(scorematchingad::evaltape(ulltape, ulltape$xtape, cbind(y_ld[, 1:p], xs, xe)))
  
  set.seed(7)
  Kstardifferent <- mclust::randomOrthogonalMatrix(p-1, p-1)
  Kstardifferent[, 1] <- det(Kstardifferent) * Kstardifferent[,1]
  badll <- sum(scorematchingad::evaltape(ulltape, S2S_constV_nota1_tovecparams(omvec = mnlink_Omega_vec(omegapar), k = k,
                               aremaining = a[-1], Kstar = Kstardifferent), cbind(y_ld[, 1:p], xs, xe)))
  expect_lt(badll, exactll)
  
  ## now try optimisation starting at true values ##
  est1 <- optim_constV(y_ld[, 1:p], x, omegapar, k, a, Gstar, xtol_rel = 1E-4)
  expect_equal(est1$solution$mean, omegapar, tolerance = 1E-1)
  expect_equal(est1$solution[c("k", "a")], list(k = k, a = a), tolerance = 1E-1, ignore_attr = TRUE)
  # check angle between estimated and true axes
    axis_distance <- function(angle1, angle2 = 0){
    diff <- abs(angle1 - angle2)
    pmin(diff, pi - diff)
  }
  expect_equal(axis_distance(acos(colSums(est1$solution$Gstar * Gstar))), rep(0, ncol(Gstar)), tolerance = 1E-2, ignore_attr = TRUE)
  
  ## now starting optimisation with bad Kstar ##
  set.seed(14)
  start <- as_mnlink_Omega(cannS2S(P = mclust::randomOrthogonalMatrix(p, p),
                               Q = mclust::randomOrthogonalMatrix(q, p),
                               B = diag(sort(runif(p-1), decreasing = TRUE))))
  set.seed(14)
  badGstar <- getHstar(omegapar$p1) %*% mclust::randomOrthogonalMatrix(p-1, p-1)[, c(2, 1)]
  est2 <- optim_constV(y_ld[, 1:p], x, omegapar, k = k, a = a, badGstar, xtol_rel = 1E-4)
  expect_equal(as_mnlink_cann(est2$solution$mean), as_mnlink_cann(omegapar), tolerance = 1E-1)
  #estimation of the concentration and scaling is quite bad without the flip of [, c(2,1)] for creating badGstar! 
  expect_equal(est2$solution[c("k", "a")], list(k = k, a = a), tolerance = 1E-1, ignore_attr = TRUE)
  # check angle between estimated and true axes
  expect_equal(axis_distance(acos(colSums(est2$solution$Gstar * Gstar))), rep(0, ncol(Gstar)), tolerance = pi/10, ignore_attr = TRUE)
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
