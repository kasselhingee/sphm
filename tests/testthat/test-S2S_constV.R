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
  rresid <- rotatedresid(y = y, ypred = ypred, base = c(1, 0, 0), path = "Absil")
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

test_that("rotationmat_amaral() rotates start to end", {
  vecs <- matrix(rnorm(6), nrow = 2)
  vecs <- vecs/sqrt(rowSums(vecs^2))
  expect_equal(drop(rotationmat_amaral(vecs[1, ], vecs[2,]) %*% vecs[1,]), vecs[2, ])
})

test_that("to and from vecparams functions correctly", {
  rmnlink_cann__place_in_env(4, 5, 4)
  omegapar <- as_mnlink_Omega(paramobj)
  k <- 30
  a <- c(1, seq(5, 0.2, length.out = p-1))
  a[-1] <- a[-1]/prod(a[-1])^(1/(p-1))
  
  set.seed(5)
  G0 <- mclust::randomOrthogonalMatrix(p, p) #axes around a random location
  G0[, p] <- det(G0) * G0[,p]
  
  set.seed(7)
  referencecoords <- mclust::randomOrthogonalMatrix(p, p)
  referencecoords[, p] <- det(referencecoords) * referencecoords[,p]
  
  # G01 identified with p1
  G01_p1 <- cbind(omegapar$p1, -JuppRmat(G0[,1], omegapar$p1) %*% G0[,-1])
  vecparams <- S2S_constV_nota1_tovecparams(omvec = mnlink_Omega_vec(omegapar), k = k,
                                            aremaining = a[-1],
                                            G0 = G01_p1,
                                            referencecoords = referencecoords,
                                            G01behaviour = "p1")
  expect_equal(S2S_constV_nota1_fromvecparamsR(vecparams, p, qs, qe, referencecoords, G01behaviour = "p1"),
               list(omvec = mnlink_Omega_vec(omegapar),
                    k = k,
                    aremaining = a[-1],
                    G0 = G01_p1), ignore_attr = TRUE)
  veclength <- length(vecparams)
  
  # G01 fixed
  vecparams <- S2S_constV_nota1_tovecparams(omvec = mnlink_Omega_vec(omegapar), k = k,
                                            aremaining = a[-1],
                                            G0 = G0,
                                            referencecoords = referencecoords,
                                            G01behaviour = "fixed")
  expect_equal(length(vecparams), veclength)
  expect_equal(S2S_constV_nota1_fromvecparamsR(vecparams, p, qs, qe, referencecoords, G01behaviour = "fixed", G01 = G0[,1]),
               list(omvec = mnlink_Omega_vec(omegapar),
                    k = k,
                    aremaining = a[-1],
                    G0 = G0), ignore_attr = TRUE)
  
  # G01 free
  vecparams <- S2S_constV_nota1_tovecparams(omvec = mnlink_Omega_vec(omegapar), k = k,
                                            aremaining = a[-1],
                                            G0 = G0,
                                            referencecoords = referencecoords,
                                            G01behaviour = "free")
  expect_equal(length(vecparams), veclength + p-1)
  expect_equal(S2S_constV_nota1_fromvecparamsR(vecparams, p, qs, qe, referencecoords, G01behaviour = "fixed", G01 = G0[,1]),
               list(omvec = mnlink_Omega_vec(omegapar),
                    k = k,
                    aremaining = a[-1],
                    G0 = G0), ignore_attr = TRUE)
  
})

test_that("maximum likelihood for parallel axes per geodesic path", {
  rmnlink_cann__place_in_env(4, 5, 4)
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
  # step 1: SvMF axes and reference location defined at P[, 1], orthogonal to P[, 1]
  set.seed(5)
  G0 <- mclust::randomOrthogonalMatrix(p, p) #axes around a random location
  G0[, p] <- det(G0) * G0[,p]
  G0star <- G0[,-1]
  # step 2: assume concentration and scales are constant
  k <- 30
  a <- c(1, seq(5, 0.2, length.out = p-1))
  a[-1] <- a[-1]/prod(a[-1])^(1/(p-1))
  #step 3: for each location define G as mean and geodesic transport of Gstar columns
  # then simulate from a SvMF, and evaluate density at the noise
  set.seed(6)
  y_ld <- t(apply(ymean, 1, function(mn){
    G <- cbind(mn, JuppRmat(G0[,1], mn) %*% G0star)
    obs <- rSvMF(1, SvMFcann(k, a, G))
    ld <- uldSvMF_cann(obs, k = k, a = a, G = G)
    return(c(obs, ld))
  }))
  
  # check ull_S2S_constV in C++
  ldCpp <- ull_S2S_constV_forR(y = y_ld[, 1:p], xs = xs, xe = xe, omvec = mnlink_Omega_vec(omegapar), k = k,
                      a1 = a[1], aremaining = a[-1], G0 = G0)
  expect_equal(ldCpp, y_ld[, p+1])
  
  # check vectorisation and reverse
  set.seed(7)
  referencecoords <- mclust::randomOrthogonalMatrix(p, p)
  vecparams <- S2S_constV_nota1_tovecparams(omvec = mnlink_Omega_vec(omegapar), k = k,
                               aremaining = a[-1], G0star = G0star, referencecoords = referencecoords)
  
  #check tape:
  expect_warning({ulltape <- tape_ull_S2S_constV_nota1(omvec = mnlink_Omega_vec(omegapar), k = k,
                            a1 = a[1], aremaining = a[-1], G0star = G0star,
                            p, qe,
                            yx = cbind(y_ld[, 1:p], xs, xe), 
                            referencecoords)}, "p!=3")
  expect_equal(ulltape$xtape, vecparams)
  expect_equal(ulltape$forward(0, ulltape$xtape), y_ld[, p+1])
  
  exactll <- sum(ulltape$forward(0, ulltape$xtape))
  
  set.seed(7)
  Kstardifferent <- mclust::randomOrthogonalMatrix(p-1, p-1)
  Kstardifferent[, 1] <- det(Kstardifferent) * Kstardifferent[,1]
  G0stardifferent <- t(rotationmat_amaral(omegapar$p1, nthpole(p))) %*% rbind(0, Kstardifferent)
  badll <- sum(ulltape$forward(0, 
                               S2S_constV_nota1_tovecparams(omvec = mnlink_Omega_vec(omegapar), k = k,
                               aremaining = a[-1], G0star = G0stardifferent, referencecoords = referencecoords)))
  expect_lt(badll, exactll)
  
  ## now try optimisation starting at true values ##
  expect_warning({est1 <- optim_constV(y_ld[, 1:p], xs, xe, omegapar, k, a, G0star, xtol_rel = 1E-4, G0reference = referencecoords)}, "p!=3")
  expect_equal(est1$solution$mean, omegapar, tolerance = 1E-1)
  expect_equal(est1$solution[c("k", "a")], list(k = k, a = a), tolerance = 1E-1)
  # check Gstar by checking angle between estimated and true axes
  axis_distance <- function(angle1, angle2 = 0){
    diff <- abs(angle1 - angle2)
    pmin(diff, pi - diff)
  }
  expect_equal(axis_distance(acos(colSums(est1$solution$G0star * G0star))), rep(0, ncol(G0star)), tolerance = 1E-1, ignore_attr = TRUE)
  
  ## now starting optimisation away from starting parameters ##
  bad_om <- as_mnlink_Omega(rmnlink_cann(p, qs, qe, preseed = 2))
  set.seed(13)
  pre <- prelim_ad(y_ld[, 1:p], xs, xe, bad_om, xtol_rel = 1E-4) #doing this preliminary estimate reduces the iterations needed by optim_constV
  badG0star <- getHstar(pre$solution$p1) %*% mclust::randomOrthogonalMatrix(p-1, p-1)
  expect_warning({est2 <- optim_constV(y_ld[, 1:p], xs, xe, pre$solution, k = 10, a = rep(1, p), G0star = -badG0star, xtol_rel = 1E-4, G0reference = referencecoords)}, "p!=3")
  expect_equal(est2$solution, est1$solution, tolerance = 1E-3)
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
