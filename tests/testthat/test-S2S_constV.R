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
                    G0 = cbind(G01_p1[,1],-G01_p1[,-1])), ignore_attr = TRUE)
  veclength <- length(vecparams)
  vecparams_p1 <- vecparams
  
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
                    G0 = cbind(G0[,1],-G0[,-1])), ignore_attr = TRUE)
  vecparams_fixed <- vecparams
  
  # G01 free
  vecparams <- S2S_constV_nota1_tovecparams(omvec = mnlink_Omega_vec(omegapar), k = k,
                                            aremaining = a[-1],
                                            G0 = G0,
                                            referencecoords = referencecoords,
                                            G01behaviour = "free")
  expect_equal(length(vecparams), veclength + p-1)
  expect_equal(S2S_constV_nota1_fromvecparamsR(vecparams, p, qs, qe, referencecoords, G01behaviour = "free"),
               list(omvec = mnlink_Omega_vec(omegapar),
                    k = k,
                    aremaining = a[-1],
                    G0 = G0), ignore_attr = TRUE)
  vecparams_free <- vecparams
  
  incommon <- 1:length((mnlink_Omega_vec(omegapar) + 1 + p-1))
  expect_equal(vecparams_fixed[incommon], vecparams_p1[incommon])
  expect_equal(vecparams_free[incommon], vecparams_p1[incommon])
})

# function for computing axial differences (ignoring sign)
axis_distance <- function(angle1, angle2 = 0){
  diff <- abs(angle1 - angle2)
  pmin(diff, pi - diff)
}
rcovars <- function(n, qs, qe){
  #generate covariates Gaussianly
  xe <- matrix(rnorm(n*qe), nrow = n)
  #generate covariates on the sphere
  xs <- matrix(rnorm(n*qs), nrow = n)
  xs <- sweep(xs, 1, apply(xs, 1, vnorm), FUN = "/")
  return(list(xe = xe, xs = xs))
}

test_that("MLE with p1 = G01", {
  rmnlink_cann__place_in_env(4, 5, 4)
  omegapar <- as_mnlink_Omega(paramobj)
  set.seed(4)
  x <- rcovars(1000, qs, qe)
  set.seed(5)
  G0_other <- mclust::randomOrthogonalMatrix(p, p) #axes around a random location
  G0_other[, p] <- det(G0_other) * G0_other[,p]
  G0 <- cbind(omegapar$p1, -JuppRmat(G0_other[,1], omegapar$p1) %*% G0_other[,-1])
  
  # simulate observations
  set.seed(7)
  k <- 30
  a <- c(1, seq(5, 0.2, length.out = p-1))
  a[-1] <- a[-1]/prod(a[-1])^(1/(p-1))
  y_ld <- rMobius_SvMF(x$xs, x$xe, mnparam = omegapar, k, a, G0)

  # check ull_S2S_constV in C++
  ldCpp <- ull_S2S_constV_forR(y = y_ld[, 1:p], xs = x$xs, xe = x$xe, omvec = mnlink_Omega_vec(omegapar), k = k,
                               a1 = a[1], aremaining = a[-1], G0 = G0)
  expect_equal(ldCpp, y_ld[, p+1])
  
  
  #check tape:
  set.seed(7)
  referencecoords <- mclust::randomOrthogonalMatrix(p, p)
  referencecoords[, p] <- det(referencecoords) * referencecoords[,p]
  vecparams <- S2S_constV_nota1_tovecparams(omvec = mnlink_Omega_vec(omegapar), k = k,
                                            aremaining = a[-1],
                                            G0 = G0, 
                                            referencecoords = referencecoords,
                                            G01behaviour = "p1")
  expect_warning({ulltape <- tape_ull_S2S_constV_nota1(omvec = mnlink_Omega_vec(omegapar), k = k,
                                                       a1 = a[1], aremaining = a[-1],
                                                       G0 = G0,
                                                       p, qe,
                                                       yx = cbind(y_ld[, 1:p], x$xs, x$xe), 
                                                       referencecoords,
                                                       G01behaviour = "p1")}, "p!=3")
  expect_equal(ulltape$xtape, vecparams)
  expect_equal(ulltape$forward(0, ulltape$xtape), y_ld[, p+1])
  
  exactll <- sum(ulltape$forward(0, ulltape$xtape))
  
  set.seed(18)
  G0_other <- mclust::randomOrthogonalMatrix(p, p) #axes around a random location
  G0_other[, p] <- det(G0_other) * G0_other[,p]
  G0stardifferent <- -JuppRmat(G0_other[,1], omegapar$p1) %*% G0_other[,-1]
  badll <- sum(ulltape$forward(0, 
                               S2S_constV_nota1_tovecparams(omvec = mnlink_Omega_vec(omegapar), k = k,
                                                            aremaining = a[-1], 
                                                            G0 = cbind(omegapar$p1, G0stardifferent),
                                                            referencecoords = referencecoords,
                                                            G01behaviour = "p1")))
  expect_lt(badll, exactll)
  
  ## now try optimisation starting at true values ##
  expect_warning({est1 <- optim_constV(y_ld[, 1:p], x$xs, x$xe, omegapar, k, a, G0, xtol_rel = 1E-4, G0reference = referencecoords, G01behaviour = "p1", 
                                       type = "Kassel", intercept = FALSE)}, "(p!=3|multimodal)")
  expect_equal(est1$mean, omegapar, tolerance = 1E-1)
  expect_equal(est1[c("k", "a")], list(k = k, a = a), tolerance = 1E-1)
  # check Gstar by checking angle between estimated and true axes
  expect_equal(axis_distance(acos(colSums(est1$G0 * G0))), rep(0, p), tolerance = 1E-1, ignore_attr = TRUE)
  # residuals should be close to MVN with covariance given by a
  # But the difference can still be detected at k = 30 when n = 1000 (see unit test of resid_SvMF_partransport below to see how large k might need to be)
  
  ## now starting optimisation away from starting parameters ##
  bad_om <- as_mnlink_Omega(rmnlink_cann(p, qs, qe, preseed = 2))
  preest <- mobius_SvMF_partransport_prelim(y_ld[, 1:p], x$xs, x$xe, mean = bad_om, type = "Kassel", G01behaviour = "p1", intercept = FALSE)
  expect_warning({est2 <- optim_constV(y_ld[, 1:p], x$xs, x$xe, 
                                       mean = preest$mean, k = preest$k, a = preest$a, G0 = preest$G0,
                                       G01behaviour = "p1",
                                       type = "Kassel", intercept = FALSE)}, "(p!=3|multimodal)")
  expect_equal(est2$mean, est1$mean, tolerance = 1E-1)
  expect_equal(est2$k, est1$k, tolerance = 0.2)
  expect_equal(est2$a, est1$a, tolerance = 1E-1, ignore_attr = "names")
  expect_equal(axis_distance(acos(colSums(est2$G0 * est1$G0))), rep(0, p), tolerance = 1E-1, ignore_attr = TRUE)
  
  # from default starts with prelim estimate
  expect_warning({est3 <- mobius_SvMF(y_ld[, 1:p], x$xs, x$xe,
                       G01behaviour = "p1",
                       type = "Kassel", intercept = FALSE)}, "(p!=3|multimodal)")
  expect_equal(est3$mean, est1$mean, tolerance = 1E-1)
  expect_equal(est3$k, est1$k, tolerance = 0.2)
  expect_equal(est3$a, est1$a, tolerance = 1E-1)
  expect_equal(axis_distance(acos(colSums(est3$G0 * est1$G0))), rep(0, p), tolerance = 1E-1, ignore_attr = TRUE)
})

test_that("MLE with G01 fixed", {
  rmnlink_cann__place_in_env(4, 5, 4)
  omegapar <- as_mnlink_Omega(paramobj)
  set.seed(4)
  x <- rcovars(1000, qs, qe)
  set.seed(5)
  G0 <- mclust::randomOrthogonalMatrix(p, p) #axes around a random location
  G0[, p] <- det(G0) * G0[,p]
  
  # simulate observations
  set.seed(6)
  k <- 30
  a <- c(1, seq(5, 0.2, length.out = p-1))
  a[-1] <- a[-1]/prod(a[-1])^(1/(p-1))
  y_ld <- rMobius_SvMF(x$xs, x$xe, mnparam = omegapar, k, a, G0)
  
  # check ull_S2S_constV in C++
  ldCpp <- ull_S2S_constV_forR(y = y_ld[, 1:p], xs = x$xs, xe = x$xe, omvec = mnlink_Omega_vec(omegapar), k = k,
                               a1 = a[1], aremaining = a[-1], G0 = G0)
  expect_equal(ldCpp, y_ld[, p+1])
  
  
  #check tape:
  set.seed(7)
  referencecoords <- mclust::randomOrthogonalMatrix(p, p)
  referencecoords[, p] <- det(referencecoords) * referencecoords[,p]
  vecparams <- S2S_constV_nota1_tovecparams(omvec = mnlink_Omega_vec(omegapar), k = k,
                                            aremaining = a[-1],
                                            G0 = G0, 
                                            referencecoords = referencecoords,
                                            G01behaviour = "fixed")
  expect_warning({ulltape <- tape_ull_S2S_constV_nota1(omvec = mnlink_Omega_vec(omegapar), k = k,
                                                       a1 = a[1], aremaining = a[-1],
                                                       G0 = G0,
                                                       p, qe,
                                                       yx = cbind(y_ld[, 1:p], x$xs, x$xe), 
                                                       referencecoords,
                                                       G01behaviour = "fixed")}, "p!=3")
  expect_equal(ulltape$xtape, vecparams)
  expect_equal(ulltape$forward(0, ulltape$xtape), y_ld[, p+1])
  
  exactll <- sum(ulltape$forward(0, ulltape$xtape))
  
  set.seed(18)
  G0_other <- mclust::randomOrthogonalMatrix(p, p) #axes around a random location
  G0_other[, p] <- det(G0_other) * G0_other[,p]
  G0stardifferent <- -JuppRmat(G0_other[,1], G0[,1]) %*% G0_other[,-1]
  badll <- sum(ulltape$forward(0, 
                               S2S_constV_nota1_tovecparams(omvec = mnlink_Omega_vec(omegapar), k = k,
                                                            aremaining = a[-1], 
                                                            G0 = cbind(G0[,1], G0stardifferent),
                                                            referencecoords = referencecoords,
                                                            G01behaviour = "fixed")))
  expect_lt(badll, exactll)
  
  ## now try optimisation starting at true values ##
  expect_warning({est1 <- optim_constV(y_ld[, 1:p], x$xs, x$xe, omegapar, k, a, G0, xtol_rel = 1E-4, G0reference = referencecoords, G01behaviour = "fixed",
                                       type = "Kassel", intercept = FALSE)}, "(p!=3|multimodal)")
  expect_equal(est1$mean, omegapar, tolerance = 1E-1)
  expect_equal(est1[c("k", "a")], list(k = k, a = a), tolerance = 1E-1)
  # check Gstar by checking angle between estimated and true axes
  expect_equal(axis_distance(acos(colSums(est1$G0 * G0)-1E-15)), rep(0, p), tolerance = 1E-1, ignore_attr = TRUE)
  
  ## now starting optimisation away from starting parameters ##
  bad_om <- as_mnlink_Omega(rmnlink_cann(p, qs, qe, preseed = 2))
  preest <- mobius_SvMF_partransport_prelim(y_ld[, 1:p], x$xs, x$xe, mean = bad_om, type = "Kassel", G0 = cbind(G0[,1], matrix(NA, p, p-1)), G01behaviour = "fixed", intercept = FALSE)
  expect_warning({est2 <- optim_constV(y_ld[, 1:p], x$xs, x$xe, 
                                       mean = preest$mean, k = preest$k, a = preest$a, G0 = preest$G0,
                                       G01behaviour = "fixed",
                                       type = "Kassel", intercept = FALSE)}, "(p!=3|multimodal)")
  expect_equal(est2$mean, est1$mean, tolerance = 1E-2)
  expect_equal(est2[c("k", "a")], est1[c("k", "a")], tolerance = 1E-1, ignore_attr = "names")
  expect_equal(est2$G0, est1$G0, tolerance = 1E-2, ignore_attr = TRUE)
  
  # from default starts with prelim estimate
  expect_warning({est3 <- mobius_SvMF(y_ld[, 1:p], x$xs, x$xe, 
                      G0 = cbind(G0[,1], matrix(NA, p, p-1)),
                      G01behaviour = "fixed",
                      type = "Kassel", intercept = FALSE)}, "(p!=3|multimodal)")
  expect_equal(est3$mean, est1$mean, tolerance = 1E-1)
  expect_equal(est3$k, est1$k, tolerance = 0.2)
  expect_equal(est3$a, est1$a, tolerance = 1E-1, ignore_attr = "names")
  expect_equal(axis_distance(acos(colSums(est3$G0 * est1$G0) - 1E-15)), rep(0, p), tolerance = 1E-1, ignore_attr = TRUE)
})

test_that("MLE with G01 free, p=5", {
  rmnlink_cann__place_in_env(5, 5, 6)
  omegapar <- as_mnlink_Omega(paramobj)
  set.seed(4)
  x <- rcovars(1000, qs, qe)
  set.seed(5)
  G0 <- mclust::randomOrthogonalMatrix(p, p) #axes around a random location
  G0[, p] <- det(G0) * G0[,p]
  
  # simulate observations
  set.seed(6)
  k <- 10
  a <- c(1, seq(0.8, 0.2, length.out = p-1))
  a[-1] <- a[-1]/prod(a[-1])^(1/(p-1))
  SvMFcann_check(SvMFcann(k, a, G0))
  y_ld <- rMobius_SvMF(x$xs, x$xe, mnparam = omegapar, k, a, G0)
  
  # check ull_S2S_constV in C++
  ldCpp <- ull_S2S_constV_forR(y = y_ld[, 1:p], xs = x$xs, xe = x$xe, omvec = mnlink_Omega_vec(omegapar), k = k,
                               a1 = a[1], aremaining = a[-1], G0 = G0)
  expect_equal(ldCpp, y_ld[, p+1])
  
  
  #check tape:
  set.seed(7)
  referencecoords <- mclust::randomOrthogonalMatrix(p, p)
  referencecoords[, p] <- det(referencecoords) * referencecoords[,p]
  vecparams <- S2S_constV_nota1_tovecparams(omvec = mnlink_Omega_vec(omegapar), k = k,
                                            aremaining = a[-1],
                                            G0 = G0, 
                                            referencecoords = referencecoords,
                                            G01behaviour = "free")
  expect_warning({ulltape <- tape_ull_S2S_constV_nota1(omvec = mnlink_Omega_vec(omegapar), k = k,
                                                       a1 = a[1], aremaining = a[-1],
                                                       G0 = G0,
                                                       p, qe,
                                                       yx = cbind(y_ld[, 1:p], x$xs, x$xe), 
                                                       referencecoords,
                                                       G01behaviour = "free")}, "p!=3")
  expect_equal(ulltape$xtape, vecparams)
  expect_equal(ulltape$forward(0, ulltape$xtape), y_ld[, p+1])
  
  exactll <- sum(ulltape$forward(0, ulltape$xtape))
  
  set.seed(18)
  G0_other <- mclust::randomOrthogonalMatrix(p, p) #axes around a random location
  G0_other[, p] <- det(G0_other) * G0_other[,p]
  badll <- sum(ulltape$forward(0, 
                               S2S_constV_nota1_tovecparams(omvec = mnlink_Omega_vec(omegapar), k = k,
                                                            aremaining = a[-1], 
                                                            G0 = G0_other,
                                                            referencecoords = referencecoords,
                                                            G01behaviour = "free")))
  expect_lt(badll, exactll)
  
  ## now try optimisation starting at true values ##
  expect_warning({est1 <- optim_constV(y_ld[, 1:p], x$xs, x$xe, omegapar, k, a, G0, xtol_rel = 1E-4, G0reference = referencecoords, G01behaviour = "free",
                                       type = "Kassel", intercept = FALSE)}, "p!=3")
  expect_equal(est1$mean, omegapar, tolerance = 1E-1)
  expect_equal(est1[c("k", "a")], list(k = k, a = a), tolerance = 1E-1)
  # check Gstar by checking angle between estimated and true axes
  expect_equal(axis_distance(acos(colSums(est1$G0 * G0))), rep(0, p), tolerance = 1E-1, ignore_attr = TRUE)
  # check likelihood returns
  externalll <- colSums(ldMobius_SvMF(y_ld[, 1:p], x$xs, x$xe, est1$mean, est1$k, est1$a, est1$G0))
  expect_equal(externalll[["R"]], est1$lLik)
  
  ## now starting optimisation away from starting parameters ##
  bad_om <- as_mnlink_Omega(rmnlink_cann(p, qs, qe, preseed = 2))
  preest <- mobius_SvMF_partransport_prelim(y_ld[, 1:p], x$xs, x$xe, mean = bad_om, type = "Kassel", intercept = FALSE, G01behaviour = "free")
  expect_warning({est2 <- optim_constV(y_ld[, 1:p], x$xs, x$xe, 
                                       mean = preest$mean, k = preest$k, a = preest$a, G0 = preest$G0,
                                       G01behaviour = "free",
                                       type = "Kassel", intercept = FALSE)}, "p!=3")
  expect_equal(est2$mean, est1$mean, tolerance = 1E-2)
  expect_equal(est2[c("k", "a")], est1[c("k", "a")], tolerance = 1E-1, ignore_attr = "names")
  expect_equal(axis_distance(acos(colSums(est2$G0 * est1$G0))), rep(0, p), tolerance = 1E-1, ignore_attr = TRUE)
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

test_that("Under high concentration, standardised residuals are the correct MVN", {
  rmnlink_cann__place_in_env(4, 5, 4)
  omegapar <- as_mnlink_Omega(paramobj)
  set.seed(4)
  x <- rcovars(1000, qs, qe)
  set.seed(5)
  G0 <- mclust::randomOrthogonalMatrix(p, p) #axes around a random location
  G0[, p] <- det(G0) * G0[,p]
  
  # simulate observations
  set.seed(6)
  k <- 500
  a <- c(1, seq(5, 0.2, length.out = p-1))
  a[-1] <- a[-1]/prod(a[-1])^(1/(p-1))
  y_ld <- rMobius_SvMF(x$xs, x$xe, mnparam = omegapar, k, a, G0)
  
  exactresids <- resid_SvMF_partransport(y_ld[, 1:p], mnlink(xs = x$xs, xe = x$xe, param = omegapar), k, a, G0)
  expect_gt(ks.test(rowSums(exactresids^2), "pchisq", df = ncol(exactresids))$p.value, 0.05)
  
  expect_warning({est1 <- optim_constV(y_ld[, 1:p], x$xs, x$xe, omegapar, k, a, G0, xtol_rel = 1E-4, G0reference = G0, G01behaviour = "fixed",
                                       type = "Kassel", intercept = FALSE)}, "p!=3")
  
  expect_gt(ks.test(rowSums(est1$rresids_std^2), "pchisq", df = ncol(est1$rresids_std))$p.value, 0.05)
})

