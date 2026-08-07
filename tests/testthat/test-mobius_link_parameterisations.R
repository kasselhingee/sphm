test_that("mobius_link_cann() objects pass check: Sph only", {
  set.seed(1)
  p <- 3
  qs <- 5
  P <- mclust::randomOrthogonalMatrix(p, p)
  set.seed(2)
  Qs <- mclust::randomOrthogonalMatrix(qs, p)
  set.seed(3)
  Bs <- diag(sort(runif(p-1), decreasing = TRUE))
  obj <- mobius_link_cann(P, Bs = Bs, Qs = Qs, check = FALSE)
  expect_silent(mobius_link_cann_check(obj))
})

test_that("mobius_link_cann() objects pass check: Euc only", {
  set.seed(1)
  p <- 3
  qe <- 5
  P <- mclust::randomOrthogonalMatrix(p, p)
  set.seed(2)
  Qe <- mclust::randomOrthogonalMatrix(qe, p)
  set.seed(3)
  Be <- diag(sort(runif(p-1), decreasing = TRUE))
  set.seed(4)
  ce <- runif(1)
  obj <- mobius_link_cann(P, Be = Be, Qe = Qe, ce = ce, check = FALSE)
  expect_silent(mobius_link_cann_check(obj))
})

test_that("mobius_link_cann() objects pass check: Sph + Euc only", {
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
  ce <- runif(1)
  obj <- mobius_link_cann(P, Bs = Bs, Qs = Qs, Be = Be, Qe = Qe, ce = ce, check = FALSE)
  expect_silent(mobius_link_cann_check(obj))
})

test_that("mobius_link_cann(): common mistakes", {
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
  ce <- runif(1)
  
  expect_error(mobius_link_cann(P, Bs = Bs, Qs = Qs, Be = Be, Qe = Qe, ce = NULL, check = TRUE))
  expect_error(mobius_link_cann(P, Bs = Bs, Qe = Qs, check = TRUE))
  expect_error(mobius_link_cann(P, Be = Be, Qs = Qe, check = TRUE))
})

test_that("Conversions work: Sph + Euc", {
  rand_mobius_link_cann__place_in_env(3, 5, 4)
  cann <- paramobj
  
  Om <- cann2Omega(cann, check = FALSE)
  

  
  # check that convert back matches
  cann2 <- Omega2cann(Om, check = FALSE)
  expect_silent(mobius_link_cann_check(cann2))
  cann2 <- P_signswitch(cann2, sign(cann2$P[1, ]) != sign(cann$P[1, ]))
  expect_equal(cann2, cann)
  
  # check properties of Om
  expect_equal(sum(diag(t(Om$Omega) %*% Om$Omega)), sum(Be^2 + Bs^2))
  expect_silent(mobius_link_Omega_check(Om))
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
  expect_gt(mobius_link_Omega_check_numerical(Om)["Omega_comm"], sqrt(.Machine$double.eps))
})

test_that("Conversions work: Sph only", {
  rand_mobius_link_cann__place_in_env(3, 5, 0)
  cann <- paramobj
  
  Om <- cann2Omega(cann, check = FALSE)
  
  # check properties of Om
  expect_equal(sum(diag(t(Om$Omega) %*% Om$Omega)), sum(Bs^2))
  expect_silent(mobius_link_Omega_check(Om))
  
  # check that convert back matches
  cann2 <- Omega2cann(Om, check = FALSE)
  expect_silent(mobius_link_cann_check(cann2))
  cann2 <- P_signswitch(cann2, sign(cann2$P[1, ]) != sign(cann$P[1, ]))
  expect_equal(cann2, cann)
})

test_that("Conversions work: Euc only", {
  rand_mobius_link_cann__place_in_env(3, 0, 4)
  cann <- paramobj
  
  Om <- cann2Omega(cann, check = FALSE)
  
  # check properties of Om
  expect_equal(sum(diag(t(Om$Omega) %*% Om$Omega)), sum(Be^2))
  expect_silent(mobius_link_Omega_check(Om))
  
  # check that convert back matches
  cann2 <- Omega2cann(Om, check = FALSE)
  expect_silent(mobius_link_cann_check(cann2))
  cann2 <- P_signswitch(cann2, sign(cann2$P[1, ]) != sign(cann$P[1, ]))
  expect_equal(cann2, cann)
})


test_that("mobius_link_Omega works directly", {
  rand_mobius_link_cann__place_in_env(3, 5, 4)
  cann <- paramobj
  Om <- cann2Omega(cann)
  
  Om2 <- mobius_link_Omega(p1 = Om$p1, qs1 = Om$qs1, qe1 = Om$qe1, Omega = Om$Omega, ce = Om$ce, check = TRUE)
  expect_equal(Om2, Om)
})
  
test_that("projections orthogonal to p1 and qe1, qs1", {
  rand_mobius_link_cann__place_in_env(3, 5, 4)
  cann <- paramobj
  Om <- cann2Omega(cann)
  
  # project Omega perpendicular to p1 and q1
  expect_equal(Omega_proj(Om), Om)
  Ommod <- Om
  Ommod$Omega <- Om$Omega * (1 + 1E-2 * matrix(runif(length(Om$Omega), -1, 1), p, qs + qe))
  expect_error(mobius_link_Omega_check(Ommod), "checks failed")
  expect_equal(Omega_proj(Ommod)$Omega, Om$Omega, tolerance = 1E-2)
  expect_equal(Omega_proj(Ommod)[c("qs1", "qe1", "p1")], Om[c("qs1", "qe1", "p1")])
  
  # just Sph
  cann <- mobius_link_cann(P, Bs = Bs, Qs = Qs)
  Om <- cann2Omega(cann)
  expect_equal(Omega_proj(Om), Om)
  Ommod <- Om
  Ommod$Omega <- Om$Omega * (1 + 1E-2 * matrix(runif(length(Om$Omega), -1, 1), p, qs))
  expect_error(mobius_link_Omega_check(Ommod), "checks failed")
  expect_equal(Omega_proj(Ommod)$Omega, Om$Omega, tolerance = 1E-2)
  expect_equal(Omega_proj(Ommod)[c("qs1", "qe1", "p1")], Om[c("qs1", "qe1", "p1")])
  
  # just Euc
  cann <- mobius_link_cann(P, Be = Be, Qe = Qe, ce = ce)
  Om <- cann2Omega(cann)
  expect_equal(Omega_proj(Om), Om)
  Ommod <- Om
  Ommod$Omega <- Om$Omega * (1 + 1E-2 * matrix(runif(length(Om$Omega), -1, 1), p, qe))
  expect_error(mobius_link_Omega_check(Ommod), "checks failed")
  expect_equal(Omega_proj(Ommod)$Omega, Om$Omega, tolerance = 1E-2)
  expect_equal(Omega_proj(Ommod)[c("qs1", "qe1", "p1")], Om[c("qs1", "qe1", "p1")])
}) 
  
test_that("vec and unvec", {
  rand_mobius_link_cann__place_in_env(3, 5, 4)
  cann <- paramobj
  Om <- cann2Omega(cann)
  
  #vec and unvec
  expect_equal(mobius_link_Omega_unvec(mobius_link_Omega_vec(Om), p, qe = qe), Om)
  
  #sph only
  Om <- cann2Omega(mobius_link_cann(P, Bs = Bs, Qs = Qs))
  expect_equal(mobius_link_Omega_unvec(mobius_link_Omega_vec(Om), p, qe = 0), Om)
  
  #Euc only
  Om <- cann2Omega(mobius_link_cann(P, Be = Be, Qe = Qe, ce = ce))
  expect_equal(mobius_link_Omega_unvec(mobius_link_Omega_vec(Om), p, qe = qe), Om)
})

test_that("check LinEuc conversion", {
  rand_mobius_link_cann__place_in_env(3, 0, 4)
  # convert to LinEuc form:
  bigQe <- rbind(0, Qe)
  bigQe[, 1] <- 0
  bigQe[1,1] <- 1
  ce <- 1
  paramobj <- mobius_link_cann(P, Be = Be, Qe = bigQe, ce = ce)
  expect_true(is_LinEuc(paramobj))
  
  # check manual
  x <- runif(qe)
  expect_equal(drop(mobius_link(xe = matrix(c(0, x), nrow = 1), param = paramobj)),
               drop(P %*% iSp(drop(Be %*% (t(bigQe[-1,-1]) %*% x)))))
  
  # check LinEuc:
  # need to check linearity of P^{-1} S(y) wrt x_e or even Qe
  # expect differences between the intermediate output to be purely due to the difference between inputs
  out <- Sp(mobius_link(xe = rbind(0, c(0,x), c(0, 2*x)), param = paramobj) %*% P)
  expect_equal(out[2, ] - out[1, ], out[3, ] - out[2, ])
})

test_that("rand_mobius_link_cann preseed gives identical results on repeated calls", {
  p1 <- rand_mobius_link_cann(3, 5, 4, preseed = 7)
  p2 <- rand_mobius_link_cann(3, 5, 4, preseed = 7)
  expect_equal(p1, p2)
})


test_that("signflip_starts() enumerates the right number of starts", {
  both    <- rand_mobius_link_cann(3, 4, 5, preseed = 11)
  xe_only <- rand_mobius_link_cann(3, 0, 5, preseed = 12)
  xs_only <- rand_mobius_link_cann(3, 4, 0, preseed = 13)
  neither <- rand_mobius_link_cann(3, 0, 0, preseed = 14)

  # 2 (p1) x 2 (qs1) x 2 (qe1) x 2 (ce regime)
  expect_length(signflip_starts(both, FALSE, FALSE), 16)
  # fix_qe1 removes both the qe1 flip and the ce variation
  expect_length(signflip_starts(both, FALSE, TRUE), 4)
  expect_length(signflip_starts(both, TRUE, TRUE), 2)
  expect_length(signflip_starts(xe_only, FALSE, TRUE), 2)
  # An absent component must drop out of the product, not collapse it: these two returned
  # zero starts when the enumeration was a nested loop with ce outermost.
  expect_length(signflip_starts(xs_only, FALSE, FALSE), 4)
  expect_length(signflip_starts(neither, FALSE, FALSE), 2)
})

test_that("signflip_starts() varies ce only when qe1 is free", {
  obj <- rand_mobius_link_cann(3, 4, 5, preseed = 15)
  ces <- function(starts) vapply(starts, function(x) x$ce, numeric(1))

  # rand_mobius_link_cann draws ce from runif(1), so ce < 10 and the alternative regime is 100
  expect_setequal(unique(ces(signflip_starts(obj, FALSE, FALSE))), c(obj$ce, 100))
  # ce and qe1 are fixed together by estprep_meanconstraints(), so a fixed qe1 must leave ce alone
  expect_identical(unique(ces(signflip_starts(obj, FALSE, TRUE))), obj$ce)

  # the opposite size regime: ce > 10 sends the alternative to 1
  big <- obj
  big$ce <- 50
  expect_setequal(unique(ces(signflip_starts(big, FALSE, FALSE))), c(50, 1))
  expect_identical(unique(ces(signflip_starts(big, FALSE, TRUE))), 50)
})

test_that("signflip_starts() preserves LinEuc form when qe1 is fixed", {
  rand_mobius_link_cann__place_in_env(3, 0, 4)
  bigQe <- rbind(0, Qe)
  bigQe[, 1] <- 0
  bigQe[1, 1] <- 1
  lin <- mobius_link_cann(P, Be = Be, Qe = bigQe, ce = 1)
  expect_true(is_LinEuc(lin))

  # LinEuc requires ce == 1 and qe1 == (1, 0, ...); every start must still satisfy that,
  # otherwise the is_LinEuc() checks in mobius_vMF()/mobius_SvMF_joint_fit() trip.
  starts <- signflip_starts(lin, FALSE, TRUE)
  expect_true(all(vapply(starts, is_LinEuc, logical(1))))
})

test_that("signflip_starts() leaves fixed poles alone and returns the original first", {
  both    <- rand_mobius_link_cann(3, 4, 5, preseed = 16)
  xs_only <- rand_mobius_link_cann(3, 4, 0, preseed = 17)

  # the unmodified starting point is always offered first
  expect_identical(signflip_starts(both, FALSE, FALSE)[[1]], both)
  expect_identical(signflip_starts(xs_only, FALSE, FALSE)[[1]], xs_only)

  fixed <- signflip_starts(both, TRUE, TRUE)
  expect_true(all(vapply(fixed, function(x) isTRUE(all.equal(x$Qs[, 1], both$Qs[, 1])), logical(1))))
  expect_true(all(vapply(fixed, function(x) isTRUE(all.equal(x$Qe[, 1], both$Qe[, 1])), logical(1))))

  free <- signflip_starts(both, FALSE, FALSE)
  # both orientations of each free pole appear, and P stays a rotation (det = +1) throughout
  expect_length(unique(lapply(free, function(x) sign(x$Qs[, 1]))), 2)
  expect_length(unique(lapply(free, function(x) sign(x$Qe[, 1]))), 2)
  expect_length(unique(lapply(free, function(x) sign(x$P[, 1]))), 2)
  expect_equal(vapply(free, function(x) det(x$P), numeric(1)), rep(det(both$P), length(free)))

  # ce is absent without Euclidean covariates and must stay that way
  expect_true(all(vapply(signflip_starts(xs_only, FALSE, FALSE), function(x) is.null(x$ce), logical(1))))
})
