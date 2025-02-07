test_that("prelim optimisation works with Euc covars", {
  rmnlink_cann__place_in_env(qs = 0)
  
  #generate covariates Gaussianly
  set.seed(4)
  x <- matrix(rnorm(1000*qe), nrow = 1000)
  ymean <- mnlink(xe = x, param = paramobj)
  
  # generate noise
  if (!requireNamespace("movMF", quietly = TRUE)){skip("Need movMF package")}
  set.seed(5)
  y <- t(apply(ymean, 1, function(mn){movMF::rmovMF(1, 30*mn)}))
  
  # optimise locally using derivative information
  # starting at the optimum
  tmp <- prelim_ad(y, xe = x, paramobj0 = as_mnlink_Omega(paramobj))
  expect_equal(tmp$solution, as_mnlink_Omega(paramobj), tolerance = 0.05)
  
  # starting away from optimum, but still within constraints
  set.seed(14)
  start <- as_mnlink_Omega(mnlink_cann(P = mclust::randomOrthogonalMatrix(p, p),
                                       Qe = mclust::randomOrthogonalMatrix(qe, p),
                                       Be = diag(sort(runif(p-1), decreasing = TRUE)),
                                       ce = rep(0, p)))
  opt2 <- prelim_ad(y, xe = x, paramobj0 = start)
  if (sign(opt2$solution$qe1[1]) != sign(as_mnlink_Omega(paramobj)$qe1[1])){
    opt2$solution <- Omega_Euc_signswitch(opt2$solution)
  }
  expect_equal(opt2$solution, as_mnlink_Omega(paramobj), tolerance = 0.05)
})

test_that("prelim optimisation works with Sph covars",{
  rmnlink_cann__place_in_env(p = 4, qs = 5, qe = 0)
  omegapar <- as_mnlink_Omega(paramobj)
  
  #generate covariates uniformly on the sphere
  set.seed(4)
  x <- matrix(rnorm(1000*qs), nrow = 1000)
  x <- sweep(x, 1, apply(x, 1, vnorm), FUN = "/")
  
  ymean <- mnlink(xs = x, param = omegapar)
  
  # generate noise
  if (!requireNamespace("movMF", quietly = TRUE)){skip("Need movMF package")}
  set.seed(5)
  y <- t(apply(ymean, 1, function(mn){movMF::rmovMF(1, 30*mn)}))
  
  # objective function in C++ and R should match when omegapar passes mnlink_Omega_check()
  objval <- prelimobj(y, xs = x, param = omegapar)
  objvalcpp <-  prelimobj_cpp(mnlink_Omega_vec(omegapar), vector(), c(p, qe), cbind(y,x))
  expect_equal(objvalcpp, objval)

  # optimise using pure R
  # opt <- prelim_R(y, x, omegapar, global = TRUE, local = TRUE)
  # expect_equal(Omega_proj(opt$solution), opt$solution, tolerance = 1E-3)
  # expect_equal(opt$solution, omegapar, tolerance = 0.05)
  
  # optimise locally using derivative information
  # starting at the optimum
  tmp <- prelim_ad(y, xs = x, paramobj0 = omegapar)
  expect_equal(tmp$solution, omegapar, tolerance = 0.05)
  
  # starting away from optimum, but still within constraints
  set.seed(14)
  start <- as_mnlink_Omega(cannS2S(P = mclust::randomOrthogonalMatrix(p, p),
          Q = mclust::randomOrthogonalMatrix(qs, p),
          B = diag(sort(runif(p-1), decreasing = TRUE))))
  opt2 <- prelim_ad(y, xs = x, paramobj0 = start)
  expect_equal(opt2$solution, omegapar, tolerance = 0.05)
  
  # check global
  opt3 <- prelim_ad(y, xs = x, paramobj0 = start, type = "Kassel", globalfirst = TRUE)
  # if (sign(opt3$solution$qe1[1]) != sign(as_mnlink_Omega(paramobj)$qe1[1])){
  #   opt3$solution <- Euc_signswitch(opt3$solution)
  # }
  # expect_equal(opt3$solution, as_mnlink_Omega(paramobj), tolerance = 0.05)
})

test_that("prelim optimisation works with Sph+Euc covars", {
  rmnlink_cann__place_in_env(3, 5, 4)
  
  #generate covariates Gaussianly
  set.seed(4)
  xe <- matrix(rnorm(1000*qe), nrow = 1000)
  #generate covariates on the sphere
  set.seed(4)
  xs <- matrix(rnorm(1000*qs), nrow = 1000)
  xs <- sweep(xs, 1, apply(xs, 1, vnorm), FUN = "/")
  
  ymean <- mnlink(xs = xs, xe = xe, param = paramobj)
  
  # generate noise
  if (!requireNamespace("movMF", quietly = TRUE)){skip("Need movMF package")}
  set.seed(5)
  y <- t(apply(ymean, 1, function(mn){movMF::rmovMF(1, 30*mn)}))
  
  # optimise locally using derivative information
  # starting at the optimum
  tmp <- prelim_ad(y, xs = xs, xe = xe, paramobj0 = as_mnlink_Omega(paramobj))
  expect_equal(tmp$solution, as_mnlink_Omega(paramobj), tolerance = 0.05)
  
  # starting away from optimum, but still within constraints
  set.seed(14)
  start <- as_mnlink_Omega(mnlink_cann(P = mclust::randomOrthogonalMatrix(p, p),
                                       Qs = mclust::randomOrthogonalMatrix(qs, p),
                                       Bs = diag(sort(runif(p-1), decreasing = TRUE)),
                                       Qe = mclust::randomOrthogonalMatrix(qe, p),
                                       Be = diag(sort(runif(p-1), decreasing = TRUE)),
                                       ce = c(1, rep(0, p-1))))
  opt2 <- prelim_ad(y, xs = xs, xe = xe, paramobj0 = start)
  if (sign(opt2$solution$qe1[1]) != sign(as_mnlink_Omega(paramobj)$qe1[1])){
    opt2$solution <- Euc_signswitch(opt2$solution)
  }
  expect_equal(opt2$solution, as_mnlink_Omega(paramobj), tolerance = 0.05)
  
  # check global
  # Global seems to get the constraints wrong instantly!!
  # (these are the 3 'h(x)' values, and they should be very close to 0)
  # opt3 <- prelim_global(y, xs = xs, xe = xe, paramobj0 = start, type = "Kassel")
  # if (sign(opt3$solution$qe1[1]) != sign(as_mnlink_Omega(paramobj)$qe1[1])){
  #   opt3$solution <- Euc_signswitch(opt3$solution)
  # }
  # expect_equal(opt3$solution, as_mnlink_Omega(paramobj), tolerance = 0.05)
})

test_that("Shogo with Sph+Euc covars", {
  rmnlink_cann__place_in_env(3, 5, 4)
  # convert to Shogo form:
  bigQe <- rbind(0, Qe)
  bigQe[, 1] <- 0
  bigQe[1,1] <- 1
  bigce <- ce
  bigce[1] <- 1
  ce <- ce[-1]
  paramobj <- mnlink_cann(P, Qs = Qs, Bs = Bs, Be = Be, Qe = bigQe, ce = bigce)
  expect_true(is_Shogo(paramobj))
  
  #generate covariates Gaussianly
  set.seed(4)
  xe <- cbind(0, matrix(rnorm(1000*qe), nrow = 1000))
  #generate covariates on the sphere
  set.seed(4)
  xs <- matrix(rnorm(1000*qs), nrow = 1000)
  xs <- sweep(xs, 1, apply(xs, 1, vnorm), FUN = "/")
  
  ymean <- mnlink(xs = xs, xe = xe, param = paramobj)
  
  # generate noise
  if (!requireNamespace("movMF", quietly = TRUE)){skip("Need movMF package")}
  set.seed(5)
  y <- t(apply(ymean, 1, function(mn){movMF::rmovMF(1, 30*mn)}))
  
  # optimise locally using derivative information
  # starting at the optimum
  tmp <- prelim_ad(y, xs = xs, xe = xe, paramobj0 = as_mnlink_Omega(paramobj), type = "Shogo")
  expect_equal(tmp$solution, as_mnlink_Omega(paramobj), tolerance = 0.05)
  
  # starting away from optimum, but still within constraints
  set.seed(14)
  start <- as_mnlink_Omega(mnlink_cann(P = mclust::randomOrthogonalMatrix(p, p),
                                       Qs = mclust::randomOrthogonalMatrix(qs, p),
                                       Bs = diag(sort(runif(p-1), decreasing = TRUE)),
                                       Qe = mclust::randomOrthogonalMatrix(qe, p),
                                       Be = diag(sort(runif(p-1), decreasing = TRUE)),
                                       ce = rep(0, p)))
  # convert to start to Shogo form:
  start$Omega <- cbind(start$Omega[,1:qs], 0, start$Omega[,qs + (1:qe)])
  start$qe1 <- c(1, rep(0, qe))
  start$ce1 <- 1
  opt2 <- prelim_ad(y, xs = xs, xe = xe, paramobj0 = start, type = "Shogo")
  if (sign(opt2$solution$qe1[1]) != sign(as_mnlink_Omega(paramobj)$qe1[1])){
    opt2$solution <- Euc_signswitch(opt2$solution)
  }
  expect_equal(opt2$solution, as_mnlink_Omega(paramobj), tolerance = 0.05)
})



test_that("C++ Omega_constraints() is zero correctly", {
  rmnlink_cann__place_in_env(3, 0, 4)
  expect_equal(Omega_constraints(mnlink_Omega_vec(cann2Omega(paramobj)), p, qe), rep(0, 1 + 2*(qs>0) + 2*(qe>0)))
  
  rmnlink_cann__place_in_env(3, 5, 0)
  expect_equal(Omega_constraints(mnlink_Omega_vec(cann2Omega(paramobj)), p, qe), rep(0, 1 + 2*(qs>0) + 2*(qe>0)))

  rmnlink_cann__place_in_env(3, 5, 4)
  expect_equal(Omega_constraints(mnlink_Omega_vec(cann2Omega(paramobj)), p, qe), rep(0, 1 + 2*(qs>0) + 2*(qe>0)))
  
  # check Jacobian - if a row is zero nlopt gives a round off error
  dims_in <- c(p, qe)
  vec_om0 <- mnlink_Omega_vec(as_mnlink_Omega(paramobj))
  constraint_tape <- tape_namedfun("Omega_constraints_wrap", vec_om0, vector(mode = "numeric"), dims_in, matrix(nrow = 0, ncol = 0), check_for_nan = FALSE)
  Jac <- matrix(constraint_tape$Jacobian(vec_om0), byrow = TRUE, ncol = length(vec_om0))
  colnames(Jac) <- names(vec_om0)
  round(Jac, 3)
  expect_true(all(apply(Jac, 1, function(x)max(abs(x))) > 0.1))
  expect_true(all(abs(svd(Jac)$d) > sqrt(.Machine$double.eps)))
})

test_that("C++ Omega_constraints() is non-zero correctly", {
  rmnlink_cann__place_in_env(3, 5, 4)
  Om <- as_mnlink_Omega(paramobj)
  Om$qe1 <- Om$qe1*2
  Om$qs1 <- Om$qs1*2
  Om$p1 <- Om$p1*2
  Om$Omega <- matrix(rnorm(p * (qs + qe)), p, qs + qe)
  Om$PBce <- Om$PBce+1
  expect_true(all(mnlink_Omega_check_numerical(Om) > 1E-3))
  expect_true(all(Omega_constraints(mnlink_Omega_vec(Om), p, qe) > 1E-3))
})

test_that("pre_est3_mod optimisation works", {
  rmnlink_cann__place_in_env(3, 3, 0)

set.seed(1)
n=100
x <- matrix(rnorm(n*qs), nrow = n)
x <- sweep(x, 1, apply(x, 1, vnorm), FUN = "/") #covariates
y <- mnlink(xs = x, param = paramobj) # assume y=mu(x) for simulation purpose

ini_value=c(0,0,0,0,0,0,0.9,0.2)  # initial value for optimization. Try multiple values.
result=nlminb(start=ini_value,objective=function(theta) pre_est3_mod(y,x,theta),lower=c(rep(-Inf,6),0,0),upper=c(rep(Inf,6),1,1))  # numerical optimization

expect_equal(cayley(result$par[1:3]), P, tolerance = 10^4*sqrt(.Machine$double.eps))
expect_equal(cayley(result$par[4:6]), Qs, tolerance = 10^4*sqrt(.Machine$double.eps))
expect_equal(result$par[7], Bs[1,1], tolerance = 1E-3)
expect_equal(result$par[8] * result$par[7], Bs[2,2], tolerance = 1E-3)
})



