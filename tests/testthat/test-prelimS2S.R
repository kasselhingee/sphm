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
  tmp <- optim_pobjS2S_parttape(y, xe = x, paramobj0 = as_mnlink_Omega(paramobj))
  expect_equal(tmp$solution, as_mnlink_Omega(paramobj), tolerance = 0.05)
  
  # starting away from optimum, but still within constraints
  set.seed(14)
  start <- as_mnlink_Omega(mnlink_cann(P = mclust::randomOrthogonalMatrix(p, p),
                                       Qe = mclust::randomOrthogonalMatrix(qe, p),
                                       Be = diag(sort(runif(p-1), decreasing = TRUE)),
                                       ce = rep(0, p)))
  opt2 <- optim_pobjS2S_parttape(y, xe = x, paramobj0 = start)
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
  objval <- pobjS2S(y, x, omegapar)
  objvalcpp <-  prelimobj_cpp(mnlink_Omega_vec(omegapar), vector(), c(p, qe), cbind(y,x))
  expect_equal(objvalcpp, objval)

  # optimise using pure R
  opt <- optim_pobjS2S_pureR(y, x, omegapar, global = TRUE, local = TRUE)
  expect_equal(Omega_proj(opt$solution), opt$solution, tolerance = 1E-3)
  expect_equal(opt$solution, omegapar, tolerance = 0.05)
  
  # optimise locally using derivative information
  # starting at the optimum
  tmp <- optim_pobjS2S_parttape(y, xs = x, paramobj0 = omegapar)
  expect_equal(tmp$solution, omegapar, tolerance = 0.05)
  
  # starting away from optimum, but still within constraints
  set.seed(14)
  start <- as_mnlink_Omega(cannS2S(P = mclust::randomOrthogonalMatrix(p, p),
          Q = mclust::randomOrthogonalMatrix(qs, p),
          B = diag(sort(runif(p-1), decreasing = TRUE))))
  opt2 <- optim_pobjS2S_parttape(y, xs = x, paramobj0 = start)
  expect_equal(opt2$solution, omegapar, tolerance = 0.05)
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
  tmp <- optim_pobjS2S_parttape(y, xs = xs, xe = xe, paramobj0 = as_mnlink_Omega(paramobj))
  expect_equal(tmp$solution, as_mnlink_Omega(paramobj), tolerance = 0.05)
  
  # starting away from optimum, but still within constraints
  set.seed(14)
  start <- as_mnlink_Omega(mnlink_cann(P = mclust::randomOrthogonalMatrix(p, p),
                                       Qs = mclust::randomOrthogonalMatrix(qs, p),
                                       Bs = diag(sort(runif(p-1), decreasing = TRUE)),
                                       Qe = mclust::randomOrthogonalMatrix(qe, p),
                                       Be = diag(sort(runif(p-1), decreasing = TRUE)),
                                       ce = rep(0, p)))
  opt2 <- optim_pobjS2S_parttape(y, xs = xs, xe = xe, paramobj0 = start)
  if (sign(opt2$solution$qe1[1]) != sign(as_mnlink_Omega(paramobj)$qe1[1])){
    opt2$solution <- Omega_Euc_signswitch(opt2$solution)
  }
  expect_equal(opt2$solution, as_mnlink_Omega(paramobj), tolerance = 0.05)
})



test_that("C++ Omega_constraints() is zero correctly", {
  rmnlink_cann__place_in_env(3, 5, 4)
  expect_equal(Omega_constraints(mnlink_Omega_vec(cann2Omega(paramobj)), p, qe), rep(0, 1 + (qs>0) + (qe>0)))
  
  rmnlink_cann__place_in_env(3, 0, 4)
  expect_equal(Omega_constraints(mnlink_Omega_vec(cann2Omega(paramobj)), p, qe), rep(0, 1 + (qs>0) + (qe>0)))
  
  rmnlink_cann__place_in_env(3, 5, 0)
  expect_equal(Omega_constraints(mnlink_Omega_vec(cann2Omega(paramobj)), p, qe), rep(0, 1 + (qs>0) + (qe>0)))
})

test_that("pre_est3_mod optimisation works", {
  rmnlink_cann__place_in_env(3, 3, 0)

set.seed(1)
n=100
x <- matrix(rnorm(n*qs), nrow = n)
x <- sweep(x, 1, apply(x, 1, vnorm), FUN = "/") #covariates
y <- mnlink(xs = u, param = paramobj) # assume y=mu(x) for simulation purpose

ini_value=c(0,0,0,0,0,0,0.9,0.2)  # initial value for optimization. Try multiple values.
result=nlminb(start=ini_value,objective=function(theta) pre_est3_mod(y,x,theta),lower=c(rep(-Inf,6),0,0),upper=c(rep(Inf,6),1,1))  # numerical optimization

expect_equal(cayley(result$par[1:3]), P, tolerance = 10^4*sqrt(.Machine$double.eps))
expect_equal(cayley(result$par[4:6]), Qs, tolerance = 10^4*sqrt(.Machine$double.eps))
expect_equal(result$par[7], Bs[1,1], tolerance = 1E-3)
expect_equal(result$par[8] * result$par[7], Bs[2,2], tolerance = 1E-3)
})



