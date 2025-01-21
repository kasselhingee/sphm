test_that("prelim optimisation works with Euc covars", {
  set.seed(1)
  p <- 3
  P <- mclust::randomOrthogonalMatrix(p, p)
  qs <- 0
  qe <- 4
  set.seed(12)
  Qe <- mclust::randomOrthogonalMatrix(qe, p)
  set.seed(13)
  Be <- diag(sort(runif(p-1), decreasing = TRUE))
  set.seed(14)
  ce <- runif(p)
  paramobj <- mnlink_cann(P, Be = Be, Qe = Qe, ce = ce, check = TRUE)
  
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
  expect_equal(opt2$solution, as_mnlink_Omega(paramobj), tolerance = 0.05)
})

test_that("prelim optimisation works with Sph covars",{
  p <- 4
  q <- 5
  # data generating parameters:
  set.seed(1)
  P <- mclust::randomOrthogonalMatrix(p, p)
  set.seed(2)
  Q <- mclust::randomOrthogonalMatrix(q, p)
  set.seed(3)
  B <- diag(sort(runif(p-1), decreasing = TRUE))
  omegapar <- as_mnlink_Omega(cannS2S(P,Q,B))
  
  #generate covariates uniformly on the sphere
  set.seed(4)
  x <- matrix(rnorm(1000*q), nrow = 1000)
  x <- sweep(x, 1, apply(x, 1, vnorm), FUN = "/")
  
  ymean <- mnlink(xs = x, param = omegapar)
  
  # generate noise
  if (!requireNamespace("movMF", quietly = TRUE)){skip("Need movMF package")}
  set.seed(5)
  y <- t(apply(ymean, 1, function(mn){movMF::rmovMF(1, 30*mn)}))
  
  # objective function in C++ and R should match when omegapar passes mnlink_Omega_check()
  objval <- pobjS2S(y, x, omegapar)
  objvalcpp <-  prelimobj_cpp(mnlink_Omega_vec(omegapar), vector(), c(p, 0), cbind(y,x))
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
          Q = mclust::randomOrthogonalMatrix(q, p),
          B = diag(sort(runif(p-1), decreasing = TRUE))))
  opt2 <- optim_pobjS2S_parttape(y, xs = x, paramobj0 = start)
  expect_equal(opt2$solution, omegapar, tolerance = 0.05)
})

test_that("prelim optimisation works with Sph+Euc covars", {
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
  ce <- runif(p)
  paramobj <- mnlink_cann(P, Bs = Bs, Qs = Qs, Be = Be, Qe = Qe, ce = ce, check = TRUE)
  
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
  expect_equal(Omega_Euc_signswitch(opt2$solution), as_mnlink_Omega(paramobj), tolerance = 0.05)
})



test_that("Omega_constraints() is zero correctly", {
  p <- 3
  q <- 5
  # data generating parameters:
  set.seed(1)
  P <- mclust::randomOrthogonalMatrix(p, p)
  set.seed(2)
  Q <- mclust::randomOrthogonalMatrix(q, p)
  set.seed(3)
  B <- diag(sort(runif(p-1), decreasing = TRUE))
  omegapar <- as_mnlink_Omega(cannS2S(P,Q,B))
  expect_equal(Omega_constraints(mnlink_Omega_vec(omegapar), p), rep(0, 1 + 1))
})

test_that("pre_est3_mod optimisation works", {

d=3  # dimension
P=Q=diag(1,d)  # d*d matrix
B=matrix(c(0.9,0,0,0.2),d-1)  # diagonal matrix
Omega=P[,2:d]%*%B%*%t(Q[,2:d])

set.seed(1)
n=100
u=rnorm(d*n)
x=y=matrix(1,d,n)
for(j in 1:n){
  x[,j]=u[d*(j-1)+1:d]/vnorm(u[d*(j-1)+1:d])  # uniform distribution on the sphere
  y[,j]=meanlinkS2S(x[,j],P,Q,B)  # assume y=mu(x) for simulation purpose
}

ini_value=c(0,0,0,0,0,0,0.9,0.2)  # initial value for optimization. Try multiple values.
result=nlminb(start=ini_value,objective=function(theta) pre_est3_mod(y,x,theta),lower=c(rep(-Inf,6),0,0),upper=c(rep(Inf,6),1,1))  # numerical optimization

expect_equal(cayley(result$par[1:3]), P, tolerance = 10^4*sqrt(.Machine$double.eps))
expect_equal(cayley(result$par[4:6]), Q, tolerance = 10^4*sqrt(.Machine$double.eps))
expect_equal(result$par[7], B[1,1], tolerance = 1E-3)
expect_equal(result$par[8] * result$par[7], B[2,2], tolerance = 1E-3)
})

test_that("taping of pobjS2S and Omega_constraints runs and evaluates", {
  p <- 3
  q <- 5
  # data generating parameters:
  set.seed(1)
  P <- mclust::randomOrthogonalMatrix(p, p)
  set.seed(2)
  Q <- mclust::randomOrthogonalMatrix(q, p)
  set.seed(3)
  B <- diag(sort(runif(p-1), decreasing = TRUE))
  omegapar <- as_mnlink_Omega(cannS2S(P,Q,B))

  
  #generate covariates uniformly on the sphere
  set.seed(4)
  x <- matrix(rnorm(1000*q), nrow = 1000)
  x <- sweep(x, 1, apply(x, 1, vnorm), FUN = "/")
  
  ymean <- meanlinkS2S(x = x, paramobj = omegapar)
  
  # generate noise
  if (!requireNamespace("movMF", quietly = TRUE)){skip("Need movMF package")}
  set.seed(5)
  y <- t(apply(ymean, 1, function(mn){movMF::rmovMF(1, 10*mn)}))

  anADFun <- tape_namedfun("prelimobj_cpp", mnlink_Omega_vec(omegapar), vector(mode = "numeric"), p, cbind(y,x), check_for_nan = FALSE)
  directeval <- prelimobj_cpp(mnlink_Omega_vec(omegapar), vector(), p, cbind(y,x))
  tapeeval <- anADFun$eval(unclass(mnlink_Omega_vec(omegapar)), vector(mode = "numeric"))
  expect_equal(tapeeval, directeval)

  # if taping worked the results should match for a different value of the parameters
  omparo <- omegapar
  omparo$Omega <- omparo$Omega + 1
  omparovec <- mnlink_Omega_vec(omparo)
  directeval <- prelimobj_cpp(omparovec, vector(), p, cbind(y,x))
  tapeeval <- anADFun$eval(unclass(omparovec), vector(mode = "numeric"))
  expect_equal(tapeeval, directeval)

  # check derivatives: pForward, via taping, and numerically
  jac_numeric <- drop(attr(numericDeriv(quote(prelimobj_cpp(omparovec, vector(), p, cbind(y,x))), c("omparovec")), "gradient"))
  jac <- anADFun$Jac(unclass(omparovec), vector(mode = "numeric"))
  expect_equal(jac, jac_numeric, tolerance = 1E-3)

  jactape <- scorematchingad:::tape_Jacobian(anADFun)
  jactapeeval <- jactape$eval(unclass(omparovec), vector(mode = "numeric"))
  expect_equal(jactapeeval, jac)

  # check constraints
  bADFun <- tape_namedfun("Omega_constraints_wrap", mnlink_Omega_vec(omegapar), vector(mode = "numeric"), p, matrix(nrow = 0, ncol = 0), check_for_nan = FALSE)
  directeval <- Omega_constraints(omparovec, p)
  tapeeval <- bADFun$eval(unclass(omparovec), vector(mode = "numeric"))
  expect_equal(tapeeval, directeval)

  jac_numeric <- attr(numericDeriv(quote(bADFun$eval(unclass(omparovec), vector(mode = "numeric"))), c("omparovec")), "gradient")
  jac <- bADFun$Jac(unclass(omparovec), vector(mode = "numeric"))
  expect_equal(matrix(jac, nrow = length(tapeeval), byrow = TRUE), jac_numeric, tolerance = 1E-3)
  jactape <- scorematchingad:::tape_Jacobian(bADFun)
  jactapeeval <- jactape$eval(unclass(omparovec), vector(mode = "numeric"))
  expect_equal(jactapeeval, jac)

  ## check inequality constraints
  cADFun <- tape_namedfun("Omega_ineqconstraints", mnlink_Omega_vec(omegapar), vector(mode = "numeric"), p, matrix(nrow = 0, ncol = 0), check_for_nan = FALSE)
  directeval <- sum(diag(t(mnlink_Omega_unvec(omparovec, p, check = FALSE)$Omega) %*% mnlink_Omega_unvec(omparovec, p, check = FALSE)$Omega)) - (p-1)
  tapeeval <- cADFun$eval(unclass(omparovec), vector(mode = "numeric"))
  expect_equal(tapeeval, directeval)

  jac_numeric <- attr(numericDeriv(quote(cADFun$eval(unclass(omparovec), vector(mode = "numeric"))), c("omparovec")), "gradient")
  jac <- cADFun$Jac(unclass(omparovec), vector(mode = "numeric"))
  expect_equal(matrix(jac, nrow = length(tapeeval), byrow = TRUE), jac_numeric, tolerance = 1E-3)
  jactape <- scorematchingad:::tape_Jacobian(cADFun)
  jactapeeval <- jactape$eval(unclass(omparovec), vector(mode = "numeric"))
  expect_equal(jactapeeval, jac)

})




