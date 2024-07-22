
test_that("optim_pobjS2S, pobjS2S() and pobjS2SCpp() works",{
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
  if (!requireNamespace("movMF", quietly = TRUE)){skip("Need movMF package")}
  set.seed(5)
  y <- t(apply(ymean, 1, function(mn){movMF::rmovMF(1, 10*mn)}))
  
  # objective function in C++ and R should match when omegapar passes OmegaS2S_check()
  objval <- pobjS2S(y, x, omegapar)
  objvalcpp <-  pobjS2Scpp(OmegaS2S_vec(omegapar), vector(), p, cbind(y,x))
  expect_equal(objvalcpp, objval)

  # optimise using pure R
  opt <- optim_pobjS2S_pureR(y, x, omegapar, global = TRUE, local = TRUE)
  expect_equal(OmegaS2S_proj(opt$solution, method = "Omega"), opt$solution, tolerance = 1E-3)
  expect_equal(opt$solution, omegapar, tolerance = 0.05)
  
  # optimise locally using derivative information
  # starting at the optimum
  tmp <- optim_pobjS2S_parttape(y, x, omegapar)
  expect_equal(tmp$solution, omegapar, tolerance = 0.05)
  
  # starting away from optimum (but also with large B)
  start <- OmegaS2S_proj(OmegaS2S_unvec(OmegaS2S_vec(omegapar) + 10, p, check = FALSE))
  opt2 <- optim_pobjS2S_parttape(y, x, start)
  if (sign(opt2$solution$p1[1]) == sign(omegapar$p1[1])){
    expect_equal(opt2$solution, omegapar, tolerance = 0.05)
  } else {
    # use an antipodal - like transform to get to another solution
    start2 <- as_cannS2S(opt2$solution)
    start2$B <- diag(1/diag(start2$B))
    start2$P <- -start2$P
    start2$Q <- -start2$Q
    opt3 <- optim_pobjS2S_parttape(y, x, as_OmegaS2S(start2))
    expect_equal(opt3$solution, omegapar, tolerance = 0.05)
    expect_equal(opt3$solution, opt$solution, tolerance = 1E-2)
    expect_equal(opt3$solution, tmp$solution, tolerance = 1E-2)
  }
})

test_that("OmegaS2S_constraints() is zero correctly", {
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
  expect_equal(OmegaS2S_constraints(OmegaS2S_vec(omegapar), p), rep(0, 1 + 1))
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

test_that("taping of pobjS2S and OmegaS2S_constraints runs and evaluates", {
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
  if (!requireNamespace("movMF", quietly = TRUE)){skip("Need movMF package")}
  set.seed(5)
  y <- t(apply(ymean, 1, function(mn){movMF::rmovMF(1, 10*mn)}))

  atapeptr <- pobjS2Stape(OmegaS2S_vec(omegapar), p, cbind(y,x))
  directeval <- pobjS2Scpp(OmegaS2S_vec(omegapar), vector(), p, cbind(y,x))
  tapeeval <- scorematchingad:::pForward0(atapeptr, unclass(OmegaS2S_vec(omegapar)), vector(mode = "numeric"))
  expect_equal(tapeeval, directeval)

  # if taping worked the results should match for a different value of the parameters
  omparo <- omegapar
  omparo$Omega <- omparo$Omega + 1
  omparovec <- OmegaS2S_vec(omparo)
  directeval <- pobjS2Scpp(omparovec, vector(), p, cbind(y,x))
  tapeeval <- scorematchingad:::pForward0(atapeptr, unclass(omparovec), vector(mode = "numeric"))
  expect_equal(tapeeval, directeval)

  # check derivatives: pForward, via taping, and numerically
  jac_numeric <- drop(attr(numericDeriv(quote(pobjS2Scpp(omparovec, vector(), p, cbind(y,x))), c("omparovec")), "gradient"))
  jac <- scorematchingad:::pJacobian(atapeptr, unclass(omparovec), vector(mode = "numeric"))
  expect_equal(jac, jac_numeric, tolerance = 1E-3)

  atape <- scorematchingad:::ADFun$new(ptr = atapeptr,
                   name = "pobjS2S",
                   xtape = unclass(OmegaS2S_vec(omegapar)), #unclass needed to pass the isa(, "numeric") check!
                   dyntape =  vector(mode = "numeric"),
                   usertheta = unclass(NA * OmegaS2S_vec(omegapar)))
  jactape <- scorematchingad::tapeJacobian(atape)
  jactapeeval <- scorematchingad:::pForward0(jactape$ptr, unclass(omparovec), vector(mode = "numeric"))
  expect_equal(jactapeeval, jac)

  # check constraints
  btapeptr <- OmegaS2S_constraintstape(OmegaS2S_vec(omegapar), p)
  directeval <- OmegaS2S_constraints(omparovec, p)
  tapeeval <- scorematchingad:::pForward0(btapeptr, unclass(omparovec), vector(mode = "numeric"))
  expect_equal(tapeeval, directeval)

  jac_numeric <- attr(numericDeriv(quote(scorematchingad:::pForward0(btapeptr, unclass(omparovec), vector(mode = "numeric"))), c("omparovec")), "gradient")
  jac <- scorematchingad:::pJacobian(btapeptr, unclass(omparovec), vector(mode = "numeric"))
  expect_equal(matrix(jac, nrow = length(tapeeval), byrow = TRUE), jac_numeric, tolerance = 1E-3)
  atape <- scorematchingad:::ADFun$new(ptr = btapeptr,
                   name = "constraint",
                   xtape = unclass(OmegaS2S_vec(omegapar)), #unclass needed to pass the isa(, "numeric") check!
                   dyntape =  vector(mode = "numeric"),
                   usertheta = unclass(NA * OmegaS2S_vec(omegapar)))
  jactape <- scorematchingad::tapeJacobian(atape)
  jactapeeval <- scorematchingad:::pForward0(jactape$ptr, unclass(omparovec), vector(mode = "numeric"))
  expect_equal(jactapeeval, jac)
})




