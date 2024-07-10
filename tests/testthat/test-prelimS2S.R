test_that("preobjS2S() works",{
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
  x <- matrix(rnorm(100*q), nrow = 100)
  x <- sweep(x, 1, apply(x, 1, vnorm), FUN = "/")
  
  y <- meanlinkS2S(x = x, paramobj = omegapar)
  
  # generate noise
  if (!requireNamespace("movMF", quietly = TRUE)){skip("Need movMF package")}
  ynoise <- t(apply(y, 1, function(mn){movMF::rmovMF(1, 10*mn)}))
  
  # evaluate objective function
  objval <- pobjS2S(ynoise, x, omegapar)
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
