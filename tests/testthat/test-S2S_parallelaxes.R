test_that("maximum likelihood for constant parallel transport axes", {
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
  
  # generate noise according to homoskedastic SvMF
  k <- 30
  a <- c(1, seq(5, 0.2, length.out = p-1))
  a[-1] <- a[-1]/prod(a[-1])^(1/(p-1))
  # make matrix V
  set.seed(6)
  Kstar <- mclust::randomOrthogonalMatrix(p-1, p-1)
  V <- Kstar %*% diag(a[-1]^2) %*% t(Kstar)
  set.seed(5)
  y <- t(apply(ymean, 1, function(mn){
    rSvMF(1, SvMFmuV(k, mn, a[1], V))
  }))
  
  # rotated resids
  rresid <- rotatedresid(y, ymean, P[, 1])
  # shift P[, 1] to north pole to reduce computational issues
  rot <- JuppRmat(P[, 1], nthpole(p))
  rresid <- rresid %*% t(rot)
  
  # return residuals to the sphere!? Seems like a should be able to do something direcly using likelihood, except likelihood is on the sphere
  rresid[, 1] <- sqrt(1-rowSums(rresid[, -1]^2))
  
  ull_SvMF_V(TFORGE::vech(V), c(rresid[1, ], k), c(a[1], nthpole(p)))
  
  
  #Scealy optimisation of kappa using getH(nth pole) = diag(c(1, -1, -1, -1...))
  # use regression form to be able to use residuals rather than data on the sphere
  c1 <- V[1,2]/sqrt(1+V[1,2]^2)  #c1 / sqrt(1-c1^2) = V[1, 2] <=> c1 = v/sqrt(1+v^2)  [square both sides, multiply by denominator solve the square, require c1 to have same sign as V[1, 2]]
  n <- nrow(rresid)
  vx <- rep(1, n) #empty covariate
  #order of bessel functions
  nu=sum(p/2,-1)
  est_V_kappa(rresid %*% diag(c(1, rep(-1, p - 1))), # rresid %*% getH(P[, 1]) should be zero in first column
              sigma3 = V[1,1] * sqrt(1-c1^2), #s / sqrt(1-c1^2) = V[1,1] s = V[1,1] * sqrt(1-c1^2)
              c1 = c1,
              delta4 = 0,
              m = k, #corresponds to 1/sigma_4 in (14) I think
              delta3 = 0,
              tol1=0.00001,
              tol2=0.000001,
              a1 = a[1])
  
})
