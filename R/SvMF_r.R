#' Simulate the SmVF Distribution
rSvMF <- function(n, param, check = TRUE)
{
  param <- as_SvMFmuV(param)
  SvMFmuV_check(param)
  browser()
  kappa <- param$k
  V <- param$V
  mu <- param$m
  a1 <- param$a1
  p <- length(mu)

	kappav=kappa
	av=a1

	betav=matrix(0,p-2,1)
	betav[1]=0*kappav
	muv=matrix(0,p,1)
	muv[1]=1
	K=diag(p)
	#simulated sample from von mises fisher distribution:
	sims=simKent(n, kappav,betav,muv,K)
	yv=sims$y

	a=eigen(V)$values
	a=rbind(av^2,t(t(a)))
	a=sqrt(a)
	
	yp=matrix(0,n,p)
	for (j in 1:p)
	{
		yp[,j]=yv[,j]*a[j]

	}

	sumsq=0
	for (j in 1:p)
	{
		sumsq=sumsq+yp[,j]^2
	}
	yp=yp/sqrt(sumsq)


	H=diag(1,p)
	H[,1]=t(t(mu))
	H[1,]=t(mu)
	mu_L=t(t(mu[2:p]))
	H[2:p,2:p]=(1/(1+H[1,1]))*mu_L%*%t(mu_L)-diag(1,sum(p,-1))

	K=diag(1,p)
	K[2:p,2:p]=eigen(V)$vectors

	Gamma=H%*%K

	ypn=matrix(0,n,p)
	fold=0
	for (i in 1:n)
	{
		ypn[i,]=yp[i,]%*%t(Gamma)
	

	}



	#response
	y=ypn



	return(list(y=y))
}

