#' Simulate the SmVF Distribution
rSvMF <- function(n, param, check = TRUE)
{
  param <- as_SvMFmuV(param)
  SvMFmuV_check(param)
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



	return(y)
}


##simKent simulates a value of y from the Kent distribution.
##kappa,beta,mu and K are the values of the parameters
simKent=function(n, kappa,beta,mu,K)
{
  p <- length(mu)
  skappa=matrix(kappa,n,1)
  
  sbeta=matrix(0,n,sum(p,-1))
  for (j in 2:sum(p,-1))
  {
    sbeta[,j]=beta[j-1]
  }
  
  #simulate sample
  
  
  cum=1
  rej=0
  rej2=0
  
  zs=matrix(0,n,p)
  zc=matrix(0,1,p)
  sig=matrix(0,n,p)
  sbetas=matrix(0,n,1)
  
  
  for (i in 2:sum(p,-1))
  {
    sbetas=sbetas+sbeta[,i]
    for (j in 1:n)
    {
      if (sbeta[j,i] > 0 ) {sig[j,i]=sqrt(skappa[j]-2*sbeta[j,i])}
      else {sig[j,i]=sqrt(skappa[j])}
      if (sbetas[j] < 0 ) {sig[j,p]=sqrt(skappa[j]+2*sbetas[j])}
      else {sig[j,p]=sqrt(skappa[j])}
      
    }
  }
  
  
  
  
  while (cum <= n)
    
  {
    
    
    for (i in 2:p)
    {
      zc[1,i]=rexp(1,rate=sig[cum,i])
    }
    
    vc=sum(zc^2)/4
    
    if (vc < 1)
    {
      
      
      r=runif(1, min=0, max=1)
      
      bz=0
      for (i in 2:sum(p,-1))
      {
        bz=bz + sbeta[cum,i]*zc[1,i]^2
      }
      ez=0
      for (i in 2:p)
      {
        ez=ez + sig[cum,i]*zc[1,i]
      }
      
      paccept=exp(((p-3)/2)*log(1-vc)-2*vc*skappa[cum]+(1-vc)*(bz-sbetas[cum]*zc[1,p]^2)+ ez-((p-1)/2))
      
      
      if (r < paccept) 
      {
        r2=runif(p, min=0, max=1)
        for (i in 2:p)
        {
          if (r2[i]< 0.5) {zc[1,i]=-1*zc[1,i]}
          zs[cum,i]=zc[1,i]
        }
        
        cum=cum+1
      }
      else
      {
        rej2=rej2+1	
      }
      
      
    }
    else
    {
      rej=rej+1	
    }
    
    
  }
  
  reject<<-rej
  ###
  H=diag(1,p)
  H[,1]=t(t(mu))
  H[1,]=t(mu)
  mu_L=t(t(mu[2:p]))
  H[2:p,2:p]=(1/(1+H[1,1]))*mu_L%*%t(mu_L)-diag(1,sum(p,-1))
  
  Gamma=H%*%K
  
  vs=matrix(0,n,1)
  ys=matrix(0,n,p)
  y=matrix(0,n,p)
  fold=0
  for (i in 1:n)
  {
    vs[i]=sum(zs[i,]^2)/4
    ys[i,1]=1-2*vs[i]
    ys[i,2:p]=((1-vs[i])^(0.5))*zs[i,2:p]
    y[i,]=ys[i,]%*%t(Gamma)
    fold_count=0
    ##count folding
    for (j in 1:p)
    {
      if (y[i,j] < 0) {fold_count=1}
    }
    fold=fold+fold_count
  }
  fold=fold/n
  
  return(list(y=y,fold=fold,ystd=ys))
}

