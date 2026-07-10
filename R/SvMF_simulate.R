# This file was largely written by Janice Scealy.

#' Simulate the SvMF Distribution
#' @param n The number of samples to simulate
#' @param param A parameter object. See [`SvMF_cann()`] or [`SvMF_muV()`].
#' @param check If TRUE, will check that the parameters are appropriate for a SvMF model.
#' @family SvMF-distribution
#' @export
rSvMF <- function(n, param, check = TRUE)
{
  # Simulation strategy (Scealy & Wood 2019): draw from a vMF distribution at the north
  # pole with concentration k, scale each component by the corresponding a_j, renormalise
  # to the sphere, then rotate into the local frame given by G (= Gamma = H %*% K below).
  param <- as_SvMF_muV(param)
  if (check){SvMF_muV_check(param)}
  kappa <- param$k
  V <- param$V
  mu <- param$m   # mean direction
  a1 <- param$a1  # first scale
  p <- length(mu)

	kappav=kappa
	av=a1

	betav=matrix(0,p-2,1)
	betav[1]=0*kappav
	muv=matrix(0,p,1)
	muv[1]=1
	K=diag(p)
	# Draw n samples from the vMF distribution at the north pole (mu = e1), which serves
	# as the unscaled base draws before the SvMF scaling and rotation are applied.
	sims=sim_kent(n, kappav,betav,muv,K)
	yv=sims$y

	# Scales: a[1] = a1, a[2:p] from sqrt of eigenvalues of V
	a=eigen(V)$values
	a=rbind(av^2,t(t(a)))
	a=sqrt(a)

	# Scale each component of yv by a[j], then renormalise to the unit sphere.
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


	# Construct the local frame Gamma = H %*% K at the mean direction mu.
	# H is the stereographic-projection-based frame at mu (same construction as get_H() in
	# SvMF_parameterisations.R). K rotates the axes according to the eigenvectors of V.
	H=diag(1,p)
	H[,1]=t(t(mu))
	H[1,]=t(mu)
	mu_L=t(t(mu[2:p]))
	H[2:p,2:p]=(1/(1+H[1,1]))*mu_L%*%t(mu_L)-diag(1,sum(p,-1))

	K=diag(1,p)
	K[2:p,2:p]=eigen(V)$vectors

	Gamma=H%*%K

	# Rotate the scaled draws into the local frame at mu.
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


# Rejection sampler for the Kent (SvMF at a fixed mean) distribution. Used internally by
# rSvMF() to draw the unscaled base samples before the SvMF scaling and rotation are applied.
# The sampler draws candidate tangential components from exponential envelopes (rates `sig`)
# and accepts with probability `paccept`, which corrects for the envelope being approximate.
# NOTE: `reject <<-` writes to the enclosing environment (side-effect counter); this is not
# re-entrant and should not be called from multiple simultaneous contexts.
sim_kent=function(n, kappa,beta,mu,K)
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
  # sig[j,i]: rate parameter of the exponential envelope for component i of observation j.
  # Chosen so the exponential tail covers the Kent density in the tangential direction.
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

    # vc = ||zc||^2 / 4 is the squared "radius" of the candidate; reject if outside the sphere
    vc=sum(zc^2)/4

    if (vc < 1)
    {


      r=runif(1, min=0, max=1)

      bz=0
      for (i in 2:sum(p,-1))
      {
        bz=bz + sbeta[cum,i]*zc[1,i]^2
      }
      # ez = sum(sig[i] * zc[i]): the log-density of the exponential envelope at zc
      ez=0
      for (i in 2:p)
      {
        ez=ez + sig[cum,i]*zc[1,i]
      }

      # paccept = target density / envelope density (unnormalised acceptance probability)
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

  # <<- writes the rejection count to the enclosing environment (side-effect; not re-entrant)
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

