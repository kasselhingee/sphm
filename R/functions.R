
########################################################################
##simKent simulates a value of y from the Kent distribution.
##kappa,beta,mu and K are the values of the parameters
############################################################################

simKent=function(kappa,beta,mu,K)
{
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
	#############################
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


########################################################################
##simProject simulates a value of y from the SvMF distribution.
##kappa,V,mu,a1 are the values of the parameters
############################################################################

simProject=function(kappa,V,mu,a1)
{

	kappav=kappa
	av=a1

	betav=matrix(0,p-2,1)
	betav[1]=0*kappav
	muv=matrix(0,p,1)
	muv[1]=1
	K=diag(p)
	#simulated sample from von mises fisher distribution:
	sims=simKent(kappav,betav,muv,K)
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



########################################################################
##momu generates the moment estimate of mu.
##y is the spherical data value.
############################################################################
momu=function(y)
{
	mu_est=matrix(0,p,1)
	S_orig=matrix(0,p,p)
	for (i in 1:n)
	{
		mu_est=mu_est+y[i,]
		S_orig=S_orig+t(t(y[i,]))%*%t(y[i,])
	}
	ybar=mu_est/n
	mag=sqrt(sum(ybar^2))
	S_orig=S_orig/n
	
	return(list(ybar=ybar,mag=mag,S=S_orig))
}

########################################################################
##momK generates the moment estimate of K given mu.
##y is the spherical sample, mu is the current value of mu.
##S and ybar are needed and they are computed in the function momu above.
############################################################################
momK=function(y,mu,S,ybar)
{
	
	H=diag(1,p)
	H[,1]=t(t(mu))
	H[1,]=t(mu)
	mu_L=t(t(mu[2:p]))
	H[2:p,2:p]=(1/(1+H[1,1]))*mu_L%*%t(mu_L)-diag(1,sum(p,-1))
	
	mag=t(H)%*%ybar
	mag=mag[1]
	B=t(H)%*%S%*%H
	spec=eigen(B[2:p,2:p])
	Ks=spec$vectors
	lan=spec$values
	K=matrix(0,p,p)
	K[1,1]=1
	K[2:p,2:p]=Ks
	return(list(K=K,lan=lan,mag=mag))
}



#matrix trace
tra=function(s)
{
	sum(diag(s))
}

########################################################################
##saddlep generates the saddlepoint approximation of kappa and beta given mu and K.
##lan are the eigenvalues and mag is the magnitude of ybar.
##kappa and beta are the initial values.
##f1 and f2 are the saddlepoint functions needed within.
############################################################################
saddlep=function(lan,mag,kappa,beta)
{
	
	f1 <- function(t)
	{
		s1=0
		for (j in 2:sum(p,-1))
		{
			s1=s1+1/(k-2*b2[j-1]-2*t)
		}
		s1=s1+1/(k+2*sum(b2)-2*t)+1/(k-2*t)+(k/(k-2*t))^2-1
		s1
	}
	
	f2 <- function(s)
	{
		kappaest=s[1]
		betaest=s[2:sum(p,-1)]
	

	
		k<<-kappaest
		b2<<-betaest
		int1=(k-2*b2[1])/2-(1/4)*p-(1/2)*(0.25*p^2+p*k^2)^(0.5)
		int2=(k-2*b2[1])/2-0.5
	
		#int1=-(1/4)*p-(1/2)*(0.25*p^2+p*k^2)^(0.5)
		#int2=(k)/2-0.5
		#int1=-1000
		#int2=-4
	
	
		t_est=uniroot(f1, lower=int1,upper=int2, maxiter = 10000000)$root
		#,tol = 1e-10
		t_sav<<-t_est
	
		gan=matrix(kappaest/2,p,1)
		for (j in 2:sum(p,-1))
		{
			gan[j]=(kappaest-2*betaest[j-1])/2
		}
		gan[p]=(kappaest+2*sum(betaest))/2
		temp2=1-t_est/gan
		temp=matrix(0,p,1)
		for (j in 1:p)
		{
			temp[j]=ifelse(temp2[j] > 0, -0.5*log(temp2[j]), Inf) 
		}
		Kt=sum(temp)+kappaest^2/(2*kappaest-4*t_est)-kappaest/2
		temp=0.5*1/(gan-t_est)^2
		K2t=sum(temp)+0.5*kappaest^2/(kappaest/2-t_est)^3
		temp=1/(gan-t_est)^3
		K3t=sum(temp)+(6/4)*kappaest^2/(kappaest/2-t_est)^4
		temp=3/(gan-t_est)^4
		K4t=sum(temp)+6*kappaest^2/(kappaest/2-t_est)^5
		TT=(1/8)*K4t/K2t^2-(5/24)*(K3t/K2t^(3/2))^2
		lf1=-0.5*log(2*pi*K2t)+Kt-t_est+TT
		temp2=gan
		temp=matrix(0,p,1)
		for (j in 1:p)
		{
			temp[j]=ifelse(temp2[j] > 0, -0.5*log(temp2[j]), Inf) 
		}
		lg=lf1+kappaest+sum(temp)+log(2)+(p/2)*log(pi)
		ss=-1*sum(betaest)*lan[p-1]
		for (j in 2:sum(p,-1))
		{
			ss=ss+betaest[j-1]*lan[j-1]
		}
		logL=ss+kappaest*mag-lg
		lg<<-lg
		-logL
	}


	v=optim(par=c(kappa,beta),fn=f2,control=list(maxit=10000000,parscale=c(kappa,beta),reltol=1e-20))
	v
}

#######################################################
##muML3 updates mu given the other parameters (ML for Kent distribution).
##mu is the initial value.
##this function assumes y_1 >>0.
#########################################################

muML3=function(ybar,mu,K,kappa,beta,S_orig)
{
	fnL=function(muL)
	{

		mu=matrix(0,p,1)
		mu[2:p]=muL
		mu[1]=sqrt(1-sum(muL^2))
	
		mu_est=mu
		H=diag(1,p)
		H[,1]=t(t(mu_est))
		H[1,]=t(mu_est)
		mu_L=t(t(mu_est[2:p]))
		H[2:p,2:p]=(1/(1+H[1,1]))*mu_L%*%t(mu_L)-diag(1,sum(p,-1))
		

	
		Dc=matrix(0,p,p)
		for (m in 2:sum(p-1))
		{
			Dc[m,m]=beta[m-1]
		}
		Dc[p,p]=-1*sum(beta)
		Sig=K%*%Dc%*%t(K)

		-kappa*t(mu)%*%ybar-tra(H%*%Sig%*%t(H)%*%S_orig)


	}

	v=optim(par=mu[2:p],fn=fnL,control=list(maxit=10000000,parscale=mu[2:p],reltol=1e-20))
	v
}



############################################################################
##MLV: update V for SvMF given current kappa and mu 
##ypn is the unstandardised SvMF response
##a1 is the tail-weight parameter
############################################################################



MLV=function(ypn,mu,kappav,a1)
{

	av=a1
	##remove mu
	H=diag(1,p)
	H[,1]=t(t(mu))
	H[1,]=t(mu)
	mu_L=t(t(mu[2:p]))
	H[2:p,2:p]=(1/(1+H[1,1]))*mu_L%*%t(mu_L)-diag(1,sum(p,-1))
	ypss=matrix(0,n,p)
	weight=matrix(0,n,1)
	V=diag(p-1)

	conv3=0
	while (conv3==0)
	{

		V_old=V

		S_orig=matrix(0,p,p)
		ss=matrix(0,n,1)
		for (i in 1:n)
		{
			ypss[i,]=ypn[i,]%*%H
			ss[i]=(1/av^2)*ypss[i,1]^2+ypss[i,2:p]%*%solve(V)%*%t(t(ypss[i,2:p]))
			weight[i]=(p-1)*(1/ss[i])+(1/av)*kappav*ypss[i,1]*(ss[i])^(-3/2)
			S_orig=S_orig+weight[i]*t(t(ypss[i,]))%*%t(ypss[i,])
		}
		S_orig=S_orig/n

		##preliminary estimate of V
		V=S_orig[2:p,2:p]*(1/det(S_orig[2:p,2:p])^(1/(p-1)))

		if (max(abs(V_old-V)) < 0.000001) {conv3=1}

	}

	return(list(V=V))
}



############################################################################
##MLmu update mu for SvMF given current kappa and V 
##ypn is the unstandardised SvMF response
###a1 is the tail-weight parameter
###mu1 is the inital value of mu
############################################################################

MLmu=function(ypn,V,kappav,a1,mu1)
{
	mu=mu1
	av=a1
	fnL=function(muL)
	{

		mu=matrix(0,p,1)
		mu[2:p]=muL
		mu[1]=sqrt(1-sum(muL^2))
	
		mu_est=mu
		H=diag(1,p)
		H[,1]=t(t(mu_est))
		H[1,]=t(mu_est)
		mu_L=t(t(mu_est[2:p]))
		H[2:p,2:p]=(1/(1+H[1,1]))*mu_L%*%t(mu_L)-diag(1,sum(p,-1))
		Vdel=matrix(0,p,p)
		Vdel[1,1]=1/av^2
		Vdel[2:p,2:p]=solve(V)
	
		lang=0
		ss=matrix(0,n,1)
		for (i in 1:n)
		{
			ss[i]=t(ypn[i,])%*%H%*%Vdel%*%t(H)%*%t(t(ypn[i,]))
			lang=lang +1*((p-1)/2)*log(ss[i])-(1/av)*kappav*t(ypn[i,])%*%mu*(ss[i])^(-1/2)
		}
	
	
		lang

	}

	v=optim(par=mu[2:p],fn=fnL,control=list(maxit=10000000,parscale=mu[2:p],reltol=1e-20))
	mu=matrix(0,p,1)
	mu[2:p]=v$par
	mu[1]=sqrt(1-sum(mu[2:p]^2))
	mu
}



############################################################################
##MLkappa update kappa for SvMF given current mu and V 
##ypn is the unstandardised SvMF response
##a1 is the tail-weight parameter
############################################################################



MLkappa=function(ypn,mu,V,a1)
{
	av=a1
	fnL=function(kappa)
	{
		mu_est=mu
		H=diag(1,p)
		H[,1]=t(t(mu_est))
		H[1,]=t(mu_est)
		mu_L=t(t(mu_est[2:p]))
		H[2:p,2:p]=(1/(1+H[1,1]))*mu_L%*%t(mu_L)-diag(1,sum(p,-1))
		Vdel=matrix(0,p,p)
		Vdel[1,1]=1/av^2
		Vdel[2:p,2:p]=solve(V)
	
		sum_weight=0
		ss=matrix(0,n,1)
		for (i in 1:n)
		{
			ss[i]=t(ypn[i,])%*%H%*%Vdel%*%t(H)%*%t(t(ypn[i,]))
			sum_weight=sum_weight+t(ypn[i,])%*%mu*(ss[i])^(-1/2)
		}
		nu=sum(p/2,-1)
		loglike=kappa*((1/av)*sum_weight-n)+n*nu*log(kappa)-n*log(besselI(kappa, nu, expon.scaled = TRUE))
		-1*loglike
	}

	v=optim(par=kappa,fn=fnL,method="Brent", lower = 0.01, upper =100000,control=list(maxit=10000000,reltol=1e-20))
	
	kappa=v$par
	kappa
}


############################################################################
##SvMFmodel fits the new SvMF IID model 
##ypn is the unstandardised SvMF response
##a1 is the tuning parameter which contols the tail-weight
############################################################################



SvMFmodel=function(ypn,a1)
{

	av=a1
	
	mu=mom_mu
	V=diag(p-1)

	kappav_old=0
	conv2=0
	while (conv2==0)
	{

		#update kappa given V and mu

		kappav=MLkappa(ypn,mu,V,a1)
		kappav

		#update V given mu and kappa
		V=MLV(ypn,mu,kappav,a1)$V
		V

		if (max(abs(kappav-kappav_old)) <  0.0001) {conv2=1}
		kappav_old=kappav

	}

	
	mu1=mu_spatial
	mu=MLmu(ypn,V,kappav,a1,mu1)
	

	return(list(mu=mu,V=V,kappav=kappav))

}

############################################################################
##Kentmodel fits the Kent IID model 
##ypn is the unstandardised Kent response
##saddle indicates whether to use the saddlepoint approximation
## for shape parameters (saddle=1 represents yes)
############################################################################



Kentmodel=function(ypn,saddle)
{


	initial=momu(ypn)
	S_orig=initial$S
	ybar=initial$ybar
	#initial mu estimate:
	mu_est=mom_mu
	####Estimate K:
	initialK=momK(ypn,mu_est,S_orig,ybar)
	K_est=initialK$K
	##eigenvalues and magnitude
	lan=initialK$lan
	mag=initialK$mag

	####Estimate kappa and beta:
	##first order approximation:
	kappaest=0
	for (j in 1:sum(p,-1))
	{
		kappaest=kappaest+(1/lan[j])
	}
	kappaest=kappaest/(p-1)
	betaest=matrix(0,sum(p,-2))
	for (j in 1:sum(p,-2))
	{
		betaest[j]=0.5*(kappaest-(1/lan[j]))
	}
	kappa_saddle=kappaest
	beta_saddle=betaest
	if (saddle == 1)
	{
		###saddlepoint approximation:
		saddle=saddlep(lan,mag,kappaest,betaest)
		kappa_saddle=saddle$par[1]
		beta_saddle=betaest
		beta_saddle[1:sum(p,-2)]=saddle$par[2:sum(p,-1)]
	}
	


	conv2=0
	while (conv2==0)
	{

		
		###re-estimate mu:
 		mu_est_update=matrix(0,p,1)
		ff=muML3(ybar,mu_est,K_est,kappa_saddle,beta_saddle,S_orig)
		mu_est_update[2:p]=ff$par
		mu_est_update[1]=sqrt(1-sum(ff$par^2))

	


		if (max(abs(mu_est_update-mu_est)) <  0.0000001) {conv2=1}
		mu_est=mu_est_update
	
		####Estimate K:
		initialK=momK(ypn,mu_est,S_orig,ybar)
		K_est=initialK$K
		
		
	
	}
	

	return(list(mu=mu_est,K=K_est,kappa=kappa_saddle,beta=beta_saddle))

}



######################Below are additional functions needed for regression


########################################################################
##mu_calc computes the current estimates of the conditional mean
##direction for each observation.
##co is the current estimate of the regression coefficients.
#######################################################################

mu_calc=function(co)
{
	mu_est=matrix(0,n,p)
	lin=matrix(0,n,sum(p,-1))
	for (k in 1:sum(p,-1))
	{
		
		for (j in qcum2[k]:qcum[k])
		{
			
			lin[,k]=co[j]*x[,j]+lin[,k]

		}

	}
	linsum=matrix(1,n,1)
	for (k in 1:sum(p,-1))
	{
		linsum=linsum+exp(lin[,k])	
	}
	mu_est[,1]=(linsum)^(-1/2)
	for (k in 1:sum(p,-1))
	{
		mu_est[,sum(k,1)]=mu_est[,1]*exp((lin[,k])/2)	
	}
	mu_est
}



#bessel functions and derivatives needed for normalising constant
g=function(kappa,nu)
{
	besselI(kappa, nu, expon.scaled = TRUE)
}

g1=function(kappa,nu)
{
	((nu/kappa))*g(kappa,nu)+ g(kappa,sum(nu,1))
}

g2=function(kappa,nu)
{
	g(kappa,nu)*((nu-1)*nu/(kappa^2))+g(kappa,sum(nu,1))*((2*nu+1)/kappa)+g(kappa,sum(nu,2))
}



############################################################################
##est_mean_projected computes SvMF ML estimates of the regression coefficients a.
##co is the current (or initial) value of the regression coefficients
##sigma3,c1 and delta4 are the current values of parameters in V
##m and delta3 are the current values of parameters in kappa
##tol1 and tol2 are tolerance values needed in the Newton-Raphson algorithm
##a1 is the tuning parameter which contols the tail-weight
#############################################################################



est_mean_projected=function(co,sigma3,c1,delta4,m,delta3,tol1,tol2,a1)
{
	av=a1

	Vdel=matrix(0,p,p)
	Vdel[1,1]=1/av^2


	conv=1
	while (conv==1)
	{
		mu_est=mu_calc(co)
		sbeta=matrix(0,Q,1)
		Jbeta=matrix(0,Q,Q)

		for (i in 1:n)
		{
			Vdel[2,2]=(vx[i]^(-2*delta4))/(sigma3*(1-c1^2)^0.5) #vx is the `depth` variable in examples, that is a covariate that isn't passed as a function argument
			Vdel[2,3]=-c1/(1-c1^2)^0.5
			Vdel[3,2]=Vdel[2,3]
			Vdel[3,3]=sigma3*(vx[i]^(2*delta4))/((1-c1^2)^(0.5))
			
			kappa=(m^(-1))*(vx[i]^(-2*delta3))

			

			H=diag(1,p)
			H[,1]=t(t(mu_est[i,]))
			H[1,]=t(mu_est[i,])
			mu_L=t(t(mu_est[i,2:p]))
			H[2:p,2:p]=(1/(1+H[1,1]))*mu_L%*%t(mu_L)-diag(1,sum(p,-1))
					

			ss=t(y[i,])%*%H%*%Vdel%*%t(H)%*%t(t(y[i,]))
			w1=(ss)^(-1/2)
			w2=-((p-1)/2)*ss^(-1)-(kappa/(2*av))*(ss^(-3/2))*t(H[,1])%*%t(t(y[i,]))
			
		

                        # prepare matrices s related to the score and J and JS related to the Jacobian of the score
			s=matrix(0,p-1,1)
			J=matrix(0,p-1,p-1)
			JS=matrix(0,p-1,p-1)
			for (j in 1:sum(p,-1))
			{
				for (mf in 1:sum(p,-1))
				{
					
				
					Hd1_a=diag(1,p)
					Hd1_b=diag(1,p)
					Hd2=diag(1,p)
				
					
					for (k in 1:p)
					{
						Hd1_a[k,1]=-1*H[k,1]*(H[sum(j,1),1])^2
						cj=sum(j,1)
						if (k==cj){Hd1_a[k,1]=Hd1_a[k,1]+H[sum(j,1),1]}
						Hd1_b[k,1]=-1*H[k,1]*(H[sum(mf,1),1])^2
						cm=sum(mf,1)
						if (k==cm){Hd1_b[k,1]=Hd1_b[k,1]+H[sum(mf,1),1]}
						delta1=0
						if (k==cj){delta1=1}
						delta2=0
						if (j==mf){delta2=1}
						sdelta3=0
						if (k==cm){sdelta3=1}
						
						Hd2[k,1]=(delta1-2*H[k,1]*H[sum(j,1),1])*(H[sum(mf,1),1]*delta2-H[sum(j,1),1]*H[sum(mf,1),1]^2)-H[sum(j,1),1]^2*(H[sum(mf,1),1]*sdelta3-H[k,1]*H[sum(mf,1),1]^2)
					
						
					}
					Hd1_a[1,1:p]=Hd1_a[1:p,1]
					Hd1_b[1,1:p]=Hd1_b[1:p,1]
					Hd2[1,1:p]=Hd2[1:p,1]
					
					for (k in 2:p)
					{
						
						for (r in 2:p)
						{
							Hd1_a[k,r]=Hd1_a[k,1]*H[r,1]/(1+H[1,1])+H[k,1]*Hd1_a[r,1]/(1+H[1,1])-H[k,1]*H[r,1]*(1+H[1,1])^(-2)*Hd1_a[1,1]
							
							Hd1_b[k,r]=Hd1_b[k,1]*H[r,1]/(1+H[1,1])+H[k,1]*Hd1_b[r,1]/(1+H[1,1])-H[k,1]*H[r,1]*(1+H[1,1])^(-2)*Hd1_b[1,1]
							
							Hd2[k,r]=Hd2[k,1]*H[r,1]*(1+H[1,1])^(-1)+Hd1_a[k,1]*Hd1_b[r,1]*(1+H[1,1])^(-1)-Hd1_a[k,1]*H[r,1]*(1+H[1,1])^(-2)*Hd1_b[1,1]+Hd1_b[k,1]*Hd1_a[r,1]*(1+H[1,1])^(-1)+H[k,1]*Hd2[r,1]*(1+H[1,1])^(-1)-H[k,1]*Hd1_a[r,1]*(1+H[1,1])^(-2)*Hd1_b[1,1]-Hd1_b[k,1]*H[r,1]*(1+H[1,1])^(-2)*Hd1_a[1,1]-H[k,1]*Hd1_b[r,1]*(1+H[1,1])^(-2)*Hd1_a[1,1]+2*H[k,1]*H[r,1]*(1+H[1,1])^(-3)*Hd1_b[1,1]*Hd1_a[1,1]-H[k,1]*H[r,1]*(1+H[1,1])^(-2)*Hd2[1,1]
								
						}	
						
					}
					
				
					J[j,mf]=(w1*(kappa/av)*t(Hd2[,1])%*%t(t(y[i,]))+w2*t(y[i,])%*%Hd2%*%Vdel%*%t(H)%*%t(t(y[i,]))+w2*t(y[i,])%*%Hd1_a%*%Vdel%*%t(Hd1_b)%*%t(t(y[i,]))+w2*t(y[i,])%*%Hd1_b%*%Vdel%*%t(Hd1_a)%*%t(t(y[i,]))+w2*t(y[i,])%*%H%*%Vdel%*%t(Hd2)%*%t(t(y[i,])))
					derb=(t(y[i,])%*%Hd1_b%*%Vdel%*%t(H)%*%t(t(y[i,]))+t(y[i,])%*%H%*%Vdel%*%t(Hd1_b)%*%t(t(y[i,])))
					w1b=-0.5*(ss^(-3/2))*derb
					w2b=((p-1)/2)*(ss^(-2))*derb+(3/4)*(kappa/av)*(ss^(-5/2))*derb*t(H[,1])%*%t(t(y[i,]))-(kappa/(2*av))*(ss^(-3/2))*t(Hd1_b[,1])%*%t(t(y[i,]))
					JS[j,mf]=(w1b*(kappa/av)*t(Hd1_a[,1])%*%t(t(y[i,])) +w2b*t(y[i,])%*%Hd1_a%*%Vdel%*%t(H)%*%t(t(y[i,])) +w2b*t(y[i,])%*%H%*%Vdel%*%t(Hd1_a)%*%t(t(y[i,])) )
				
				}
				s[j]=(w1*(kappa/av)*t(Hd1_a[,1])%*%t(t(y[i,])) +w2*t(y[i,])%*%Hd1_a%*%Vdel%*%t(H)%*%t(t(y[i,])) +w2*t(y[i,])%*%H%*%Vdel%*%t(Hd1_a)%*%t(t(y[i,])) )
				

			}
		
                        # compute the relevant score function sbeta at measurement i
			ds=matrix(0,Q,Q)
			for (k in 1:sum(p-1))
			{
				for (j in qcum2[k]:qcum[k])
				{
					ds[j,j]=s[k]
					
				}
			}
			#sbeta=sbeta+ds%*%t(t(x[i,]))
			sbeta=sbeta+0.5*ds%*%t(t(x[i,]))

                        # compute the derivative of the score function Jbeta
			comp=t(t(x[i,]))%*%t(x[i,])
			Jtemp=matrix(0,Q,Q)
			for (k in 1:sum(p,-1))
			{
				for (j in 1:sum(p,-1))
				{
					Jtemp[qcum2[k]:qcum[k],qcum2[j]:qcum[j]]=J[k,j]+JS[k,j]

				}

			}
			Jbeta=Jbeta+0.25*comp*Jtemp
			#Jbeta=Jbeta+comp*Jtemp

		}
		#print(co)
		cprev=co
		co=co-solve(Jbeta)%*%sbeta #this is the Newton-Raphson step
		dec=max(abs(co-cprev)/(abs(co)+tol1))
		#conv=0
		if (dec < tol2){conv=0}
	
	}
	info=-Jbeta
	return(list(a=co,info=info,sbeta=sbeta))
	
}



############################################################################
##est_V_kappa computes SvMF ML estimates of V and kappa.
##res is the current residual vector res=H^ty
##sigma3,c1 and delta4 are the initial values of parameters in V
##m and delta3 are the initial values of parameters in kappa
##tol1 and tol2 are tolerance values needed in the Newton-Raphson algorithm
##a1 is the tuning parameter which contols the tail-weight
#############################################################################



est_V_kappa=function(res,sigma3,c1,delta4,m,delta3,tol1,tol2,a1)
{

	av=a1

	par_old=matrix(0,5,1)
	par_old[1]=sigma3
	par_old[2]=c1
	par_old[3]=delta4
	par_old[4]=m
	par_old[5]=delta3

	
	conv=1
	while (conv==1)
	{
		
		w=0

	
		Vdel=matrix(0,p,p)
		Vdel[1,1]=1/av^2
		
		ss=matrix(0,n,1)
		f_sigma3=0
		f_c1=0
		f_delta4=0
		J_sigma3_c1=0
		J_sigma3_delta4=0
		J_c1_delta4=0
		J_c1_c1=0
		J_sigma3_sigma3=0
		J_delta4_delta4=0

	

		f_m=0
		f_delta3=0
		J_m_delta3=0
		J_m_m=0
		J_delta3_delta3=0

		J_m_sigma3=0
		J_m_c1=0
		J_m_delta4=0

		J_delta3_sigma3=0
		J_delta3_c1=0
		J_delta3_delta4=0
	

		for (i in 1:n)
		{
			Vdel[2,2]=(vx[i]^(-2*delta4))/(sigma3*(1-c1^2)^0.5)
			Vdel[2,3]=-c1/(1-c1^2)^0.5
			Vdel[3,2]=Vdel[2,3]
			Vdel[3,3]=sigma3*(vx[i]^(2*delta4))/((1-c1^2)^(0.5))
			ss[i]=t(res[i,])%*%Vdel%*%t(t(res[i,]))
			w=(p-1)*(1/ss[i])+(1/av)*(m^(-1))*(vx[i]^(-2*delta3))*res[i,1]*(ss[i])^(-3/2)
			#w=1
		
			Vin_sigma3=matrix(0,2,2)
			Vin_sigma3[1,1]=-1*(vx[i]^(-2*delta4))/((sigma3^2)*(1-c1^2)^0.5)
			Vin_sigma3[2,2]=(vx[i]^(2*delta4))/((1-c1^2)^0.5)

			Vin_c1=matrix(0,2,2)
			Vin_c1[1,1]=((1-c1^2)^(-1.5))*c1*(vx[i]^(-2*delta4))/(sigma3)
			Vin_c1[2,2]=((1-c1^2)^(-1.5))*c1*sigma3*(vx[i]^(2*delta4))
			Vin_c1[1,2]=-1*((1-c1^2)^(-1.5))*c1^2-((1-c1^2)^(-0.5))
			Vin_c1[2,1]=Vin_c1[1,2]
		
			Vin_delta4=matrix(0,2,2)
			Vin_delta4[1,1]=-2*log(vx[i])*(vx[i]^(-2*delta4))/(sigma3*(1-c1^2)^0.5)
			Vin_delta4[2,2]=2*log(vx[i])*sigma3*(vx[i]^(2*delta4))/((1-c1^2)^(0.5))

			f_sigma3=f_sigma3+w*t(res[i,2:3])%*%Vin_sigma3%*%t(t(res[i,2:3]))
			f_c1=f_c1+w*t(res[i,2:3])%*%Vin_c1%*%t(t(res[i,2:3]))
			f_delta4=f_delta4+w*t(res[i,2:3])%*%Vin_delta4%*%t(t(res[i,2:3]))


			w2=-1*(p-1)*(1/ss[i]^2)-1.5*(1/av)*(m^(-1))*(vx[i]^(-2*delta3))*res[i,1]*(ss[i])^(-5/2)
			#w2=0

			Vin_sigma3_c1=matrix(0,2,2)
			Vin_sigma3_c1[1,1]=-c1*(vx[i]^(-2*delta4))/((sigma3^2)*(1-c1^2)^1.5)
			Vin_sigma3_c1[2,2]=c1*(vx[i]^(2*delta4))/((1-c1^2)^1.5)

			J_sigma3_c1=J_sigma3_c1+w2*t(res[i,2:3])%*%Vin_sigma3%*%t(t(res[i,2:3]))*t(res[i,2:3])%*%Vin_c1%*%t(t(res[i,2:3]))+w*t(res[i,2:3])%*%Vin_sigma3_c1%*%t(t(res[i,2:3]))		
	
			Vin_sigma3_delta4=matrix(0,2,2)
			Vin_sigma3_delta4[1,1]=2*log(vx[i])*(vx[i]^(-2*delta4))/((sigma3^2)*(1-c1^2)^0.5)
			Vin_sigma3_delta4[2,2]=2*log(vx[i])*(vx[i]^(2*delta4))/((1-c1^2)^0.5)

			J_sigma3_delta4=J_sigma3_delta4+w2*t(res[i,2:3])%*%Vin_sigma3%*%t(t(res[i,2:3]))*t(res[i,2:3])%*%Vin_delta4%*%t(t(res[i,2:3]))+w*t(res[i,2:3])%*%Vin_sigma3_delta4%*%t(t(res[i,2:3]))		
	
			Vin_c1_delta4=matrix(0,2,2)
			Vin_c1_delta4[1,1]=-2*log(vx[i])*c1*(vx[i]^(-2*delta4))/((sigma3)*(1-c1^2)^1.5)
			Vin_c1_delta4[2,2]=2*log(vx[i])*sigma3*c1*(vx[i]^(2*delta4))/((1-c1^2)^1.5)

			J_c1_delta4=J_c1_delta4+w2*t(res[i,2:3])%*%Vin_c1%*%t(t(res[i,2:3]))*t(res[i,2:3])%*%Vin_delta4%*%t(t(res[i,2:3]))+w*t(res[i,2:3])%*%Vin_c1_delta4%*%t(t(res[i,2:3]))		
	
			#need other 3 derivatives

			Vin_c1_c1=matrix(0,2,2)
			Vin_c1_c1[1,1]=((1-c1^2)^(-1.5))*(vx[i]^(-2*delta4))/(sigma3)+3*(c1^2)*((1-c1^2)^(-2.5))*(vx[i]^(-2*delta4))/(sigma3)
			Vin_c1_c1[2,2]=((1-c1^2)^(-1.5))*(vx[i]^(2*delta4))*(sigma3)+3*(c1^2)*((1-c1^2)^(-2.5))*(vx[i]^(2*delta4))*(sigma3)
			Vin_c1_c1[1,2]=-3*((1-c1^2)^(-1.5))*c1-3*(c1^3)*((1-c1^2)^(-2.5))
			Vin_c1_c1[2,1]=Vin_c1_c1[1,2]

			J_c1_c1=J_c1_c1+w2*t(res[i,2:3])%*%Vin_c1%*%t(t(res[i,2:3]))*t(res[i,2:3])%*%Vin_c1%*%t(t(res[i,2:3]))+w*t(res[i,2:3])%*%Vin_c1_c1%*%t(t(res[i,2:3]))		
	
		
			Vin_sigma3_sigma3=matrix(0,2,2)
			Vin_sigma3_sigma3[1,1]=2*(vx[i]^(-2*delta4))/((sigma3^3)*(1-c1^2)^0.5)
		
			J_sigma3_sigma3=J_sigma3_sigma3+w2*t(res[i,2:3])%*%Vin_sigma3%*%t(t(res[i,2:3]))*t(res[i,2:3])%*%Vin_sigma3%*%t(t(res[i,2:3]))+w*t(res[i,2:3])%*%Vin_sigma3_sigma3%*%t(t(res[i,2:3]))		
		
			Vin_delta4_delta4=matrix(0,2,2)
			Vin_delta4_delta4[1,1]=((2*log(vx[i]))^2)*(vx[i]^(-2*delta4))/(sigma3*(1-c1^2)^0.5)
			Vin_delta4_delta4[2,2]=((2*log(vx[i]))^2)*sigma3*(vx[i]^(2*delta4))/((1-c1^2)^(0.5))

			J_delta4_delta4=J_delta4_delta4+w2*t(res[i,2:3])%*%Vin_delta4%*%t(t(res[i,2:3]))*t(res[i,2:3])%*%Vin_delta4%*%t(t(res[i,2:3]))+w*t(res[i,2:3])%*%Vin_delta4_delta4%*%t(t(res[i,2:3]))		
		
			w3=(1/av)*res[i,1]*(ss[i])^(-1/2)
		
			kappa=(m^(-1))*(vx[i]^(-2*delta3))
			kappa1=-(m^(-2))*(vx[i]^(-2*delta3))
			f_m=f_m+kappa1*w3-kappa1*g1(kappa,nu)/g(kappa,nu)+kappa1*nu/kappa

			kappa=(m^(-1))*(vx[i]^(-2*delta3))
			kappa1=-2*log(vx[i])*(m^(-1))*(vx[i]^(-2*delta3))
			f_delta3=f_delta3+kappa1*w3-kappa1*g1(kappa,nu)/g(kappa,nu)+kappa1*nu/kappa

			kappa=(m^(-1))*(vx[i]^(-2*delta3))
			kappa1=-(m^(-2))*(vx[i]^(-2*delta3))
			kappa2=2*(m^(-3))*(vx[i]^(-2*delta3))
			J_m_m=J_m_m+kappa2*w3-kappa1*kappa1*g2(kappa,nu)/g(kappa,nu)+kappa1*kappa1*g1(kappa,nu)*g1(kappa,nu)/(g(kappa,nu))^2-kappa2*g1(kappa,nu)/g(kappa,nu)+kappa2*nu/kappa-kappa1*kappa1*nu/(kappa^2)
	
			kappa=(m^(-1))*(vx[i]^(-2*delta3))
			kappa1=-2*log(vx[i])*(m^(-1))*(vx[i]^(-2*delta3))
			kappa2=((2*log(vx[i]))^2)*(m^(-1))*(vx[i]^(-2*delta3))
			J_delta3_delta3=J_delta3_delta3+kappa2*w3-kappa1*kappa1*g2(kappa,nu)/g(kappa,nu)+kappa1*kappa1*g1(kappa,nu)*g1(kappa,nu)/(g(kappa,nu))^2-kappa2*g1(kappa,nu)/g(kappa,nu)+kappa2*nu/kappa-kappa1*kappa1*nu/(kappa^2)
	

			kappa=(m^(-1))*(vx[i]^(-2*delta3))
			kappa1a=-(m^(-2))*(vx[i]^(-2*delta3))
			kappa1b=-2*log(vx[i])*(m^(-1))*(vx[i]^(-2*delta3))
			kappa2=2*log(vx[i])*(m^(-2))*(vx[i]^(-2*delta3))
			J_m_delta3=J_m_delta3+kappa2*w3-kappa1a*kappa1b*g2(kappa,nu)/g(kappa,nu)+kappa1a*kappa1b*g1(kappa,nu)*g1(kappa,nu)/(g(kappa,nu))^2-kappa2*g1(kappa,nu)/g(kappa,nu)+kappa2*nu/kappa-kappa1a*kappa1b*nu/(kappa^2)
	
			kappa=(m^(-1))*(vx[i]^(-2*delta3))
			kappa1=-(m^(-2))*(vx[i]^(-2*delta3))
			J_m_sigma3=J_m_sigma3-0.5*(1/av)*kappa1*res[i,1]*t(res[i,2:3])%*%Vin_sigma3%*%t(t(res[i,2:3]))*(ss[i])^(-3/2)
			J_m_c1=J_m_c1-0.5*(1/av)*kappa1*res[i,1]*t(res[i,2:3])%*%Vin_c1%*%t(t(res[i,2:3]))*(ss[i])^(-3/2)
			J_m_delta4=J_m_delta4-0.5*(1/av)*kappa1*res[i,1]*t(res[i,2:3])%*%Vin_delta4%*%t(t(res[i,2:3]))*(ss[i])^(-3/2)
	
			kappa=(m^(-1))*(vx[i]^(-2*delta3))
			kappa1=-2*log(vx[i])*(m^(-1))*(vx[i]^(-2*delta3))
			J_delta3_sigma3=J_delta3_sigma3-0.5*(1/av)*kappa1*res[i,1]*t(res[i,2:3])%*%Vin_sigma3%*%t(t(res[i,2:3]))*(ss[i])^(-3/2)
			J_delta3_c1=J_delta3_c1-0.5*(1/av)*kappa1*res[i,1]*t(res[i,2:3])%*%Vin_c1%*%t(t(res[i,2:3]))*(ss[i])^(-3/2)
			J_delta3_delta4=J_delta3_delta4-0.5*(1/av)*kappa1*res[i,1]*t(res[i,2:3])%*%Vin_delta4%*%t(t(res[i,2:3]))*(ss[i])^(-3/2)
	
		}
		
	

		J=matrix(0,5,5)
		J[1,1]=J_sigma3_sigma3
		J[1,2]=J_sigma3_c1
		J[1,3]=J_sigma3_delta4
		J[2,2]=J_c1_c1
		J[2,3]=J_c1_delta4
		J[3,3]=J_delta4_delta4
		J[2,1]=J[1,2]
		J[3,1]=J[1,3]
		J[3,2]=J[2,3]
		J[4,4]=J_m_m
		J[4,5]=J_m_delta3
		J[5,5]=J_delta3_delta3
		J[5,4]=J[4,5]
		J[1,4]=J_m_sigma3
		J[4,1]=J[1,4]
		J[2,4]=J_m_c1
		J[4,2]=J[2,4]
		J[3,4]=J_m_delta4
		J[4,3]=J[3,4]
		J[1,5]=J_delta3_sigma3
		J[5,1]=J[1,5]
		J[2,5]=J_delta3_c1
		J[5,2]=J[2,5]
		J[3,5]=J_delta3_delta4
		J[5,3]=J[3,5]


		Jin=solve(J)

		f=matrix(0,5,1)
		f[1]=f_sigma3
		f[2]=f_c1
		f[3]=f_delta4
		f[4]=f_m
		f[5]=f_delta3
	

		par_new=par_old-Jin%*%f
		#print(par_old)
	
		dec=max(abs(par_new-par_old)/(abs(par_old)+tol1))
		#conv=0
		if (dec < tol2){conv=0}


		par_old=par_new
		sigma3=par_old[1]
		c1=par_old[2]
		delta4=par_old[3]
		m=par_old[4]
		delta3=par_old[5]




	}

	info=-J
	return(list(sigma3=sigma3,c1=c1,delta4=delta4,m=m,delta3=delta3,info=info,f=f))

}


############################################################################
##SvMFreg fits the SvMF regression model.
##co is the initial value of the regression coefficients
##sigma3,c1 and delta4 are the initial values of parameters in V
##m and delta3 are the initial values of parameters in kappa
##tol1 and tol2 are tolerance values
##a1 is the tuning parameter which contols the tail-weight
############################################################################



SvMFreg=function(co,sigma3,c1,delta4,m,delta3,tol1,tol2,a1)
{

	conv=1
	while (conv==1)
	{
		cprev=co
		###update regression coefficients:
		mm=est_mean_projected(co,sigma3,c1,delta4,m,delta3,0.00001,0.000001,a1)
		mu_est=mu_calc(mm$a)
		#print(mm$a)
		info=mm$info


		###update residuals:
		ynew2=matrix(0,n,p)
		for (i in 1:n)
		{

			#Calculating Gamma for the given observation
			H=diag(1,p)
			H[,1]=t(t(mu_est[i,]))
			H[1,]=t(mu_est[i,])
			mu_L=t(t(mu_est[i,2:p]))
			H[2:p,2:p]=(1/(1+H[1,1]))*mu_L%*%t(mu_L)-diag(1,sum(p,-1))

			ynew2[i,]=y[i,]%*%H


		}
		res=ynew2


		##update V and kappa
		mm3=est_V_kappa(res,sigma3,c1,delta4,m,delta3,0.00001,0.000001,a1)
		sigma3=mm3$sigma3
		c1=mm3$c1
		delta4=mm3$delta4
		m=mm3$m
		delta3=mm3$delta3
		#print(sigma3)

		co=mm$a
		dec=max(abs(co-cprev)/(abs(co)+tol1))
		#conv=0
		if (dec < tol2){conv=0}



	}

		

	return(list(a=co,sigma3=sigma3,c1=c1,delta4=delta4,m=m,delta3=delta3,info=info,res=res))

}


############################################################################
##est_mean_kent computes Kent ML estimates of the regression coefficients.
##co is the current (or initial) value of the regression coefficients
##sigma3,c1 and delta4 are the current values of parameters in V
##m and delta3 are the current values of parameters in kappa
##tol1 and tol2 are tolerance values needed in the Newton-Raphson algorithm
#############################################################################



est_mean_kent=function(co,sigma3,c1,delta4,m,delta3,tol1,tol2)
{

	sigma4=m

	conv=1
	while (conv==1)
	{
		mu_est=mu_calc(co)
		sbeta=matrix(0,Q,1)
		Jbeta=matrix(0,Q,Q)

		for (i in 1:n)
		{
			V=matrix(0,2,2)
			V[1,1]=sigma3*(vx[i]^(2*delta4))/((1-c1^2)^0.5)
			V[1,2]=c1/(1-c1^2)^0.5
			V[2,1]=V[1,2]
			V[2,2]=(vx[i]^(-2*delta4))/(sigma3*(1-c1^2)^0.5)
			V=V*(vx[i]^(2*delta3))*sigma4
			e=eigen(V)
			lan=e$values
			
			Ke=e$vectors
			K=diag(1,p)
			K[2:p,2:p]=Ke


			kappa=0
			for (j in 1:sum(p,-1))
			{
				kappa=kappa+(1/lan[j])
			}
			kappa=kappa/(p-1)
			beta=matrix(0,sum(p,-2))
			for (j in 1:sum(p,-2))
			{
				beta[j]=0.5*(kappa-(1/lan[j]))
			}

			

			C=matrix(0,p,p)
			for(j in 2:sum(p,-1))
			{
				C[j,j]=beta[j-1]
			}
			C[p,p]=-sum(beta)
			
			s=matrix(0,p-1,1)
			J=matrix(0,p-1,p-1)
			for (j in 1:sum(p,-1))
			{
				for (mf in 1:sum(p,-1))
				{
					
					H=diag(1,p)
					Hd1_a=diag(1,p)
					Hd1_b=diag(1,p)
					Hd2=diag(1,p)
				
					H[,1]=t(t(mu_est[i,]))
					H[1,]=t(mu_est[i,])
					mu_L=t(t(mu_est[i,2:p]))
					H[2:p,2:p]=(1/(1+H[1,1]))*mu_L%*%t(mu_L)-diag(1,sum(p,-1))
					
					for (k in 1:p)
					{
						Hd1_a[k,1]=-1*H[k,1]*(H[sum(j,1),1])^2
						cj=sum(j,1)
						if (k==cj){Hd1_a[k,1]=Hd1_a[k,1]+H[sum(j,1),1]}
						Hd1_b[k,1]=-1*H[k,1]*(H[sum(mf,1),1])^2
						cm=sum(mf,1)
						if (k==cm){Hd1_b[k,1]=Hd1_b[k,1]+H[sum(mf,1),1]}
						delta1=0
						if (k==cj){delta1=1}
						delta2=0
						if (j==mf){delta2=1}
						delta3=0
						if (k==cm){delta3=1}
						
						Hd2[k,1]=(delta1-2*H[k,1]*H[sum(j,1),1])*(H[sum(mf,1),1]*delta2-H[sum(j,1),1]*H[sum(mf,1),1]^2)-H[sum(j,1),1]^2*(H[sum(mf,1),1]*delta3-H[k,1]*H[sum(mf,1),1]^2)
					
						
					}
					Hd1_a[1,1:p]=Hd1_a[1:p,1]
					Hd1_b[1,1:p]=Hd1_b[1:p,1]
					Hd2[1,1:p]=Hd2[1:p,1]
					
					for (k in 2:p)
					{
						
						for (r in 2:p)
						{
							Hd1_a[k,r]=Hd1_a[k,1]*H[r,1]/(1+H[1,1])+H[k,1]*Hd1_a[r,1]/(1+H[1,1])-H[k,1]*H[r,1]*(1+H[1,1])^(-2)*Hd1_a[1,1]
							
							Hd1_b[k,r]=Hd1_b[k,1]*H[r,1]/(1+H[1,1])+H[k,1]*Hd1_b[r,1]/(1+H[1,1])-H[k,1]*H[r,1]*(1+H[1,1])^(-2)*Hd1_b[1,1]
							
							Hd2[k,r]=Hd2[k,1]*H[r,1]*(1+H[1,1])^(-1)+Hd1_a[k,1]*Hd1_b[r,1]*(1+H[1,1])^(-1)-Hd1_a[k,1]*H[r,1]*(1+H[1,1])^(-2)*Hd1_b[1,1]+Hd1_b[k,1]*Hd1_a[r,1]*(1+H[1,1])^(-1)+H[k,1]*Hd2[r,1]*(1+H[1,1])^(-1)-H[k,1]*Hd1_a[r,1]*(1+H[1,1])^(-2)*Hd1_b[1,1]-Hd1_b[k,1]*H[r,1]*(1+H[1,1])^(-2)*Hd1_a[1,1]-H[k,1]*Hd1_b[r,1]*(1+H[1,1])^(-2)*Hd1_a[1,1]+2*H[k,1]*H[r,1]*(1+H[1,1])^(-3)*Hd1_b[1,1]*Hd1_a[1,1]-H[k,1]*H[r,1]*(1+H[1,1])^(-2)*Hd2[1,1]
								
						}	
						
					}
					
				
					J[j,mf]=(kappa*t(Hd2[,1])%*%t(t(y[i,]))+t(y[i,])%*%Hd2%*%K%*%C%*%t(K)%*%t(H)%*%t(t(y[i,]))+t(y[i,])%*%Hd1_a%*%K%*%C%*%t(K)%*%t(Hd1_b)%*%t(t(y[i,]))+t(y[i,])%*%Hd1_b%*%K%*%C%*%t(K)%*%t(Hd1_a)%*%t(t(y[i,]))+t(y[i,])%*%H%*%K%*%C%*%t(K)%*%t(Hd2)%*%t(t(y[i,])))

				
				}
				s[j]=(kappa*t(Hd1_a[,1])%*%t(t(y[i,])) +t(y[i,])%*%Hd1_a%*%K%*%C%*%t(K)%*%t(H)%*%t(t(y[i,])) +t(y[i,])%*%H%*%K%*%C%*%t(K)%*%t(Hd1_a)%*%t(t(y[i,])) )
				

			}
			
			ds=matrix(0,Q,Q)
			for (k in 1:sum(p-1))
			{
				for (j in qcum2[k]:qcum[k])
				{
					ds[j,j]=s[k]
					
				}
			}
			sbeta=sbeta+0.5*ds%*%t(t(x[i,]))
			comp=t(t(x[i,]))%*%t(x[i,])
			Jtemp=matrix(0,Q,Q)
			for (k in 1:sum(p,-1))
			{
				for (j in 1:sum(p,-1))
				{
					Jtemp[qcum2[k]:qcum[k],qcum2[j]:qcum[j]]=J[k,j]

				}

			}
			Jbeta=Jbeta+0.25*comp*Jtemp

		}
		#print(co)
		cprev=co
		co=co-solve(Jbeta)%*%sbeta
		dec=max(abs(co-cprev)/(abs(co)+tol1))
		#conv=0
		if (dec < tol2){conv=0}
	
	}
	info=-Jbeta
	return(list(a=co,info=info))
	

	
}


############################################################################
##est_V_kappa_kent computes approximate ML estimates of V and kappa.
##res is the current residual vector res=H^ty
##sigma3,c1 and delta4 are the initial values of parameters in V
##m and delta3 are the initial values of parameters in kappa
##tol1 and tol2 are tolerance values needed in the Newton-Raphson algorithm
#############################################################################



est_V_kappa_kent=function(res,sigma3,c1,delta4,m,delta3,tol1,tol2)
{



	par_old=matrix(0,5,1)
	par_old[1]=sigma3
	par_old[2]=c1
	par_old[3]=delta4
	par_old[4]=m
	par_old[5]=delta3

	
	conv=1
	while (conv==1)
	{
		
		w=0

	
		Vdel=matrix(0,p,p)
		Vdel[1,1]=1
		
		ss=matrix(0,n,1)
		f_sigma3=0
		f_c1=0
		f_delta4=0
		J_sigma3_c1=0
		J_sigma3_delta4=0
		J_c1_delta4=0
		J_c1_c1=0
		J_sigma3_sigma3=0
		J_delta4_delta4=0

	

		f_m=0
		f_delta3=0
		J_m_delta3=0
		J_m_m=0
		J_delta3_delta3=0

		J_m_sigma3=0
		J_m_c1=0
		J_m_delta4=0

		J_delta3_sigma3=0
		J_delta3_c1=0
		J_delta3_delta4=0
	

		for (i in 1:n)
		{
			Vdel[2,2]=(vx[i]^(-2*delta4))/(sigma3*(1-c1^2)^0.5)
			Vdel[2,3]=-c1/(1-c1^2)^0.5
			Vdel[3,2]=Vdel[2,3]
			Vdel[3,3]=sigma3*(vx[i]^(2*delta4))/((1-c1^2)^(0.5))
			ss[i]=t(res[i,])%*%Vdel%*%t(t(res[i,]))
			w=(m^(-1))*(vx[i]^(-2*delta3))
			
		
			Vin_sigma3=matrix(0,2,2)
			Vin_sigma3[1,1]=-1*(vx[i]^(-2*delta4))/((sigma3^2)*(1-c1^2)^0.5)
			Vin_sigma3[2,2]=(vx[i]^(2*delta4))/((1-c1^2)^0.5)

			Vin_c1=matrix(0,2,2)
			Vin_c1[1,1]=((1-c1^2)^(-1.5))*c1*(vx[i]^(-2*delta4))/(sigma3)
			Vin_c1[2,2]=((1-c1^2)^(-1.5))*c1*sigma3*(vx[i]^(2*delta4))
			Vin_c1[1,2]=-1*((1-c1^2)^(-1.5))*c1^2-((1-c1^2)^(-0.5))
			Vin_c1[2,1]=Vin_c1[1,2]
		
			Vin_delta4=matrix(0,2,2)
			Vin_delta4[1,1]=-2*log(vx[i])*(vx[i]^(-2*delta4))/(sigma3*(1-c1^2)^0.5)
			Vin_delta4[2,2]=2*log(vx[i])*sigma3*(vx[i]^(2*delta4))/((1-c1^2)^(0.5))

			f_sigma3=f_sigma3+w*t(res[i,2:3])%*%Vin_sigma3%*%t(t(res[i,2:3]))
			f_c1=f_c1+w*t(res[i,2:3])%*%Vin_c1%*%t(t(res[i,2:3]))
			f_delta4=f_delta4+w*t(res[i,2:3])%*%Vin_delta4%*%t(t(res[i,2:3]))


			w2=0

			Vin_sigma3_c1=matrix(0,2,2)
			Vin_sigma3_c1[1,1]=-c1*(vx[i]^(-2*delta4))/((sigma3^2)*(1-c1^2)^1.5)
			Vin_sigma3_c1[2,2]=c1*(vx[i]^(2*delta4))/((1-c1^2)^1.5)

			J_sigma3_c1=J_sigma3_c1+w2*t(res[i,2:3])%*%Vin_sigma3%*%t(t(res[i,2:3]))*t(res[i,2:3])%*%Vin_c1%*%t(t(res[i,2:3]))+w*t(res[i,2:3])%*%Vin_sigma3_c1%*%t(t(res[i,2:3]))		
	
			Vin_sigma3_delta4=matrix(0,2,2)
			Vin_sigma3_delta4[1,1]=2*log(vx[i])*(vx[i]^(-2*delta4))/((sigma3^2)*(1-c1^2)^0.5)
			Vin_sigma3_delta4[2,2]=2*log(vx[i])*(vx[i]^(2*delta4))/((1-c1^2)^0.5)

			J_sigma3_delta4=J_sigma3_delta4+w2*t(res[i,2:3])%*%Vin_sigma3%*%t(t(res[i,2:3]))*t(res[i,2:3])%*%Vin_delta4%*%t(t(res[i,2:3]))+w*t(res[i,2:3])%*%Vin_sigma3_delta4%*%t(t(res[i,2:3]))		
	
			Vin_c1_delta4=matrix(0,2,2)
			Vin_c1_delta4[1,1]=-2*log(vx[i])*c1*(vx[i]^(-2*delta4))/((sigma3)*(1-c1^2)^1.5)
			Vin_c1_delta4[2,2]=2*log(vx[i])*sigma3*c1*(vx[i]^(2*delta4))/((1-c1^2)^1.5)

			J_c1_delta4=J_c1_delta4+w2*t(res[i,2:3])%*%Vin_c1%*%t(t(res[i,2:3]))*t(res[i,2:3])%*%Vin_delta4%*%t(t(res[i,2:3]))+w*t(res[i,2:3])%*%Vin_c1_delta4%*%t(t(res[i,2:3]))		
	
			#need other 3 derivatives

			Vin_c1_c1=matrix(0,2,2)
			Vin_c1_c1[1,1]=((1-c1^2)^(-1.5))*(vx[i]^(-2*delta4))/(sigma3)+3*(c1^2)*((1-c1^2)^(-2.5))*(vx[i]^(-2*delta4))/(sigma3)
			Vin_c1_c1[2,2]=((1-c1^2)^(-1.5))*(vx[i]^(2*delta4))*(sigma3)+3*(c1^2)*((1-c1^2)^(-2.5))*(vx[i]^(2*delta4))*(sigma3)
			Vin_c1_c1[1,2]=-3*((1-c1^2)^(-1.5))*c1-3*(c1^3)*((1-c1^2)^(-2.5))
			Vin_c1_c1[2,1]=Vin_c1_c1[1,2]

			J_c1_c1=J_c1_c1+w2*t(res[i,2:3])%*%Vin_c1%*%t(t(res[i,2:3]))*t(res[i,2:3])%*%Vin_c1%*%t(t(res[i,2:3]))+w*t(res[i,2:3])%*%Vin_c1_c1%*%t(t(res[i,2:3]))		
	
		
			Vin_sigma3_sigma3=matrix(0,2,2)
			Vin_sigma3_sigma3[1,1]=2*(vx[i]^(-2*delta4))/((sigma3^3)*(1-c1^2)^0.5)
		
			J_sigma3_sigma3=J_sigma3_sigma3+w2*t(res[i,2:3])%*%Vin_sigma3%*%t(t(res[i,2:3]))*t(res[i,2:3])%*%Vin_sigma3%*%t(t(res[i,2:3]))+w*t(res[i,2:3])%*%Vin_sigma3_sigma3%*%t(t(res[i,2:3]))		
		
			Vin_delta4_delta4=matrix(0,2,2)
			Vin_delta4_delta4[1,1]=((2*log(vx[i]))^2)*(vx[i]^(-2*delta4))/(sigma3*(1-c1^2)^0.5)
			Vin_delta4_delta4[2,2]=((2*log(vx[i]))^2)*sigma3*(vx[i]^(2*delta4))/((1-c1^2)^(0.5))

			J_delta4_delta4=J_delta4_delta4+w2*t(res[i,2:3])%*%Vin_delta4%*%t(t(res[i,2:3]))*t(res[i,2:3])%*%Vin_delta4%*%t(t(res[i,2:3]))+w*t(res[i,2:3])%*%Vin_delta4_delta4%*%t(t(res[i,2:3]))		
		
			w3=t(res[i,2:p])%*%Vdel[2:p,2:p]%*%t(t(res[i,2:p]))
		
			kappa=(m^(-1))*(vx[i]^(-2*delta3))
			kappa1=-(m^(-2))*(vx[i]^(-2*delta3))
			f_m=f_m+kappa1*w3+(2/m)

			kappa=(m^(-1))*(vx[i]^(-2*delta3))
			kappa1=-2*log(vx[i])*(m^(-1))*(vx[i]^(-2*delta3))
			f_delta3=f_delta3+kappa1*w3+4*log(vx[i])

			kappa=(m^(-1))*(vx[i]^(-2*delta3))
			kappa1=-(m^(-2))*(vx[i]^(-2*delta3))
			kappa2=2*(m^(-3))*(vx[i]^(-2*delta3))
			J_m_m=J_m_m+kappa2*w3-2/(m^2)
	
			kappa=(m^(-1))*(vx[i]^(-2*delta3))
			kappa1=-2*log(vx[i])*(m^(-1))*(vx[i]^(-2*delta3))
			kappa2=((2*log(vx[i]))^2)*(m^(-1))*(vx[i]^(-2*delta3))
			J_delta3_delta3=J_delta3_delta3+kappa2*w3
	

			kappa=(m^(-1))*(vx[i]^(-2*delta3))
			kappa1a=-(m^(-2))*(vx[i]^(-2*delta3))
			kappa1b=-2*log(vx[i])*(m^(-1))*(vx[i]^(-2*delta3))
			kappa2=2*log(vx[i])*(m^(-2))*(vx[i]^(-2*delta3))
			J_m_delta3=J_m_delta3+kappa2*w3
	
			kappa=(m^(-1))*(vx[i]^(-2*delta3))
			kappa1=-(m^(-2))*(vx[i]^(-2*delta3))
			J_m_sigma3=J_m_sigma3+kappa1*t(res[i,2:3])%*%Vin_sigma3%*%t(t(res[i,2:3]))
			J_m_c1=J_m_c1+kappa1*t(res[i,2:3])%*%Vin_c1%*%t(t(res[i,2:3]))
			J_m_delta4=J_m_delta4+kappa1*t(res[i,2:3])%*%Vin_delta4%*%t(t(res[i,2:3]))
	
			kappa=(m^(-1))*(vx[i]^(-2*delta3))
			kappa1=-2*log(vx[i])*(m^(-1))*(vx[i]^(-2*delta3))
			J_delta3_sigma3=J_delta3_sigma3+kappa1*t(res[i,2:3])%*%Vin_sigma3%*%t(t(res[i,2:3]))
			J_delta3_c1=J_delta3_c1+kappa1*t(res[i,2:3])%*%Vin_c1%*%t(t(res[i,2:3]))
			J_delta3_delta4=J_delta3_delta4+kappa1*t(res[i,2:3])%*%Vin_delta4%*%t(t(res[i,2:3]))
	
		}
		
	

		J=matrix(0,5,5)
		J[1,1]=J_sigma3_sigma3
		J[1,2]=J_sigma3_c1
		J[1,3]=J_sigma3_delta4
		J[2,2]=J_c1_c1
		J[2,3]=J_c1_delta4
		J[3,3]=J_delta4_delta4
		J[2,1]=J[1,2]
		J[3,1]=J[1,3]
		J[3,2]=J[2,3]
		J[4,4]=J_m_m
		J[4,5]=J_m_delta3
		J[5,5]=J_delta3_delta3
		J[5,4]=J[4,5]
		J[1,4]=J_m_sigma3
		J[4,1]=J[1,4]
		J[2,4]=J_m_c1
		J[4,2]=J[2,4]
		J[3,4]=J_m_delta4
		J[4,3]=J[3,4]
		J[1,5]=J_delta3_sigma3
		J[5,1]=J[1,5]
		J[2,5]=J_delta3_c1
		J[5,2]=J[2,5]
		J[3,5]=J_delta3_delta4
		J[5,3]=J[3,5]


		Jin=solve(J)

		f=matrix(0,5,1)
		f[1]=f_sigma3
		f[2]=f_c1
		f[3]=f_delta4
		f[4]=f_m
		f[5]=f_delta3
	

		par_new=par_old-Jin%*%f
		#print(par_old)
	
		dec=max(abs(par_new-par_old)/(abs(par_old)+tol1))
		#conv=0
		if (dec < tol2){conv=0}


		par_old=par_new
		sigma3=par_old[1]
		c1=par_old[2]
		delta4=par_old[3]
		m=par_old[4]
		delta3=par_old[5]




	}

	info=-J
	return(list(sigma3=sigma3,c1=c1,delta4=delta4,m=m,delta3=delta3,info=info,f=f))

}


############################################################################
##Kentreg fits the Kent regression model.
##co is the initial value of the regression coefficients
##sigma3,c1 and delta4 are the initial values of parameters in V
##m and delta3 are the initial values of parameters in kappa
##tol1 and tol2 are tolerance values
############################################################################



Kentreg=function(co,sigma3,c1,delta4,m,delta3,tol1,tol2)
{

	conv=1
	while (conv==1)
	{
		cprev=co
		###update regression coefficients:
		mm=est_mean_kent(co,sigma3,c1,delta4,m,delta3,0.00001,0.000001)
		mu_est=mu_calc(mm$a)
		#print(mm$a)
		info=mm$info


		###update residuals:
		ynew2=matrix(0,n,p)
		for (i in 1:n)
		{

			#Calculating Gamma for the given observation
			H=diag(1,p)
			H[,1]=t(t(mu_est[i,]))
			H[1,]=t(mu_est[i,])
			mu_L=t(t(mu_est[i,2:p]))
			H[2:p,2:p]=(1/(1+H[1,1]))*mu_L%*%t(mu_L)-diag(1,sum(p,-1))

			ynew2[i,]=y[i,]%*%H


		}
		res=ynew2


		##update V and kappa
		mm3=est_V_kappa_kent(res,sigma3,c1,delta4,m,delta3,0.00001,0.000001)
		sigma3=mm3$sigma3
		c1=mm3$c1
		delta4=mm3$delta4
		m=mm3$m
		delta3=mm3$delta3
		#print(sigma3)

		co=mm$a
		dec=max(abs(co-cprev)/(abs(co)+tol1))
		#conv=0
		if (dec < tol2){conv=0}



	}

		

	return(list(a=co,sigma3=sigma3,c1=c1,delta4=delta4,m=m,delta3=delta3,info=info,res=res))

}
