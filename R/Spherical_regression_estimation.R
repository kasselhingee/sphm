#-----------------------------------------------------------------------
# FLEXIBLE SPHERICAL REGRESSION USING A SCALED MOBIUS TRANSFORMATION
#    BY SHOGO KATO, ANDREW T.A. WOOD AND JANICE L. SCEALY
#
# R CODE FOR:
#    PARAMETER ESTIMATION OF SPHERICAL REGRESSION MODEL 
#    USING A PRELIMINARY ESTIMATOR IN SECTION 4
#-----------------------------------------------------------------------

# We consider the case p=q=3 of our regression model

library(scatterplot3d)

# objective function
pre_est3_mod=function(y,x,theta){
  
  p12=theta[1]
  p13=theta[2]
  p23=theta[3]
  q12=theta[4]
  q13=theta[5]
  q23=theta[6]
  b1=theta[7]
  b2=theta[8]
  
  p=matrix(c(0,p12,p13,-p12,0,p23,-p13,-p23,0),3,3)
  q=matrix(c(0,q12,q13,-q12,0,q23,-q13,-q23,0),3,3)
  
  P=(diag(1,3)-p)%*%solve(diag(1,3)+p)
  Q=(diag(1,3)-q)%*%solve(diag(1,3)+q)
  B=b1*diag(c(1,b2))
  
  value=0
  for(j in 1:(dim(y)[2])) value=value+sum(y[,j]*mu(x[,j],P,Q,B))
  return(-value)
}


d=3  # dimension
P=Q=diag(1,d)  # d*d matrix
B=matrix(c(0.9,0,0,0.2),d-1)  # diagonal matrix
Omega=P[,2:d]%*%B%*%t(Q[,2:d])
e1=c(1,rep(0,d-1))

n=100
u=rnorm(d*n)
x=y=matrix(1,d,n)
for(j in 1:n){
  x[,j]=u[d*(j-1)+1:d]/vnorm(u[d*(j-1)+1:d])  # uniform distribution on the sphere
  y[,j]=mu(x[,j],P,Q,B)  # assume y=mu(x) for simulation purpose
}

ini_value=c(0,0,0,0,0,0,0.9,0.2)  # initial value for optimization. Try multiple values.
result=nlminb(start=ini_value,objective=function(theta) pre_est3_mod(y,x,theta),lower=c(rep(-Inf,6),0,0),upper=c(rep(Inf,6),1,1))  # numerical optimization

theta=result[[1]]

p12=theta[1]
p13=theta[2]
p23=theta[3]
q12=theta[4]
q13=theta[5]
q23=theta[6]
b1=theta[7]
b2=theta[8]

p=matrix(c(0,p12,p13,-p12,0,p23,-p13,-p23,0),3,3)
q=matrix(c(0,q12,q13,-q12,0,q23,-q13,-q23,0),3,3)

Ph=(diag(1,3)-p)%*%solve(diag(1,3)+p)  # Cayley transform of hat{P}
Qh=(diag(1,3)-q)%*%solve(diag(1,3)+q)  # Cayley transform of hat{Q}
Bh=b1*diag(c(1,b2))
Omegah=Ph[,2:d]%*%Bh%*%t(Qh[,2:d])

# plot of x, y and residuals.
res=y
for(j in 1:n) res[,j]=y[,j]-mu(x[,j],Ph,Qh,Bh)  # residual: y-mu(x)

res_max1=max(res[1,]); res_min1=min(res[1,])
res_max2=max(res[2,]); res_min2=min(res[2,])
res_max3=max(res[3,]); res_min3=min(res[3,])

par(mfrow=c(1,3))
scatterplot3d(x[1,],x[2,],x[3,],xlim=c(-1,1),ylim=c(-1,1),zlim=c(-1,1),color="black")
scatterplot3d(y[1,],y[2,],y[3,],xlim=c(-1,1),ylim=c(-1,1),zlim=c(-1,1),color="blue")
scatterplot3d(res[1,],res[2,],res[3,],xlim=c(res_min1,res_max1),ylim=c(res_min2,res_max2),zlim=c(res_min3,res_max3),color="red")


ini_value  # initial value for optimization
P 
Ph  # hat{P}
Q
Qh  # hat{Q}
B
Bh  # hat{B}
Omega-Omegah  # Omega-hat{Omega}
sum(res[,]^2)  # squared sum of residuals
