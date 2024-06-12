#' Preliminary objective function for S2S Link with p=q
#' @details Uses Cayley transform to parameterise P and Q. 
#' + Could be more accurate the closer P and Q are to the identity (check notes with Andy).
#' + Might be missing sign stuff to get negative determinants for Cayley transform
#' @export
pre_est3_mod=function(y,x,theta){
  
  b1=theta[7]
  b2=theta[8]
  
  P=cayley(theta[1:3])
  Q=cayley(theta[4:6])
  B=b1*diag(c(1,b2))
  
  value=0
  for(j in 1:(dim(y)[2])) value=value+sum(y[,j]*meanlinkS2S(x[,j],P,Q,B))
  return(-value)
}


