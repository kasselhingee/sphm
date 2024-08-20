#' Log-likelihood of SvMF
#' @param y is row-vectors of data
SvMF_ll_cann <- function(y, param){
  list2env(param, envir = environment())
  p <- nrow(G)
  if (p == 3){lconst <- -log(2*pi*(exp(k) - exp(-k))/k) - log(a[1])}
  if (p != 3){stop("vMF normalising constant for p != 3 implemented yet")}
  Gscal <- t(t(G)/a) #scale columns of Gamma by a
  denom <- sqrt(rowSums((y %*% Gscal)^2)) #the denominator in the density exponent
  
  ll <- lconst - (p-1) * log(denom) + drop(k * y %*% Gscal[,1])/denom
  return(ll)
}

#' How to calculate the vMF normalisation `c_p(kappa)`? When p=3, it simple:
#' `c_p(kappa) = 2*pi * (exp(kappa) - exp(-kappa))/kappa`. When p!=3 need to approximate c_p or use Score Matching.


