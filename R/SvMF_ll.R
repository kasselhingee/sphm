#' Log-likelihood of SvMF
#' @param y is row-vectors of data
#' @param param is a set of parameters from either [`SvMFcann()`] or [`SvMFmuV()`].
#' @param log Return log density?
#' @export
#' @return A vector of values.
dSvMF <- function(y, param, log = FALSE){
  if (inherits(param, "SvMFmuV")){ll <- SvMF_ll_muV(y, param)}
  if (inherits(param, "SvMFcann")){ll <- SvMF_ll_cann(y, param)}
  if (log){return(ll)}
  else {return(exp(ll))}
}


SvMF_ll_cann <- function(y, param){
  list2env(param, envir = environment())
  p <- nrow(G)
  lconst <- -log(vMFnormconst(k, p)) - log(a[1]) #from Scealy and Wood 2019, this nice and simple for p = 3
  Gscal <- t(t(G)/a) #scale columns of Gamma by a
  denom <- sqrt(rowSums((y %*% Gscal)^2)) #the denominator in the density exponent
  
  ll <- lconst - (p-1) * log(denom) + drop(k * y %*% Gscal[,1])/denom
  return(ll)
}

SvMF_ll_muV <- function(y, param){
  list2env(param, envir = environment())
  p <- length(m)
  lconst <- -log(vMFnormconst(k, p)) - log(a1) #from Scealy and Wood 2019, this nice and simple for p = 3
  Hstar <- getHstar(m)
  ystarstarL <- y %*% Hstar
  denom <- sqrt(drop((y %*% m/a1)^2) + rowSums((ystarstarL %*% solve(V))*ystarstarL))

  ll <- lconst - (p-1) * log(denom) + drop(k * y %*% m/a1)/denom
  return(ll)
}

vMFnormconst <- function(k, p){
  if (p == 3){return(2*pi*(exp(k) - exp(-k))/k)} #from Scealy and Wood 2019, this nice and simple for p = 3
  if (p != 3){stop("vMF normalising constant for p != 3 implemented yet")}
}
