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

# This is the normalising constant w.r.t. the Lebesgue measure on the sphere.
vMFnormconst <- function(k, p, method = 'base'){
  if (p == 3){return(2*pi*(exp(k) - exp(-k))/k)} #from Scealy and Wood 2019, this nice and simple for p = 3
    if (method == 'base'){
      # formula in Scealy and Wood 2019
      return((2*pi)^(p/2) * besselI(k, p/2 - 1) /k^(p/2 - 1))  #p/2-1 kind
    }
    if (method == 'Bessel'){
      requireNamespace("Bessel")
      return((2*pi)^(p/2) * Bessel::BesselI(k, p/2 - 1) /k^(p/2 - 1))  #p/2-1 kind
    }
    if (method == "movMF"){
      requireNamespace("movMF")
      # the normalising constant pops out from evaluating the vMF density at any point orthogonal to the mean
      # From Scealy and Wood vs Hornik and Grun, movMF seems to be calculating density wrt uniform distribution on the sphere
      # there is a constant ration between the two of gamma(p/2)/((2*pi)^(p/2))
      # Another 2^(p/2 - 1) is needed, but I'm not sure why - perhaps something to do with multiple clusters in movMF?
      const_wrtunif <- 
        movMF::dmovMF(matrix(c(0, 1, rep(0, p-2)), nrow = 1),
                      matrix(c(1, rep(0, p-1)), nrow = 1) * k,
                      log = FALSE)
      const_wrtLeb <- const_wrtunif * 2^(p/2 - 1) * gamma(p/2)/((2*pi)^(p/2))
      return(1/const_wrtLeb)
    }
}
