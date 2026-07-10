#' Log-likelihood of SvMF
#' @param y is row-vectors of data
#' @param param is a set of parameters from either [`SvMF_cann()`] or [`SvMF_muV()`].
#' @param log Return log density?
#' @family SvMF-distribution
#' @export
#' @return A vector of values.
dSvMF <- function(y, param, log = FALSE){
  if (inherits(param, "SvMF_muV")){ll <- SvMF_log_lik_muV(y, param)}
  if (inherits(param, "SvMF_cann")){ll <- SvMF_log_lik_cann(y, param)}
  if (log){return(ll)}
  else {return(exp(ll))}
}


SvMF_log_lik_cann <- function(y, param){
  list2env(param, envir = environment())
  p <- nrow(G)
  lconst <- -vMF_log_norm_const_exact(k, p) - log(a[1]) #from Scealy and Wood 2019, this nice and simple for p = 3
  # Gscal divides each column of G by the corresponding scale a_j, normalising each axis.
  Gscal <- t(t(G)/a) #scale columns of Gamma by a
  # denom = ||Gscal^T y|| is the effective radius of y in the scaled frame; it appears in
  # the log-likelihood formula of Scealy & Wood (2019) as the normalising factor.
  denom <- sqrt(rowSums((y %*% Gscal)^2)) #the denominator in the density exponent

  ll <- lconst - (p-1) * log(denom) + drop(k * y %*% Gscal[,1])/denom
  return(ll)
}

SvMF_log_lik_muV <- function(y, param){
  list2env(param, envir = environment())
  p <- length(m)
  lconst <- -vMF_log_norm_const_exact(k, p) - log(a1) #from Scealy and Wood 2019, this nice and simple for p = 3
  Hstar <- get_Hstar(m)
  ystarstarL <- y %*% Hstar
  denom <- sqrt(drop((y %*% m/a1)^2) + rowSums((ystarstarL %*% solve(V))*ystarstarL))

  ll <- lconst - (p-1) * log(denom) + drop(k * y %*% m/a1)/denom
  return(ll)
}

# This is the __log__ of the normalising constant w.r.t. the Lebesgue measure on the sphere.
# Three methods are available:
#   'base'  — uses base R besselI() with exponential scaling; numerically stable at large k.
#   'Bessel' — uses the Bessel package (unscaled); may overflow at large k.
#   'movMF'  — extracts the constant from movMF::dmovMF(); requires a convention correction
#              since movMF uses the uniform-on-sphere measure rather than Lebesgue measure.
vMF_log_norm_const_exact <- function(k, p, method = 'base'){
  if (p == 3){return(log(2*pi) + log(1 - exp(-2*k)) + k - log(k))} #from Scealy and Wood 2019, this nice and simple for p = 3, with exponential scaling to avoid Inf at large k: log(exp(k) - exp(-k)) = log(1-exp(-2k)) + k
    if (method == 'base'){
      # formula in Scealy and Wood 2019
      # using expon.scaled so that logarithm behaves well at large k
      return(p/2 * log(2*pi) + log(besselI(k, p/2 - 1, expon.scaled = TRUE)) + k - (p/2 - 1) * log(k))  #p/2-1 kind
    }
    if (method == 'Bessel'){
      requireNamespace("Bessel")
      return(p/2 * log(2*pi) + log(Bessel::BesselI(k, p/2 - 1)) - (p/2 - 1) * log(k))  #p/2-1 kind
    }
    if (method == "movMF"){
      requireNamespace("movMF")
      # the normalising constant pops out from evaluating the vMF density at any point orthogonal to the mean
      # movMF uses the uniform-on-sphere measure; convert to Lebesgue measure via the factor
      # gamma(p/2) / (2*pi)^{p/2}. An additional 2^{p/2-1} correction is also needed
      # (reason unclear — possibly a movMF-internal convention for multiple-cluster models).
      const_wrtunif <-
        movMF::dmovMF(matrix(c(0, 1, rep(0, p-2)), nrow = 1),
                      matrix(c(1, rep(0, p-1)), nrow = 1) * k,
                      log = TRUE)
      const_wrtLeb <- const_wrtunif + (p/2 - 1)*log(2) + log(gamma(p/2)) - (p/2) * log(2*pi)
      return(-const_wrtLeb)
    }
}
