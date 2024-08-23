#' The Scealy and Wood 2019 Link for V and kappa
#' Map from covariates to V and kappa with parameters using the method in Scealy and Wood 2019.

#' @param x is a single covariate, equivalent to the v_i in Scealy and Wood 2019
Vonk <- function(x, sigma, delta, c1){
  p <- length(sigma) + 1
  stopifnot(p == 3)
  diagmat <- diag(x^delta[1], sigma[2]*x^delta[2])
  sigma[1]^2 * diagmat %*% 
    matrix(c(1, c1, c1, 1), p-1,p-1, byrow = FALSE) %*%
    diagmat
}



