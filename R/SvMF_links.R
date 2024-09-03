#' Links for SvMF Variance Distribution
#' @rdname SvMF_varlinks
#' These functions map from covariates to the variance distribution of the SvMF.
#' Different links are: 
#' + Fixed: independent of covariates (a useful first step)
#' + SW19: from Scealy and Wood 2019. Need params for the function g().
#' + Paine20: from Paine et al 2020 for p=3 only, use two free values gamma_1 and gamma_2, to specify the variance. The gamma_1 and gamma_2 dependence on covariates needs specifying. (Looks involved to program initially).
#' + aligned: axes given by removing the mean from the directions given by the columns of P (looks easiest to implement)
NULL

#' @describeIn SvMF_varlinks Aligns the columns of the Mobius-link rotation matrix `P` for the mean to the columns of G. Note that the first column of the returned G is the given mean. Returns the matrix G.
alignedG_ <- function(m, P){
  mproj <- m %*% t(m)
  
  #first remove the mean direction from all directions in P. P 'no m'
  Pnom <- (diag(1, length(m)) - mproj) %*% P
  
  # then progressively remove the j and earlier directions from P1
  for (j in 2:(length(m)-1)){
    Pnom[, j] <- Pnom[, j]/vnorm(Pnom[, j])
    Pnom[, (j+1):length(m)] <-  (diag(1, length(m)) - Pnom[, j] %*% t(Pnom[, j])) %*% Pnom[, (j+1):length(m)]
  }
  Pnom[, length(m)] <- Pnom[, length(m)]/vnorm(Pnom[, length(m)])
  G <- Pnom
  G[,1] <- m
  return(G)
}
