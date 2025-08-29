#' @title Gram-Schmidt Orthogonalisation to get Orientation Axes
#' @description Given a rotation matrix `P`, uses Gram-Schmidt orthogonalisation to obtain a new set of orthogonal axes that include the chosen direction `m` as the first axis.
#' @param m a unit vector
#' @param P an orthogonal matrix
alignedG <- function(m, P){
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
