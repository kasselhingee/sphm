#' Jupp's Rotated Residuals

#' @export
rotatedresid <- function(y, ypred, base){
  #crude residuals
  cresids <- cruderesid(y, ypred)

  # rotate them
  rresids <- t(sapply(1:nrow(y), function(i){
    JuppRmat(ypred[i, ], base) %*%  cresids[i, ]
  }))  
  return(rresids)
}

cruderesid <- function(y, ypred){
  t(sapply(1:nrow(y), function(i){
    projmat <- ypred[i,] %*% t(ypred[i, ])
    y[i, ] - projmat %*% y[i, ]
  }))  
}

JuppRmat <- function(y, base){
  (y+base) %*% t(y+base)/(1+drop(y%*%base)) - diag(1, length(y))
}


rotationmat_amaral  <- function(a, b){ #assumes a and b are unit vectors
  ab <- (a %*% b)[[1]]
  alpha <- acos(ab)
  c <- b - a*ab
  c <- c/sqrt(c %*% c)[[1]]
  A <- a%o%c - c%o%a
  Q = diag(length(a)) + sin(alpha)*A + (cos(alpha) - 1)*(a%o%a + c%o%c)
  return(Q)
}
