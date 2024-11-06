#' Jupp's Rotated Residuals

#' @export
rotatedresid <- function(y, ypred, base, path = "geo"){
  #crude residuals
  cresids <- cruderesid(y, ypred)

  # rotate them
  transportmat <- if (path == "Jupp"){JuppRmat} else {rotationmat_amaral}
  rresids <- t(sapply(1:nrow(y), function(i){
    transportmat(ypred[i, ], base) %*%  cresids[i, ]
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

# Rotation of tangent vectors from base to y
rotationmat_amaral  <- function(base, a){ #assumes a and b are unit vectors
  ab <- (a %*% base)[[1]]
  alpha <- acos(ab)
  c <- base - a*ab
  c <- c/sqrt(c %*% c)[[1]]
  A <- a%o%c - c%o%a
  Q = diag(length(a)) + sin(alpha)*A + (cos(alpha) - 1)*(a%o%a + c%o%c)
  return(Q)
}
