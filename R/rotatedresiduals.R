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
