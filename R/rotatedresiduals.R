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
  attr(rresids, "samehemisphere") <- attr(cresids, "samehemisphere")
  return(rresids)
}

cruderesid <- function(y, ypred){
  out <- t(sapply(1:nrow(y), function(i){
    projmat <- ypred[i,] %*% t(ypred[i, ])
    y[i, ] - projmat %*% y[i, ]
  }))
  samehemisphere <- rowSums(y * ypred) >= 0
  attr(out, "samehemisphere") <- samehemisphere
  return(out)
}

JuppRmat <- function(y, base){
  (y+base) %*% t(y+base)/(1+drop(y%*%base)) - diag(1, length(y))
}

# Rotation of tangent vectors from base to y
rotationmat_amaral  <- function(start, end){ #assumes a and b are unit vectors
  ab <- (end %*% start)[[1]]
  alpha <- acos(ab)
  c <- start - end*ab
  c <- c/sqrt(c %*% c)[[1]]
  A <- end%o%c - c%o%end
  Q = diag(length(end)) + sin(alpha)*A + (cos(alpha) - 1)*(end%o%end + c%o%c)
  return(Q)
}
