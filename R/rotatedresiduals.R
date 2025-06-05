#' Jupp's Rotated Residuals

#' @export
rotatedresid <- function(y, ypred, base, path = "geo"){
  #crude residuals
  cresids <- cruderesid(y, ypred)

  # rotate them
  transportmat <- switch(path, Jupp = JuppRmat, geo = rotationmat_amaral, Amaral = rotationmat_amaral, Absil = partransportmat)
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

# Parallel transport matrix from Absil, Mahony and Sepulchre equation (8.4)
partransportmat <- function(start, end){
  alpha <- drop(acos(start %*% end))
  u <- (end - cos(alpha) * start)/sin(alpha)
  out <- diag(length(start)) - sin(alpha) * start %*% t(u) + (cos(alpha) - 1) * u %*% t(u)
}

resid_SvMF_partransport <- function(y, ypred, k, a, G0){
  rresids_tmp <- rotatedresid(y, ypred, G0[,1])
  rresids_tmp <- rresids_tmp %*% G0
  rresids <- rresids_tmp[, -1]
  rresids_std <- sqrt(k) * rresids %*% diag(a[1]/a[-1])
  attr(rresids_std, "samehemisphere") <-  attr(rresids_tmp, "samehemisphere")
  colnames(rresids_std) <- paste0("r", 1:ncol(rresids_std))
  return(rresids_std)
}