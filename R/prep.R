#  preplist <- list(y = y, xs = xs, xe = xe, start = start)
  ### if needed, add Euclidean covariates and update start accordingly ###
addEuccovars <- function(preplist, type, intercept){
  # if Shogo add a zeros covariate
  if ((type == "Shogo") && (!is.null(preplist$xe))){
    if (!is.null(preplist$xe)){
      if (any(preplist$xe[,1]^2 > .Machine$double.eps)){
        preplist$xe <- cbind("dummyzero" = 0, preplist$xe)
      }
    }
    if (!is.null(preplist$start)){stopifnot(is_Shogo(preplist$start))}
  }
 
  # if intercept==TRUE add a ones 
  onescovaridx <- 0
  if (intercept && is.null(preplist$xe)){warning("Intercept==TRUE will be ignored because no Euclidean covariates")}
  if (intercept && !is.null(preplist$xe)){
    # search for 1s covariate, otherwise add it to end
    constxe <- (apply(preplist$xe, 2, sd) < sqrt(.Machine$double.eps))
    onescovaridx <- which(constxe)[which(abs(colMeans(preplist$xe[,constxe, drop = FALSE]) - 1) < .Machine$double.eps)]
    if (length(onescovaridx) == 0){
      preplist$xe <- cbind(preplist$xe, "ones" = 1)
      onescovaridx <- ncol(preplist$xe)
      if (!is.null(preplist$start)){
        preplist$start <- as_mnlink_cann(preplist$start)
        # if start is wrong dimension add a row of zeros
        if (dim(preplist$start)["qe"] != ncol(preplist$xe)){
          stopifnot(dim(preplist$start)["qe"] == ncol(preplist$xe) - 1)
          preplist$start$Qe <- rbind(preplist$start$Qe, ones = 0)
        }
      }
    }
    if (length(onescovaridx) > 1){onescovaridx <- onescovaridx[1]}
  }
  return(c(
    preplist,
    list(onescovaridx = onescovaridx)
  ))
}

### If start not supplied, choose start close to identities since data standardised ###
defaultstart <- function(preplist, type){
  if (is.null(preplist$start)){
    p <- ncol(preplist$y)
    if (!is.null(preplist$xe)){stopifnot(ncol(preplist$xe) >= p)}
    if (!is.null(preplist$xs)){stopifnot(ncol(preplist$xs) >= p)}
    preplist$start <- mnlink_cann(
                P = diag(p),
                Bs = if (!is.null(preplist$xs)){diag(0.9, p-1)},
                Qs = if (!is.null(preplist$xs)){diag(1, ncol(preplist$xs), p)},
                Be = if (!is.null(preplist$xe)){diag(0.9, p-1)},
                Qe = if (!is.null(preplist$xe)){diag(1, ncol(preplist$xe), p)},
                ce = if (!is.null(preplist$xe)){1}
    )
    if ((type == "Kassel") && !is.null(preplist$xe)){ #default to be larger an all values of -xe
      preplist$start$ce[1] <- max(-preplist$xe)  +  0.1*IQR(preplist$xe)
    }
  }
  return(preplist)
}
