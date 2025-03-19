#' Preliminary objective function for S2S Link
#' @importClassesFrom scorematchingad Rcpp_ADFun
#' @details Assumes that the distribution is isotropic around the mean with constant concentration, thus maximising
#' \deqn{\sum_i=1^n y_i^T \mu(x_i).}
#' @param y A set of unit vectors in embedded coordinates, each row corresponds to a single unit vector.
#' @param x A set of covariate vectors (also unit vectors), each row corresponds to the same row in `y`.
#' @param paramobj A set of link parameters. See [`cannS2S()`] and [`OmegaS2S()`].
#' @inheritParams mnlink
#' @export
prelimobj <- function(y, xs = NULL, xe = NULL, param){
  predictedmeans <- mnlink(xs = xs, xe = xe, param = param, check = FALSE)
  stopifnot(nrow(y) == nrow(predictedmeans))
  stopifnot(ncol(y) == ncol(predictedmeans))
  return(-mean(rowSums(y * predictedmeans)))
}

#' Optimisation of the Preliminary Objective Function
#' @details Uses `nloptr`. `NLopt` doesn't have any algorithms for global optimisation with non-linear equality constraints that use provided gradients. So `_parttape` only does local optimisation and uses `NLOPT_LD_SLSQP` which is the only algorithm that takes advantage of derivatives and can handle non-linear equality constraints.
#' 
#' Before fitting, standardises y, xs and xe (*the latter needs implementing*). If supplied, `start`, is updated accordingly.
#' Note that if standardised y has a vMF distribution with the given means, the unstandardised y *does not* because of the second-moment standardisation (I would expect is to not be isotropic).
#' 
#' If `type == "Shogo"` a column of zeros called `'dummyzero'` is added to the front of `xe`.
#' 
#' Default scaling of 0.9 avoids being on the inequality boundary at the start of the search.
#' @param start is a starting parameter object. For Shogo mean link the Qe matrix must have an extra row and column that at the front/top, with 1 in the first entry (and zero elsewhere).
#' @param ... Passed as options to [`nloptr()`]. 
#' @export
prelim <- function(y, xs = NULL, xe = NULL, type = "Kassel", method = "local", start = NULL, ssqOmbuffer = 2, ...){
  # standardise inputs
  y_stdmat <- standardise_mat(y)
  y <- standardise(y, y_stdmat)
  if (!is.null(xs)){
    xs_stdmat <- standardise_mat(xs)
    xs <- standardise(xs, xs_stdmat)
  }
  
  # Center and rotate xe if it exists. 
  # Link isn't equivariant to scaling so dont do it
  if (!is.null(xe)){
    xe_names <- colnames(xe)
    xe_centers <- colMeans(xe)
    xe <- t(t(xe) - xe_centers)
    xe_pcares <- princomp(xe)
    xe <- xe_pcares$scores #ie xe <- xe %*% xe_pcares$loadings
  }
  
  if (!is.null(start)){
    start <- as_mnlink_cann(start)
    start$P <- t(y_stdmat) %*% start$P
    if (!is.null(xs)){
      start$Qs <- t(xs_stdmat) %*% start$Qs
    }
    if (!is.null(xe)){
      start$ce <- drop(t(start$Qe) %*% c(if(type=="Shogo"){0}, xe_centers) + start$ce)
      tmploadings <- xe_pcares$loadings
      if (type=="Shogo"){
        tmploadings <- diag(1, nrow(start$Qe))
        tmploadings[-1, -1] <- xe_pcares$loadings
        }
      start$Qe <- t(tmploadings) %*% start$Qe
    }
  }
  
  
  if (is.null(start)){
    p <- ncol(y)
    if (!is.null(xe)){stopifnot(ncol(xe) + (type == "Shogo") >= p)}
    if (!is.null(xs)){stopifnot(ncol(xs) >= p)}
    start <- mnlink_cann(
                P = diag(p),
                Bs = if (!is.null(xs)){diag(0.9, p-1)},
                Qs = if (!is.null(xs)){diag(1, ncol(xs), p)},
                Be = if (!is.null(xe)){diag(0.9, p-1)},
                Qe = if (!is.null(xe)){diag(1, ncol(xe) + (type == "Shogo"), p)},
                ce = if (!is.null(xe)){c(1, rep(0, p-1))}
    )
    if ((type == "Shogo") && !is.null(xe)){
      start$ce[1] <- 1
    }
    if ((type == "Kassel") && !is.null(xe)){ #default to be larger an all values of -xe
      start$ce[1] <- max(-xe)  +  0.1*IQR(xe)
    }
  }
  
  # RUN
  if (type == "Shogo"){
    stopifnot(is_Shogo(start))
    if (!is.null(xe)){
      xe <- cbind(0, xe)
      xe_names <- c("dummyzero", xe_names)
    }
  }
  if (method == "local"){
    out <- prelim_ad(y = y, xs = xs, xe = xe, paramobj0 = start, type = type, ssqOmbuffer = ssqOmbuffer, ...)
  }
  if (method == "global"){
    out <- prelim_R(y = y, xs = xs, xe = xe, paramobj0 = start, type = type, ...)
  }
  
  # Aspect of the fit using standardised coordinates
  pred <- mnlink(xs = xs, xe = xe, param = out$solution)
  dists <- acos(rowSums(pred * y))
  rresids_tmp <- rotatedresid(y, pred, nthpole(ncol(y)))
  rresids <- rresids_tmp[, -1]
  attr(rresids, "samehemisphere") <-  attr(rresids_tmp, "samehemisphere")
  colnames(rresids) <- paste0("r", 1:ncol(rresids))
  
  # revert estimated parameters and pred to pre-standardisation coordinates
  est <- as_mnlink_cann(out$solution)
  est$P <- y_stdmat %*% est$P
  if (!is.null(xs)){
    est$Qs <- xs_stdmat %*% est$Qs #xs_stdmat has colnames included
  }
  if (!is.null(xe)){
    if (type == "Shogo"){
      est$ce[-1] <- est$ce[-1] - t(est$Qe[-1,-1]) %*% t(xe_pcares$loadings) %*% xe_centers
      est$Qe[-1, -1] <- xe_pcares$loadings %*% est$Qe[-1, -1]
      rownames(est$Qe) <- xe_names #could also get all but the dummyzero name from rownames(xe_pcares$loadings)
    } else {
      warning("destandardisation not implemented")
    }
  }
  est <- as_mnlink_Omega(est)
  pred <- pred %*% t(y_stdmat) #y_stdmat has colnames included
  
  # destandardise xe
  if (type == "Shogo"){
    xe <- cbind(xe[,1], t(t(xe[, -1] %*% t(xe_pcares$loadings)) + xe_centers))
  }
  niceout <- list(
    est = est,
    obj = out$loc_nloptr$objective,
    solution = out$solution, #non-standardised solution
    opt = out$loc_nloptr,
    y = destandardise(y, y_stdmat),
    xs = destandardise(xs, xs_stdmat),
    xe = xe,
    pred = pred,
    rresids = rresids,
    dists = dists
  )
  return(niceout)
}
