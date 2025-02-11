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
#' @param paramobj0 is a starting parameter object.
#' @param ... Passed as options to [`nloptr()`]. 
#' @export
prelim <- function(y, xs = NULL, xe = NULL, type = "Kassel", method = "local", start = NULL, ...){
  if (is.null(start)){
    p <- ncol(y)
    if ((type == "Shogo") && !is.null(xe)){xe <- cbind(0, xe)}
    if (!is.null(xe)){stopifnot(ncol(xe) >= p)}
    if (!is.null(xs)){stopifnot(ncol(xs) >= p)}
    start <- mnlink_cann(
                P = diag(p),
                Bs = if (!is.null(xs)){diag(p-1)},
                Qs = if (!is.null(xs)){diag(1, ncol(xs), p)},
                Be = if (!is.null(xe)){diag(p-1)},
                Qe = if (!is.null(xe)){diag(1, ncol(xe), p)},
                ce = if (!is.null(xe)){c(1, rep(0, p-1))}
    )
    if ((type == "Shogo") && !is.null(xe)){
      start$ce[1] <- 1
    }
  }
  
  # check inputs:
  if (type == "Shogo"){stopifnot(is_Shogo(start))}
  if (method == "local"){
    out <- prelim_ad(y = y, xs = xs, xe = xe, paramobj0 = start, type = type, ...)
  }
  if (method == "global"){
    out <- prelim_R(y = y, xs = xs, xe = xe, paramobj0 = start, type = type, ...)
  }
  
  # some aspects of the fit:
  pred <- mnlink(xs = xs, xe = xe, param = out$solution)
  colnames(pred) <- colnames(y)
  rresids <- rotatedresid(y, pred, nthpole(ncol(y)))[, -1]
  colnames(rresids) <- paste0("r", 1:ncol(rresids))
  dists <- rowSums(pred * y)
  
  colnames(out$solution$Omega) <- c(colnames(xs), colnames(xe))
  names(out$solution$qe1) <- colnames(xe)
  names(out$solution$qs1) <- colnames(xs)
  rownames(out$solution$Omega) <- colnames(y)
  names(out$solution$p1) <- colnames(y)
  
  niceout <- list(
    est = out$solution,
    obj = out$loc_nloptr$objective,
    solution = out$solution,
    opt = out$loc_nloptr,
    y = y,
    xs = xs,
    xe = xe,
    pred = pred,
    rresids = rresids,
    dists = dists
  )
  return(niceout)
}
