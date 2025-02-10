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
    start <- mnlink_cann(
                P = diag(p),
                Bs = if (!is.null(xs)){diag(p-1)},
                Qs = diag(p, ncol(xs)),
                Be = if (!is.null(xe)){diag(p-1)},
                Qe = diag(p, ncol(xe)),
                ce = c(1, rep(0, ncol(xe)-1))
    )
    if ((type == "Shogo") && !is.null(xe)){
      start$ce[1] <- 1
    }
  }
  
  # check inputs:
  if (type == "Shogo"){stopifnot(is_Shogo(start))}
  if (method == "local"){
    out <- prelim_ad(y = y, xs = xs, xe = xe, paramobj0 = start, ...)
  }
  if (method == "global"){
    out <- prelim_R(y = y, xs = xs, xe = xe, paramobj0 = start, ...)
  }
  return(out)
}
