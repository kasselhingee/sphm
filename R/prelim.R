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
  y <- standardise_sph(y)
  if (!is.null(xs)){xs <- standardise_sph(xs)}
  if ((type == "Shogo") && (!is.null(xe))){
    if (!is.null(xe)){xe <- cbind("dummyzero" = 0, xe)}
    if (!is.null(start)){stopifnot(is_Shogo(start))}
  }
  if (!is.null(xe)){xe <- standardise_Euc(xe)}
  
  # update starting coordinates if provided
  if (!is.null(start)){
    start <- recoordinate_Omega(as_mnlink_Omega(start), 
                                yrot = attr(y, "std_rotation"), 
                                xsrot = attr(xs, "std_rotation"), #if xs/xe is NULL then attr(xs/xe, ..) is NULL too
                                xerot = attr(xe, "std_rotation"), 
                                xecenter = attr(xe, "std_center"))
  }

  # If start not supplied, choose start close to identities
  if (is.null(start)){
    p <- ncol(y)
    if (!is.null(xe)){stopifnot(ncol(xe) >= p)}
    if (!is.null(xs)){stopifnot(ncol(xs) >= p)}
    start <- mnlink_cann(
                P = diag(p),
                Bs = if (!is.null(xs)){diag(0.9, p-1)},
                Qs = if (!is.null(xs)){diag(1, ncol(xs), p)},
                Be = if (!is.null(xe)){diag(0.9, p-1)},
                Qe = if (!is.null(xe)){diag(1, ncol(xe), p)},
                ce = if (!is.null(xe)){c(1, rep(0, p-1))}
    )
    if ((type == "Kassel") && !is.null(xe)){ #default to be larger an all values of -xe
      start$ce[1] <- max(-xe)  +  0.1*IQR(xe)
    }
  }
  
  # RUN
  if (method == "local"){
    out <- prelim_ad(y = y, xs = xs, xe = xe, paramobj0 = start, type = type, ssqOmbuffer = ssqOmbuffer, ...)
  }
  if (method == "global"){
    out <- prelim_R(y = y, xs = xs, xe = xe, paramobj0 = start, type = type, ...)
  }
  
  # Aspects of the fit that are invariant to coordinates used
  # distances in response space
  pred <- mnlink(xs = xs, xe = xe, param = out$solution)
  dists <- acos(rowSums(pred * y))
  rresids_tmp <- rotatedresid(y, pred, nthpole(ncol(y)))
  rresids <- rresids_tmp[, -1]
  attr(rresids, "samehemisphere") <-  attr(rresids_tmp, "samehemisphere")
  colnames(rresids) <- paste0("r", 1:ncol(rresids))
  
  # revert estimated parameters and pred to pre-standardisation coordinates
  est <- undo_recoordinate_Omega(as_mnlink_Omega(out$solution), 
                          yrot = attr(y, "std_rotation"), 
                          xsrot = attr(xs, "std_rotation"), #if xs/xe is NULL then attr(xs/xe, ..) is NULL too
                          xerot = attr(xe, "std_rotation"), 
                          xecenter = attr(xe, "std_center"))
  
  niceout <- list(
    est = est,
    obj = out$loc_nloptr$objective,
    solution = out$solution, #non-standardised solution
    opt = out$loc_nloptr,
    y = destandardise_sph(y, attr(y, "std_rotation")),
    xs = destandardise_sph(xs, attr(xs, "std_rotation")),
    xe = destandardise_Euc(xe, attr(xe, "std_center"), attr(xe, "std_rotation")),
    pred = destandardise_sph(pred, tG = attr(y, "std_rotation")),
    rresids = rresids,
    dists = dists
  )
  return(niceout)
}
