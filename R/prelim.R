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
#' @param intercept `TRUE` to include a Euclidean intercept term using a covariate that is always `1`. This is needed for centering of Euclidean covariates, which is part of standardising the covariates. If `intercept = FALSE` then the Euclidean covariates will not be standardised.
#' @export
prelim <- function(y, xs = NULL, xe = NULL, type = "Kassel", fix_qs1 = FALSE, start = NULL, intercept = TRUE, ssqOmbuffer = 2, ...){
  preplist <- list(y = y, xs = xs, xe = xe, start = start)
  # if needed, add Euclidean covariates and update start accordingly
  preplist <- addEuccovars(preplist, type = type, intercept = intercept)
  # standardise y, xe and xe and update start accordingly. Dont standardise xe if intercept = FALSE
  preplist <- standardise_data(preplist, intercept)
  # If start not supplied, choose start close to identities since data standardised
  preplist <- defaultstart(preplist, type)

  # Check Shogo link initiated properly
  if ((type == "Shogo") && (!is.null(preplist$xe))){
    stopifnot(is_Shogo(preplist$start))
    stopifnot(all(preplist$xe[, 1]^2 < sqrt(.Machine$double.eps)))
  }
  
  # RUN
  out <- mobius_vMF(y = preplist$y, xs = preplist$xs, xe = preplist$xe, paramobj0 = preplist$start, fix_qs1 = fix_qs1, fix_qe1 = (type == "Shogo"), ssqOmbuffer = ssqOmbuffer, ...)
  
  # Aspects of the fit that are invariant to coordinates used
  # distances in response space
  pred <- mnlink(xs = preplist$xs, xe = preplist$xe, param = out$solution)
  dists <- acos(rowSums(pred * preplist$y))
  rresids_tmp <- rotatedresid(preplist$y, pred, nthpole(ncol(preplist$y)))
  rresids <- rresids_tmp[, -1]
  attr(rresids, "samehemisphere") <-  attr(rresids_tmp, "samehemisphere")
  colnames(rresids) <- paste0("r", 1:ncol(rresids))
  
  ### revert estimated parameters and pred to pre-standardisation coordinates ###
  est <- undo_recoordinate_Omega(as_mnlink_Omega(out$solution), 
                          yrot = attr(preplist$y, "std_rotation"), 
                          xsrot = attr(preplist$xs, "std_rotation"), #if xs/xe is NULL then attr(xs/xe, ..) is NULL too
                          xerot = attr(preplist$xe, "std_rotation"), 
                          xecenter = attr(preplist$xe, "std_center"),
                          onescovaridx = preplist$onescovaridx)
  
  niceout <- list(
    est = est,
    obj = out$nlopt$objective,
    solution = out$solution, #non-standardised solution
    opt = out$nlopt,
    y = y,
    xs = xs,
    xe = if (!is.null(xe)){if (intercept){destandardise_Euc(preplist$xe, attr(preplist$xe, "std_center"), attr(preplist$xe, "std_rotation"))} else {xe}}, #this recovers any added covariates too
    pred = destandardise_sph(pred, tG = attr(preplist$y, "std_rotation")),
    rresids = rresids,
    dists = dists
  )
  return(niceout)
}
