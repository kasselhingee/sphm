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
