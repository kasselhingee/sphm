#' Links for SvMF Variance Distribution
#' These objects classes are useful for computing log-likelihood given covariates.
#' Different links are: 
#' + Fixed: independent of covariates (a useful first step)
#' + SW19: from Scealy and Wood 2019. Need params for the function g().
#' + Paine20: from Paine et al 2020 for p=3 only, use two free values gamma_1 and gamma_2, to specify the variance. The gamma_1 and gamma_2 dependence on covariates needs specifying. (Looks involved to program initially).
#' + aligned: axes given by removing the mean from the directions given by the columns of P (looks easiest to implement)

