#' Homosckedastic SvMF Regression
#' @details
#' The mean is assumed to follow the usual mean link.
#' The concentration and scaling in the SvMF is assumed constant across observations. By constant 'axes' we mean that parallel transport of the axes to some base point according to Jupp's rotated residuals method yields the same set of axes regardless of the location of the mean. The base point will be the first column of the matrix `P` from the mean link so that the base point is not antipodal to any mean location.
#' 
#' KLH: I haven't yet nailed down whether Jupp's parallel transport is parallel transport along the geodesic. I would have thought parallel transport along the geodesic would look more like the matrix `Q` in Amaral et al 2007.


