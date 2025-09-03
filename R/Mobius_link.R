#' Calculate the Mean Given Covariates
#' @description
#' Implements mean link:
#' \deqn{\mu(x) = P\mathcal{S}^{-1}\left(B_s \mathcal{S}(Q_s^\top x_s)  +  \frac{B_e(Q_e[,-1]^\top x_e\right)}{Qe[,1]^\top x_e + c_e}.}
#' @param xs A matrix of row-vectors of the spherical covariate.
#' @param xe A matrix of row-vectors of the Euclidean covariates.
#' @param param Parameters of the mean link. As an object of class "mnlink_Omega" or "mnlink_cann". See [`mnlink_params`]. 
#' @details
#' This general form of the mean link encompases the primary form of the mean link in "Regression for spherical responses with linear and spherical covariates using a scaled link function" and a more general form that uses the a stereographic-like projection of the Euclidean covariates.
#' See [`mnlink_params`] for further details.
#'
#' If `param` is of class "mnlink_Omega" then means are computed as
#' \deqn{\mu(x) = \frac{(1-\|\tilde{y}(x)\|^2) P[,1] + 2 \tilde{y}(x)}{1+\|\tilde{y}(x)\|^2}}
#' where
#' \deqn{\tilde{y}(x) = \frac{\Omega_s x_s}{1+Q_s[,1]^\top x_s} + \frac{\Omega_e x_e}{c_e+{Q_e[,1]}^\top x_e}}
#' and \eqn{x_s} and \eqn{x_e} are the spherical and Euclidean covariate.
#' @export
mnlink <- function(xs = NULL, xe = NULL, param = NULL, check = TRUE){
  if (!is.null(xs)){
    if (inherits(xs, "mnlink_Omega") | inherits(xs, "mnlink_cann")){
      stop("xs is a parameter object (mnlink_Omega or mnlink_cann), but should be a matrix of covariate values.")
    }
    stopifnot(inherits(xs, "matrix"))}
  if (!is.null(xe)){
    if (inherits(xe, "mnlink_Omega") | inherits(xe, "mnlink_cann")){
      stop("xe is a parameter object (mnlink_Omega or mnlink_cann), but should be a matrix of covariate values.")
    }
    stopifnot(inherits(xe, "matrix"))}
  if (inherits(param, "mnlink_Omega")){
    # by C++
    if (is.null(xs)){xs <- matrix(ncol = 0, nrow = nrow(xe))}
    if (is.null(xe)){xe <- matrix(ncol = 0, nrow = nrow(xs))}
    # Checks that xs and xe compatible with param
    stopifnot(ncol(xs) == length(param$qs1))
    stopifnot(ncol(xe) == length(param$qe1))
    # Evaluate
    out <- mnlink_cpp(xs, xe, mnlink_Omega_vec(param), length(param$p1))
  } else if (inherits(param, "mnlink_cann")){
    out <- mnlink_pred_cann(xs, xe, param)
  } else {
    stop("param is not of the correct class")
  }
  return(out)
}

# following Remark 1, not the main one.
mnlink_pred_cann <- function(xs = NULL, xe = NULL, paramobj){
  stopifnot(inherits(paramobj, "mnlink_cann"))
  if (!is.null(xs) && !is.null(xe)){stopifnot(nrow(xs) == nrow(xe))}
  y <- matrix(0, nrow = max(nrow(xs), nrow(xe)), ncol = ncol(paramobj$P) - 1)
  
  if (!is.null(paramobj$Qs)){
    y <- y + Sp(xs %*% paramobj$Qs) %*% paramobj$Bs
  }
  if (!is.null(paramobj$Qe)){
    xetilde <- xe %*% paramobj$Qe #first column is used in denominator
    numerator <- (xetilde[, -1, drop = FALSE]) %*% paramobj$Be
    denominator <- xetilde[,1] + paramobj$ce
    y <- y + (numerator/denominator)
  }
  out <- iSp(y) %*% t(paramobj$P)
  return(out)
}
