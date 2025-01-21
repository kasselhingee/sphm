#' Calculate the Mean Given Covariates
#' @description
#' Implements mean following Remark 1:
#' \deqn{\mu(x) = P\mathcal{S}^{-1}\left(B_s \mathcal{S}(Q_s^\top x_s)  +  \frac{B_e(Q_e[,-1]^\top x_e + c_e[-1])}{Q_e[,1]^\top x_e + c_e[1]}\right)}
#' @param xs A matrix of row-vectors of the spherical covariate.
#' @param xe A matrix of row-vectors of the Euclidean covariates.
#' @param paramobj Parameters of the mean link. As an object of class "mnlink_Omega" or "mnlink_cann". See [`mnlink_params`]. 
#' @details
#' If `paramobj` is of class "mnlink_Omega" then means are computed as
#' \deqn{\mu(x) = \frac{(1-\|\tilde{y}(x)\|^2) P[,1] + 2 \tilde{y}(x)}{1+\|\tilde{y}(x)\|^2}}
#' where
#' \deqn{\tilde{y}(x) = \frac{\Omega_s x_s}{1+Q_s[,1]^\top x_s} + \frac{\Omega_e x_e + \tilde{c}_e}{c_e[1]+{Q_e[,1]}^\top x_e}}
#' and \eqn{x_s} and \eqn{x_e} are the spherical and Euclidean covariate and \eqn{\tilde{c}_e = P[,-1]B_e c_e[-1]} is the `PBce` element of `paramobj`.
#' @export
mnlink <- function(xs = NULL, xe = NULL, paramobj = NULL, check = TRUE){
  if (inherits(paramobj, "mnlink_Omega")){
    # by C++
    meanlinkS2Scpp(xs, xe, mnlink_Omega_vec(paramobj), length(paramobj$p1))
  }
  if (inherits(paramobj, "mnlink_cann")){
    mnlink_pred_cann(xs, xe, paramobj)
  }
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
    xetilde <- t(t(xe %*% paramobj$Qe) + paramobj$ce) #first column is used in denominator
    numerator <- (xetilde[, -1, drop = FALSE]) %*% paramobj$Be
    denominator <- xetilde[,1]
    y <- y + (numerator/denominator)
  }
  out <- iSp(y) %*% t(paramobj$P)
  return(out)
}