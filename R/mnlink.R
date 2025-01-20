#' Calculate the Mean Given Covariates
#' Implements mean following Remark 1:
#' \deqn{\mu(x) = PS^{-1}\left(B_s S(Q_s^\top x_s)  +  \frac{B_e(Q_e[,-1]^\top x_e + c_e[-1])}{Q_e[,1]^\top x_e + c_e[1]}\right)}
#' @export
mnlink <- function(xs = NULL, xe = NULL, paramobj = NULL, check = TRUE){
  if (inherits(paramobj, "mnlink_Omega")){
    # by C++
    meanlinkS2Scpp(xs, xe, mnlink_Omega_vec(paramobj), length(paramobj$p1))
  }
  if (inherits(paramobj, "mnlink_cann")){
    # implement here by R
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