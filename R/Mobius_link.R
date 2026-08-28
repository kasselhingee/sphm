#' Calculate the Mean Given Covariates
#' @description
#' Implements mean link:
#' \deqn{\mu(x) = P\mathcal{S}^{-1}\left(B_s \mathcal{S}(Q_s^\top x_s)  +  \frac{B_e\left(Q_e[,-1]^\top x_e\right)}{Qe[,1]^\top x_e + c_e}\right).}
#' @param xs A matrix of row-vectors of the spherical covariate.
#' @param xe A matrix of row-vectors of the Euclidean covariates.
#' @param param Parameters of the mean link. An object of class "mobius_link_Omega" or
#'   "mobius_link_cann" for the general link, or "mobius_link_prop4i" or
#'   "mobius_link_prop4ii" for the constrained links of Proposition 4.
#'   See [`mobius_link_params`].
#' @details
#' This general form of the mean link encompases the primary form of the mean link in "Regression for spherical responses with linear and spherical covariates using a scaled link function" and a more general form that uses the a stereographic-like projection of the Euclidean covariates.
#' See [`mobius_link_params`] for further details.
#'
#' If `param` is of class "mobius_link_Omega" then means are computed as
#' \deqn{\mu(x) = \frac{(1-\|\tilde{y}(x)\|^2) P[,1] + 2 \tilde{y}(x)}{1+\|\tilde{y}(x)\|^2}}
#' where
#' \deqn{\tilde{y}(x) = \frac{\Omega_s x_s}{1+Q_s[,1]^\top x_s} + \frac{\Omega_e x_e}{c_e+{Q_e[,1]}^\top x_e}}
#' and \eqn{x_s} and \eqn{x_e} are the spherical and Euclidean covariate.
#' @family link-function
#' @export
mobius_link <- function(xs = NULL, xe = NULL, param = NULL, check = TRUE){
  paramclasses <- c("mobius_link_Omega", "mobius_link_cann",
                    "mobius_link_prop4i", "mobius_link_prop4ii")
  if (!is.null(xs)){
    if (inherits(xs, paramclasses)){
      stop("xs is a parameter object (", paste(paramclasses, collapse = ", "), "), but should be a matrix of covariate values.")
    }
    stopifnot(inherits(xs, "matrix"))}
  if (!is.null(xe)){
    if (inherits(xe, paramclasses)){
      stop("xe is a parameter object (", paste(paramclasses, collapse = ", "), "), but should be a matrix of covariate values.")
    }
    stopifnot(inherits(xe, "matrix"))}
  if (inherits(param, "mobius_link_Omega")){
    # by C++
    if (is.null(xs)){xs <- matrix(ncol = 0, nrow = nrow(xe))}
    if (is.null(xe)){xe <- matrix(ncol = 0, nrow = nrow(xs))}
    # Checks that xs and xe compatible with param
    stopifnot(ncol(xs) == length(param$qs1))
    stopifnot(ncol(xe) == length(param$qe1))
    # Evaluate
    out <- mobius_link_cpp(xs, xe, mobius_link_Omega_vec(param), length(param$p1))
  } else if (inherits(param, "mobius_link_cann")){
    out <- mobius_link_pred_cann(xs, xe, param)
  } else if (inherits(param, "mobius_link_prop4ii")){
    # Proposition 4(ii): Bs fixed at I_(p-1). See Mobius_prop4ii.R.
    out <- mobius_link_pred_prop4ii(xs, xe, param, check = check)
  } else if (inherits(param, "mobius_link_prop4i")){
    # Proposition 4(i): Bs fixed at beta_s I_(p-1). See Mobius_prop4i.R.
    out <- mobius_link_pred_prop4i(xs, xe, param, check = check)
  } else {
    stop("param is not of the correct class")
  }
  return(out)
}

# Evaluate the mean link using the canonical parameterisation. Follows the SpEuc form
# (Remark 1 of the manuscript), which is slightly more general than the primary form.
mobius_link_pred_cann <- function(xs = NULL, xe = NULL, paramobj){
  stopifnot(inherits(paramobj, "mobius_link_cann"))
  if (!is.null(xs) && !is.null(xe)){stopifnot(nrow(xs) == nrow(xe))}
  y <- matrix(0, nrow = max(nrow(xs), nrow(xe)), ncol = ncol(paramobj$P) - 1)
  
  if (!is.null(paramobj$Qs)){
    # xs stored as row vectors, so xs %*% Qs right-multiplies; Sp is applied row-wise.
    y <- y + Sp(xs %*% paramobj$Qs) %*% paramobj$Bs
  }
  if (!is.null(paramobj$Qe)){
    xetilde <- xe %*% paramobj$Qe #first column is used in denominator
    # xetilde[,1] = Qe[,1]^T xe acts as the denominator in the stereographic-like Euclidean
    # projection; Qe[,1] is the "north-pole" reference direction in covariate space (the
    # singular direction is -Qe[,1]), and ce is an offset ensuring the denominator stays
    # positive for all training points.
    numerator <- (xetilde[, -1, drop = FALSE]) %*% paramobj$Be
    denominator <- xetilde[,1] + paramobj$ce
    y <- y + (numerator/denominator)
  }
  # iSp maps the accumulated displacement y back to S^{p-1}; t(P) rotates the result
  # to the target orientation (transpose because the output is stored as a row vector).
  out <- iSp(y) %*% t(paramobj$P)
  return(out)
}
