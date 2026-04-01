#' @keywords internal
"_PACKAGE"

#' @section Overview:
#' `sphm` implements regression for data where each observation is a unit vector
#' (a point on a sphere \eqn{S^{p-1}}). Covariates may be Euclidean vectors,
#' spherical (unit vector) covariates, or both.
#'
#' The regression model uses a **scaled Möbius link function** to map covariates
#' to predicted mean directions, combined with an elliptically symmetric error
#' distribution. See the manuscript "Regression for spherical responses with
#' linear and spherical covariates using a scaled link function" (Kato, Hingee,
#' Scealy, Wood) for full details.
#'
#' @section Typical workflow:
#'
#' **Step 1 — preliminary vMF regression** with [`mobius_vMF()`]
#'
#' Fits a fast preliminary regression using the von Mises–Fisher (vMF) error
#' distribution. The optimiser maximises \eqn{\sum_i y_i^\top \mu(x_i)},
#' which is equivalent to maximising the vMF log-likelihood over the mean link
#' parameters only. The result provides good starting values for Step 2.
#'
#' **Step 2 — full SvMF regression** with [`mobius_SvMF()`]
#'
#' Fits the full model using the scaled von Mises–Fisher (SvMF) error
#' distribution, which allows anisotropic (elliptical) spread around the
#' predicted mean. All parameters are estimated jointly by maximum likelihood
#' using automatic differentiation. The axes of symmetry of the error
#' distribution are constructed via parallel transport (see [`parallel_transport_mat()`]).
#'
#' ```r
#' # Example: simulate then fit on S^2 with a spherical covariate
#' set.seed(1)
#' n <- 50; p <- 3
#' xs <- matrix(rnorm(n * p), n, p)
#' xs <- xs / sqrt(rowSums(xs^2))
#' mean_params <- rand_mobius_link_cann(p = p, qs = p, qe = 0)
#' sim <- rMobius_SvMF(xs = xs, xe = NULL, mean = mean_params,
#'                     k = 5, a = c(1, 2, 0.5), G0 = diag(p))
#' y <- sim[, 1:p]
#'
#' fit_vMF  <- mobius_vMF(y = y, xs = xs)           # Step 1
#' fit_SvMF <- mobius_SvMF(y = y, xs = xs,           # Step 2
#'                          start = fit_vMF$mean)
#' ```
#'
#' @section Key functions:
#'
#' | Task | Function |
#' |---|---|
#' | **Fit regression** | [`mobius_SvMF()`], [`mobius_vMF()`], [`mobius_vMF_refit()`] |
#' | **Simulate data** | [`rMobius_SvMF()`], [`rSvMF()`] |
#' | **Evaluate link** | [`mobius_link()`] |
#' | **Residuals** | [`rotated_resid()`], [`parallel_transport_mat()`] |
#' | **SvMF density** | [`dSvMF()`], [`mobius_SvMF_log_lik()`] |
#' | **SvMF prelim. estimates** | [`SvMF_moment_axes()`], [`SvMF_prelim_scales()`] |
#' | **Link parameterisations** | [`mobius_link_cann()`], [`mobius_link_Omega()`], [`cann2Omega()`] |
#' | **SvMF parameterisations** | [`SvMF_cann()`], [`SvMF_muV()`] |
#' | **Data standardisation** | [`standardise_sph()`], [`second_moment_mat()`] |
#' | **Degrees of freedom** | [`mobius_dof()`] |
#'
#' @section Link function parameterisations:
#'
#' The mean link has two equivalent parameterisations:
#'
#' - **Canonical** ([`mobius_link_cann()`]): parameters `P`, `Bs`, `Qs`, `Be`, `Qe`
#'   directly match the notation in the manuscript (Equation 1). Use this for
#'   interpreting fitted results.
#' - **Omega** ([`mobius_link_Omega()`]): a reparameterisation that merges the rotation
#'   and scale matrices into a single matrix `Omega`. Used internally during
#'   optimisation to avoid working on Stiefel manifolds.
#'
#' Convert between them with [`as_mobius_link_cann()`], [`as_mobius_link_Omega()`], and
#' [`cann2Omega()`].
#'
#' @section C++ and automatic differentiation:
#' The log-likelihood (and its gradient) are implemented
#' in C++ using `Eigen` and `CppAD` (accessed via the `scorematchingad` package) to allow for complicated gradient-based optimisation.
#' The C++ code is in `src/`.
#'
#' @section Reproducing manuscript results:
#' Two Quarto documents in `vignettes/` reproduce all numerical results from the
#' manuscript:
#' - `reproduce_midatlantic.qmd` — midatlantic ridge data on \eqn{S^2}
#' - `reproduce_earthquakes.qmd` — earthquake moment tensors on \eqn{S^4}
#'
#' @name sphm-package
#' @aliases sphm
NULL
