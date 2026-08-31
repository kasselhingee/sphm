# Simulation fixtures shared across the regression test files.
#
# These were defined inside test-mobius_SvMF.R, which made them invisible to every other
# test file; testthat sources helper-*.R before each file, so they belong here.

# function for computing axial differences (ignoring sign)
axis_distance <- function(angle1, angle2 = 0){
  diff <- abs(angle1 - angle2)
  pmin(diff, pi - diff)
}
rcovars <- function(n, qs, qe){
  #generate covariates Gaussianly
  xe <- matrix(rnorm(n*qe), nrow = n)
  #generate covariates on the sphere
  xs <- matrix(rnorm(n*qs), nrow = n)
  xs <- sweep(xs, 1, apply(xs, 1, vnorm), FUN = "/")
  return(list(xe = xe, xs = xs))
}
