
ull_SvMF_V_R <- function(Vvec, yk, a1, m = m){
  p <- 1+(-1 + sqrt(8 * length(Vvec) + 1))/2
  V <- matrix(NA, nrow = p-1, ncol = p-1)
  V[lower.tri(V, diag = TRUE)] <- Vvec
  V[upper.tri(V)] <- t(V)[upper.tri(V)]
  y <- matrix(yk[1:p], nrow = 1)
  k <- yk[p+1]
  uldSvMF_muV(y, k =k, m = m, a1 = a1, V = V)
}
