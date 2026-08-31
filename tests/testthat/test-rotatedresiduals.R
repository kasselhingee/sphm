test_that("rotatedresiduals() rotates residuals to the northpole along a geodesic correctly", {
  eps <- 0.1
  ypred <- rbind(c(0, 1, 0),
             c(0, 1, 0),
             c(0, 1, 0),
             c(0, 0, 1),
             c(0, 1/sqrt(2), 1/sqrt(2)))
  y <- rbind(c(0, 1, 0) + eps * c(1, 0, 0), #tangent to geodesic
                 c(0, 1, 0) + eps * c(0, 0, 1), #orthogonal to geodesic
                 c(0, 1, 0) + eps * c(1/sqrt(2), 0, 1/sqrt(2)), #pi/4 off from geodesic
                 c(0, 0, 1) + eps * c(1, 0, 0), #tangent to geodesic
                 c(0, 1/sqrt(2), 1/sqrt(2)) + eps * c(1, 0, 0)) #tangent to geodesic
  rresid <- rotated_resid(y = y, ypred = ypred, base = c(1, 0, 0), path = "Absil")
  expect_equal(rresid[1, ], eps * c(0, -1, 0))
  expect_equal(rresid[2, ], eps * c(0, 0, 1))
  c(1/sqrt(2), 0, 1/sqrt(2)) %*% c(0, 0, 1)
  expect_equal(rresid[2, ] %*% rresid[3, ], eps^2 * c(1/sqrt(2), 0, 1/sqrt(2)) %*% c(0, 0, 1))
  expect_equal(rresid[4, ], eps * c(0, 0, -1))
  expect_equal(rresid[5, ], eps * c(0, -1, -1)/sqrt(2))
  
  #and actually Jupp's transport seems to give the reflected residual of geodesic transport
  rresid2 <- rotated_resid(y = y, ypred = ypred, base = c(1, 0, 0), path = "Jupp")
  expect_equal(rresid2, -rresid)
})

# These two tests draw their vectors randomly. Without a seed they depend on the RNG state
# left behind by whichever test files ran earlier, so adding or removing a test elsewhere
# silently changes what they check -- and the antipodal assertion below is not true for
# every vector (see its comment). Seeded here so the file is reproducible on its own.
test_that("parallel_transport_mat() rotates start to end", {
  set.seed(20)
  vecs <- matrix(rnorm(6), nrow = 2)
  vecs <- vecs/sqrt(rowSums(vecs^2))
  expect_equal(drop(parallel_transport_mat(vecs[1, ], vecs[2,]) %*% vecs[1,]), vecs[2, ])
})

# THE SECOND ASSERTION BELOW IS EXPECTED TO FAIL: it documents an open defect in
# parallel_transport_mat() at the antipodal case, end = -start.
#
# (-v) %*% v evaluates to -0.9999999999999998 rather than exactly -1, so
# cvec <- start - end*ab comes out a rounding-level nonzero multiple of v instead of the
# exact zero the algebra gives. The guard at R/rotatedresiduals.R:76 tests
# (cvec %*% cvec)[[1]] > 0 rather than a tolerance, so cvec is normalised to a unit vector,
# the end%o%end + cvec%o%cvec term doubles, and Q comes back as I - 4 v v^T instead of
# I - 2 v v^T. The result is not orthogonal (max|Q^T Q - I| = 2.9), does not map v to -v
# (error 1.1), and has an eigenvalue of -2.96 where a reflection must have -1.
#
# This affects 658 of 2000 random unit vectors -- whichever ones make the dot product land
# just above -1. The test previously had no set.seed() and so drew its vector from whatever
# RNG state earlier test files left behind, which is why it passed before a test was moved
# out of test-mobius_prop4.R. Seeding it makes the failure deterministic rather than
# dependent on file ordering. Widening the guard at R/rotatedresiduals.R:76 to a tolerance
# is a one-line fix.
test_that("parallel_transport_mat() is identity when start=end", {
  set.seed(21)
  myvec <- rnorm(6)
  myvec <- myvec/sqrt(sum(myvec^2))
  expect_equal(parallel_transport_mat(myvec,myvec), diag(6))
  expect_equal(parallel_transport_mat(myvec,-myvec), diag(6) - 2 * myvec %*% t(myvec))
})
