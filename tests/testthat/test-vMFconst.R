test_that("Bessel::besselIasym() matches my C++ implementation of it", {
  expect_equal(besselIasym(1000, 3, 5, log_result = TRUE),
    Bessel::besselIasym(1000, 3, 5, log = TRUE))
  expect_equal(besselIasym(1000, 3, 10, log_result = TRUE), log(Bessel::BesselI(1000, 3, expon.scaled = TRUE)) + 1000)
  # and expect it to match at low values, but be a poor approximation of the true value
  expect_equal(besselIasym(pi/2 + 0.1, 3, 10, log_result = FALSE), Bessel::besselIasym(pi/2 + 0.1, 3, 10))
  expect_error(expect_equal(besselIasym(pi/2 + 0.1, 3, 10, log_result = FALSE), Bessel::BesselI(pi/2 + 0.1, 3)))
  
  # middle sizes of x:
  expect_equal(besselIasym(15, 5, 10, log_result = TRUE), log(Bessel::BesselI(15, 5)))
})

test_that("my implementation of the series for besselI works", {
  expect_equal(Bessel::BesselI(pi/2, 3), besselItrunc(pi/2, 3, 6, log_result = FALSE))
  expect_equal(Bessel::BesselI(pi/2, 5), besselItrunc(pi/2, 5, 6, log_result = FALSE))
  expect_equal(Bessel::BesselI(10, 5), besselItrunc(10, 5, 15, log_result = FALSE))
  expect_equal(Bessel::BesselI(15, 5), besselItrunc(15, 5, 15, log_result = FALSE), tolerance = 1E-5)
  expect_equal(Bessel::BesselI(17, 5), besselItrunc(17, 5, 15, log_result = FALSE), tolerance = 1E-4)
})

test_that("besselImixed() evaluation works", {
  expect_equal(besselImixed(12, 5, 15, 10, log_result = FALSE), besselItrunc(12, 5, 10, log_result = FALSE))
  expect_equal(besselImixed(18, 5, 15, 10, log_result = TRUE), besselIasym(18, 5, 10, log_result = TRUE))
})