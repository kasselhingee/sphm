test_that("Bessel::besselIasym() matches my C++ implementation of it", {
  expect_equal(besselIasym(1000, 3, 5, log_result = TRUE),
    Bessel::besselIasym(1000, 3, 5, log = TRUE))
  # and expect it to match at low values, but be a poor approximation of the true value
  expect_equal(besselIasym(pi/2 + 0.1, 3, 10, log_result = FALSE), Bessel::besselIasym(pi/2 + 0.1, 3, 10))
  expect_error(expect_equal(besselIasym(pi/2 + 0.1, 3, 10, log_result = FALSE), Bessel::BesselI(pi/2 + 0.1, 3)))
})

test_that("My implementation of the series for besselI works", {
  expect_equal(Bessel::BesselI(pi/2, 3), besselItrunc(pi/2, 3, 6, log_result = FALSE))
  expect_equal(Bessel::BesselI(pi/2, 5), besselItrunc(pi/2, 5, 6, log_result = FALSE))
  expect_equal(Bessel::BesselI(10, 5), besselItrunc(10, 5, 15, log_result = FALSE))
})