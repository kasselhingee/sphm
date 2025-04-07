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
  
  # for d = 5, nu = 1.5
  expect_equal(besselItrunc(15, 1.5, 15, log_result = FALSE), Bessel::BesselI(15, 1.5), tolerance = 1E-5)
  # for d = 3, nu = 3/2 - 1 = 0.5
  expect_equal(besselItrunc(14, 0.5, 15, log_result = FALSE), Bessel::BesselI(14, 0.5), tolerance = 1E-5)
})

test_that("besselImixed() evaluation works", {
  expect_equal(besselImixed(12, 5, 15, 10, log_result = FALSE), besselItrunc(12, 5, 10, log_result = FALSE))
  expect_equal(besselImixed(18, 5, 15, 10, log_result = TRUE), besselIasym(18, 5, 10, log_result = TRUE))
})

test_that("differentiation of besselImixed() works and allows x below and above threshold", {
  mytape <- tape_besselImixed(9, 0.5, 10, 15) #10 is the threshold
  # evaluation works:
  expect_equal(mytape$forward(0, 9.99), besselImixed(9.99, 0.5, 10, 15))
  expect_equal(mytape$forward(0, 10.1), besselImixed(10.1, 0.5, 10, 15))
  expect_equal(mytape$forward(0, 20), besselImixed(20, 0.5, 10, 15))
  
  # differentiation at 9
  x <- 9
  expect_equal(mytape$Jacobian(x), drop(attr(numericDeriv(quote(besselImixed(x, 0.5, 10, 15)), "x"), "gradient")))
  # differentiation at 11, above the threshold
  x <- 11
  expect_equal(mytape$Jacobian(x), drop(attr(numericDeriv(quote(besselImixed(x, 0.5, 10, 15)), "x"), "gradient")))
})