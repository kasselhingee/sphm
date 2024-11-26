testthat("Bessel::besselIasym() matches my C++ implementation of it", {
  expect_equal(besselIasym(1000, 3, 5, log_result = TRUE),
    Bessel::besselIasym(1000, 3, 5, log = TRUE))
  Bessel::BesselI(pi/2 + 0.1, 3)
  Bessel::besselIasym(pi/2 + 0.1, 3, 10)
})