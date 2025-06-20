set.seed(1234)
savednoisemat <- matrix(runif(100*100), 100, 100)
usethis::use_data(savednoisemat, internal = TRUE, overwrite = TRUE)