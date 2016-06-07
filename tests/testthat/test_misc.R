library(rutilstimflutre)
context("Quantgen")

test_that("removeFileExtension", {
  f <- "data.txt.gz"

  expected <- "data.txt"

  observed <- removeFileExtension(file=f, fileext="gz")

  expect_equal(observed, expected)
})
