library(rutilstimflutre)
context("Misc")

test_that("removeFileExtension", {
  f <- "data.txt.gz"

  expected <- "data.txt"
  observed <- removeFileExtension(file=f, fileext=".gz")
  expect_equal(observed, expected)

  expected <- "data"
  observed <- removeFileExtension(file=f, fileext=".txt.gz")
  expect_equal(observed, expected)
})
