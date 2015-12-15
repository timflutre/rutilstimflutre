library(rutilstimflutre)
context("Quantgen")

test_that("estimMaf", {
  N <- 2 # individuals
  P <- 4 # SNPs
  X <- matrix(c(1,1, 1,0, 2,1, 1,NA), nrow=N, ncol=P)
  expected <- c(2/4, 1/4, 1/4, 1/2)
  observed <- estimMaf(X=X)
  expect_equal(observed, expected)
})
