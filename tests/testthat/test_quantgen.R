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

test_that("estimGenRel", {
  N <- 2 # individuals
  P <- 4 # SNPs
  X <- matrix(c(1,1, 1,0, 2,1, 1,0), nrow=N, ncol=P)

  mafs <- estimMaf(X=X)

  expected <- matrix(data=0, nrow=N, ncol=N)
  for(p in 1:P)
    expected <- expected +
      (matrix(X[,p] - 2 * mafs[p]) %*% t(X[,p] - 2 * mafs[p])) /
      (4 * mafs[p] * (1 - mafs[p]))
  expected <- expected / P

  observed <- estimGenRel(X=X, mafs=NULL, thresh=0, relationships="additive",
                          method="astle-balding", verbose=0)

  expect_equal(observed, expected)
})
