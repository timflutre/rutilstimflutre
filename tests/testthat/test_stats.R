library(rutilstimflutre)
context("Stats")

test_that("isSingular", {
  ## example of an ANOVA with 2 factors
  L1 <- 3 # nb of levels of the first factor
  L2 <- 4 # nb of levels of the second factor
  R <- 3 # nb of replicates per treatment
  dat <- data.frame(f1=gl(n=L1, k=L2*R, labels=c("low", "med", "high")),
                    f2=gl(n=L2, k=R, length=L1*L2*R, labels=letters[1:L2]))
  X <- model.matrix(~ 1 + f1 + f2, dat)
  expected <- FALSE
  observed <- isSingular(X)
  expect_equal(observed, expected)

  X2 <- cbind(X, 2*X[,1])
  expected <- TRUE
  observed <- isSingular(X2)
  expect_equal(observed, expected)
})

test_that("corrMatAR1", {
  rho <- 0.7

  n <- 0
  expected <- matrix(data=0, nrow=n, ncol=n)
  observed <- corrMatAR1(n, rho)
  expect_equal(observed, expected)

  n <- 1
  expected <- matrix(data=1, nrow=n, ncol=n)
  observed <- corrMatAR1(n, rho)
  expect_equal(observed, expected)

  n <- 2
  expected <- matrix(data=c(1, rho, rho, 1), nrow=n, ncol=n)
  observed <- corrMatAR1(n, rho)
  expect_equal(observed, expected)

  n <- 3
  expected <- matrix(data=c(1,rho,rho^2, rho,1,rho, rho^2,rho,1),
                     nrow=n, ncol=n)
  observed <- corrMatAR1(n, rho)
  expect_equal(observed, expected)

  n <- 4
  expected <- matrix(data=c(1,rho,rho^2,rho^3,
                            rho,1,rho,rho^2,
                            rho^2,rho,1,rho,
                            rho^3,rho^2,rho,1), nrow=n, ncol=n)
  observed <- corrMatAR1(n, rho)
  expect_equal(observed, expected)
})
