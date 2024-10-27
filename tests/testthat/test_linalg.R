library(rutilstimflutre)
context("LinAlg")

test_that("LDLT", {
  ## example matrix from Misztal and Perez-Enciso (1993):
  W <- matrix(NA, nrow=5, ncol=5)
  W[1,] <- c(5, 0, 0, 3, 0)
  W[2,-1] <- c(3, 0, 0, 1)
  W[3,-c(1:2)] <- c(4, 1, 1)
  W[4,-c(1:3)] <- c(4, 0)
  W[5,-c(1:4)] <- 2
  W[lower.tri(W)] <- t(W)[lower.tri(W)]
  W

  ## expected outputs after rounding with 2 digits (from the same paper):
  round_vecD <- c(5.0, 3.0, 4.0, 1.95, 1.38)
  round_U <- matrix(0, nrow=5, ncol=5)
  round_U[1,] <- c(1, 0, 0, 0.6, 0)
  round_U[2,-1] <- c(1, 0, 0, 0.33)
  round_U[3,-c(1:2)] <- c(1, 0.25, 0.25)
  round_U[4,-c(1:3)] <- c(1, -0.13)
  round_U[5,-c(1:4)] <- 1

  observed <- LDLT(x=W, returnFormat="list")
  expect_equal(round(observed$vecD, 2),
               round_vecD)
  expect_equal(round(t(observed$L), 2),
               round_U)

  ## expected outputs without rounding:
  L_c <- t(chol(W))
  diag_L_c <- diag(L_c)
  L <- L_c %*% diag(1 / diag_L_c)
  vecD <- diag_L_c^2

  expected <- L
  diag(expected) <- vecD
  observed <- LDLT(x=W, returnFormat="matrix")
  expect_equal(observed, expected)

  expected <- list(L=L, vecD=vecD)
  observed <- LDLT(x=W, returnFormat="list")
  expect_equal(observed, expected)
})
