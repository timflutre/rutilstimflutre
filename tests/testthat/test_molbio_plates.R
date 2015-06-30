library(rutilstimflutre)
context("Plates")

test_that("init.plates", {
  expected <- list("a"=matrix(data=NA, nrow=3, ncol=5,
                       dimnames=list(LETTERS[1:3], as.character(1:5))),
                   "b"=matrix(data=NA, nrow=8, ncol=12,
                       dimnames=list(LETTERS[1:8], as.character(1:12))))
  observed <- init.plates(n=2, nrow=c(3,8), ncol=c(5,12), names=letters[1:2])
  expect_equal(observed, expected)
})
