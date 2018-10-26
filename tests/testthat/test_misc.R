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

test_that("inlineFctForm", {
  form <- "log(y) ~ 1 + x"
  expected <- list("log(y)"=c("y", "log"))
  observed <- inlineFctForm(form=form, only.resp=TRUE)
  expect_equal(observed, expected)

  form <- "y ~ 1 + x"
  expected <- list("y"=NA)
  observed <- inlineFctForm(form=form, only.resp=TRUE)
  expect_equal(observed, expected)

  form <- "sqrt(log(y)) ~ 1 + x"
  expected <- list("sqrt(log(y))"=c("y", "sqrt", "log"))
  observed <- inlineFctForm(form=form, only.resp=TRUE)
  expect_equal(observed, expected)

  form <- "y"
  expected <- list("y"=NA)
  observed <- inlineFctForm(form=form, only.resp=TRUE)
  expect_equal(observed, expected)

  form <- "log(y)"
  expected <- list("log(y)"=c("y", "log"))
  observed <- inlineFctForm(form=form, only.resp=TRUE)
  expect_equal(observed, expected)
})
