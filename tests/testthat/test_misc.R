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

test_that("matWide2Long", {
  mat <- matrix(c(1,2,3, 2,3,4, 5,6,7), byrow=TRUE,
                nrow=3, ncol=3,
                dimnames=list(c("row1","row2","row3"),
                              c("col1","col2","col3")))

  expected <- data.frame(row=c("row1",
                               "row1","row2",
                               "row1","row2","row3"),
                         col=c("col1",
                               "col2","col2",
                               "col3","col3","col3"),
                         val=c(1,
                               2,3,
                               3,4,7))

  observed <- matWide2Long(mat)

  expect_equal(observed, expected)
})
