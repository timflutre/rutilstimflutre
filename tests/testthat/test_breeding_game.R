library(rutilstimflutre)
context("Breeding_game")

test_that("makeDfPhenos", {
  ind.ids <- c("ind08", "ind05")
  nb.plots.per.ind <- c(2, 3)
  year <- 2005
  pathogen <- TRUE

  expected <- data.frame(ind=c("ind05","ind05","ind05","ind08","ind08"),
                         year=factor(2005),
                         plot=factor(c(1, 2, 3, 4, 5)),
                         pathogen=TRUE,
                         trait1.raw=NA,
                         trait1=NA,
                         trait2=NA,
                         trait3=NA)

  observed <- makeDfPhenos(ind.ids=ind.ids,
                           nb.plots.per.ind=nb.plots.per.ind,
                           year=year,
                           pathogen=pathogen)

  expect_equal(observed, expected)
})
