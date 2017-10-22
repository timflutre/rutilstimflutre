library(rutilstimflutre)
context("VarCompEstimReml")

## example from Mrode (2005), 11.6
## concerning the pre-weaning gain (WWG) of beef calves
.makeDataMrode116 <- function(){
  ## responses
  y <- c(2.6,0.1,1.0,3.0,1.0)

  ## pedigree
  ped <- data.frame(sire=c(0,0,0,1,3,1,4,3),
                 dam=c(0,0,0,0,2,2,5,6))
  Ainv <- quass(ped$sire, ped$dam)

  ## incidence matrix X
  X <- matrix(c(1,0,0,1,1,0,1,1,0,0), nrow=5, ncol=2)

  ## incidence matrix Z
  z1 <- matrix(0, ncol=3, nrow=5)
  z2 <- matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1), ncol=5)
  Z <- cbind(z1,z2)

  return(list(y=y, ped=ped, Ainv=Ainv, X=X, Z=Z))
}

test_that("emreml", {
  inputs <- .makeDataMrode116()

  expected <- c(varE=0.485, varU=0.549)

  observed <- emreml(Ainv=inputs$Ainv, y=inputs$y,
                     X=inputs$X, Z=inputs$Z,
                     initE=0.4, initU=0.2, verbose=0)

  expect_equal(observed, expected, tolerance=0.01, scale=1)
})
