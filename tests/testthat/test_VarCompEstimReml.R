library(rutilstimflutre)
context("VarCompEstimReml")

## example from Mrode (2005), 11.6
## concerning the pre-weaning gain (WWG) of beef calves
.makeDataMrode116 <- function(){
  ## responses (5x1)
  y <- c(2.6, 0.1, 1.0, 3.0, 1.0)

  ## pedigree
  ped <- data.frame(sire=c(0,0,0,1,3,1,4,3),
                    dam=c(0,0,0,0,2,2,5,6))
  A <- createA(ped$sire, ped$dam)
  Ainv <- quass(ped$sire, ped$dam)

  ## incidence matrix X (5x2)
  X <- matrix(c(1,0,0,1,1,0,1,1,0,0), nrow=5, ncol=2)

  ## incidence matrix Z (5x8)
  z1 <- matrix(0, ncol=3, nrow=5)
  z2 <- matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1), ncol=5)
  Z <- cbind(z1,z2)

  ## in this particular design, the Hessian Matrix is going to be singular
  ## thus, we will add one data point
  y.ns <- c(y, 2.0)
  X.ns <- rbind(X,
                matrix(c(0,1), ncol=2))
  Z.ns <- rbind(Z,
                matrix(c(0,0,0,0,0,1,0,0), ncol=8))

  return(list(ped=ped, A=A, Ainv=Ainv,
              y=y, X=X, Z=Z,
              y.ns=y.ns, X.ns=X.ns, Z.ns=Z.ns))
}

test_that("emreml", {
  inputs <- .makeDataMrode116()

  expected <- c(varE=0.485, varU=0.549)

  observed <- emreml(Ainv=inputs$Ainv, y=inputs$y,
                     X=inputs$X, Z=inputs$Z,
                     initE=0.4, initU=0.2, verbose=0)

  expect_equal(observed, expected, tolerance=0.01, scale=1)
})

test_that("aireml", {
  inputs <- .makeDataMrode116()

  expected <- c(varE=0.489, varU=1.128)

  observed <- aireml(A=inputs$A, y=inputs$y.ns,
                     X=inputs$X.ns, Z=inputs$Z.ns,
                     initE=0.4, initU=0.2, verbose=0)

  expect_equal(observed, expected, tolerance=0.01, scale=1)
})
