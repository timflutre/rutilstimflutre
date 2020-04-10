library(rutilstimflutre)
context("Stats")

test_that("binaryClassif", {
  known.nulls <- c(TRUE, TRUE, FALSE, FALSE, TRUE)
  called.nulls <- c(TRUE, TRUE, FALSE, TRUE, FALSE)

  n <- 5; n0 <- 3; n1 <- 2
  a <- 3; tn <- 2; fn <- 1
  r <- 2; tp <- 1; fp <- 1
  tpp <- tp / n1; fnp <- fn / n1
  tnp <- tn / n0; fpp <- fp / n0
  fdp <- fp / r; ppv <- tp / r
  acc <- (tp + tn) / n
  mcc <- (tp * tn - fp * fn) /
    sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  expected <- c(n=n, n0=n0, n1=n1,
                r=r, tp=tp, fp=fp,
                a=a, tn=tn, fn=fn,
                tpp=tpp, fnp=fnp, tnp=tnp, fpp=fpp, fdp=fdp, ppv=ppv, acc=acc,
                mcc=mcc)

  observed <- binaryClassif(known.nulls=known.nulls, called.nulls=called.nulls)

  expect_equal(observed, expected)
})

test_that("log10WeightedSum_without_weights", {
  x <- c(0.3, 0.01, 0.7)

  expected <- log10(mean(10^x))

  observed <- log10WeightedSum(x=x)

  expect_equal(observed, expected)
})

test_that("log10WeightedSum_with_weights", {
  x <- c(0.3, 0.01, 0.7)
  weights <- c(0.2, 0.5, 0.3)

  expected <- log10(sum(weights * (10^x)))

  observed <- log10WeightedSum(x=x, weights=weights)

  expect_equal(observed, expected)
})

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

test_that("significantTests", {
  P <- 6 # tests
  pvalues <- setNames(c(0.01, 0.2, 0.05, NA, 0.001, 0.8),
                      paste0("snp", 1:P))

  expected <- list(pvalues=pvalues)
  expected[["fwer.bonf"]] <-
    data.frame(pv.adj=p.adjust(pvalues, method="bonferroni"))
  expected[["fwer.bonf"]][["0.01"]] <-
    expected[["fwer.bonf"]][["pv.adj"]] <= 0.01
  expected[["fwer.bonf"]][["0.05"]] <-
    expected[["fwer.bonf"]][["pv.adj"]] <= 0.05
  expected[["fdr.bh"]] <-
    data.frame(pv.adj=p.adjust(pvalues, method="BH"))
  expected[["fdr.bh"]][["0.01"]] <-
    expected[["fdr.bh"]][["pv.adj"]] <= 0.01
  expected[["fdr.bh"]][["0.05"]] <-
    expected[["fdr.bh"]][["pv.adj"]] <= 0.05

  observed <- significantTests(pvalues=pvalues,
                               thresh.fwer.bonf=c(0.01, 0.05),
                               thresh.fdr.bh=c(0.01, 0.05),
                               thresh.fdr.storey=NULL)

  expect_equal(observed, expected)
})
