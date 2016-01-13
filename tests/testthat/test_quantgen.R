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

test_that("estimGenRel_astle-balding", {
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

test_that("getIndNamesFromHaplos", {
  N <- 2 # individuals
  P <- 3 # SNPs
  haplos <- list(chr1=matrix(data=c(1,1,1,0, 1,0,0,1), nrow=2*N, ncol=P-1,
                             dimnames=list(c("ind1_h1", "ind1_h2", "ind2_h1", "ind2_h2"),
                                           c("snp1", "snp2"))),
                 chr2=matrix(data=c(0,0,1,0), nrow=2*N, ncol=1,
                             dimnames=list(c("ind1_h1", "ind1_h2", "ind2_h1", "ind2_h2"),
                                           c("snp3"))))

  expected <- c("ind1", "ind2")

  observed <- getIndNamesFromHaplos(haplos)

  expect_equal(observed, expected)
})

test_that("getHaplosInd", {
  N <- 2 # individuals
  P <- 3 # SNPs

  haplos <- list(chr1=matrix(data=c(1,1,1,0, 1,0,0,1), nrow=2*N, ncol=P-1,
                             dimnames=list(c("ind1_h1", "ind1_h2", "ind2_h1", "ind2_h2"),
                                           c("snp1", "snp2"))),
                 chr2=matrix(data=c(0,0,1,0), nrow=2*N, ncol=1,
                             dimnames=list(c("ind1_h1", "ind1_h2", "ind2_h1", "ind2_h2"),
                                           c("snp3"))))

  expected <- list(chr1=matrix(data=c(1,0, 0,1), nrow=2, ncol=P-1,
                               dimnames=list(c("ind2_h1", "ind2_h2"),
                                             c("snp1", "snp2"))),
                   chr2=matrix(data=c(1,0), nrow=2, ncol=1,
                               dimnames=list(c("ind2_h1", "ind2_h2"),
                                             c("snp3"))))

  observed <- getHaplosInd(haplos, "ind2")

  expect_equal(observed, expected)
})

test_that("getHaplosInds", {
  N <- 3 # individuals
  P <- 3 # SNPs

  haplos <- list(chr1=matrix(data=c(1,1,1,0,0,1, 1,0,0,1,1,0), nrow=2*N, ncol=P-1,
                             dimnames=list(c("ind1_h1", "ind1_h2", "ind2_h1",
                                             "ind2_h2", "ind3_h1", "ind3_h2"),
                                           c("snp1", "snp2"))),
                 chr2=matrix(data=c(0,0,1,0,1,1), nrow=2*N, ncol=1,
                             dimnames=list(c("ind1_h1", "ind1_h2", "ind2_h1",
                                             "ind2_h2", "ind3_h1", "ind3_h2"),
                                           c("snp3"))))

  expected <- list(chr1=matrix(data=c(1,1,1,0, 1,0,0,1), nrow=2*(N-1), ncol=P-1,
                               dimnames=list(c("ind1_h1", "ind1_h2", "ind2_h1",
                                               "ind2_h2"),
                                             c("snp1", "snp2"))),
                   chr2=matrix(data=c(0,0,1,0), nrow=2*(N-1), ncol=1,
                               dimnames=list(c("ind1_h1", "ind1_h2", "ind2_h1",
                                               "ind2_h2"),
                                             c("snp3"))))

  observed <- getHaplosInds(haplos, c("ind1", "ind2"))

  expect_equal(observed, expected)
})

test_that("makeGameteSingleIndSingleChrom", {
  P <- 4 # SNPs
  haplos.par.chr <- matrix(data=c(1,0, 0,1, 1,0, 1,0), nrow=2, ncol=P,
                           dimnames=list(c("ind2_h1", "ind2_h2"),
                                         paste0("snp", 1:P)))
  loc.crossovers <- c(2)

  expected <- matrix(c(1,0, 0,0), nrow=1,
                     dimnames=list("ind2", paste0("snp", 1:P)))

  observed <- makeGameteSingleIndSingleChrom(haplos.par.chr, loc.crossovers)

  expect_equal(observed, expected)
})

test_that("makeGameteSingleInd", {
  P <- 8 # SNPs
  haplos.par <- list(chr1=matrix(data=c(1,0, 0,1, 1,0, 1,0), nrow=2, ncol=P/2,
                                 dimnames=list(c("ind2_h1", "ind2_h2"),
                                               paste0("snp", 1:(P/2)))),
                     chr2=matrix(data=c(1,0, 0,1, 1,0, 1,0), nrow=2, ncol=P/2,
                                 dimnames=list(c("ind2_h1", "ind2_h2"),
                                               paste0("snp", (P/2+1):P))))
  loc.crossovers <- list(chr1=c(2), chr2=c(3))

  expected <- list(chr1=matrix(c(1,0, 0,0), nrow=1,
                               dimnames=list("ind2", paste0("snp", 1:(P/2)))),
                   chr2=matrix(c(1,0,1,0), nrow=1,
                               dimnames=list("ind2", paste0("snp", (P/2+1):P))))

  observed <- makeGameteSingleInd(haplos.par, loc.crossovers)

  expect_equal(observed, expected)
})

test_that("fecundation", {
  P <- 8 # SNPs
  gam1 <- list(chr1=matrix(c(1,0,0,0), nrow=1,
                           dimnames=list("ind27", paste0("snp", 1:(P/2)))),
               chr2=matrix(c(1,0,1,0), nrow=1,
                           dimnames=list("ind27", paste0("snp", (P/2+1):P))))
  gam2 <- list(chr1=matrix(c(1,0,1,1), nrow=1,
                           dimnames=list("ind10", paste0("snp", 1:(P/2)))),
               chr2=matrix(c(1,0,1,0), nrow=1,
                           dimnames=list("ind10", paste0("snp", (P/2+1):P))))

  expected <- list(chr1=matrix(c(1,0,0,0, 1,0,1,1), nrow=2, byrow=TRUE,
                               dimnames=list(c("ind27_x_ind10_h1",
                                               "ind27_x_ind10_h2"),
                                             paste0("snp", 1:(P/2)))),
                   chr2=matrix(c(1,0,1,0, 1,0,1,0), nrow=2, byrow=TRUE,
                               dimnames=list(c("ind27_x_ind10_h1",
                                               "ind27_x_ind10_h2"),
                                             paste0("snp", (P/2+1):P))))

  observed <- fecundation(gam1, gam2)

  expect_equal(observed, expected)
})

test_that("makeDoubledHaploids", {
  N <- 2 # individuals
  P <- 8 # SNPs
  haplos <- list(chr1=matrix(data=c(1,1,1,0, 1,1,0,1, 1,1,1,0, 1,1,1,0), nrow=2*N, ncol=P/2,
                             dimnames=list(c("ind1_h1", "ind1_h2", "ind2_h1", "ind2_h2"),
                                           paste0("snp", 1:(P/2)))),
                 chr2=matrix(data=c(0,0,1,0, 0,0,0,1, 0,0,1,1, 0,0,1,0), nrow=2*N, ncol=P/2,
                             dimnames=list(c("ind1_h1", "ind1_h2", "ind2_h1", "ind2_h2"),
                                           paste0("snp", (P/2+1):P))))
  ind.names <- c("ind2")
  loc.crossovers <- list(ind2=list(chr1=c(2), chr2=c(3)))
  child.names <- c("ind3")

  expected <- list(chr1=matrix(rep(c(1,0,0,0), 2), nrow=2, ncol=P/2, byrow=TRUE,
                               dimnames=list(c("ind3_h1", "ind3_h2"),
                                             paste0("snp", 1:(P/2)))),
                   chr2=matrix(rep(c(1,0,1,0), 2), nrow=2, ncol=P/2, byrow=TRUE,
                               dimnames=list(c("ind3_h1", "ind3_h2"),
                                             paste0("snp", (P/2+1):P))))

  observed <- makeDoubledHaploids(haplos, ind.names, loc.crossovers,
                                  child.names)

  expect_equal(observed, expected)
})
