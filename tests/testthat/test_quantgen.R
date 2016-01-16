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
                               dimnames=list(c("ind27-x-ind10_h1",
                                               "ind27-x-ind10_h2"),
                                             paste0("snp", 1:(P/2)))),
                   chr2=matrix(c(1,0,1,0, 1,0,1,0), nrow=2, byrow=TRUE,
                               dimnames=list(c("ind27-x-ind10_h1",
                                               "ind27-x-ind10_h2"),
                                             paste0("snp", (P/2+1):P))))

  observed <- fecundation(gam1, gam2, "ind27-x-ind10")

  expect_equal(observed, expected)
})

test_that("drawLocCrossovers", {
  set.seed(1)

  crosses <- data.frame(parent1=c("ind2", "ind1"),
                        parent2=c("ind1", NA),
                        child=c("ind3", "ind1-hd"),
                        stringsAsFactors=FALSE)
  nb.snps <- setNames(c(1437,1502), c("chr1", "chr2"))

  ## nb.crossovers <- c(1, 1, 2, 4, 1, 4)
  expected <- list(ind3=list(ind2=list(chr1=c(1357),
                                       chr2=c(992)),
                             ind1=list(chr1=c(89, 904),
                                       chr2=c(265, 310, 576, 1030))),
                   "ind1-hd"=list(ind1=list(chr1=c(1106),
                                            chr2=c(570, 748, 1077, 1487))))

  observed <- drawLocCrossovers(crosses, nb.snps)

  expect_equal(observed, expected)
})

test_that("makeCrosses_dh", {
  N <- 2 # individuals
  P <- 8 # SNPs
  haplos <- list(chr1=matrix(data=c(1,1,1,0, 1,1,0,1, 1,1,1,0, 1,1,1,0),
                             nrow=2*N, ncol=P/2,
                             dimnames=list(c("ind1_h1", "ind1_h2",
                                             "ind2_h1", "ind2_h2"),
                                           paste0("snp", 1:(P/2)))),
                 chr2=matrix(data=c(0,0,1,0, 0,0,0,1, 0,0,1,1, 0,0,1,0),
                             nrow=2*N, ncol=P/2,
                             dimnames=list(c("ind1_h1", "ind1_h2",
                                             "ind2_h1", "ind2_h2"),
                                           paste0("snp", (P/2+1):P))))
  crosses <- data.frame(parent1=c("ind2"),
                        parent2=NA,
                        child=c("ind2-hd"),
                        stringsAsFactors=FALSE)
  loc.crossovers <- list("ind2-hd"=list(
                             ind2=list(
                                 chr1=c(2),
                                 chr2=c(3))))

  expected <- list(chr1=matrix(rep(c(1,0,0,0), 2), nrow=2, ncol=P/2, byrow=TRUE,
                               dimnames=list(c("ind2-hd_h1", "ind2-hd_h2"),
                                             paste0("snp", 1:(P/2)))),
                   chr2=matrix(rep(c(1,0,1,0), 2), nrow=2, ncol=P/2, byrow=TRUE,
                               dimnames=list(c("ind2-hd_h1", "ind2-hd_h2"),
                                             paste0("snp", (P/2+1):P))))

  observed <- makeCrosses(haplos, crosses, loc.crossovers, verbose=0)

  expect_equal(observed, expected)
})

test_that("makeCrosses", {
  N <- 2 # individuals
  P <- 8 # SNPs
  haplos <- list(chr1=matrix(data=c(1,1,1,0, 1,1,0,1, 1,1,1,0, 1,0,1,0),
                             nrow=2*N, ncol=P/2,
                             dimnames=list(c("ind1_h1", "ind1_h2",
                                             "ind2_h1", "ind2_h2"),
                                           paste0("snp", 1:(P/2)))),
                 chr2=matrix(data=c(0,0,1,0, 1,0,0,1, 0,0,1,1, 0,0,1,0),
                             nrow=2*N, ncol=P/2,
                             dimnames=list(c("ind1_h1", "ind1_h2",
                                             "ind2_h1", "ind2_h2"),
                                           paste0("snp", (P/2+1):P))))
  crosses <- data.frame(parent1=c("ind2"),
                        parent2=c("ind1"),
                        child=c("ind3"),
                        stringsAsFactors=FALSE)
  loc.crossovers <- list("ind3"=list(
                             ind2=list(
                                 chr1=c(2),
                                 chr2=c(3)),
                             ind1=list(
                                 chr1=c(1),
                                 chr2=c(2))))

  expected <- list(chr1=matrix(c(1,0,0,0, 1,1,1,0), nrow=2, ncol=P/2, byrow=TRUE,
                               dimnames=list(c("ind3_h1", "ind3_h2"),
                                             paste0("snp", 1:(P/2)))),
                   chr2=matrix(c(1,0,1,0, 0,1,0,0), nrow=2, ncol=P/2, byrow=TRUE,
                               dimnames=list(c("ind3_h1", "ind3_h2"),
                                             paste0("snp", (P/2+1):P))))

  observed <- makeCrosses(haplos, crosses, loc.crossovers, verbose=0)

  expect_equal(observed, expected)
})

test_that("distSnpPairs", {
  pair.snps <- data.frame(loc1=c("snp1", "snp1", "snp2"),
                          loc2=c("snp2", "snp3", "snp3"))
  snp.coords <- data.frame(chr=rep("chr1", 3),
                           pos=c(1, 9, 13))
  rownames(snp.coords) <- paste0("snp", 1:3)

  expected <- c(7, 11, 3)

  observed <- distSnpPairs(pair.snps, snp.coords)

  expect_equal(observed, expected)
})
