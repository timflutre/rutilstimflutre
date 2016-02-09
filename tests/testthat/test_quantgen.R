library(rutilstimflutre)
context("Quantgen")

test_that("estimAf", {
  N <- 2 # individuals
  P <- 4 # SNPs
  X <- matrix(c(1,1, 1,0, 2,1, 1,NA), nrow=N, ncol=P,
              dimnames=list(paste0("ind", 1:N), paste0("snp", 1:P)))

  expected <- setNames(c(2/4, 1/4, 3/4, 1/2), colnames(X))

  observed <- estimAf(X=X)

  expect_equal(observed, expected)
})

test_that("estimMaf", {
  N <- 2 # individuals
  P <- 4 # SNPs
  X <- matrix(c(1,1, 1,0, 2,1, 1,NA), nrow=N, ncol=P,
              dimnames=list(paste0("ind", 1:N), paste0("snp", 1:P)))

  expected <- setNames(c(2/4, 1/4, 1/4, 1/2), colnames(X))

  observed <- estimMaf(X=X)

  expect_equal(observed, expected)

  afs <- setNames(c(2/4, 1/4, 3/4, 1/2), colnames(X))
  observed <- estimMaf(afs=afs)

  expect_equal(observed, expected)
})

test_that("estimGenRel_vanraden1", {
  N <- 3 # individuals
  P <- 4 # SNPs
  X <- matrix(c(0,0,2, 1,1,0, 0,1,0, 1,0,0), nrow=N, ncol=P,
              dimnames=list(paste0("ind", 1:N), paste0("snp", 1:P)))

  afs <- setNames(c(0.383, 0.244, 0.167, 0.067), colnames(X))
  tmp <- matrix(rep(2*(afs-0.5), N), nrow=N, ncol=P, byrow=TRUE)
  Z <- X - 1 - tmp
  denom <- 0
  for(p in 1:P)
    denom <- denom + afs[p] * (1 - afs[p])
  denom <- 2 * denom
  expected <- (Z %*% t(Z)) / denom

  observed <- estimGenRel(X=X, afs=afs, thresh=0, relationships="additive",
                          method="vanraden1", verbose=0)

  expect_equal(observed, expected)
})

test_that("estimGenRel_astle-balding", {
  N <- 2 # individuals
  P <- 4 # SNPs
  X <- matrix(c(1,1, 1,0, 2,1, 1,0), nrow=N, ncol=P,
              dimnames=list(paste0("ind", 1:N), paste0("snp", 1:P)))

  afs <- estimAf(X=X)

  expected <- matrix(data=0, nrow=N, ncol=N,
                     dimnames=list(rownames(X), rownames(X)))
  for(p in 1:P)
    expected <- expected +
      (matrix(X[,p] - 2 * afs[p]) %*% t(X[,p] - 2 * afs[p])) /
      (4 * afs[p] * (1 - afs[p]))
  expected <- 2 * (1/P) * expected

  observed <- estimGenRel(X=X, afs=afs, thresh=0, relationships="additive",
                          method="astle-balding", verbose=0)

  expect_equal(observed, expected)
})

test_that("estimGenRel_yang", {
  N <- 2 # individuals
  P <- 4 # SNPs
  X <- matrix(c(1,1, 1,0, 2,1, 1,1), nrow=N, ncol=P,
              dimnames=list(paste0("ind", 1:N), paste0("snp", 1:P)))

  afs <- estimAf(X=X)

  expected <- matrix(data=0, nrow=N, ncol=N,
                     dimnames=list(rownames(X), rownames(X)))
  for(i in 1:N){
    summands <- rep(NA, P)
    for(p in 1:P)
      summands[p] <- (X[i,p]^2 - (1 + 2 * afs[p]) * X[i,p] + 2 * afs[p]^2) /
        (2 * afs[p] * (1 - afs[p]))
    expected[i,i] <- 1 + (1/P) * sum(summands)
  }
  summands <- rep(NA, P)
  for(p in 1:P)
    summands[p] <- ((X[1,p] - 2 * afs[p]) * (X[2,p] - 2 * afs[p])) /
      (2 * afs[p] * (1 - afs[p]))
  expected[1,2] <- (1/P) * sum(summands)
  expected[2,1] <- expected[1,2]

  observed <- estimGenRel(X=X, afs=afs, thresh=0, relationships="additive",
                          method="yang", verbose=0)

  expect_equal(observed, expected)
})

test_that("estimGenRel_su", {
  N <- 3 # individuals
  P <- 4 # SNPs
  X <- matrix(c(0,0,2, 1,1,0, 0,1,0, 1,0,0), nrow=N, ncol=P,
              dimnames=list(paste0("ind", 1:N), paste0("snp", 1:P)))

  afs <- setNames(c(0.383, 0.244, 0.167, 0.067), colnames(X))
  H <- matrix(c(0 - 2 * afs[1] * (1 - afs[1]),
                0 - 2 * afs[1] * (1 - afs[1]),
                0 - 2 * afs[1] * (1 - afs[1]),
                1 - 2 * afs[2] * (1 - afs[2]),
                1 - 2 * afs[2] * (1 - afs[2]),
                0 - 2 * afs[2] * (1 - afs[2]),
                0 - 2 * afs[3] * (1 - afs[3]),
                1 - 2 * afs[3] * (1 - afs[3]),
                0 - 2 * afs[3] * (1 - afs[3]),
                1 - 2 * afs[4] * (1 - afs[4]),
                0 - 2 * afs[4] * (1 - afs[4]),
                0 - 2 * afs[4] * (1 - afs[4])),
              nrow=N, ncol=P,
              dimnames=list(rownames(X), colnames(X)))
  denom <- 0
  for(p in 1:P)
    denom <- denom + 2 * afs[p] * (1 - afs[p]) *
      (1 - 2 * afs[p] * (1 - afs[p]))
  expected <- (H %*% t(H)) / denom

  observed <- estimGenRel(X=X, afs=afs, thresh=0, relationships="dominance",
                          method="su", verbose=0)

  expect_equal(observed, expected)
})

test_that("estimGenRel_vitezica", {
  N <- 3 # individuals
  P <- 4 # SNPs
  X <- matrix(c(0,0,2, 1,1,0, 0,1,0, 1,0,0), nrow=N, ncol=P,
              dimnames=list(paste0("ind", 1:N), paste0("snp", 1:P)))

  afs <- setNames(c(0.383, 0.244, 0.167, 0.067), colnames(X))
  W <- matrix(c(- 2 * (1 - afs[1])^2,
                - 2 * (1 - afs[1])^2,
                - 2 * afs[1]^2,
                2 * afs[2] * (1 - afs[2]),
                2 * afs[2] * (1 - afs[2]),
                - 2 * (1 - afs[2])^2,
                - 2 * (1 - afs[3])^2,
                2 * afs[3] * (1 - afs[3]),
                - 2 * (1 - afs[3])^2,
                2 * afs[4] * (1 - afs[4]),
                - 2 * (1 - afs[4])^2,
                - 2 * (1 - afs[4])^2),
              nrow=N, ncol=P,
              dimnames=list(rownames(X), colnames(X)))
  denom <- 0
  for(p in 1:P)
    denom <- denom + (2 * afs[p] * (1 - afs[p]))^2
  expected <- (W %*% t(W)) / denom

  observed <- estimGenRel(X=X, afs=afs, thresh=0, relationships="dominance",
                          method="vitezica", verbose=0)

  expect_equal(observed, expected)
})

test_that("getIndNamesFromHaplos", {
  N <- 2 # individuals
  P <- 3 # SNPs
  haplos <- list(chr1=matrix(data=c(1,1,1,0, 1,0,0,1), nrow=2*N, ncol=P-1,
                             dimnames=list(c("ind1_h1", "ind1_h2", "ind2_h1",
                                             "ind2_h2"),
                                           c("snp1", "snp2"))),
                 chr2=matrix(data=c(0,0,1,0), nrow=2*N, ncol=1,
                             dimnames=list(c("ind1_h1", "ind1_h2", "ind2_h1",
                                             "ind2_h2"),
                                           c("snp3"))))

  expected <- c("ind1", "ind2")

  observed <- getIndNamesFromHaplos(haplos)

  expect_equal(observed, expected)
})

test_that("getHaplosInd", {
  N <- 2 # individuals
  P <- 3 # SNPs

  haplos <- list(chr1=matrix(data=c(1,1,1,0, 1,0,0,1), nrow=2*N, ncol=P-1,
                             dimnames=list(c("ind1_h1", "ind1_h2", "ind2_h1",
                                             "ind2_h2"),
                                           c("snp1", "snp2"))),
                 chr2=matrix(data=c(0,0,1,0), nrow=2*N, ncol=1,
                             dimnames=list(c("ind1_h1", "ind1_h2", "ind2_h1",
                                             "ind2_h2"),
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

test_that("splitGenomesTrainTest", {
  nb.inds.train <- 1
  nb.inds.test <- 1
  nb.inds <- nb.inds.train + nb.inds.test
  nb.snps <- 2
  snp.names <- c("snp1", "snp2")

  chr1 <- matrix(c(1,0,1,0, 1,1,0,1), nrow=2 * nb.inds, ncol=nb.snps,
                 dimnames=list(c("ind1_h1", "ind1_h2", "ind2_h1", "ind2_h2"),
                               snp.names))
  genomes <- list(haplos=list(chr1=chr1),
                  genos=matrix(c(1,1,2,1), nrow=nb.inds, ncol=nb.snps,
                               dimnames=list(c("ind1", "ind2"),
                                             snp.names)))

  expected <- list(genomes=list(haplos=list(
                                    chr1=matrix(c(1,0, 1,1),
                                                nrow=2 * nb.inds.train,
                                                ncol=nb.snps,
                                                dimnames=list(c("ind1_h1",
                                                                "ind1_h2"),
                                                              snp.names))),
                                genos=matrix(c(1,2),
                                             nrow=nb.inds.train,
                                             ncol=nb.snps,
                                             dimnames=list(c("ind1"),
                                                           snp.names))),
                   genomes.pred=list(haplos=list(
                                         chr1=matrix(c(1,0, 0,1),
                                                     nrow=2 * nb.inds.test,
                                                     ncol=nb.snps,
                                                     dimnames=list(c("ind2_h1",
                                                                     "ind2_h2"),
                                                                   snp.names))),
                                     genos=matrix(c(1,1),
                                                  nrow=nb.inds.test,
                                                  ncol=nb.snps,
                                                  dimnames=list(c("ind2"),
                                                                snp.names))))

  observed <- splitGenomesTrainTest(genomes, nb.inds.pred=1)

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
  snp.pairs <- data.frame(loc1=c("snp1", "snp1", "snp2"),
                          loc2=c("snp2", "snp3", "snp3"),
                          stringsAsFactors=FALSE)
  snp.coords <- data.frame(chr=rep("chr1", 3),
                           pos=c(1, 9, 13))
  rownames(snp.coords) <- paste0("snp", 1:3)

  expected <- c(7, 11, 3)

  observed <- distSnpPairs(snp.pairs, snp.coords)

  expect_equal(observed, expected)
})

test_that("calcAvgPwDiffBtwHaplos", {
  ## figure 1.4 from Wakeley (2008)
  haplos.chr <- matrix(data=c(1,1,1,0,0, 1,1,0,1,1, 1,0,0,0,0, 1,1,0,0,0),
                       nrow=5, ncol=4)

  expected <- 2

  observed <- calcAvgPwDiffBtwHaplos(haplos.chr)

  expect_equal(observed, expected)
})

test_that("estimLd_cor-r2", {
  N <- 4 # individuals
  P <- 5 # SNPs
  X <- matrix(c(2,1,0,1, 2,0,0,1, 1,0,0,1, 0,0,0,1, 1,0,1,2), nrow=N, ncol=P,
              dimnames=list(inds=paste0("ind", 1:N), snps=paste0("snp", 1:P)))
  snp.coords <- data.frame(chr=c(rep("chr1", P-1), "chr2"),
                           pos=c(2, 17, 25, 33, 5),
                           stringsAsFactors=FALSE)
  rownames(snp.coords) <- colnames(X)

  expected <- data.frame(loc1=c(rep("snp1",3), rep("snp2",2), "snp3"),
                         loc2=c("snp2","snp3","snp4","snp3","snp4","snp4"),
                         cor2=c(cor(X[,"snp1"], X[,"snp2"])^2,
                                cor(X[,"snp1"], X[,"snp3"])^2,
                                cor(X[,"snp1"], X[,"snp4"])^2,
                                cor(X[,"snp2"], X[,"snp3"])^2,
                                cor(X[,"snp2"], X[,"snp4"])^2,
                                cor(X[,"snp3"], X[,"snp4"])^2),
                         stringsAsFactors=TRUE)

  observed <- estimLd(X=X, snp.coords=snp.coords, only.chr="chr1")

  expect_equal(observed, expected)

  ## check that LDcorSV returns the same results
  colnames(expected)[3] <- "r2"
  observed <- estimLd(X=X, snp.coords=snp.coords, only.chr="chr1",
                      use.ldcorsv=TRUE)
  expect_equal(observed, expected)
})

test_that("distConsecutiveSnps", {
  P <- 5 # SNPs
  snp.coords <- data.frame(chr=c(rep("chr1", P-2), rep("chr2", 2)),
                           pos=c(2, 17, 28, 5, 100),
                           stringsAsFactors=FALSE)
  rownames(snp.coords) <- paste0("snp", 1:P)

  expected <- list(chr1=setNames(c(14, 10), c("snp2-snp1", "snp3-snp2")),
                   chr2=setNames(c(94), c("snp5-snp4")))

  observed <- distConsecutiveSnps(snp.coords=snp.coords)

  expect_equal(observed, expected)
})

test_that("thinSnps_index", {
  P <- 5 # SNPs
  snp.coords <- data.frame(chr=c(rep("chr1", P-2), rep("chr2", 2)),
                           coord=c(17, 2, 28, 5, 100))
  rownames(snp.coords) <- paste0("snp", 1:P)

  threshold <- 2

  expected <- c("snp2", "snp3", "snp4")

  observed <- thinSnps(method="index", threshold=threshold,
                       snp.coords=snp.coords)

  expect_equal(observed, expected)
})

test_that("thinSnps_coord", {
  P <- 5 # SNPs
  snp.coords <- data.frame(chr=c(rep("chr1", P-2), rep("chr2", 2)),
                           coord=c(17, 2, 28, 5, 100))
  rownames(snp.coords) <- paste0("snp", 1:P)

  threshold <- 15

  expected <- c("snp2", "snp1", "snp4", "snp5")

  observed <- thinSnps(method="coord", threshold=threshold,
                       snp.coords=snp.coords)

  expect_equal(observed, expected)
})
