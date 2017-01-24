library(rutilstimflutre)
context("Quantgen")

test_that("reformatGenoClasses", {
  N <- 3 # individuals
  P <- 4 # SNPs
  genoClasses <- data.frame(ind1=c("AA", "GC", "AA", "GC"),
                            ind2=c("AA", "GC", "AT", "??"),
                            ind3=c("TA", "UU", "AA", "CG"),
                            row.names=paste0("snp", 1:P))

  expected <- matrix(c("AA", "CG", "AA", "CG",
                       "AA", "CG", "AT", NA,
                       "AT", NA, "AA", "CG"),
                     byrow=TRUE,
                     nrow=N,
                     ncol=P,
                     dimnames=list(colnames(genoClasses),
                                   rownames(genoClasses)))

  observed <- reformatGenoClasses(x=genoClasses, na.string="??", verbose=0)

  expect_equal(observed, expected)
})

test_that("genoClasses2genoDoses", {
  N <- 3 # individuals
  P <- 4 # SNPs
  genoClasses <- data.frame(snp=paste0("snp", 1:P),
                            ind1=c("AA", "GC", "AA", "GC"),
                            ind2=c("AA", "GC", "AT", "??"),
                            ind3=c("TT", "GG", "AA", "CC"),
                            stringsAsFactors=FALSE)

  X <- matrix(c(0,0,2,
                1,1,0,
                0,1,0,
                1,NA,0),
              nrow=N, ncol=P,
              dimnames=list(colnames(genoClasses)[-1], genoClasses[,1]))
  alleles <- matrix(c("A","T", "G","C", "A","T", "C","G"),
                    byrow=TRUE, nrow=P, ncol=2,
                    dimnames=list(colnames(X), c("major", "minor")))
  expected <- list(geno.doses=X,
                   alleles=alleles)

  observed <- genoClasses2genoDoses(x=genoClasses, na.string="??", verbose=0)

  expect_equal(observed, expected)
})

test_that("updateJoinMap", {
  N <- 4 # individuals
  P <- 5 # loci
  x <- data.frame(locus=paste0("snp", 1:P),
                  seg=c("<abxcd>", "<abxac>", "<abxab>", "<abxaa>", "<aaxab>"),
                  phase=c("{00}", "{01}", "{10}", "{0-}", "{-1}"),
                  clas=NA,
                  ind1=c("ac", "aa", "aa", "aa", "aa"),
                  ind2=c("ad", "ac", "ab", "aa", "ab"),
                  ind3=c("bc", "ba", "ba", "ba", "aa"),
                  ind4=c("bd", "bc", "bb", "ba", "--"),
                  stringsAsFactors=FALSE)

  expected <- data.frame(locus=paste0("snp", 1:P),
                         seg=c("<abxcd>", "<efxeg>", "<hkxhk>", "<lmxll>", "<nnxnp>"),
                         phase=x$phase,
                         clas=NA,
                         ind1=c("ac", "ee", "hh", "ll", "nn"),
                         ind2=c("ad", "eg", "hk", "ll", "np"),
                         ind3=c("bc", "fe", "kh", "ml", "nn"),
                         ind4=c("bd", "fg", "kk", "ml", "--"),
                         stringsAsFactors=FALSE)

  observed <- updateJoinMap(x=x, verbose=0)

  expect_equal(observed, expected)
})

test_that("genoClasses2JoinMap", {
  nb.offs <- 4 # offsprings
  N <- 2 + nb.offs
  P <- 5 # SNPs
  x <- data.frame(par1=c("AA", "GC", "CG", "AT", NA),
                  par2=c("AT", "GC", "GG", "AT", "AT"),
                  off1=c("AA", "GG", "CG", "AA", "AA"),
                  off2=c("AT", "GG", "CG", "AT", "AT"),
                  off3=c("AT", "GG", "GG", "TT", "TT"),
                  off4=c(NA, NA, NA, NA, NA),
                  row.names=paste0("snp", 1:P),
                  stringsAsFactors=FALSE)

  expected <- data.frame(par1=c("nn", "CG", "lm", "hk", NA),
                         par2=c("np", "CG", "ll", "hk", "AT"),
                         p1.A=c("A", "C", "C", "A", NA),
                         p1.B=c("A", "G", "G", "T", NA),
                         p2.C=c("A", "C", "G", "A", "A"),
                         p2.D=c("T", "G", "G", "T", "T"),
                         seg=c("<nnxnp>", NA, "<lmxll>", "<hkxhk>", NA),
                         off1=c("nn", "GG", "lm", "hh", "AA"),
                         off2=c("np", "GG", "lm", "hk", "AT"),
                         off3=c("np", "GG", "ll", "kk", "TT"),
                         off4=as.character(c(NA, NA, NA, NA, NA)),
                         row.names=rownames(x),
                         stringsAsFactors=FALSE)

  observed <- genoClasses2JoinMap(x=x, reformat.input=TRUE, verbose=0)

  expect_equal(observed, expected)
})

test_that("filterSegreg", {
  nb.offs <- 6 # offsprings
  N <- 2 + nb.offs
  P <- 4 # SNPs
  x <- data.frame(seg=c("<nnxnp>", NA, "<lmxll>", "<hkxhk>"),
                  off1=c("nn", "GG", "lm", "hh"),
                  off2=c("np", "GG", "lm", "hk"),
                  off3=c("np", "GG", "ll", "kk"),
                  off4=c("np", "GG", "ll", "kk"),
                  off5=c("nn", "GG", "ll", "hk"),
                  off6=c("np", NA, NA, "hk"),
                  row.names=paste0("snp", 1:P),
                  stringsAsFactors=FALSE)

  tmp <- data.frame(nb.classes=c(2, NA, 2, 3),
                    class1=c("nn", NA, "ll", "hh"),
                    class2=c("np", NA, "lm", "hk"),
                    class3=c(NA, NA, NA, "kk"),
                    class4=c(NA, NA, NA, NA),
                    obs1=c(2, NA, 3, 1),
                    obs2=c(4, NA, 2, 3),
                    obs3=c(NA, NA, NA, 2),
                    obs4=c(NA, NA, NA, NA),
                    exp1=c(0.5*6, NA, 0.5*5, 0.25*6),
                    exp2=c(0.5*6, NA, 0.5*5, 0.5*6),
                    exp3=c(NA, NA, NA, 0.25*6),
                    exp4=c(NA, NA, NA, NA),
                    row.names=rownames(x))
  tmp$chi2 <- c((tmp$obs1[1] - tmp$exp1[1])^2 / tmp$exp1[1] +
                (tmp$obs2[1] - tmp$exp2[1])^2 / tmp$exp2[1],
                NA,
                (tmp$obs1[3] - tmp$exp1[3])^2 / tmp$exp1[3] +
                (tmp$obs2[3] - tmp$exp2[3])^2 / tmp$exp2[3],
                (tmp$obs1[4] - tmp$exp1[4])^2 / tmp$exp1[4] +
                (tmp$obs2[4] - tmp$exp2[4])^2 / tmp$exp2[4] +
                (tmp$obs3[4] - tmp$exp3[4])^2 / tmp$exp3[4])
  tmp$pvalue <- pchisq(q=tmp$chi2, df=tmp$nb.classes - 1, lower.tail=FALSE)
  expected <- as.matrix(tmp[, c("chi2", "pvalue")])

  observed <- filterSegreg(x=x, verbose=0)

  expect_equal(observed, expected)
})

test_that("genoDoses2genoClasses", {
  N <- 2 # individuals
  P <- 4 # SNPs
  X <- matrix(c(1,1, NA,0, 2,1, 1,NA), nrow=N, ncol=P,
              dimnames=list(paste0("ind", 1:N), paste0("snp", 1:P)))
  alleles <- data.frame(first=c("T","T","A","C"),
                        second=c("A","A","T","G"),
                        stringsAsFactors=FALSE)
  rownames(alleles) <- colnames(X)

  expected <- data.frame(ind1=c("TA", "??", "TT", "CG"),
                         ind2=c("TA", "TT", "AT", "??"),
                         stringsAsFactors=FALSE)
  rownames(expected) <- colnames(X)

  observed <- genoDoses2genoClasses(tX=t(X), alleles=alleles, na.string="??",
                                    verbose=0)

  expect_equal(observed, expected)
})

test_that("haplosAlleles2num", {
  nb.inds <- 2
  nb.snps <- 3
  haplos <- matrix(data=c("T","T","T","T", "T","A","T","A", "A","T","A","A"),
                   nrow=2*nb.inds, ncol=nb.snps,
                   dimnames=list(c("ind1_h1","ind1_h2","ind2_h1","ind2_h2"),
                                 paste0("snp", 1:nb.snps)))
  alleles <- as.data.frame(matrix(data=c(rep("A", nb.snps), rep("T", nb.snps)),
                                  nrow=nb.snps, ncol=2,
                                  dimnames=list(paste0("snp", 1:nb.snps),
                                                c("minor", "major"))))

  expected <- matrix(data=c(0,0,0,0, 0,1,0,1, 1,0,1,1),
                     nrow=2*nb.inds, ncol=nb.snps,
                     dimnames=dimnames(haplos))

  observed <- haplosAlleles2num(haplos=haplos, alleles=alleles)

  expect_equal(observed, expected)
})

test_that("calcFreqMissSnpGenosPerSnp", {
  N <- 2 # individuals
  P <- 4 # SNPs
  X <- matrix(c(1,1, NA,NA, 2,1, 1,NA), nrow=N, ncol=P,
              dimnames=list(paste0("ind", 1:N), paste0("snp", 1:P)))

  expected <- setNames(c(0/N, 2/N, 0/N, 1/N), colnames(X))

  observed <- calcFreqMissSnpGenosPerSnp(X)

  expect_equal(observed, expected)
})

test_that("calcFreqMissSnpGenosPerGeno", {
  N <- 2 # individuals
  P <- 4 # SNPs
  X <- matrix(c(1,1, NA,NA, 2,1, 1,NA), nrow=N, ncol=P,
              dimnames=list(paste0("ind", 1:N), paste0("snp", 1:P)))

  expected <- setNames(c(1/P, 2/P), rownames(X))

  observed <- calcFreqMissSnpGenosPerGeno(X)

  expect_equal(observed, expected)
})

test_that("discardMarkersMissGenos", {
  N <- 2 # individuals
  P <- 4 # SNPs
  X <- matrix(c(1,1, NA,NA, 2,1, 1,NA), nrow=N, ncol=P,
              dimnames=list(paste0("ind", 1:N), paste0("snp", 1:P)))

  expected <- X[, -c(2,4)]

  observed <- discardMarkersMissGenos(X=X, verbose=0)

  expect_equal(observed, expected)
})

test_that("estimSnpAf", {
  N <- 2 # individuals
  P <- 4 # SNPs
  X <- matrix(c(1,1, 1,0, 2,1, 1,NA), nrow=N, ncol=P,
              dimnames=list(paste0("ind", 1:N), paste0("snp", 1:P)))

  expected <- setNames(c(2/4, 1/4, 3/4, 1/2), colnames(X))

  observed <- estimSnpAf(X=X)

  expect_equal(observed, expected)
})

test_that("estimSnpMaf", {
  N <- 2 # individuals
  P <- 4 # SNPs
  X <- matrix(c(1,1, 1,0, 2,1, 1,NA), nrow=N, ncol=P,
              dimnames=list(paste0("ind", 1:N), paste0("snp", 1:P)))

  expected <- setNames(c(2/4, 1/4, 1/4, 1/2), colnames(X))

  observed <- estimSnpMaf(X=X)

  expect_equal(observed, expected)

  afs <- setNames(c(2/4, 1/4, 3/4, 1/2), colnames(X))
  observed <- estimSnpMaf(afs=afs)

  expect_equal(observed, expected)
})

test_that("discardSnpsLowMaf", {
  N <- 3 # individuals
  P <- 4 # SNPs
  X <- matrix(c(0,0,2, 1,1,0, 0,1,0, 1,0,0), nrow=N, ncol=P,
              dimnames=list(paste0("ind", 1:N), paste0("snp", 1:P)))

  observed <- discardSnpsLowMaf(X=X, mafs=NULL, thresh=0, verbose=0)
  expected <- X
  expect_equal(observed, expected)

  observed <- discardSnpsLowMaf(X=X, mafs=NULL, thresh=0.2, verbose=0)
  expected <- X[, c("snp1", "snp2")]
  expect_equal(observed, expected)

  mafs <- estimSnpMaf(X)
  observed <- discardSnpsLowMaf(X=X, mafs=mafs, thresh=0.2, verbose=0)
  expected <- X[, c("snp1", "snp2")]
  expect_equal(observed, expected)
})

test_that("recodeGenosMinorSnpAllele", {
  N <- 4 # individuals
  P <- 3 # SNPs
  X <- matrix(c(2,0,2,1, 1,1,0,1, 1,NA,0,2), nrow=N, ncol=P,
              dimnames=list(paste0("ind", 1:N), paste0("snp", 1:P)))
  alleles <- data.frame(A=c("A", "G", "T"),
                        B=c("T", "C", "A"),
                        row.names=colnames(X),
                        stringsAsFactors=FALSE)

  expected <- list(X=matrix(c(0,2,0,1, 1,1,0,1, 1,NA,0,2), nrow=N, ncol=P,
                            dimnames=dimnames(X)),
                   alleles=data.frame(major=c("T", "G", "T"),
                                      minor=c("A", "C", "A"),
                                      row.names=colnames(X),
                                      stringsAsFactors=FALSE))

  observed <- recodeGenosMinorSnpAllele(X=X, alleles=alleles, verbose=0)

  expect_equal(observed, expected)
})

test_that("countGenotypicClasses", {
  N <- 3 # individuals
  P <- 4 # SNPs
  X <- matrix(c(0,0,2, 1,1,0, 0,1,0, 1,NA,0), nrow=N, ncol=P,
              dimnames=list(paste0("ind", 1:N), paste0("snp", 1:P)))

  expected <- matrix(c(2,0,1,0, 1,2,0,0, 2,1,0,0, 1,1,0,1), byrow=TRUE,
                     nrow=P, ncol=4,
                     dimnames=list(colnames(X), c("0", "1", "2", "NA")))

  observed <- countGenotypicClasses(X=X)

  expect_equal(observed, expected)
})

test_that("chiSqSnpGenos", {
  N <- 3 # individuals
  P <- 5 # SNPs
  X <- matrix(c(0,0,2, 1,1,0, 0,1,0, 1,NA,0, 0,0,0), nrow=N, ncol=P,
              dimnames=list(paste0("ind", 1:N), paste0("snp", 1:P)))
  alleles <- data.frame(first=c("A", "G", "T", "T", "C"),
                        second=c("T", "C", "A", "A", "G"),
                        row.names=colnames(X),
                        stringsAsFactors=FALSE)

  ## external check
  if(FALSE){
    library(HardyWeinberg) # available on CRAN
    cts <- countGenotypicClasses(X=X)[, -4]
    colnames(cts) <- c("AA","AB","BB")
    suppressWarnings(HWChisqMat(X=cts, cc=0, verbose=FALSE))
    suppressWarnings(HWChisqMat(X=cts, cc=0.5, verbose=FALSE))
  }

  n.AB <- c(0, 2, 1, 1, 0)
  n <- c(N, N, N, N-1, N)
  p <- c(2, 2, 1, 1, 0) / (2 * n)
  q <- 1 - p
  e.AB <- 2 * n * p * q
  D <- 0.5 * (n.AB - e.AB)
  expected <- setNames(D^2 / (p^2 * q^2 * n),
                       colnames(X))

  observed <- chiSqSnpGenos(X=X, c=0, calc.with.D=TRUE)
  expect_equal(observed[,"chi2"], expected)

  observed <- chiSqSnpGenos(X=X, c=0, calc.with.D=FALSE)
  expect_equal(observed[,"chi2"], expected)

  c <- 0.5
  expected <- setNames((D^2 - 2*c*abs(D)*(1-p*q) + c^2*(1-(3/2)*p*q)) /
                       (p^2 * q^2 * n),
                       colnames(X))
  observed <- chiSqSnpGenos(X=X, c=c, thresh.c=0, calc.with.D=TRUE)
  expect_equal(observed[,"chi2"], expected)
  observed <- chiSqSnpGenos(X=X, c=c, thresh.c=0, calc.with.D=FALSE)
  expect_equal(observed[,"chi2"], expected)
})

test_that("imputeGenosWithMean", {
  N <- 4 # individuals
  P <- 4 # SNPs
  X <- matrix(c(0,1,2,NA,
                1,1,0,NA,
                0,0,0,NA,
                0,NA,NA,NA),
              nrow=N, ncol=P,
              dimnames=list(paste0("ind", 1:N), paste0("snp", 1:P)))

  expected <- matrix(c(0,1,2,mean(c(0,1,2)),
                       1,1,0,mean(c(1,1,0)),
                       0,0,0,NA,
                       0,NA,NA,NA),
                     nrow=N, ncol=P,
                     dimnames=list(paste0("ind", 1:N), paste0("snp", 1:P)))
  observed <- imputeGenosWithMean(X=X, min.maf=0.1, max.miss=0.3,
                                  rm.still.miss=FALSE)
  expect_equal(observed, expected)

  expected <- expected[, -c(3,4)]
  observed <- imputeGenosWithMean(X=X, min.maf=0.1, max.miss=0.3,
                                  rm.still.miss=TRUE)
  expect_equal(observed, expected)
})

test_that("recodeIntoDominant", {
  N <- 3 # individuals
  P <- 4 # SNPs
  X <- matrix(c(0,0,2, 1,1,0, 0,1,0, 1,0,2), nrow=N, ncol=P,
              dimnames=list(paste0("ind", 1:N), paste0("snp", 1:P)))

  expected <- matrix(c(0,0,0, 1,1,0, 0,1,0, 1,0,0), nrow=N, ncol=P,
                     dimnames=dimnames(X))

  observed <- recodeIntoDominant(X=X)

  expect_equal(observed, expected)
})

test_that("recodeIntoDominant_imputed", {
  N <- 3 # individuals
  P <- 4 # SNPs
  X <- matrix(c(0,0,2, 1.1,1,0, 0.2,1,0, 1,0,2), nrow=N, ncol=P,
              dimnames=list(paste0("ind", 1:N), paste0("snp", 1:P)))

  expected <- matrix(c(0,0,0, 1,1,0, 0,1,0, 1,0,0), nrow=N, ncol=P,
                     dimnames=dimnames(X))

  observed <- recodeIntoDominant(X=X, simplify.imputed=TRUE)

  expect_equal(observed, expected)
})

test_that("estimGenRel_vanraden1", {
  N <- 3 # individuals
  P <- 4 # SNPs
  X <- matrix(c(0,0,2, 1,1,0, 0,1,0, 1,0,0), nrow=N, ncol=P,
              dimnames=list(paste0("ind", 1:N), paste0("snp", 1:P)))

  afs <- setNames(c(0.383, 0.244, 0.567, 0.067), colnames(X))

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

test_that("estimGenRel_vanraden1_wNA", {
  N <- 3 # individuals
  P <- 5 # SNPs
  X <- matrix(c(0,0,2, 1,1,0, 0,1,0, 1,0,0, 1,NA,0), nrow=N, ncol=P,
              dimnames=list(paste0("ind", 1:N), paste0("snp", 1:P)))

  afs <- setNames(c(0.383, 0.244, 0.567, 0.067, 0.17), colnames(X))

  X2 <- discardMarkersMissGenos(X=X, verbose=0)
  P2 <- ncol(X2)
  afs2 <- afs[colnames(X2)]
  tmp <- matrix(rep(2*(afs2 - 0.5), N), nrow=N, ncol=P2, byrow=TRUE)
  Z <- X2 - 1 - tmp
  denom <- 0
  for(p in 1:P2)
    denom <- denom + afs2[p] * (1 - afs2[p])
  denom <- 2 * denom
  expected <- (Z %*% t(Z)) / denom

  observed <- estimGenRel(X=X, afs=afs, thresh=0, relationships="additive",
                          method="vanraden1", verbose=0)

  expect_equal(observed, expected)
})

test_that("estimGenRel_vanraden1_wMAF", {
  N <- 3 # individuals
  P <- 4 # SNPs
  X <- matrix(c(0,0,2, 1,1,0, 0,1,0, 1,0,0), nrow=N, ncol=P,
              dimnames=list(paste0("ind", 1:N), paste0("snp", 1:P)))

  afs <- setNames(c(0.383, 0.244, 0.167, 0.067), colnames(X))

  thresh <- 0.1
  mafs <- estimSnpMaf(afs=afs)
  afs2 <- afs[afs >= thresh]
  X2 <- X[, names(afs2)]
  P <- ncol(X2)
  tmp <- matrix(rep(2*(afs2-0.5), N), nrow=N, ncol=P, byrow=TRUE)
  Z <- X2 - 1 - tmp
  denom <- 0
  for(p in 1:P)
    denom <- denom + afs2[p] * (1 - afs2[p])
  denom <- 2 * denom
  expected <- (Z %*% t(Z)) / denom

  observed <- estimGenRel(X=X, afs=afs, thresh=thresh, relationships="additive",
                          method="vanraden1", verbose=0)

  expect_equal(observed, expected)
})

test_that("estimGenRel_astle-balding", {
  N <- 2 # individuals
  P <- 4 # SNPs
  X <- matrix(c(1,1, 1,0, 2,1, 1,0), nrow=N, ncol=P,
              dimnames=list(paste0("ind", 1:N), paste0("snp", 1:P)))

  afs <- estimSnpAf(X=X)

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

  afs <- estimSnpAf(X=X)

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

  observed <- estimGenRel(X=X, afs=afs, thresh=0, relationships="dominant",
                          method="su", verbose=0)

  expect_equal(observed, expected)
})

test_that("estimGenRel_vitezica", {
  N <- 3 # individuals
  P <- 4 # SNPs
  X <- matrix(c(0,0,2, 1,1,0, 0,1,0, 1,0,0), nrow=N, ncol=P,
              dimnames=list(paste0("ind", 1:N), paste0("snp", 1:P)))

  afs <- setNames(c(0.383, 0.244, 0.167, 0.067), colnames(X))

  W <- matrix(c(- 2 * afs[1]^2,
                - 2 * afs[1]^2,
                - 2 * (1 - afs[1])^2,

                2 * afs[2] * (1 - afs[2]),
                2 * afs[2] * (1 - afs[2]),
                - 2 * afs[2]^2,

                - 2 * afs[3]^2,
                2 * afs[3] * (1 - afs[3]),
                - 2 * afs[3]^2,

                2 * afs[4] * (1 - afs[4]),
                - 2 * afs[4]^2,
                - 2 * afs[4]^2),
              nrow=N, ncol=P,
              dimnames=list(rownames(X), colnames(X)))
  denom <- 0
  for(p in 1:P)
    denom <- denom + (2 * afs[p] * (1 - afs[p]))^2
  expected <- (W %*% t(W)) / denom

  observed <- estimGenRel(X=X, afs=afs, thresh=0, relationships="dominant",
                          method="vitezica", verbose=0)

  expect_equal(observed, expected)
})

test_that("haplosList2Matrix", {
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

  expected <- matrix(c(1,1,1,0, 1,0,0,1, 0,0,1,0),
                     nrow=2*N, ncol=P,
                     dimnames=list(rownames(haplos[[1]]),
                                   paste0("snp", 1:P)))

  observed <- haplosList2Matrix(haplos)

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
                             dimnames=list(c("ind1_h1", "ind1_h2", "ind10_h1",
                                             "ind10_h2"),
                                           c("snp1", "snp2"))),
                 chr2=matrix(data=c(0,0,1,0), nrow=2*N, ncol=1,
                             dimnames=list(c("ind1_h1", "ind1_h2", "ind10_h1",
                                             "ind10_h2"),
                                           c("snp3"))))

  expected <- list(chr1=matrix(data=c(1,1, 1,0), nrow=2, ncol=P-1,
                               dimnames=list(c("ind1_h1", "ind1_h2"),
                                             c("snp1", "snp2"))),
                   chr2=matrix(data=c(0,0), nrow=2, ncol=1,
                               dimnames=list(c("ind1_h1", "ind1_h2"),
                                             c("snp3"))))

  observed <- getHaplosInd(haplos, "ind1")

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
  if(all(requireNamespace("LDcorSV"))){
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

    observed <- estimLd(X=X, snp.coords=snp.coords, only.chr="chr1",
                        verbose=0)

    expect_equal(observed, expected)

    ## check that LDcorSV returns the same results
    colnames(expected)[3] <- "r2"
    observed <- estimLd(X=X, snp.coords=snp.coords, only.chr="chr1",
                        use.ldcorsv=TRUE, verbose=0)
    expect_equal(observed, expected)
  }
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

test_that("mme", {
  ## see Mrode (2005), section 3.2, example 3.1 page 42
  X <- matrix(c(1,0,0,1,1,
                0,1,1,0,0),
              ncol=2)
  Z <- matrix(c(0,0,0,1,0,0,0,0,
                0,0,0,0,1,0,0,0,
                0,0,0,0,0,1,0,0,
                0,0,0,0,0,0,1,0,
                0,0,0,0,0,0,0,1),
              ncol=8, byrow=TRUE)
  y <- matrix(c(4.5, 2.9, 3.9, 3.5, 5.0), ncol=1)

  sigma.a.2 <- 20
  sigma.e.2 <- 40
  alpha <- sigma.e.2 / sigma.a.2
  Ainv <- matrix(c(1.833,0.5,0,-0.667,0,-1,0,0,
                   0.5,2,0.5,0,-1,-1,0,0,
                   0,0.5,2,0,-1,0.5,0,-1,
                   -0.667,0,0,1.833,0.5,0,-1,0,
                   0,-1,-1,0.5,2.5,0,-1,0,
                   -1,-1,0.5,0,0,2.5,0,-1,
                   0,0,0,-1,-1,0,2,0,
                   0,0,-1,0,0,-1,0,2),
                 ncol=8, byrow=TRUE)
  G <- sigma.a.2 * solve(Ainv)
  R <- sigma.e.2 * diag(length(y))
  V <- Z %*% G %*% t(Z) + R
  Vinv <- solve(V)

  b.hat <- solve(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv %*% y
  a.hat <- G %*% t(Z) %*% Vinv %*% (y - X %*% b.hat)

  expected <- c(b.hat, a.hat)

  observed <- mme(y=y, W=X, Z=Z, sigma.A2=sigma.a.2, Ainv=Ainv, V.E=sigma.e.2)

  expect_equal(format(observed, digits=1), format(expected, digits=1))
})


test_that("rearrangeInputsForAssoGenet", {
  ids <- data.frame(cultivar.code=c(34, 150, 19),
                    cultivar.name=c("Grenache", "Syrah", "Carigan"),
                    accession.code=c("34Mtp6", "150Mtp11", "18Mtp17"),
                    stringsAsFactors=FALSE)
  y <- setNames(c(1.72, 0.98), c("18Mtp17", "34Mtp6"))
  X <- matrix(c(1,0, 2,0), nrow=2, ncol=2,
              dimnames=list(c("34Mtp6", "150Mtp11"), c("snp2", "snp5")))
  snp.coords <- data.frame(chr=c("chr7", "chr2", "chr2", "chr3"),
                           coord=c(726354, 12536, 18700, 763542),
                           row.names=c("snp1", "snp2", "snp5", "snp8"))
  alleles <- data.frame(minor=c("A", "G", "C"),
                        major=c("T", "C", "G"),
                        row.names=c("snp2", "snp5", "snp8"),
                        stringsAsFactors=FALSE)

  cultivars.tokeep <- c("34")
  snps.tokeep <- c("snp2", "snp5")
  exp <- list(ids=ids[1, , drop=FALSE],
              y=data.frame(y=y["34Mtp6"],
                           row.names=cultivars.tokeep,
                           stringsAsFactors=FALSE),
              X=X["34Mtp6", snps.tokeep, drop=FALSE],
              snp.coords=droplevels(snp.coords[snps.tokeep,]),
              alleles=droplevels(alleles[snps.tokeep,]))
  exp$ids$cultivar.code <- as.character(exp$ids$cultivar.code)
  rownames(exp$X) <- "34"

  obs <- rearrangeInputsForAssoGenet(ids=ids, y=y, X=X, snp.coords=snp.coords,
                                     alleles=alleles, verbose=0)
  expect_equal(obs$ids, exp$ids)
  expect_equal(obs$y, exp$y)
  expect_equal(obs$X, exp$X)
  expect_equal(obs$snp.coords, exp$snp.coords)
  expect_equal(obs$alleles, exp$alleles)
})
