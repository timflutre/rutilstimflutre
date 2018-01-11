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
                     nrow=P,
                     ncol=N,
                     dimnames=dimnames(genoClasses))

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
  P <- 8 # SNPs
  x <- data.frame(par1=c("AA", "GC", "CG", "AT", NA, "AA", "AA", "GG"),
                  par2=c("AT", "GC", "GG", "AT", "AT", "AT", "TT", "GG"),
                  off1=c("AA", "GG", "CG", "AA", "AA", "AT", "AT", "GG"),
                  off2=c("AT", "GG", "CG", "AT", "AT", "AA", "AT", "GG"),
                  off3=c("AT", "GG", "GG", "TT", "TT", NA, "AT", "GG"),
                  off4=c(NA, NA, NA, NA, NA, NA, NA, NA),
                  row.names=paste0("snp", 1:P),
                  stringsAsFactors=FALSE)

  expected <- data.frame(par1=c("nn", "CG", "lm", "hk", NA, "AA", "AA", "GG"),
                         par2=c("np", "CG", "ll", "hk", "AT", "AT", "TT", "GG"),
                         p1.A=c("A", "C", "C", "A", NA, "A", "A", "G"),
                         p1.B=c("A", "G", "G", "T", NA, "A", "A", "G"),
                         p2.C=c("A", "C", "G", "A", "A", "A", "T", "G"),
                         p2.D=c("T", "G", "G", "T", "T", "T", "T", "G"),
                         seg.pars=c("<nnxnp>", "<hkxhk>", "<lmxll>", "<hkxhk>",
                                    NA, "<nnxnp>", NA, NA),
                         seg.offs=c("<lmxll>_<nnxnp>", NA, "<lmxll>_<nnxnp>",
                                    "<hkxhk>", "<hkxhk>", NA, NA, NA),
                         seg=c("<nnxnp>", NA, "<lmxll>", "<hkxhk>", NA, NA, NA, NA),
                         off1=c("nn", "GG", "lm", "hh", "AA", "AT", "AT", "GG"),
                         off2=c("np", "GG", "lm", "hk", "AT", "AA", "AT", "GG"),
                         off3=c("np", "GG", "ll", "kk", "TT", NA, "AT", "GG"),
                         off4=as.character(c(NA, NA, NA, NA, NA, NA, NA, NA)),
                         row.names=rownames(x),
                         stringsAsFactors=FALSE)

  observed <- genoClasses2JoinMap(x=x, reformat.input=TRUE, thresh.na=2,
                                  verbose=0)

  expect_equal(observed, expected)
})

test_that("genoClasses2JoinMap_F2", {
  nb.offs <- 4 # offsprings
  N <- 2 + 1 + nb.offs
  P <- 11 # SNPs
  x <- data.frame(par1=c("AA", "AA", # 1. AAxAT=AA; AAxAT=AT
                         "GC", "GC", "GC", # 2. GCxGC=GG; GCxGC=GC; GCxGC=CC
                         "CG", "GG", # 3. GCxGG=GG; GGxCG=GC
                         NA,         # 4. NAxAT=AT
                         "AA",       # 5. AAxAT=AT (and NA in off3)
                         "AA",       # 6. AAxTT=AT
                         "GG"),      # 7. GGxGG=GG
                  par2=c("AT", "AT",
                         "GC", "GC", "GC",
                         "GG", "CG",
                         "AT",
                         "AT",
                         "TT",
                         "GG"),
                  f1=c("AA", "AT",
                       "GG", "GC", "CC",
                       "GG", "GC",
                       "AT",
                       "AT",
                       "AT",
                       "GG"),
                  off1=c("AA", "AA",
                         "GG", "GG", "CC",
                         "GG", "GG",
                         "AA",
                         "AA",
                         "AA",
                         "GG"),
                  off2=c("AA", "AT",
                         "GG", "GC", "CC",
                         "GG", "GC",
                         "AT",
                         "AT",
                         "AT",
                         "GG"),
                  off3=c("AA", "TT",
                         "GG", "CC", "CC",
                         "GG", "CC",
                         "TT",
                         NA,
                         "TT",
                         "GG"),
                  off4=c(NA, NA,
                         NA, NA, NA,
                         NA, NA,
                         NA,
                         NA,
                         NA,
                         NA),
                  row.names=paste0("snp", 1:P),
                  stringsAsFactors=FALSE)

  expected <- data.frame(par1=c("AA", "AA",
                                "CG", "CG", "CG",
                                "CG", "GG",
                                NA,
                                "AA",
                                "AA",
                                "GG"),
                         par2=c("AT", "AT",
                                "CG", "CG", "CG",
                                "GG", "CG",
                                "AT",
                                "AT",
                                "TT",
                                "GG"),
                         p1.A=c("A", "A",
                                "C", "C", "C",
                                "C", "G",
                                NA,
                                "A",
                                "A",
                                "G"),
                         p1.B=c("A", "A",
                                "G", "G", "G",
                                "G", "G",
                                NA,
                                "A",
                                "A",
                                "G"),
                         p2.C=c("A", "A",
                                "C", "C", "C",
                                "G", "C",
                                "A",
                                "A",
                                "T",
                                "G"),
                         p2.D=c("T", "T",
                                "G", "G", "G",
                                "G", "G",
                                "T",
                                "T",
                                "T",
                                "G"),
                         seg.pars=c("A=A", "A=A",
                                    "A=G", "A=C", "A=C",
                                    "A=G", "A=G",
                                    NA,
                                    "A=A",
                                    "A=A",
                                    NA),
                         seg.offs=rep(NA, P),
                         seg=c("F2", "F2",
                               "F2", "F2", "F2",
                               "F2", "F2",
                               NA,
                               "F2",
                               "F2",
                               NA),
                         off1=c("A", "A",
                                "A", "B", "A",
                                "A", "A",
                                NA,
                                "A",
                                "A",
                                NA),
                         off2=c("A", "H",
                                "A", "H", "A",
                                "A", "H",
                                NA,
                                "H",
                                "H",
                                NA),
                         off3=c("A", "B",
                                "A", "A", "A",
                                "A", "B",
                                NA,
                                NA,
                                "B",
                                NA),
                         off4=as.character(rep(NA, P)),
                         row.names=rownames(x),
                         stringsAsFactors=FALSE)

  observed <- genoClasses2JoinMap(x=x[,-3], reformat.input=TRUE,
                                  thresh.na=NULL, thresh.counts=NULL,
                                  is.F2=TRUE, verbose=0)

  expect_equal(observed, expected)
})

test_that("readSegregJoinMap", {
  tmpd <- tempdir()

  jm.file <- paste0(tmpd, "/genos.loc")

  nb.locus <- 6
  nb.genos <- 4
  jm <- data.frame(seg=c("<abxcd>", "<efxeg>", "<hkxhk>", "<lmxll>", "<nnxnp>",
                         NA),
                   phase=rep("{??}", nb.locus),
                   geno1=c("ac", "ee", "hh", "ll", "nn", NA),
                   geno2=c("ad", "eg", "hk", "ll", "np", NA),
                   geno3=c("bc", "ef", "hk", "lm", "nn", NA),
                   geno4=c("bd", "fg", "kk", "lm", NA, NA),
                   row.names=paste0("loc", 1:nb.locus),
                   stringsAsFactors=FALSE)
  writeSegregJoinMap(pop.name="test", pop.type="CP",
                     locus=rownames(jm), segregs=jm$seg,
                     genos=jm[,-c(1:2)], phases=jm$phase,
                     file=jm.file, save.ind.names=TRUE,
                     na.string="--", verbose=0)

  expected <- jm

  observed <- readSegregJoinMap(file=jm.file, na.string="--", verbose=0)

  expect_equal(observed, expected)

  if(file.exists(jm.file))
    file.remove(jm.file)
})

test_that("joinMap2backcross_qtl", {
  nb.locus <- 6
  nb.genos <- 4
  jm <- data.frame(seg=c("<abxcd>", "<efxeg>", "<hkxhk>", "<lmxll>", "<nnxnp>",
                         NA),
                   phase=rep(NA, nb.locus),
                   clas=rep(NA, nb.locus),
                   geno1=c("ac", "ee", "hh", "ll", "nn", NA),
                   geno2=c("ad", "eg", "hk", "ll", "np", NA),
                   geno3=c("bc", "ef", "hk", "lm", "nn", NA),
                   geno4=c("bd", "fg", "kk", "lm", NA, NA),
                   row.names=paste0("loc", 1:nb.locus),
                   stringsAsFactors=FALSE)

  expected <- list(par1=matrix(data=NA, nrow=2 * (nb.locus - 2),
                               ncol=nb.genos),
                   par2=matrix(data=NA, nrow=2 * (nb.locus - 2),
                               ncol=nb.genos))
  rownames(expected$par1) <- c("loc1", "loc2", "loc3", "loc4",
                               "loc1_m", "loc2_m", "loc3_m", "loc4_m")
  colnames(expected$par1) <- paste0("geno", 1:nb.genos)
  expected$par1["loc1",] <- c(1, 1, 2, 2)
  expected$par1["loc2",] <- c(1, 1, 2, 2)
  expected$par1["loc3",] <- c(1, NA, NA, 2)
  expected$par1["loc4",] <- c(1, 1, 2, 2)
  expected$par1["loc1_m",] <- c(2, 2, 1, 1)
  expected$par1["loc2_m",] <- c(2, 2, 1, 1)
  expected$par1["loc3_m",] <- c(2, NA, NA, 1)
  expected$par1["loc4_m",] <- c(2, 2, 1, 1)
  rownames(expected$par2) <- c("loc1", "loc2", "loc3", "loc5",
                               "loc1_m", "loc2_m", "loc3_m", "loc5_m")
  colnames(expected$par2) <- paste0("geno", 1:nb.genos)
  expected$par2["loc1",] <- c(1, 2, 1, 2)
  expected$par2["loc2",] <- c(1, 2, 1, 2)
  expected$par2["loc3",] <- c(1, NA, NA, 2)
  expected$par2["loc5",] <- c(1, 2, 1, NA)
  expected$par2["loc1_m",] <- c(2, 1, 2, 1)
  expected$par2["loc2_m",] <- c(2, 1, 2, 1)
  expected$par2["loc3_m",] <- c(2, NA, NA, 1)
  expected$par2["loc5_m",] <- c(2, 1, 2, NA)

  observed <- joinMap2backcross(x=jm, alias.hom=1, alias.het=2, alias.dup=0,
                                alias.miss=NA, parent.names=c("par1", "par2"),
                                verbose=0)

  expect_equal(observed, expected)
})

test_that("joinMap2backcross_CarthaGene", {
  nb.locus <- 6
  nb.genos <- 4
  jm <- data.frame(locus=paste0("loc", 1:nb.locus),
                   seg=c("<abxcd>", "<efxeg>", "<hkxhk>", "<lmxll>", "<nnxnp>",
                         NA),
                   phase=rep(NA, nb.locus),
                   clas=rep(NA, nb.locus),
                   geno1=c("ac", "ee", "hh", "ll", "nn", ""),
                   geno2=c("ad", "eg", "hk", "ll", "np", ""),
                   geno3=c("bc", "ef", "hk", "lm", "nn", ""),
                   geno4=c("bd", "fg", "kk", "lm", NA, ""),
                   stringsAsFactors=FALSE)

  expected <- list(ind088=data.frame(geno1=c("A", "A", "A", "A"),
                                     geno2=c("A", "A", "-", "A"),
                                     geno3=c("H", "H", "-", "H"),
                                     geno4=c("H", "H", "H", "H"),
                                     row.names=paste0("loc", c(1,2,3,4)),
                                     stringsAsFactors=FALSE),
                   ind284=data.frame(geno1=c("A", "A", "A", "A"),
                                     geno2=c("H", "H", "-", "H"),
                                     geno3=c("A", "A", "-", "A"),
                                     geno4=c("H", "H", "H", "-"),
                                     row.names=paste0("loc", c(1,2,3,5)),
                                     stringsAsFactors=FALSE))
  expected$ind088 <- rbind(expected$ind088,
                            data.frame(geno1=c("H", "H", "H", "H"),
                                       geno2=c("H", "H", "-", "H"),
                                       geno3=c("A", "A", "-", "A"),
                                       geno4=c("A", "A", "A", "A"),
                                       row.names=paste0(rownames(expected$ind088),
                                                        "_m"),
                                       stringsAsFactors=FALSE))
  expected$ind284 <- rbind(expected$ind284,
                            data.frame(geno1=c("H", "H", "H", "H"),
                                       geno2=c("A", "A", "-", "A"),
                                       geno3=c("H", "H", "-", "H"),
                                       geno4=c("A", "A", "A", "-"),
                                       row.names=paste0(rownames(expected$ind284),
                                                        "_m"),
                                       stringsAsFactors=FALSE))
  expected$ind088 <- as.matrix(expected$ind088)
  expected$ind284 <- as.matrix(expected$ind284)

  observed <- joinMap2backcross(x=jm, alias.hom="A", alias.het="H",
                                parent.names=c("ind088", "ind284"),
                                verbose=0)

  expect_equal(observed, expected)
})

test_that("setJoinMapPhasesFromParentalLinkGroups", {
  nb.locus <- 9 # one of them (loc8) is NA
  ## 2 chromosomes:
  ## chr1: 4 locus, genetically ordered as:
  ##  loc5 (lmxll) <- only informative for parent 1
  ##  loc6 (nnxnp) <- only informative for parent 2
  ##  loc4 (hkxhk) <- informative for both parents
  ##  loc9 (lmxll) <- only informative for parent 1
  ## chr2: 4 locus, genetically ordered as:
  ##  loc1 (nnxnp) <- only informative for parent 2
  ##  loc7 (hkxhk) <- informative for both parents
  ##  loc3 (lmxll) <- only informative for parent 1
  ##  loc2 (nnxnp) <- only informative for parent 2
  ## no inversion between parent 1 and parent 2
  nb.genos <- 3
  jm <- data.frame(seg=c("<nnxnp>", "<nnxnp>", "<lmxll>", "<hkxhk>", "<lmxll>",
                         "<nnxnp>", "<hkxhk>", NA, "<lmxll>"),
                   geno1=c("nn", "np", "lm", "hk", "ll", "np", "kk", NA, "lm"),
                   geno2=c("np", "np", "ll", "hk", "ll", "np", "hk", NA, "ll"),
                   geno3=c("np", "nn", "ll", "hh", "lm", "np", "hk", NA, "lm"))
  rownames(jm) <- paste0("loc", 1:nb.locus)
  gm1 <- data.frame(linkage.group=c(rep("chr2", 2),
                                    rep("chr1", 3)),
                    locus=c("loc7", "loc3_m",
                            "loc5", "loc4_m", "loc9"),
                    genetic.distance=c(0.0, 8.1,
                                       0.0, 1.5, 3.9))
  gm2 <- data.frame(linkage.group=c(rep("chr1", 2),
                                    rep("chr2", 3)),
                    locus=c("loc6", "loc4_m",
                            "loc1_m", "loc7_m", "loc2"),
                    genetic.distance=c(0.0, 7.3,
                                       0.0, 4.2, 5.1))

  expected <- data.frame(seg=jm$seg,
                         phase=c("{-1}", "{-0}", "{1-}", "{11}", "{0-}",
                                 "{-0}", "{01}", NA, "{0-}"),
                         geno1=jm$geno1,
                         geno2=jm$geno2,
                         geno3=jm$geno3)
  rownames(expected) <- rownames(jm)
  expected <- convertFactorColumnsToCharacter(expected)

  observed <- setJoinMapPhasesFromParentalLinkGroups(x=jm,
                                                     lg.par1=gm1,
                                                     lg.par2=gm2)

  expect_equal(observed, expected)
})

test_that("phasedJoinMapCP2qtl", {
  nb.locus <- 17
  nb.genos <- 4
  jm <- data.frame(locus=paste0("loc", 1:nb.locus),
                   seg=c("<abxcd>", "<abxcd>", "<abxcd>", "<abxcd>",
                         "<efxeg>", "<efxeg>", "<efxeg>", "<efxeg>",
                         "<hkxhk>", "<hkxhk>", "<hkxhk>", "<hkxhk>",
                         "<lmxll>", "<lmxll>",
                         "<nnxnp>", "<nnxnp>",
                         NA),
                   phase=c("{00}", "{01}","{10}","{11}",
                           "{00}", "{01}","{10}","{11}",
                           "{00}", "{01}","{10}","{11}",
                           "{0-}", "{1-}",
                           "{-0}", "{-1}",
                           NA),
                   clas=rep(NA, nb.locus),
                   geno1=c("ac", "ad", "bc", "bd",
                           "ee", "eg", "ef", "fg",
                           "hh", "hk", "hk", "kk",
                           "ll", "lm",
                           "nn", "np",
                           NA),
                   geno2=c("ad", "ac", "bd", "bc",
                           "eg", "ee", "fg", "ef",
                           "hk", "hh", "kk", "hk",
                           "ll", "lm",
                           "np", "nn",
                           NA),
                   geno3=c("bc", "bd", "ac", "ad",
                           "ef", "fg", "ee", "eg",
                           "hk", "kk", "hh", "hk",
                           "lm", "ll",
                           "nn", "np",
                           NA),
                   geno4=c("bd", "bc", "ad", "ac",
                           "fg", "ef", "eg", "ee",
                           "kk", "hk", "hk", "hh",
                           "lm", "ll",
                           "np", "nn",
                           NA),
                   stringsAsFactors=FALSE)

  expected <- matrix(data=NA, nrow=nb.locus, ncol=nb.genos,
                     dimnames=list(paste0("loc", 1:nb.locus),
                                   paste0("geno", 1:nb.genos)))
  ## abxcd
  expected[1,] <- c(1, 3, 2, 4)
  expected[2,] <- c(1, 3, 2, 4)
  expected[3,] <- c(1, 3, 2, 4)
  expected[4,] <- c(1, 3, 2, 4)
  ## efxeg
  expected[5,] <- c(1, 3, 2, 4)
  expected[6,] <- c(1, 3, 2, 4)
  expected[7,] <- c(1, 3, 2, 4)
  expected[8,] <- c(1, 3, 2, 4)
  ## hkxhk
  expected[9,] <- c(1, 10, 10, 4)
  expected[10,] <- c(9, 3, 2, 9)
  expected[11,] <- c(9, 3, 2, 9)
  expected[12,] <- c(1, 10, 10, 4)
  ## lmxll
  expected[13,] <- c(5, 5, 6, 6)
  expected[14,] <- c(5, 5, 6, 6)
  ## nnxnp
  expected[15,] <- c(7, 8, 7, 8)
  expected[16,] <- c(7, 8, 7, 8)
  ## missing
  expected[17,] <- rep(NA, nb.genos)

  observed <- phasedJoinMapCP2qtl(x=jm, verbose=0)

  expect_equal(observed, expected)
})

test_that("joinMap2designMatrix_without-phase", {
  nb.locus <- 5
  nb.genos <- 4
  jm <- data.frame(locus=paste0("loc", 1:nb.locus),
                   seg=c("<abxcd>", "<efxeg>", "<hkxhk>", "<lmxll>", "<nnxnp>"),
                   phase=rep(NA, nb.locus),
                   clas=rep(NA, nb.locus),
                   geno1=c("ac", "ee", "hh", "ll", "nn"),
                   geno2=c("ad", "eg", "hk", "ll", "np"),
                   geno3=c("bc", "ef", "hk", "lm", "nn"),
                   geno4=c("bd", "fg", "kk", "lm", "np"),
                   stringsAsFactors=FALSE)

  expected <- matrix(data=0, nrow=nb.genos,
                     ncol=1+(4+4)+(3+4)+(2+3)+(2+2)+(2+2))
  rownames(expected) <- paste0("geno", 1:nb.genos)
  colnames(expected) <- c("intercept",
                          "loc1.a", "loc1.b", "loc1.c", "loc1.d",
                          "loc1.ac", "loc1.ad", "loc1.bc", "loc1.bd",
                          "loc2.e", "loc2.f", "loc2.g",
                          "loc2.ee", "loc2.ef", "loc2.eg", "loc2.fg",
                          "loc3.h", "loc3.k",
                          "loc3.hh", "loc3.hk", "loc3.kk",
                          "loc4.l", "loc4.m",
                          "loc4.ll", "loc4.lm",
                          "loc5.n", "loc5.p",
                          "loc5.nn", "loc5.np")
  expected[, "intercept"] <- 1
  expected["geno1", c(2,4,6, 10,13, 17,19, 22,24, 26,28)] <-
    c(1,1,1, 2,1, 2,1, 2,1, 2,1)
  expected["geno2", c(2,5,7, 10,12,15, 17,18,20, 22,24, 26,27,29)] <-
    c(1,1,1, 1,1,1, 1,1,1, 2,1, 1,1,1)
  expected["geno3", c(3,4,8, 10,11,14, 17,18,20, 22,23,25, 26,28)] <-
    c(1,1,1, 1,1,1, 1,1,1, 1,1,1, 2,1)
  expected["geno4", c(3,5,9, 11,12,16, 18,21, 22,23,25, 26,27,29)] <-
    c(1,1,1, 1,1,1, 2,1, 1,1,1, 1,1,1)

  observed <- joinMap2designMatrix(jm=jm, parameterization="allele",
                                   constraints=NULL, rm.col.zeros=TRUE,
                                   verbose=0)

  expect_equal(observed, expected)

  ## test rm.dom=TRUE
  expected <- expected[, -c(6:9, 13:16, 19:21, 24:25, 28:29)]
  observed <- joinMap2designMatrix(jm=jm, parameterization="allele",
                                   constraints=NULL, rm.col.zeros=TRUE,
                                   rm.dom=TRUE, verbose=0)
  expect_equal(observed, expected)
})

test_that("joinMap2designMatrix_with-phase", {
  nb.locus <- 10
  nb.genos <- 4
  jm <- data.frame(locus=paste0("loc", 1:nb.locus),
                   seg=c("<abxcd>", "<abxcd>", "<abxcd>", "<abxcd>",
                         "<efxeg>", "<hkxhk>",
                         "<lmxll>", "<lmxll>", "<nnxnp>", "<nnxnp>"),
                   phase=c("{00}", "{01}", "{10}", "{11}",
                           "{00}", "{11}",
                           "{0-}", "{1-}", "{-1}", "{-0}"),
                   clas=rep(NA, nb.locus),
                   geno1=c("ac", "ac", "ac", "ac",
                           "ee", "hh",
                           "ll", "ll", "nn", "nn"),
                   geno2=c("ad", "ad", "ad", "ad",
                           "eg", "hk",
                           "ll", "ll", "np", "np"),
                   geno3=c("bc", "bc", "bc", "bc",
                           "ef", "hk",
                           "lm", "lm", "nn", "nn"),
                   geno4=c("bd", "bd", "bd", "bd",
                           "fg", "kk",
                           "lm", "lm", "np", "np"),
                   stringsAsFactors=FALSE)

  expected <- matrix(data=0, nrow=nb.genos, ncol=1 + (4 * nb.locus))
  rownames(expected) <- paste0("geno", 1:nb.genos)
  colnames(expected) <- c("intercept",
                          paste0(rep(paste0("loc", 1:nb.locus), each=4),
                                 rep(paste0(rep(paste0(".par", 1:2), each=2),
                                            rep(paste0(".haplo", 1:2), 2)), nb.locus)))
  expected[, "intercept"] <- 1
  expected[, grep("loc1.par[12].haplo1", colnames(expected))] <- 1
  expected[, grep("loc2.par1.haplo1|loc2.par2.haplo2", colnames(expected))] <- 1
  expected[, grep("loc3.par1.haplo2|loc3.par2.haplo1", colnames(expected))] <- 1
  expected[, grep("loc4.par[12].haplo2", colnames(expected))] <- 1
  expected[, grep("loc5.par[12].haplo1", colnames(expected))] <- 1
  expected[, grep("loc6.par[12].haplo2", colnames(expected))] <- 1
  expected[, grep("loc7.par1.haplo1", colnames(expected))] <- 1
  expected[, grep("loc8.par1.haplo2", colnames(expected))] <- 1
  expected[, grep("loc9.par2.haplo2", colnames(expected))] <- 1
  expected[, grep("loc10.par2.haplo1", colnames(expected))] <- 1

  observed <- joinMap2designMatrix(jm=jm, use.phase=TRUE, rm.col.zeros=FALSE,
                                   verbose=0)

  expect_equal(observed, expected)
})

test_that("getSegregatingLocusPerParent", {
  N <- 2
  P <- 10
  tX <- matrix(data=c(rep(0:2, each=3),NA, rep(0:2, times=3),1),
               nrow=P, ncol=N,
               dimnames=list(paste0("snp", 1:P), c("par1", "par2")))

  expected <- list(parent1=setNames(4:6, paste0("snp", 4:6)),
                   parent2=setNames(c(seq(2,10,3),10),
                                    paste0("snp", c(seq(2,10,3),10))))

  observed <- getSegregatingLocusPerParent(tX)

  expect_equal(observed, expected)
})

test_that("filterSegreg", {
  nb.offs <- 6 # offsprings
  N <- 2 + nb.offs
  P <- 5 # SNPs
  x <- data.frame(seg=c("<nnxnp>", NA, "<lmxll>", "<hkxhk>", "<nnxnp>"),
                  off1=c("nn", "GG", "lm", "hh", "np"),
                  off2=c("np", "GG", "lm", "hk", "np"),
                  off3=c("np", "GG", "ll", "kk", "np"),
                  off4=c("np", "GG", "ll", "kk", "np"),
                  off5=c("nn", "GG", "ll", "hk", NA),
                  off6=c("np", NA, NA, "hk", "np"),
                  row.names=paste0("snp", 1:P),
                  stringsAsFactors=FALSE)

  tmp <- data.frame(seg=x$seg,
                    nb.classes=c(2, NA, 2, 3, 2),
                    class1=c("nn", NA, "ll", "hh", "nn"),
                    class2=c("np", NA, "lm", "hk", "np"),
                    class3=c(NA, NA, NA, "kk", NA),
                    class4=c(NA, NA, NA, NA, NA),
                    obs1=c(2, NA, 3, 1, 0),
                    obs2=c(4, NA, 2, 3, 5),
                    obs3=c(NA, NA, NA, 2, NA),
                    obs4=c(NA, NA, NA, NA, NA),
                    exp1=c(0.5*6, NA, 0.5*5, 0.25*6, 0.5*5),
                    exp2=c(0.5*6, NA, 0.5*5, 0.5*6, 0.5*5),
                    exp3=c(NA, NA, NA, 0.25*6, NA),
                    exp4=c(NA, NA, NA, NA, NA),
                    row.names=rownames(x),
                    stringsAsFactors=FALSE)
  tmp$chi2 <- c((tmp$obs1[1] - tmp$exp1[1])^2 / tmp$exp1[1] +
                (tmp$obs2[1] - tmp$exp2[1])^2 / tmp$exp2[1],
                NA,
                (tmp$obs1[3] - tmp$exp1[3])^2 / tmp$exp1[3] +
                (tmp$obs2[3] - tmp$exp2[3])^2 / tmp$exp2[3],
                (tmp$obs1[4] - tmp$exp1[4])^2 / tmp$exp1[4] +
                (tmp$obs2[4] - tmp$exp2[4])^2 / tmp$exp2[4] +
                (tmp$obs3[4] - tmp$exp3[4])^2 / tmp$exp3[4],
                (tmp$obs1[5] - tmp$exp1[5])^2 / tmp$exp1[5] +
                (tmp$obs2[5] - tmp$exp2[5])^2 / tmp$exp2[5])
  tmp$pvalue <- pchisq(q=tmp$chi2, df=tmp$nb.classes - 1, lower.tail=FALSE)
  tmp$pvalue.bonf <- stats::p.adjust(p=tmp$pvalue, method="bonferroni")
  tmp$pvalue.bh <- stats::p.adjust(p=tmp$pvalue, method="BH")

  expected <- tmp[, c("chi2", "pvalue", "pvalue.bonf", "pvalue.bh")]
  observed <- filterSegreg(x=x, verbose=0)
  expect_equal(observed, expected)

  expected <- tmp
  observed <- filterSegreg(x=x, return.counts=TRUE, verbose=0)
  expect_equal(observed, expected)

  expected <- tmp
  observed <- filterSegreg(x=x, return.counts=TRUE, nb.cores=2, verbose=0)
  expect_equal(observed, expected)
})

test_that("filterSegreg_F2", {
  nb.offs <- 6 # offsprings
  N <- 2 + nb.offs
  P <- 4 # SNPs
  x <- data.frame(seg=c("F2", NA, "F2", "F2"),
                  off1=c("A", NA, "A", "A"),
                  off2=c("A", NA, "A", "A"),
                  off3=c("H", NA, "H", "A"),
                  off4=c("H", NA, "H", NA),
                  off5=c("B", NA, "A", NA),
                  off6=c("B", NA, NA, "A"),
                  row.names=paste0("snp", 1:P),
                  stringsAsFactors=FALSE)

  tmp <- data.frame(seg=x$seg,
                    nb.classes=c(3, NA, 3, 3),
                    class1=c("A", NA, "A", "A"),
                    class2=c("H", NA, "H", "H"),
                    class3=c("B", NA, "B", "B"),
                    class4=c(NA, NA, NA, NA),
                    obs1=c(2, NA, 3, 4),
                    obs2=c(2, NA, 2, 0),
                    obs3=c(2, NA, 0, 0),
                    obs4=c(NA, NA, NA, NA),
                    exp1=c(0.25*6, NA, 0.25*5, 0.25*4),
                    exp2=c(0.50*6, NA, 0.50*5, 0.50*4),
                    exp3=c(0.25*6, NA, 0.25*5, 0.25*4),
                    exp4=c(NA, NA, NA, NA),
                    row.names=rownames(x),
                    stringsAsFactors=FALSE)
  tmp$chi2 <- c((tmp$obs1[1] - tmp$exp1[1])^2 / tmp$exp1[1] +
                (tmp$obs2[1] - tmp$exp2[1])^2 / tmp$exp2[1] +
                (tmp$obs3[1] - tmp$exp3[1])^2 / tmp$exp3[1],
                NA,
                (tmp$obs1[3] - tmp$exp1[3])^2 / tmp$exp1[3] +
                (tmp$obs2[3] - tmp$exp2[3])^2 / tmp$exp2[3] +
                (tmp$obs3[3] - tmp$exp3[3])^2 / tmp$exp3[3],
                (tmp$obs1[4] - tmp$exp1[4])^2 / tmp$exp1[4] +
                (tmp$obs2[4] - tmp$exp2[4])^2 / tmp$exp2[4] +
                (tmp$obs3[4] - tmp$exp3[4])^2 / tmp$exp3[4])
  tmp$pvalue <- pchisq(q=tmp$chi2, df=tmp$nb.classes - 1, lower.tail=FALSE)
  tmp$pvalue.bonf <- stats::p.adjust(p=tmp$pvalue, method="bonferroni")
  tmp$pvalue.bh <- stats::p.adjust(p=tmp$pvalue, method="BH")

  expected <- tmp[, c("chi2", "pvalue", "pvalue.bonf", "pvalue.bh")]
  observed <- filterSegreg(x=x, is.F2=TRUE, verbose=0)
  expect_equal(observed, expected)

  expected <- tmp
  observed <- filterSegreg(x=x, is.F2=TRUE, return.counts=TRUE, verbose=0)
  expect_equal(observed, expected)

  expected <- tmp
  observed <- filterSegreg(x=x, is.F2=TRUE, return.counts=TRUE, nb.cores=2,
                           verbose=0)
  expect_equal(observed, expected)
})

test_that("genoDoses2genoClasses", {
  N <- 2 # individuals
  P <- 4 # SNPs
  X <- matrix(c(1,1, NA,0, 2,1, 1,NA), nrow=N, ncol=P,
              dimnames=list(paste0("ind", 1:N), paste0("snp", 1:P)))
  alleles <- data.frame(first=c("T","T","A","C","C"),
                        second=c("A","A","T","G","G"),
                        stringsAsFactors=FALSE)
  rownames(alleles) <- c(colnames(X), paste0("snp", P+1))

  expected <- data.frame(ind1=c("TA", "??", "TT", "CG"),
                         ind2=c("TA", "TT", "AT", "??"),
                         stringsAsFactors=FALSE)
  rownames(expected) <- colnames(X)
  expected <- as.matrix(expected)

  observed <- genoDoses2genoClasses(X=X, alleles=alleles, na.string="??",
                                    verbose=0)
  expect_equal(observed, expected)

  observed <- genoDoses2genoClasses(tX=t(X), alleles=alleles, na.string="??",
                                    nb.cores=2, verbose=0)
  expect_equal(observed, expected)
})

.expect_equal_VCFfile <- function(observed, expected){
  ## see https://support.bioconductor.org/p/74013/
  expect_equal(VariantAnnotation::info(observed),
               VariantAnnotation::info(expected))
  expect_equal(as.character(VariantAnnotation::ref(observed)),
               as.character(VariantAnnotation::ref(expected)))
  expect_equal(VariantAnnotation::alt(observed),
               VariantAnnotation::alt(expected))
  expect_equal(VariantAnnotation::samples(VariantAnnotation::header(observed)),
               VariantAnnotation::samples(VariantAnnotation::header(expected)))
  expect_equal(VariantAnnotation::geno(observed)$GT,
               VariantAnnotation::geno(expected)$GT)
}

test_that("genoDoses2vcf", {
  if(all(requireNamespace("Biostrings"),
         requireNamespace("VariantAnnotation"),
         requireNamespace("S4Vectors"))){
    nb.snps <- 4
    snp.ids <- paste0("snp", 1:nb.snps)
    snp.coords <- data.frame(chr=c(rep("chr1", 2), rep("chr2", 2)),
                             pos=c(3, 7, 635, 789),
                             row.names=snp.ids,
                             stringsAsFactors=FALSE)
    alleles <- simulRefAltSnpAlleles(snp.ids=snp.ids, verbose=0)
    nb.inds <- 3
    ind.ids <- paste0("ind", 1:nb.inds)
    X <- matrix(c(0,0,1, 1,2,0, NA,0,1, 1,2,NA),
                nrow=nb.inds, ncol=nb.snps,
                dimnames=list(ind.ids, snp.ids))

    gr <- seqIdStartEnd2GRanges(seq.id=snp.coords$chr,
                                seq.start=snp.coords$pos,
                                seq.end=snp.coords$pos,
                                subseq.name=snp.ids)
    Df <- S4Vectors::DataFrame(Samples=1:nb.inds,
                               row.names=ind.ids)
    expected <- VariantAnnotation::VCF(rowRanges=gr, colData=Df)
    VariantAnnotation::header(expected) <-
      VariantAnnotation::VCFHeader(samples=ind.ids)
    VariantAnnotation::geno(VariantAnnotation::header(expected)) <-
      S4Vectors::DataFrame(Number="1", Type="String",
                           Description="Genotype",
                           row.names="GT")
    VariantAnnotation::ref(expected) <- Biostrings::DNAStringSet(alleles$ref)
    VariantAnnotation::alt(expected) <- Biostrings::DNAStringSetList(
                                                        as.list(alleles$alt))
    VariantAnnotation::fixed(expected)[c("REF", "ALT")]
    VariantAnnotation::geno(expected) <-
      S4Vectors::SimpleList(
                     GT=matrix(c("0/0","0/1","./.","0/1",  # 1st ind
                                 "0/0","1/1","0/0","1/1",  # 2nd ind
                                 "0/1","0/0","0/1","./."), # ...
                               nrow=nb.snps, ncol=nb.inds,
                               dimnames=list(snp.ids, ind.ids)))

    observed <- genoDoses2Vcf(X=X, snp.coords=snp.coords, alleles=alleles,
                              verbose=0)

    .expect_equal_VCFfile(observed, expected)
  }
})

test_that("simulRefAltSnpAlleles", {
  nb.snps <- 4
  snp.ids <- paste0("snp", 1:nb.snps)

  expected <- data.frame(first=rep(NA, nb.snps),
                         second=NA,
                         stringsAsFactors=FALSE,
                         row.names=snp.ids)

  observed <- simulRefAltSnpAlleles(snp.ids=snp.ids,
                                    colnames=colnames(expected),
                                    verbose=0)

  expect_equal(is.data.frame(observed), is.data.frame(expected))
  expect_equal(dim(observed), dim(expected))
  expect_equal(dimnames(observed), dimnames(expected))
})

test_that("segSites2allDoses", {
  nb.inds <- 2
  nb.chrs <- 2
  nb.snps <- c(3, 2)
  ind.ids <- paste0("ind", 1:nb.inds)
  snp.ids <- paste0("snp", 1:sum(nb.snps))
  haplos <- list(chr1=matrix(c(0,0,0,0, 0,1,1,0, 1,1,0,1),
                             nrow=2 * nb.inds, ncol=nb.snps[1]),
                 chr2=matrix(c(1,0,0,1, 0,0,1,1), nrow=2 * nb.inds, ncol=nb.snps[2]))

  expected <- matrix(c(0,0, 1,1, 2,1, 1,1, 0,2),
                     nrow=nb.inds, ncol=sum(nb.snps),
                     dimnames=list(ind.ids, snp.ids))

  observed <- segSites2allDoses(seg.sites=haplos, ind.ids=ind.ids,
                                snp.ids=snp.ids)

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
                                                c("first", "second"))))

  expected <- matrix(data=c(1,1,1,1, 1,0,1,0, 0,1,0,0),
                     nrow=2*nb.inds, ncol=nb.snps,
                     dimnames=dimnames(haplos))

  observed <- haplosAlleles2num(haplos=haplos, alleles=alleles)

  expect_equal(observed, expected)
})

test_that("permuteAllelesInHaplosNum", {
  nb.inds <- 2
  nb.snps <- c(3, 2)
  ind.ids <- c("ind1_h1","ind1_h2","ind2_h1","ind2_h2")
  snp.ids <- list(paste0("snp", 1:nb.snps[1]),
                  paste0("snp", (nb.snps[1]+1):sum(nb.snps)))
  haplos <- list(chr1=matrix(data=c(0,0,0,0, 0,1,0,1, 1,0,1,1),
                             nrow=2*nb.inds, ncol=nb.snps[1],
                             dimnames=list(ind.ids, snp.ids[[1]])),
                 chr2=matrix(data=c(1,0,0,0, 0,1,0,1),
                             nrow=2*nb.inds, ncol=nb.snps[2],
                             dimnames=list(ind.ids, snp.ids[[2]])))

  snps.toperm <- c("snp2", "snp4")
  expected <- list(chr1=matrix(data=c(0,0,0,0, 1,0,1,0, 1,0,1,1),
                               nrow=2*nb.inds, ncol=nb.snps[1],
                               dimnames=dimnames(haplos[[1]])),
                   chr2=matrix(data=c(0,1,1,1, 0,1,0,1),
                               nrow=2*nb.inds, ncol=nb.snps[2],
                               dimnames=dimnames(haplos[[2]])))

  observed <- permuteAllelesInHaplosNum(haplos=haplos,
                                        snps.toperm=snps.toperm,
                                        verbose=0)

  expect_equal(observed, expected)
})

test_that("permuteAllelesInGenosDose", {
  nb.inds <- 2
  nb.snps <- 5
  ind.ids <- paste0("ind", 1:nb.inds)
  snp.ids <- paste0("snp", 1:nb.snps)
  X <- matrix(c(0,0, 1,1, 1,2, 1,0, 1,1), nrow=nb.inds, ncol=nb.snps,
              dimnames=list(ind.ids, snp.ids))

  snps.toperm <- c("snp3", "snp4")
  expected <- matrix(c(0,0, 1,1, 1,0, 1,2, 1,1), nrow=nb.inds, ncol=nb.snps,
                     dimnames=dimnames(X))

  observed <- permuteAllelesInGenosDose(X=X, snps.toperm=snps.toperm,
                                        verbose=0)

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

test_that("calcFreqMissSnpGenosPerSnp_vcf", {
  if(all(requireNamespace("Rsamtools"),
         requireNamespace("VariantAnnotation"))){

    tmpd <- tempdir()

    vcf.file <- system.file("extdata", "example.vcf",
                            package="rutilstimflutre")
    N <- 3
    snp.ids <- c("snp1", "snp2", "indel1")

    expected <- setNames(c(0/N, 2/N, 0/N), snp.ids)

    if(! file.exists(paste0(tmpd, "/", basename(vcf.file))))
      file.copy(from=vcf.file, to=tmpd)
    bgz.file <- Rsamtools::bgzip(file=paste0(tmpd, "/", basename(vcf.file)),
                                 overwrite=TRUE)
    Rsamtools::indexTabix(bgz.file, "vcf")
    observed <- calcFreqMissSnpGenosPerSnp(vcf.file=bgz.file)

    expect_equal(observed, expected)
  }
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

  ## as in VanRaden (2008)
  tmp <- matrix(rep(2*(afs-0.5), N), nrow=N, ncol=P, byrow=TRUE)
  Z <- X - 1 - tmp
  denom <- 0
  for(p in 1:P)
    denom <- denom + afs[p] * (1 - afs[p])
  denom <- 2 * denom
  expected <- (Z %*% t(Z)) / denom

  ## as in Toro et al (2011), equation 14
  if(FALSE){ # for debugging purposes
    expected <- matrix(data=NA, nrow=N, ncol=N,
                       dimnames=list(rownames(X), rownames(X)))
    for(i in 1:N){
      for(j in i:N){
        numerator <- 0
        denominator <- 0
        for(k in 1:P){
          numerator <- numerator + (X[i,k]/2 - afs[k]) * (X[j,k]/2 - afs[k])
          denominator <- denominator + afs[k] * (1 - afs[k])
        }
        expected[i,j] <- 2 * numerator / denominator
      }
    }
    expected[lower.tri(expected)] <- expected[upper.tri(expected)]
  }

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

test_that("estimGenRel_toro2011_eq10", {
  N <- 3 # individuals
  P <- 4 # SNPs
  X <- matrix(c(0,0,2, 1,1,0, 0,1,0, 1,0,0), nrow=N, ncol=P,
              dimnames=list(paste0("ind", 1:N), paste0("snp", 1:P)))

  afs <- setNames(c(0.383, 0.244, 0.567, 0.067), colnames(X))

  expected <- matrix(data=NA, nrow=N, ncol=N,
                     dimnames=list(rownames(X), rownames(X)))
  Cov_M <- matrix(data=0, nrow=N, ncol=N,
                  dimnames=list(rownames(X), rownames(X)))
  p.bar <- mean(afs)
  q.bar <- mean(1 - afs)
  var.p <- var(afs)
  for(i in 1:N){
    for(j in i:N){
      g.bar.i <- (1 / P) * sum(X[i,] / 2)
      g.bar.j <- (1 / P) * sum(X[j,] / 2)
      for(k in 1:P)
        Cov_M[i,j] <- Cov_M[i,j] + (X[i,k]/2 - g.bar.i) * (X[j,k]/2 - g.bar.j)
      Cov_M[i,j] <- Cov_M[i,j] / (P - 1) # use P-1 instead of P
      expected[i,j] <- (1 / (p.bar * q.bar - var.p)) * Cov_M[i,j] -
        var.p / (p.bar * q.bar - var.p) # coancestry
      expected[i,j] <- 2 * expected[i,j] # relationship
    }
  }
  expected[lower.tri(expected)] <- expected[upper.tri(expected)]

  observed <- estimGenRel(X=X, afs=afs, relationships="additive",
                          method="toro2011_eq10", verbose=0)

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
  P <- 5 # SNPs
  haplos.par.chr <- matrix(data=c(1,0, 0,1, 1,0, 1,0, 1,0), nrow=2, ncol=P,
                           dimnames=list(c("ind2_h1", "ind2_h2"),
                                         paste0("snp", 1:P)))

  loc.crossovers <- c(2, 4)
  expected <- matrix(c(1, 0, 0, 0, 1), nrow=1,
                     dimnames=list("ind2", paste0("snp", 1:P)))
  observed <- makeGameteSingleIndSingleChrom(haplos.par.chr, loc.crossovers,
                                             start.haplo=1)
  expect_equal(observed, expected)

  loc.crossovers <- c(1)
  expected <- matrix(c(1, 1, 0, 0, 0), nrow=1,
                     dimnames=list("ind2", paste0("snp", 1:P)))
  observed <- makeGameteSingleIndSingleChrom(haplos.par.chr, loc.crossovers,
                                             start.haplo=1)
  expect_equal(observed, expected)

  loc.crossovers <- c(2, 4)
  expected <- matrix(c(0, 1, 1, 1, 0), nrow=1,
                     dimnames=list("ind2", paste0("snp", 1:P)))
  observed <- makeGameteSingleIndSingleChrom(haplos.par.chr, loc.crossovers,
                                             start.haplo=2)
  expect_equal(observed, expected)

  loc.crossovers <- c(1)
  expected <- matrix(c(0, 0, 1, 1, 1), nrow=1,
                     dimnames=list("ind2", paste0("snp", 1:P)))
  observed <- makeGameteSingleIndSingleChrom(haplos.par.chr, loc.crossovers,
                                             start.haplo=2)
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
  observed <- makeGameteSingleInd(haplos.par=haplos.par,
                                  loc.crossovers=loc.crossovers)
  expect_equal(observed, expected)

  start.haplos <- as.list(c(1, 2))
  expected <- list(chr1=matrix(c(1,0, 0,0), nrow=1,
                               dimnames=list("ind2", paste0("snp", 1:(P/2)))),
                   chr2=matrix(c(0,1,0,1), nrow=1,
                               dimnames=list("ind2", paste0("snp", (P/2+1):P))))
  observed <- makeGameteSingleInd(haplos.par=haplos.par,
                                  loc.crossovers=loc.crossovers,
                                  start.haplos=start.haplos)
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

  observed <- makeCrosses(haplos, crosses, loc.crossovers, 1, verbose=0)

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

  observed <- makeCrosses(haplos, crosses, loc.crossovers, 1, verbose=0)

  expect_equal(observed, expected)
})

test_that("subsetDiffHaplosWithinParent", {
  P <- 10 # SNPs
  haplos.chr.par1 <- matrix(data=c(1,1, 0,1, 1,1, 0,0, 1,0,
                                   0,1, 0,1, 1,0, 0,0, 1,1),
                            nrow=2, ncol=P,
                            dimnames=list(c("ind2_h1", "ind2_h2"),
                                          paste0("snp", 1:P)))

  snps.co <- c("snp3", "snp7")
  expected <- haplos.chr.par1[, c(2,3,5,6,7,8)]

  observed <- subsetDiffHaplosWithinParent(haplos.chr=haplos.chr.par1,
                                           snps.tokeep=snps.co)

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

test_that("solveMme", {
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

  sigma.u.2 <- 20
  sigma.e.2 <- 40
  lambda <- sigma.e.2 / sigma.u.2
  Ainv <- matrix(c(1.833,0.5,0,-0.667,0,-1,0,0,
                   0.5,2,0.5,0,-1,-1,0,0,
                   0,0.5,2,0,-1,0.5,0,-1,
                   -0.667,0,0,1.833,0.5,0,-1,0,
                   0,-1,-1,0.5,2.5,0,-1,0,
                   -1,-1,0.5,0,0,2.5,0,-1,
                   0,0,0,-1,-1,0,2,0,
                   0,0,-1,0,0,-1,0,2),
                 ncol=8, byrow=TRUE)
  G <- sigma.u.2 * solve(Ainv)
  R <- sigma.e.2 * diag(length(y))
  V <- Z %*% G %*% t(Z) + R
  Vinv <- solve(V)

  b.hat <- solve(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv %*% y
  u.hat <- G %*% t(Z) %*% Vinv %*% (y - X %*% b.hat)

  expected <- c(b.hat, u.hat)

  observed <- solveMme(y=y, X=X, Z=Z, sigma.u2=sigma.u.2, Ainv=Ainv,
                       sigma2=sigma.e.2)

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

test_that("calcAsymptoticBayesFactorWakefield", {
  beta <- c(0.7, 0.7, -0.1)
  se <- c(0.2, 0.1, 0.4)
  grid = c(0.1, 0.2, 0.4, 0.8, 1.6)

  z2 <- (beta/se)^2
  v2 <- se^2
  expected <- log10(sapply(1:length(z2), function(i){
    mean(sapply(grid, function(phi){
      phi2 <- phi^2
      sqrt(v2[i] / (v2[i] + phi2)) * exp(0.5 * phi2 * z2[i] / (phi2 + v2[i]))
    }))
  }))

  observed <- calcAsymptoticBayesFactorWakefield(theta.hat=beta, V=se^2,
                                                 W=grid^2)

  expect_equal(observed, expected)
})

