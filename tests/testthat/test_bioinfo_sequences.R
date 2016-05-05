library(rutilstimflutre)
context("VCF")

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

test_that("setGt2Na", {
  genome <- "fakeGenomeV0"
  yieldSize <- 100

  vcf.init.file <- system.file("extdata", "example.vcf",
                               package="rutilstimflutre")
  vcf.init.file.bgz <- Rsamtools::bgzip(file=vcf.init.file,
                                        overwrite=TRUE)
  vcf.init.file.bgz.idx <- Rsamtools::indexTabix(file=vcf.init.file.bgz,
                                                 format="vcf")
  vcf.init <- VariantAnnotation::readVcf(file=vcf.init.file.bgz,
                                         genome=genome)
  expected <- vcf.init
  VariantAnnotation::geno(expected)[["GT"]][3, c(2,3)] <- "."

  vcf.obs.file <- sprintf("%s.vcf", tempfile())
  library(IRanges) # see https://support.bioconductor.org/p/73957/#73962
  vcf.obs.file.bgz <- setGt2Na(vcf.file=vcf.init.file.bgz,
                               genome=genome,
                               out.file=vcf.obs.file,
                               yieldSize=yieldSize,
                               min.gq=90,
                               verbose=0)
  observed <- VariantAnnotation::readVcf(file=vcf.obs.file.bgz,
                                         genome=genome)

  .expect_equal_VCFfile(observed, expected)
})

test_that("filterVariantCalls", {
  genome <- "fakeGenomeV0"
  yieldSize <- 100

  vcf.init.file <- system.file("extdata", "example.vcf",
                               package="rutilstimflutre")
  vcf.init.file.bgz <- Rsamtools::bgzip(file=vcf.init.file,
                                        overwrite=TRUE)
  vcf.init.file.bgz.idx <- Rsamtools::indexTabix(file=vcf.init.file.bgz,
                                                 format="vcf")
  vcf.init <- VariantAnnotation::readVcf(file=vcf.init.file.bgz,
                                         genome=genome)
  expected <- vcf.init[-c(2,3)] # remove indel 1 and snp2

  vcf.obs.file <- sprintf("%s.vcf", tempfile())
  library(IRanges) # see https://support.bioconductor.org/p/73957/#73962
  vcf.obs.file.bgz <- filterVariantCalls(vcf.file=vcf.init.file.bgz,
                                         genome=genome,
                                         out.file=vcf.obs.file,
                                         yieldSize=yieldSize,
                                         is.snv=TRUE,
                                         max.var.prop.gt.na=0.5,
                                         verbose=0)
  observed <- VariantAnnotation::readVcf(file=vcf.obs.file.bgz,
                                         genome=genome)

  .expect_equal_VCFfile(observed, expected)
})

test_that("summaryVariant", {
  genome <- "fakeGenomeV0"
  yieldSize <- 100

  vcf.init.file <- system.file("extdata", "example.vcf",
                               package="rutilstimflutre")
  vcf.init.file.bgz <- Rsamtools::bgzip(file=vcf.init.file,
                                        overwrite=TRUE)
  vcf.init.file.bgz.idx <- Rsamtools::indexTabix(file=vcf.init.file.bgz,
                                                 format="vcf")
  vcf.init <- VariantAnnotation::readVcf(file=vcf.init.file.bgz,
                                         genome=genome)
  expected <- matrix(data=NA,
                     nrow=nrow(vcf.init),
                     ncol=9,
                     dimnames=list(rownames(vcf.init),
                                   c("n","na","mean","sd","min","q1","med",
                                     "q3","max")))
  gq <- VariantAnnotation::geno(vcf.init)[["GQ"]]
  for(i in 1:nrow(vcf.init))
    expected[i,] <- c(ncol(gq),
                      sum(is.na(gq[i,])),
                      mean(gq[i,]),
                      sd(gq[i,]),
                      min(gq[i,]),
                      quantile(gq[i,], 0.25),
                      median(gq[i,]),
                      quantile(gq[i,], 0.75),
                      max(gq[i,]))

  observed <- summaryVariant(vcf.file=vcf.init.file.bgz,
                             genome=genome,
                             yieldSize=yieldSize,
                             field="GQ",
                             verbose=0)

  expect_equal(observed, expected)
})

test_that("vcf2dosage", {
  genome <- "fakeGenomeV0"
  yieldSize <- 100

  vcf.init.file <- system.file("extdata", "example.vcf",
                               package="rutilstimflutre")
  vcf.init.file.bgz <- Rsamtools::bgzip(file=vcf.init.file,
                                        overwrite=TRUE)
  vcf.init.file.bgz.idx <- Rsamtools::indexTabix(file=vcf.init.file.bgz,
                                                 format="vcf")
  vcf.init <- VariantAnnotation::readVcf(file=vcf.init.file.bgz,
                                         genome=genome)
  vcf.init <- vcf.init[S4Vectors::elementLengths(VariantAnnotation::alt(vcf.init)) == 1L]
  expected <- list(genotypes=gtVcf2dose(vcf.init),
                   map=rngVcf2df(vcf.init))

  pre.obs.files <- tempfile()
  gdose.file <- paste0(pre.obs.files, "_genos-dose.txt.gz")
  amap.file <- paste0(pre.obs.files, "_alleles-map.txt.gz")
  library(IRanges) # see https://support.bioconductor.org/p/73957/#73962
  obs.files <- vcf2dosage(vcf.file=vcf.init.file.bgz,
                          genome=genome,
                          yieldSize=yieldSize,
                          gdose.file=gdose.file,
                          amap.file=amap.file,
                          verbose=0)
  gtmp <- read.table(file=gdose.file, row.names=1)
  atmp <- read.table(file=amap.file, header=TRUE, stringsAsFactors=FALSE)
  observed <- list(genotypes=as.matrix(gtmp),
                   map=atmp)

  expect_equal(observed, expected)
})
