library(rutilstimflutre)
context("VCF")

.expected_equal_VCFfile <- function(observed, expected){
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

  .expected_equal_VCFfile(observed, expected)
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

  .expected_equal_VCFfile(observed, expected)
})
