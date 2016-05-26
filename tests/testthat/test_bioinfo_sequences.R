library(rutilstimflutre)
context("Bioinfo")

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

test_that("coverageBams", {
  if(all(file.exists(Sys.which("bwa")), file.exists(Sys.which("samtools")))){
    tmpd <- tempdir()
    set.seed(1)

    ## reference genome: 2 100-bp long, identical chromosomes
    faFile <- paste0(tmpd, "/refgenome.fa")
    chrs <- Biostrings::DNAStringSet(c(chr1=paste(sample(c("A","T","G","C"),
                                                         100, replace=TRUE),
                                                  collapse=""),
                                       chr2=paste(sample(c("A","T","G","C"),
                                                         100, replace=TRUE),
                                                  collapse="")))
    Biostrings::writeXStringSet(x=chrs, filepath=faFile, format="fasta")

    ## reads: 1 perfectly matching, 1 matching with internal indel
    fqFile <- paste0(tmpd, "/reads.fq")
    r1 <- chrs[[1]][11:80]
    r2 <- c(r1[1:30], r1[38:length(r1)])
    r3 <- chrs[[2]][11:80]
    reads <- Biostrings::DNAStringSet(c(read1=r1, read2=r2, read3=r3))
    Biostrings::writeXStringSet(x=reads, filepath=fqFile, format="fastq")

    ## alignment
    bamFile <- paste0(tmpd, "/align.bam")
    args <- paste0("index -p ", tmpd, "/", strsplit(basename(faFile),
                                                    "\\.")[[1]][1],
                   " ", faFile)
    system2("bwa", args, stdout=NULL, stderr=NULL)
    cmd <- paste0("bwa mem -R \'@RG\tID:ind1' -M ", tmpd, "/",
                  strsplit(basename(faFile), "\\.")[[1]][1],
                  " ", fqFile,
                  " | samtools fixmate -O bam - - ",
                  " | samtools sort -o ", bamFile, " -O bam -")
    system(cmd, ignore.stdout=TRUE, ignore.stderr=TRUE)
    args <- paste0("index ", bamFile, " ", bamFile, ".bai")
    system2("samtools", args)

    ## view alignments
    args <- paste0("view ", bamFile)
    system2("samtools", args)

    ## coverage via samtools
    cmd <- paste0("samtools depth -Q 5 ", bamFile)
    ## cmd <- paste0("samtools depth -Q 5 -r chr1:1-1000000 ", bamFile)
    cvg <- data.table::fread(input=cmd,
                             col.names=c("chr", "pos", "count"))
    nrow(cvg)
    sum(cvg$count)

    ## expected
    expected <- list()
    expected[[basename(bamFile)]] <-
      matrix(c(Biostrings::width(chrs), c(70, 70), c(133, 70)),
             nrow=length(chrs), ncol=3,
             dimnames=list(names(chrs),
                           c("nb.bases", "mapped.bases", "count.sum")))
    expected[[1]] <- cbind(expected[[1]],
                           breadth=expected[[1]][,"mapped.bases"] /
                             expected[[1]][,"nb.bases"],
                           depth=expected[[1]][,"count.sum"] /
                             expected[[1]][,"nb.bases"])

    ## observed
    observed <- coverageBams(bamFiles=bamFile)

    expect_equal(observed, expected)
  }
})

test_that("setGt2Na", {
  if(all(requireNamespace("Rsamtools"),
         requireNamespace("VariantAnnotation"),
         requireNamespace("IRanges"))){
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
    vcf.obs.file.bgz <- setGt2Na(vcf.file=vcf.init.file.bgz,
                                 genome=genome,
                                 out.file=vcf.obs.file,
                                 yieldSize=yieldSize,
                                 min.gq=90,
                                 verbose=0)
    observed <- VariantAnnotation::readVcf(file=vcf.obs.file.bgz,
                                           genome=genome)

    .expect_equal_VCFfile(observed, expected)
  }
})

test_that("filterVariantCalls", {
  if(all(requireNamespace("Rsamtools"),
         requireNamespace("VariantAnnotation"),
         requireNamespace("IRanges"))){
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
  }
})

test_that("summaryVariant", {
  if(all(requireNamespace("Rsamtools"),
         requireNamespace("VariantAnnotation"))){
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
  }
})

test_that("vcf2dosage", {
  if(all(requireNamespace("Rsamtools"),
         requireNamespace("VariantAnnotation"),
         requireNamespace("IRanges"),
         requireNamespace("S4Vectors"))){
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
  }
})
