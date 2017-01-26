library(rutilstimflutre)
context("Bioinfo")

test_that("chromNames2integers", {
  x <- c("chr1_random", "chr10", "chrUn", "chr1", "chr10_random", "chr2")

  expected <- data.frame(original=x,
                         renamed=c(11L, 10L, 21L, 1L, 20L, 2L),
                         stringsAsFactors=FALSE)

  observed <- chromNames2integers(x=x)

  expect_equal(observed, expected)
})

test_that("extractFasta", {
  if(all(requireNamespace("Biostrings"),
         requireNamespace("GenomicRanges"),
         requireNamespace("S4Vectors"),
         requireNamespace("IRanges"),
         requireNamespace("GenomeInfoDb"),
         requireNamespace("BSgenome"))){

    tmpd <- tempdir()

    records <- Biostrings::DNAStringSet(x=c("chr1 1872"="AAATTA",
                                            "chr2 5339"="CCCGGC"))
    in.fa <- paste0(tmpd, "/chromosomes.fa.gz")
    Biostrings::writeXStringSet(x=records, filepath=in.fa, compress=TRUE)

    sub.info <- data.frame(seq="chr2", start=4, end=5, name="loc")

    out.fa <- paste0(tmpd, "/subsequences.fa")
    if(file.exists(out.fa))
      file.remove(out.fa)

    expected <- Biostrings::DNAStringSet(x=c("loc"="GG"))

    extractFasta(in.fa=in.fa, sub.info=sub.info, out.fa=out.fa,
                 split.names=" ", verbose=0)

    observed <- Biostrings::readDNAStringSet(filepath=out.fa,
                                             format="fasta")
    expect_equal(observed, expected)

    if(file.exists(in.fa))
      file.remove(in.fa)
    if(file.exists(out.fa))
      file.remove(out.fa)
  }
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

test_that("coverageBams", {
  if(all(file.exists(Sys.which("bwa")),
         file.exists(Sys.which("samtools")))){
    tmpd <- tempdir()
    set.seed(1)
    all.files <- c()

    ## reference genome: 2 100-bp chromosomes
    faFile <- paste0(tmpd, "/refgenome.fa")
    chrs <- Biostrings::DNAStringSet(c(chr1=paste(sample(c("A","T","G","C"),
                                                         100, replace=TRUE),
                                                  collapse=""),
                                       chr2=paste(sample(c("A","T","G","C"),
                                                         100, replace=TRUE),
                                                  collapse="")))
    Biostrings::writeXStringSet(x=chrs, filepath=faFile, format="fasta")
    all.files <- append(all.files, faFile)

    ## reads: 1 perfectly matching, 1 matching with internal indel
    fqFile <- paste0(tmpd, "/reads.fq")
    r1 <- chrs[[1]][11:80] # 70-bp read on chr1
    r2 <- c(r1[1:30], r1[38:length(r1)]) # same as read 1 with 7-bp indel
    r3 <- chrs[[2]][11:80] # 70-bp read on chr2
    reads <- Biostrings::DNAStringSet(c(read1=r1, read2=r2, read3=r3))
    Biostrings::writeXStringSet(x=reads, filepath=fqFile, format="fastq")
    all.files <- append(all.files, fqFile)

    ## alignment
    bamFile <- paste0(tmpd, "/align.bam")
    args <- paste0("index -p ", tmpd, "/", strsplit(basename(faFile),
                                                    "\\.")[[1]][1],
                   " ", faFile)
    system2("bwa", args, stdout=NULL, stderr=NULL)
    all.files <- append(all.files, paste0(tmpd, "/refgenome.",
                                          c("amb", "ann", "bwt",
                                            "pac", "sa")))
    cmd <- paste0("bwa mem -R \'@RG\tID:ind1' -M ", tmpd, "/",
                  strsplit(basename(faFile), "\\.")[[1]][1],
                  " ", fqFile,
                  " | samtools fixmate -O bam - - ",
                  " | samtools sort -o ", bamFile, " -O bam -")
    system(cmd, ignore.stdout=TRUE, ignore.stderr=TRUE)
    all.files <- append(all.files, bamFile)
    args <- paste0("index ", bamFile, " ", bamFile, ".bai")
    system2("samtools", args)
    all.files <- append(all.files, paste0(bamFile, ".bai"))

    ## view alignments
    args <- paste0("view ", bamFile)
    ## system2("samtools", args)

    ## coverage via samtools
    depthFile <- paste0(tmpd, "/depth.txt")
    args <- paste0("depth -Q 5 ", bamFile, " > ", depthFile)
    ## args <- paste0("depth -Q 5 -r chr1:1-1000000 ", bamFile, " > ", depthFile)
    all.files <- append(all.files, depthFile)
    system2("samtools", args)
    cvg <- utils::read.table(file=depthFile, header=FALSE,
                             col.names=c("chr", "pos", "count"))
    nrow(cvg) # 140
    sum(cvg$count) # 203

    ## expected
    expected <- list()
    expected[[basename(bamFile)]] <-
      matrix(c(Biostrings::width(chrs),
               c(sum(cvg$chr == "chr1"),
                 sum(cvg$chr == "chr2")),
               c(sum(cvg$count[cvg$chr == "chr1"]),
                 sum(cvg$count[cvg$chr == "chr2"]))),
             nrow=length(chrs),
             ncol=3,
             dimnames=list(names(chrs),
                           c("nb.bases", "mapped.bases", "count.sum")))
    expected[[1]] <- cbind(expected[[1]],
                           breadth=expected[[1]][,"mapped.bases"] /
                             expected[[1]][,"nb.bases"],
                           depth=expected[[1]][,"count.sum"] /
                             expected[[1]][,"nb.bases"],
                           depth.map=expected[[1]][,"count.sum"] /
                             expected[[1]][,"mapped.bases"])

    ## observed
    observed <- coverageBams(bamFiles=bamFile)

    expect_equal(observed, expected)

    for(f in all.files)
      if(file.exists(f))
        file.remove(f)
  }
})

test_that("setGt2Na", {
  if(all(requireNamespace("Rsamtools"),
         requireNamespace("VariantAnnotation"),
         requireNamespace("IRanges"))){
    genome <- "fakeGenomeV0"
    yieldSize <- 100
    all.files <- c()

    vcf.init.file <- system.file("extdata", "example.vcf",
                                 package="rutilstimflutre")
    vcf.init.file.bgz <- Rsamtools::bgzip(file=vcf.init.file,
                                          overwrite=TRUE)
    all.files <- c(all.files, vcf.init.file.bgz)
    vcf.init.file.bgz.idx <- Rsamtools::indexTabix(file=vcf.init.file.bgz,
                                                   format="vcf")
    all.files <- c(all.files, vcf.init.file.bgz.idx)
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
    all.files <- c(all.files, vcf.obs.file.bgz,
                   paste0(vcf.obs.file.bgz, ".tbi"))
    observed <- VariantAnnotation::readVcf(file=vcf.obs.file.bgz,
                                           genome=genome)

    .expect_equal_VCFfile(observed, expected)

    for(f in all.files)
      if(file.exists(f))
        file.remove(f)
  }
})

test_that("filterVariantCalls", {
  if(all(requireNamespace("Rsamtools"),
         requireNamespace("VariantAnnotation"),
         requireNamespace("IRanges"))){
    genome <- "fakeGenomeV0"
    yieldSize <- 100
    all.files <- c()

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
                                           max.var.perc.gt.na=50,
                                           verbose=0)
    all.files <- c(all.files, vcf.obs.file.bgz,
                   paste0(vcf.obs.file.bgz, ".tbi"))
    observed <- VariantAnnotation::readVcf(file=vcf.obs.file.bgz,
                                           genome=genome)

    .expect_equal_VCFfile(observed, expected)

    for(f in all.files)
      if(file.exists(f))
        file.remove(f)
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
         requireNamespace("S4Vectors"),
         requireNamespace("BiocInstaller"))){
    genome <- "fakeGenomeV0"
    yieldSize <- 100
    all.files <- c()

    vcf.init.file <- system.file("extdata", "example.vcf",
                                 package="rutilstimflutre")
    vcf.init.file.bgz <- Rsamtools::bgzip(file=vcf.init.file,
                                          overwrite=TRUE)
    vcf.init.file.bgz.idx <- Rsamtools::indexTabix(file=vcf.init.file.bgz,
                                                   format="vcf")
    vcf.init <- VariantAnnotation::readVcf(file=vcf.init.file.bgz,
                                           genome=genome)
    if(utils::compareVersion(as.character(BiocInstaller::biocVersion()),
                             "3.4") < 0){
      vcf.init <- vcf.init[S4Vectors::elementLengths(VariantAnnotation::alt(vcf.init)) == 1L]
    } else
      vcf.init <- vcf.init[S4Vectors::elementNROWS(VariantAnnotation::alt(vcf.init)) == 1L]
    expected <- list(genotypes=gtVcf2dose(vcf.init),
                     map=rngVcf2df(vcf.init))

    pre.obs.files <- tempfile()
    gdose.file <- paste0(pre.obs.files, "_genos-dose.txt.gz")
    amap.file <- paste0(pre.obs.files, "_alleles-map.txt.gz")
    all.files <- c(all.files, c(gdose.file, amap.file))
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

    for(f in all.files)
      if(file.exists(f))
        file.remove(f)
  }
})

test_that("vcf2genoClasses", {
  if(all(requireNamespace("Rsamtools"),
         requireNamespace("VariantAnnotation"),
         requireNamespace("IRanges"),
         requireNamespace("S4Vectors"))){
    genome <- "fakeGenomeV0"
    yieldSize <- 100
    all.files <- c()

    vcf.init.file <- system.file("extdata", "example.vcf",
                                 package="rutilstimflutre")
    vcf.init.file.bgz <- Rsamtools::bgzip(file=vcf.init.file,
                                          overwrite=TRUE)
    vcf.init.file.bgz.idx <- Rsamtools::indexTabix(file=vcf.init.file.bgz,
                                                   format="vcf")
    expected <- matrix(c("AA", "AC", "CC",
                         "AC", "NN", "NN"),
                       byrow=TRUE, nrow=2, ncol=3,
                       dimnames=list(c("snp1", "snp2"),
                                     c("ind1", "ind2", "ind3")))

    prefix.obs <- tempfile()
    gclasses.file <- paste0(prefix.obs, "_genoClasses.txt.gz")
    all.files <- c(all.files, gclasses.file)
    obs.file <- vcf2genoClasses(vcf.file=vcf.init.file.bgz,
                                genome=genome,
                                yieldSize=yieldSize,
                                gclasses.file=gclasses.file,
                                na.string="NN",
                                verbose=0)
    observed <- as.matrix(read.table(file=gclasses.file, header=TRUE,
                                     sep="\t", row.names=1))

    expect_equal(observed, expected)

    for(f in all.files)
      if(file.exists(f))
        file.remove(f)
  }
})

test_that("invertGRanges", {
  if(all(requireNamespace("S4Vectors"),
         requireNamespace("BiocGenerics"),
         requireNamespace("IRanges"),
         requireNamespace("GenomicRanges"),
         requireNamespace("GenomeInfoDb"))){

    ## 3 alignments, with a mix of +/- strands on both references and queries
    in.gr <- GenomicRanges::GRanges(
        seqnames=S4Vectors::Rle(c("chr1_v1", "chr2_v1", "chr2_v1")),
        ranges=IRanges::IRanges(start=c(1, 11, 31),
                                end=c(10, 20, 40)),
        strand=c("+", "+", "-"))
    names(in.gr) <- c("chr1_v2", "chr2_v2", "chr2_v2")
    S4Vectors::mcols(in.gr) <- data.frame(
        qry.start=c(101, 120, 131),
        qry.end=c(110, 111, 140))

    expected <- GenomicRanges::GRanges(
        seqnames=S4Vectors::Rle(c("chr1_v2", "chr2_v2", "chr2_v2")),
        ranges=IRanges::IRanges(start=c(101, 111, 131),
                                end=c(110, 120, 140)),
        strand=c("+", "-", "+"))
    names(expected) <- c("chr1_v1", "chr2_v1", "chr2_v1")
    S4Vectors::mcols(expected) <- data.frame(
        qry.start=c(1, 11, 40),
        qry.end=c(10, 20, 31))

    observed <- invertGRanges(in.gr)

    expect_equal(observed, expected)
  }
})
