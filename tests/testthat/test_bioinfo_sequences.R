library(rutilstimflutre)
context("Bioinfo")

test_that("chromNames2integers_grape", {
  x <- c("chr1_random", "chr10", "chrUn", "chr1", "chr10_random", "chr2")

  expected <- data.frame(original=x,
                         renamed=c(10+1, 10, 2*10+1, 1, 2*10, 2),
                         stringsAsFactors=FALSE)

  observed <- chromNames2integers(x=x)

  expect_equal(observed, expected)
})

test_that("chromNames2integers_apple", {
  x <- c("Chr15", "Chr01", "Chr02", "Chr00", "Chr02")

  expected <- data.frame(original=x,
                         renamed=c(15, 1, 2, ((2*15)+1), 2),
                         stringsAsFactors=FALSE)

  observed <- chromNames2integers(x=x)

  expect_equal(observed, expected)
})

test_that("chromNames2integers_cherry", {
  x <- c("Super-Scaffold_14", "Super-Scaffold_4374", "Super-Scaffold_27")

  expected <- data.frame(original=x,
                         renamed=c(1, 3, 2),
                         stringsAsFactors=FALSE)

  observed <- chromNames2integers(x=x, prefix="Super-Scaffold_",
                                  thresh.max.chr.int=500)

  expect_equal(observed, expected)
})

test_that("chromNames2integers_apricot", {
  x <- c("Pp08", "scaffold_51", "Pp02", "Pp01", "scaffold_23")

  expected <- data.frame(original=x,
                         renamed=c(8, 8+51, 2, 1, 8+23),
                         stringsAsFactors=FALSE)

  observed <- chromNames2integers(x=x, prefix="Pp", toreplace2=NULL,
                                  prefix2="scaffold_")

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

test_that("formatReadCountsPerLane_simple", {
  reads <- data.frame(ind=c("ind1","ind2","ind3"),
                      lane=c("fc1_lane1", "fc1_lane1",
                             "fc1_lane2"),
                      assigned=c(1, 1, 1))

  expected <- matrix(data=c(1, 1, 0,
                            0, 0, 1),
                     nrow=2, ncol=3, byrow=TRUE,
                     dimnames=list(c("fc1_lane1","fc1_lane2"),
                                   c("ind1","ind2","ind3")))

  observed <- formatReadCountsPerLane(x=reads)

  expect_equal(observed, expected)
})

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
                  " 2> /dev/null",
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
    observed <- coverageBams(bamFiles=bamFile, verbose=0)

    expect_equal(observed, expected)

    for(f in all.files)
      if(file.exists(f))
        file.remove(f)
  }
})

test_that("seqIdStartEnd2GRanges", {
  if(all(requireNamespace("IRanges"),
         requireNamespace("GenomicRanges"))){
    ir <- IRanges::IRanges(start=c(13, 875),
                           end=c(13, 1000))
    expected <- GenomicRanges::GRanges(seqnames=c("chr1", "chr2"),
                                       ranges=ir)
    names(expected) <- c("snp099", "gene13")

    observed <- seqIdStartEnd2GRanges(seq.id=c("chr1", "chr2"),
                                      seq.start=c(13, 875),
                                      seq.end=c(13, 1000),
                                      subseq.name=c("snp099", "gene13"))

    expect_equal(observed, expected)
  }
})

test_that("readBcftoolsCounts", {

  tmpd <- tempdir()

  input.files <- c(paste0(tmpd, "/stdout_bcftools_counts_1.txt"),
                   paste0(tmpd, "/stdout_bcftools_counts_2.txt"))
  txt <- "Number of samples: 100"
  txt <- paste0(txt, "\nNumber of SNPs:    1000")
  txt <- paste0(txt, "\nNumber of INDELs:  300")
  txt <- paste0(txt, "\nNumber of MNPs:    0")
  txt <- paste0(txt, "\nNumber of others:  13")
  txt <- paste0(txt, "\nNumber of sites:   1230")
  txt <- paste0(txt, "\n")
  cat(txt, file=input.files[1])
  txt <- "WARNING: bcftools version mismatch .. "
  txt <- paste0(txt, "\nNumber of samples: 50")
  txt <- paste0(txt, "\nNumber of SNPs:    500")
  txt <- paste0(txt, "\nNumber of INDELs:  170")
  txt <- paste0(txt, "\nNumber of MNPs:    0")
  txt <- paste0(txt, "\nNumber of others:  3")
  txt <- paste0(txt, "\nNumber of sites:   638")
  txt <- paste0(txt, "\n")
  cat(txt, file=input.files[2])

  expected <- data.frame(samples=c(100L, 50L),
                         SNPs=c(1000L, 500L),
                         INDELs=c(300L, 170L),
                         MNPs=c(0L, 0L),
                         others=c(13L, 3L),
                         sites=c(1230L, 638L),
                         stringsAsFactors=FALSE)
  rownames(expected) <- input.files
  expected <- as.matrix(expected)

  observed <- readBcftoolsCounts(files=input.files)

  expect_equal(observed, expected)

  for(f in input.files)
    if(file.exists(f))
      file.remove(f)
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

.makeDummyVcf <- function(){
  gr <- GenomicRanges::GRanges("chr1",
                               IRanges::IRanges(1:3*3, width=c(1, 2, 1)))
  nb.variants <- length(gr)
  names(gr) <- paste0("var", 1:nb.variants)
  nb.samples <- 5
  Df <- S4Vectors::DataFrame(Samples=1:nb.samples,
                             row.names=paste0("ind", 1:nb.samples))
  vcf <- VariantAnnotation::VCF(rowRanges=gr, colData=Df)
  VariantAnnotation::header(vcf) <-
    VariantAnnotation::VCFHeader(samples=paste0("ind", 1:nb.samples))
  VariantAnnotation::geno(VariantAnnotation::header(vcf)) <-
    S4Vectors::DataFrame(Number="1", Type="String",
                         Description="Genotype",
                         row.names="GT")
  VariantAnnotation::ref(vcf) <- Biostrings::DNAStringSet(c("G", c("AA"), "T"))
  VariantAnnotation::alt(vcf) <- Biostrings::DNAStringSetList("A", c("TT"), c("G", "A"))
  VariantAnnotation::fixed(vcf)[c("REF", "ALT")]
  VariantAnnotation::geno(vcf) <-
    S4Vectors::SimpleList(
                   GT=matrix(c("0/0","0/0","0/0",  # 1st sample
                               "0/1","0/1","./.",  # 2nd sample
                               "./.","0/0","1/1",  # ...
                               "1/1","1/1","1/2",
                               "1/1","./.","0/2"),
                             nrow=nb.variants, ncol=nb.samples,
                             dimnames=list(names(gr),
                                           rownames(Df))))
  return(vcf)
}

test_that("dimVcf", {
  if(all(requireNamespace("Rsamtools"),
         requireNamespace("VariantAnnotation"))){
    tmpd <- tempdir()

    vcf.file <- system.file("extdata", "example.vcf",
                            package="rutilstimflutre")
    bgz.file <- Rsamtools::bgzip(vcf.file,
                                 paste0(tmpd, "/", basename(vcf.file), ".gz"),
                                 overwrite=TRUE)

    expected <- stats::setNames(object=c(3, 3),
                                nm=c("sites", "samples"))

    observed <- dimVcf(bgz.file)

    expect_equal(observed, expected)

    file.remove(bgz.file)
  }
})

test_that("tableVcfAlt", {
  if(all(requireNamespace("Biostrings"),
         requireNamespace("VariantAnnotation"),
         requireNamespace("S4Vectors"))){
    vcf.file <- system.file("extdata", "example.vcf",
                            package="rutilstimflutre")

    expected <- Biostrings::DNAStringSetList(c("C"), c("C"), c("T"))
    names(expected) <- c("snp1", "snp2", "indel1")

    observed <- tableVcfAlt(vcf.file=vcf.file, verbose=0)

    expect_equal(observed, expected)
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
    vcf.init.file.bgz.idx <- paste0(vcf.init.file.bgz, ".tbi")
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
    all.files <- c(all.files, vcf.init.file.bgz)
    vcf.init.file.bgz.idx <- paste0(vcf.init.file.bgz, ".tbi")
    all.files <- c(all.files, vcf.init.file.bgz.idx)
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
    all.files <- c()

    vcf.init.file <- system.file("extdata", "example.vcf",
                                 package="rutilstimflutre")
    vcf.init.file.bgz <- Rsamtools::bgzip(file=vcf.init.file,
                                          overwrite=TRUE)
    all.files <- c(all.files, vcf.init.file.bgz)
    vcf.init.file.bgz.idx <- paste0(vcf.init.file.bgz, ".tbi")
    all.files <- c(all.files, vcf.init.file.bgz.idx)
    vcf.init <- VariantAnnotation::readVcf(file=vcf.init.file.bgz,
                                           genome=genome)
    expected <-
      list(GQ=matrix(data=NA,
                     nrow=nrow(vcf.init),
                     ncol=9,
                     dimnames=list(rownames(vcf.init),
                                   c("n","na","mean","sd","min","q1","med",
                                     "q3","max"))))
    gq <- VariantAnnotation::geno(vcf.init)[["GQ"]]
    for(i in 1:nrow(vcf.init))
      expected$GQ[i,] <- c(ncol(gq),
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
                               fields="GQ",
                               verbose=0)

    expect_equal(observed, expected)

    for(f in all.files)
      if(file.exists(f))
        file.remove(f)
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
    all.files <- c(all.files, vcf.init.file.bgz)
    vcf.init.file.bgz.idx <- paste0(vcf.init.file.bgz, ".tbi")
    all.files <- c(all.files, vcf.init.file.bgz.idx)
    vcf.init <- VariantAnnotation::readVcf(file=vcf.init.file.bgz,
                                           genome=genome)
    if(utils::compareVersion(as.character(BiocInstaller::biocVersion()),
                             "3.4") < 0){
      vcf.init <- vcf.init[S4Vectors::elementLengths(VariantAnnotation::alt(vcf.init)) == 1L]
    } else
      vcf.init <- vcf.init[S4Vectors::elementNROWS(VariantAnnotation::alt(vcf.init)) == 1L]
    expected <- list(genos=matrix(c(0, 1, 2,
                                    1, NA, NA),
                                  byrow=TRUE, nrow=2, ncol=3,
                                  dimnames=list(c("snp1", "snp2"),
                                                c("ind1", "ind2", "ind3"))),
                     ca=data.frame(chr=c("chr1", "chr1"),
                                   pos=c(3L, 7L),
                                   allele.ref=c("A", "A"),
                                   allele.alt=c("C", "C"),
                                   row.names=c("snp1", "snp2"),
                                   stringsAsFactors=FALSE))

    pre.obs.files <- tempfile()
    gdose.file <- paste0(pre.obs.files, "_genos-dose.txt.gz")
    ca.file <- paste0(pre.obs.files, "_coords-alleles.txt.gz")
    all.files <- c(all.files, c(gdose.file, ca.file))
    vcf2dosage(vcf.file=vcf.init.file.bgz,
               genome=genome,
               yieldSize=yieldSize,
               gdose.file=gdose.file,
               ca.file=ca.file,
               verbose=0)
    observed <- list(genos=as.matrix(read.table(file=gdose.file,
                                                header=TRUE, sep="\t",
                                                row.names=1)),
                     ca=read.table(file=ca.file,
                                   header=TRUE, sep="\t",
                                   stringsAsFactors=FALSE))

    expect_equal(observed, expected)

    for(f in all.files)
      if(file.exists(f))
        file.remove(f)
  }
})

test_that("gtVcf2genoClasses", {
  if(all(requireNamespace("VariantAnnotation"),
         requireNamespace("GenomicRanges"),
         requireNamespace("IRanges"),
         requireNamespace("S4Vectors"),
         requireNamespace("Biostrings"))){
    vcf.init <- .makeDummyVcf()
    sample.names <- VariantAnnotation::samples(VariantAnnotation::header(vcf.init))
    vcf.init <- vcf.init[VariantAnnotation::isSNV(vcf.init, singleAltOnly=FALSE),]

    ## case: SNV with 2 alleles
    expected <- matrix(c("GG",  # 1st sample
                         "GA",  # 2nd sample
                         NA,    # ...
                         "AA",
                         "AA"),
                       nrow=1, ncol=ncol(vcf.init),
                       dimnames=list(names(vcf.init)[1], sample.names))
    observed <- gtVcf2genoClasses(vcf=vcf.init, na.string=NA, single.alt=TRUE)
    expect_equal(observed, expected)

    ## case: SNV with 2 alleles or more
    expected <- matrix(c("GG","TT",  # 1st sample
                         "GA",NA,    # 2nd sample
                         NA,"GG",    # ...
                         "AA","GA",
                         "AA","TA"),
                       nrow=2, ncol=ncol(vcf.init),
                       dimnames=list(names(vcf.init),
                                     sample.names))
    observed <- gtVcf2genoClasses(vcf=vcf.init, na.string=NA, single.alt=FALSE)
    expect_equal(observed, expected)
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
    all.files <- c(all.files, vcf.init.file.bgz)
    vcf.init.file.bgz.idx <- paste0(vcf.init.file.bgz, ".tbi")
    all.files <- c(all.files, vcf.init.file.bgz.idx)
    expected <- list(genos=matrix(c("AA", "AC", "CC",
                                    "AC", "NN", "NN"),
                                  byrow=TRUE, nrow=2, ncol=3,
                                  dimnames=list(c("snp1", "snp2"),
                                                c("ind1", "ind2", "ind3"))),
                     ca=data.frame(chr=c("chr1", "chr1"),
                                   pos=c(3L, 7L),
                                   allele.ref=c("A", "A"),
                                   allele.alt=c("C", "C"),
                                   row.names=c("snp1", "snp2"),
                                   stringsAsFactors=FALSE))

    prefix.obs <- tempfile()
    gclasses.file <- paste0(prefix.obs, "_geno-classes.txt.gz")
    ca.file <- paste0(prefix.obs, "_coords-alleles.txt.gz")
    all.files <- c(all.files, gclasses.file, ca.file)
    vcf2genoClasses(vcf.file=vcf.init.file.bgz,
                    genome=genome,
                    yieldSize=yieldSize,
                    gclasses.file=gclasses.file,
                    ca.file=ca.file,
                    na.string="NN",
                    single.alt=TRUE,
                    verbose=0)
    observed <- list(genos=as.matrix(read.table(file=gclasses.file,
                                                header=TRUE, sep="\t",
                                                row.names=1)),
                     ca=read.table(file=ca.file,
                                   header=TRUE, sep="\t",
                                   stringsAsFactors=FALSE))

    expect_equal(observed, expected)

    for(f in all.files)
      if(file.exists(f))
        file.remove(f)
  }
})

test_that("calcFreqNaVcf", {
  if(all(requireNamespace("S4Vectors"),
         requireNamespace("VariantAnnotation"),
         requireNamespace("SummarizedExperiment"))){
    vcf.file <- system.file("extdata", "example.vcf",
                            package="rutilstimflutre")
    vcf <- VariantAnnotation::readVcf(vcf.file, genome="")
    N <- 3
    snp.ids <- c("snp1", "snp2", "indel1")

    expected <- setNames(c(0/N, 2/N, 0/N), snp.ids)

    observed <- calcFreqNaVcf(vcf=vcf, with.coords=FALSE)

    expect_equal(observed, expected)
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
