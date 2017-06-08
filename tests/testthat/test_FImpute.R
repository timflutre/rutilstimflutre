library(rutilstimflutre)
context("FImpute")

test_that("readGenosFimpute_haplos2genos", {
  if(file.exists(Sys.which("FImpute"))){

    tmpd <- tempdir()

    N <- 3 # nb of genotypes
    P <- 4 # nb of SNPs
    input <- data.frame(ID=paste0("ind", 1:N),
                        Chip=rep(1, N),
                        "Call..."=c("0234", "1567", "0289"))
    input.file <- paste0(tmpd, "/fimpute_genotypes_imp.txt")
    utils::write.table(x=input, file=input.file, quote=FALSE, sep="\t",
                       row.names=FALSE, col.names=TRUE)

    expected <- matrix(data=c(0,2,1,1, 1,NA,NA,NA, 0,2,NA,NA), byrow=TRUE,
                       nrow=N, ncol=P,
                       dimnames=list(paste0("ind", 1:N), paste0("snp", 1:P)))

    observed <- readGenosFimpute(file=input.file, snp.ids=paste0("snp", 1:P),
                                 input.haplos=TRUE, output.genos=TRUE)

    expect_equal(observed, expected)

    if(file.exists(input.file))
      file.remove(input.file)
  }
})
