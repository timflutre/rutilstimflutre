## Contains functions useful to wrap fastPHASE.

##' Write inputs for fastPHASE
##'
##' Use the external software fcGENE from \href{dx.plos.org/10.1371/journal.pone.0097589}{Roshyara and Scholz (2014)}.
##' @param X matrix of bi-allelic SNP genotypes encoded, for each SNP, in number of copies of its second allele, i.e. as allele doses in {0,1,2}, with genotypes in rows and SNPs in columns; the "second" allele is arbitrary, it corresponds to the second column of \code{alleles}, which can be the minor or the major allele
##' @param snp.coords data.frame with SNP identifiers as row names, and two columns, "chr" and "coord" or "pos"
##' @param alleles data.frame with SNPs in rows (names as row names) and alleles in columns (exactly 2 columns are required); the second column should correspond to the allele which number of copies is counted at each SNP in \code{X}
##' @param file write the genotype data to this file (for instance 'genotypes_fastphase.txt' or 'genotypes_fastphase.txt.gz', but don't use \code{\link[base]{gzfile}})
##' @param verbose verbosity level (0/1/2)
##' @return nothing
##' @author Timothee Flutre
##' @export
writeInputsFastphase <- function(X, snp.coords, alleles, file, verbose=0){
  stopifnot(file.exists(Sys.which("fcgene")))
  stopIfNotValidGenosDose(X, check.noNA=FALSE)
  stopifnot(.isValidSnpCoords(snp.coords),
            is.data.frame(alleles),
            ncol(alleles) == 2)

  if(file.exists(file))
    file.remove(file)

  alleles <- convertFactorColumnsToCharacter(alleles)

  ## put in same order
  snp.coords <- snp.coords[colnames(X),]
  alleles <- alleles[colnames(X),]

  if(verbose > 0)
    write("save SNP genotypes as 'r-formatted' file for fcGENE...", stdout())
  tmp.geno.file <- tempfile(pattern="genotypes_", fileext=".txt")
  sep <- "\t"
  txt <- paste0("SAMPLE_ID", sep, paste0(colnames(X), collapse=sep), "\n")
  cat(txt, file=tmp.geno.file, sep="")
  tmp <- do.call(rbind, lapply(1:nrow(X), function(i){
    paste0(rownames(X)[i], sep, paste0(X[i,], collapse=sep), "\n")
  }))
  cat(tmp, file=tmp.geno.file, sep="", append=TRUE)

  if(verbose > 0)
    write("save 'snpinfo'...", stdout())
  tmp.snpinfo.file <- tempfile(pattern="snpinfo_", fileext=".txt")
  snp.info <- data.frame(snpid=rownames(snp.coords),
                         rsid=rownames(snp.coords),
                         position=snp.coords$pos,
                         allele1=alleles[,2],
                         allele2=alleles[,1],
                         stringsAsFactors=FALSE)
  utils::write.table(x=snp.info, file=tmp.snpinfo.file, quote=FALSE, sep="\t",
                     row.names=FALSE)

  if(verbose > 0)
    write("convert with fcGENE...", stdout())
  tmp.out.prefix <- tempfile(pattern="output")
  args <- paste0("--rgeno ", tmp.geno.file,
                 " --snpinfo ", tmp.snpinfo.file,
                 " --oformat fastphase",
                 " --out ", tmp.out.prefix)
  retVal <- system2(command="fcgene", args=args, stdout=TRUE, stderr=TRUE)
  if(! is.null(attributes(retVal)) && attr(retVal, "status") != 0){
    print(retVal)
    stop()
  }
  if(verbose > 1)
    print(retVal)

  if(verbose > 0)
    write("rename (and compress if necessary)...", stdout())
  use.gzip <- grepl(".gz", file)
  if(use.gzip)
    file <- removeFileExtension(file, "gz")
  file.copy(from=paste0(tmp.out.prefix, ".inp"), to=file)
  if(use.gzip){
    if(file.exists(paste0(file, ".gz")))
      file.remove(paste0(file, ".gz"))
    system2(command="gzip", args=file)
  }

  if(verbose > 0)
    write("clean...", stdout())
  tmp <- file.remove(c(tmp.geno.file, tmp.snpinfo.file,
                       paste0(tmp.out.prefix, ".inp"),
                       paste0(tmp.out.prefix, "_fcgene.log"),
                       paste0(tmp.out.prefix, "_pedinfo.txt"),
                       paste0(tmp.out.prefix, "_snpinfo.txt")))
}

##' Read genotypes from fastPHASE
##'
##' Extract the genotypes from the output file of fastPHASE.
##' @param file character
##' @param snp.ids vector of SNP identifiers
##' @return matrix with genotypes in rows (2 per genotype) and SNPs in columns
##' @author Timothee Flutre
##' @export
readGenosFastphase <- function(file, snp.ids){
  stopifnot(file.exists(file),
            is.vector(snp.ids))

  lines <- readLines(file, warn=FALSE)
  idx <- which(grepl("SSS", lines))
  lines <- lines[(idx+1):length(lines)]
  stopifnot(length(lines) %% 3 == 0)

  genos <- matrix(data=NA, nrow=2 * (length(lines) / 3), ncol=nchar(lines[2]))
  rownames(genos) <- paste0(rep(lines[seq(1,length(lines),3)], each=2),
                             c("_g1", "_g2"))
  colnames(genos) <- snp.ids
  for(i in seq(1, length(lines), 3)){
    ind.id <- paste0(lines[i], "_g1")
    ind.idx <- which(rownames(genos) == ind.id)
    genos[ind.idx:(ind.idx+1),] <- do.call(rbind, strsplit(lines[(i+1):(i+2)], ""))
  }

  return(genos)
}

##' Read haplotypes from fastPHASE
##'
##' Extract the haplotypes from fastPHASE's output files ("hapguess_switch" and "hapguess_indiv").
##' @param file character
##' @param snp.ids vector of SNP identifiers
##' @return matrix with haplotypes in rows (2 per genotype) and SNPs in columns
##' @author Timothee Flutre
##' @export
readHaplosFastphase <- function(file, snp.ids){
  stopifnot(file.exists(file),
            is.vector(snp.ids))

  lines <- readLines(file, warn=FALSE)
  idx <- (which(grepl("BEGIN GENOTYPES", lines))+1):(which(grepl("END GENOTYPES", lines))-1)
  lines <- lines[idx]
  lines <- gsub(" ", "", lines)
  lines <- lines[! lines == ""]
  stopifnot(length(lines) %% 3 == 0)

  haplos <- matrix(data=NA, nrow=2 * (length(lines) / 3), ncol=nchar(lines[2]))
  rownames(haplos) <- paste0(rep(lines[seq(1,length(lines),3)], each=2),
                             c("_h1", "_h2"))
  colnames(haplos) <- snp.ids
  for(i in seq(1, length(lines), 3)){
    ind.id <- paste0(lines[i], "_h1")
    ind.idx <- which(rownames(haplos) == ind.id)
    haplos[ind.idx:(ind.idx+1),] <- do.call(rbind, strsplit(lines[(i+1):(i+2)], ""))
  }

  return(haplos)
}

##' Genotype imputation via fastPHASE
##'
##' Impute SNP genotypes via fastPHASE (Scheet and Stephens, 2006).
##' @param X matrix of bi-allelic SNP genotypes encoded in number of copies of the 2nd allele, i.e. as allele doses in {0,1,2}, with genotypes in rows and SNPs in columns
##' @param snp.coords data.frame with SNP identifiers as row names, and two columns, "chr" and "coord" or "pos"
##' @param alleles data.frame with SNPs in rows (names as row names) and alleles in columns, named "minor" and "major", in whatever order as long as the second column corresponds to the allele which number of copies is counted at each SNP in \code{X}
##' @param out.dir directory in which the output files will be saved
##' @param task.id identifier of the task (used in temporary and output file names)
##' @param nb.starts number of random starts of the EM algorithm
##' @param nb.iters number of iterations of the EM algorithm
##' @param nb.samp.haplos number of haplotypes sampled from the posterior distribution obtained from a particular random start of the EM algorithm
##' @param estim.haplos estimate haplotypes by minimizing genotype error
##' @param nb.clusters number of clusters
##' @param seed seed for the pseudo-random number generator for fastPHASE
##' @param clean remove files
##' @param verbose verbosity level (0/1)
##' @return list with matrix of genotypes before imputation, matrix of haplotypes (2 per genotype, in rows) after imputation, stdout/stderr
##' @author Timothee Flutre
##' @seealso \code{\link{writeInputsFastphase}}, \code{\link{readGenosFastphase}}, \code{\link{readHaplosFastphase}}, \code{\link{runFimpute}}
##' @examples
##' \dontrun{## simulate haplotypes and genotypes in a single population
##' set.seed(1859)
##' nb.inds <- 5*10^2
##' nb.chrs <- 1
##' Ne <- 10^4
##' chrom.len <- 1*10^6
##' mu <- 10^(-8)
##' c.rec <- 10^(-7)
##' genomes <- simulCoalescent(nb.inds=nb.inds,
##'                            nb.reps=nb.chrs,
##'                            pop.mut.rate=4 * Ne * mu * chrom.len,
##'                            pop.recomb.rate=4 * Ne * c.rec * chrom.len,
##'                            chrom.len=chrom.len,
##'                            get.alleles=TRUE)
##' nb.snps <- nrow(genomes$snp.coords)
##' plotHaplosMatrix(genomes$haplos[[1]]) # quick view of the amount of LD
##'
##' ## discard some genotypes according to a "microarray" design:
##' ## some inds with high density of genotyped SNPs, and the others with
##' ## low density of SNPs, these being on both microarrays
##' ind.names <- rownames(genomes$genos)
##' inds.high <- sample(x=ind.names, size=floor(0.4 * nb.inds))
##' inds.low <- setdiff(ind.names, inds.high)
##' snp.names <- colnames(genomes$genos)
##' mafs <- estimSnpMaf(X=genomes$genos)
##' length(idx <- which(mafs >= 0.05))
##' snps.high <- names(idx)
##' snps.high.only <- sample(x=snps.high, size=floor(0.7*length(snps.high)))
##' dim(X <- genomes$genos[ind.names, snps.high])
##' X.na <- X
##' X.na[inds.low, snps.high.only] <- NA
##' sum(is.na(X.na)) / length(X.na)
##' plotGridMissGenos(X=X.na)
##'
##' ## perform imputation
##' snp.coords <- genomes$snp.coords[colnames(X.na),]
##' alleles <- genomes$alleles[colnames(X.na),]
##' out.fP <- runFastphase(X=X.na, snp.coords=snp.coords, alleles=alleles,
##'                        nb.starts=3, clean=TRUE)
##' out.fP$stdouterr
##' tmp <- haplosAlleles2num(haplos=out.fP$haplos.switch, alleles=alleles)
##' X.imp.fP <- segSites2allDoses(seg.sites=list(tmp), ind.ids=ind.names)
##'
##' ## assess imputation accuracy
##' genomes$haplos[[1]][1:4, 1:7]
##' genomes$genos[1:2, 1:7]
##' X.na[1:2, 1:6]
##' head(alleles)
##' out.fP$genos[1:4, 1:6]
##' out.fP$haplos.switch[1:4, 1:6]
##' X.imp.fP[1:2, 1:6]
##' 100 * sum(X.imp.fP != X) / sum(is.na(X.na))
##' }
##' @export
runFastphase <- function(X, snp.coords, alleles, out.dir=getwd(),
                         task.id="fastphase", nb.starts=20, nb.iters=25,
                         nb.samp.haplos=50, estim.haplos=TRUE,
                         nb.clusters=10, seed=1859, clean=FALSE,
                         verbose=1){
  exe.name <- "fastPHASE"
  stopifnot(file.exists(Sys.which(exe.name)))
  stopIfNotValidGenosDose(X, check.noNA=FALSE)

  out <- list()

  out.dir <- path.expand(out.dir)
  task.id <- paste0(out.dir, "/", task.id)

  if(verbose > 0){
    msg <- "prepare input files..."
    write(msg, stdout())
  }
  tmp.files <- c()
  tmp.files["genos"] <- paste0(task.id, "_input-genos.txt")
  writeInputsFastphase(X=X, snp.coords=snp.coords, alleles=alleles,
                       file=tmp.files["genos"], verbose=verbose-1)

  if(verbose > 0){
    msg <- paste0("run ", exe.name, "...")
    write(msg, stdout())
  }
  args <- paste0("-o", task.id,
                 " -T", nb.starts,
                 " -C", nb.iters,
                 " -H", nb.samp.haplos,
                 " -K", nb.clusters,
                 " -S", seed)
  if(estim.haplos)
    args <- paste0(args, " -i")
  args <- paste0(args, " ", tmp.files["genos"])
  if(verbose > 0){
    msg <- paste(exe.name, args)
    write(msg, stdout())
  }
  tmp.files["stdouterr"] <- paste0(task.id, "_stdouterr.txt")
  ret <- system2(command=exe.name, args=args,
                 stdout=tmp.files["stdouterr"], stderr=tmp.files["stdouterr"],
                 wait=TRUE)

  if(file.exists(tmp.files["stdouterr"]))
    out$stdouterr <- readLines(tmp.files["stdouterr"])
  if(ret != 0){
    msg <- paste0("command '", exe.name, "' returned '", ret, "'")
    warning(msg)
  } else{
    if(verbose > 0){
      msg <- "load files..."
      write(msg, stdout())
    }
    out$genos <- readGenosFastphase(tmp.files["genos"],
                                    rownames(snp.coords))
    if(estim.haplos){
      tmp.files["haplos_switch"] <- paste0(task.id, "_hapguess_switch.out")
      tmp.files["haplos_indiv"] <- paste0(task.id, "_hapguess_indiv.out")
      tmp.files["finallik"] <- paste0(task.id, "_finallikelihoods")
      tmp.files["origchars"] <- paste0(task.id, "_origchars")
      out$haplos.switch <- readHaplosFastphase(tmp.files["haplos_switch"],
                                               rownames(snp.coords))
      out$haplos.indiv <- readHaplosFastphase(tmp.files["haplos_indiv"],
                                              rownames(snp.coords))
    }
    if(clean)
      for(f in tmp.files)
        file.remove(f)
  }

  return(out)
}
