## Contains functions useful to wrap findhap.
## https://aipl.arsusda.gov/software/findhap/

##' Write pedigree for findhap
##'
##' Write a pedigree into a file formatted for findhap.
##' @param ped data frame of pedigree with four columns in that order, genotype identifiers (same as in \code{X}), parent1 identifiers (considered as father/sire), parent2 identifiers (considered as mother/dam), and sex (as M or F); use \code{NA} for founders' parents (will be replaced by \code{0} in the output file)
##' @param ped.file path to the file in which the pedigree will be saved
##' @return nothing
##' @author Timothee Flutre
##' @export
writePedFileFindhap <- function(ped, ped.file){
  stopifnot(is.data.frame(ped),
            ncol(ped) == 4,
            is.character(ped.file))

  colnames(ped) <- c("ID", "Sire", "Dam", "Sex")

  utils::write.table(x=ped, file=ped.file, quote=FALSE, sep="\t",
                     row.names=FALSE, col.names=TRUE, na="0")
}

##' Write SNP information for findhap
##'
##' Write SNP information into a file formatted for findhap.
##' @param snp.coords data.frame with SNP identifiers as row names, and two compulsory columns in that order, chromosome identifiers and coordinate/position identifiers; chromosome identifiers should be numeric (the highest number for X-specific SNPs, Y-specific SNPs being unsupported); the "within" and "overall" columns will be added, as well as "n_chips" and "chip1" if there is a single microarray
##' @param snp.info.file path to the file in which SNP information will be saved
##' @return nothing
##' @author Timothee Flutre
##' @export
writeSnpInfoFindhap <- function(snp.coords, snp.info.file){
  stopifnot(is.data.frame(snp.coords),
            ! is.null(rownames(snp.coords)),
            ncol(snp.coords) >= 2,
            is.numeric(snp.coords[,1]),
            is.numeric(snp.coords[,2]),
            ! is.null(snp.info.file))
  if(! is.null(colnames(snp.coords)))
    stopifnot(! "within" %in% colnames(snp.coords),
              ! "overall" %in% colnames(snp.coords))

  snp.coords <- cbind(rownames(snp.coords), snp.coords)
  rownames(snp.coords) <- NULL
  colnames(snp.coords)[1:3] <- c("SNPname", "chrome", "location")
  last.colnames <- colnames(snp.coords)[seq(3,ncol(snp.coords))]
  tmp <- data.frame(overall=1:nrow(snp.coords),
                    within=NA)
  tmp$within <- tapply(snp.coords$location,
                       list(as.factor(snp.coords$chrome)),
                       function(x){1:length(x)})[[1]]
  snp.coords <- cbind(snp.coords[, c("SNPname","chrome")],
                      tmp[, c("within","overall")],
                      snp.coords[, seq(3,ncol(snp.coords))])
  colnames(snp.coords)[seq(5,ncol(snp.coords))] <- last.colnames
  if(ncol(snp.coords) == 5){
    snp.coords$n_chips <- rep(1, nrow(snp.coords))
    snp.coords$chip1 <- 1:nrow(snp.coords)
  }

  snp.coords <- snp.coords[order(snp.coords$chrome,
                                 snp.coords$location),]

  utils::write.table(x=snp.coords, file=snp.info.file, quote=FALSE, sep="\t",
                     row.names=FALSE, col.names=TRUE)
}

##' Write SNP genotypes for findhap
##'
##' Write SNP genotypes into a file formatted for findhap.
##' @param X matrix of bi-allelic SNP genotypes encoded in number of copies of the 2nd allele, i.e. as allele doses in \{0,1,2\}, with genotypes in rows and SNPs in columns; missing values should be encoded as NA; the maximum length of genotypes identifiers is 10 characters
##' @param genos.file path to the file in which genotypes will be saved
##' @return invisible data frame which is written in \code{genos.file}
##' @author Timothee Flutre
##' @examples
##' \dontrun{## fake genotypes
##' set.seed(1)
##' nb.inds <- 3
##' nb.snps <- 5
##' X <- matrix(data=sample(c(0,1,2,NA), size=nb.inds * nb.snps, replace=TRUE),
##'             nrow=nb.inds, ncol=nb.snps,
##'             dimnames=list(paste0("ind", 1:nb.inds),
##'                           paste0("snp", 1:nb.snps)))
##' f <- tempfile()
##' out <- writeGenosFindhap(X, f)
##' }
##' @export
writeGenosFindhap <- function(X, genos.file){
  stopIfNotValidGenosDose(X, check.noNA=FALSE)
  stopifnot(all(nchar(rownames(X)) <= 10),
            ! is.null(genos.file))

  X[is.na(X)] <- 5

  tmp <- cbind(sprintf("%10s", rownames(X)),
               sprintf("%10i", rep(1, nrow(X))),
               sprintf("%10i ", rep(ncol(X), nrow(X))),
               apply(X, 1, paste, collapse=""))
  rownames(tmp) <- NULL
  colnames(tmp) <- c("animal.id", "chip.id", "nb.snps", "snp.genotypes")

  tmp <- tmp[order(tmp[,"animal.id"]),]

  utils::write.table(x=tmp, file=genos.file, quote=FALSE, sep="",
                     row.names=FALSE, col.names=FALSE)

  invisible(tmp)
}

##' Write inputs for findhap
##'
##' Write inputs into files formatted for findhap
##' @param opt.file path to the options file in which findhap's configuration will be saved
##' @param X matrix of bi-allelic SNP genotypes encoded in number of copies of the 2nd allele, i.e. as allele doses in \{0,1,2\}, with genotypes in rows and SNPs in columns; missing values should be encoded as NA; if not NULL, will be used in priority even if \code{vcf.file} is not NULL
##' @param genos.file path to the file in which genotypes will be saved
##' @param chips if several microarrays were used, provide a vector with chip numbers which names are genotype identifiers (same as in \code{X})
##' @param snp.coords data.frame with SNP identifiers as row names, and two compulsory columns in that order, chromosome identifiers and coordinate/position identifiers; chromosome identifiers should be numeric; other optional column(s) contain the SNP order on the microarray(s) (maximum 10); compulsory if \code{X} is specified
##' @param snp.info.file path to the file in which SNP information will be saved
##' @param vcf.file path to the VCF file (if the bgzip index doesn't exist in the same directory, it will be created); used only if \code{X=NULL}
##' @param yieldSize number of records to yield each time the file is read from  (see \code{?TabixFile})
##' @param ped data frame of pedigree with four columns in that order, genotype identifiers (same as in \code{X}), parent1 identifiers (considered as father/sire), parent2 identifiers (considered as mother/dam), and sex (as M or F)
##' @param ped.file path to the file in which the pedigree will be saved
##' @param iters maximum iterations for haplotype imputation
##' @param Xchrom final chromosome number for X, or 0 if no X
##' @param maxlen maximum length of segments to check when haplotyping
##' @param minlen minimum length of segments to check when haplotyping
##' @param steps number of steps to get from maxlen to minlen
##' @param maxhap maximum different haplotypes within any segment
##' @param hapout output full haplotypes or not
##' @param genout output filled genotypes or not
##' @param damout output imputed genotypes of parents or not
##' @param listout output haplotypes list or not
##' @param errate fraction of miscalled genotypes allowed when matching
##' @param verbose verbosity level (0/1)
##' @return nothing
##' @author Timothee Flutre
##' @export
writeInputsFindhap <- function(opt.file,
															 X=NULL, genos.file, chips=NULL,
															 snp.coords=NULL, snp.info.file,
															 vcf.file=NULL, yieldSize=10000,
															 ped=NULL, ped.file=NULL,
                               iters=4, Xchrom=0,
                               maxlen=1000, minlen=100, steps=3,
                               maxhap=5000, hapout=1, genout=1,
                               damout=1, listout=-1, errate=0.003,
															 verbose=1){
  stopifnot(! all(is.null(X), is.null(vcf.file)))
	if(is.null(X)){
		file.prefix <- paste0(dirname(vcf.file), "/tmpfile-for-findhap")
		vcf2dosage(vcf.file=vcf.file, yieldSize=yieldSize,
							 gdose.file=paste0(file.prefix, "_gdose.txt.gz"),
							 ca.file=paste0(file.prefix, "_ca.txt.gz"),
							 verbose=verbose - 1)
		X <- t(as.matrix(
				utils::read.table(file=paste0(file.prefix, "_gdose.txt.gz"),
													check.names=FALSE)))
		snp.coords <- utils::read.table(file=paste0(file.prefix, "_ca.txt.gz"),
																		header=TRUE, stringsAsFactors=FALSE)
		snp.coords <- snp.coords[,1:2]
		if(! is.numeric(snp.coords$chr))
			snp.coords$chr <- chromNames2integers(x=snp.coords$chr)$renamed
		## file.remove(paste0(file.prefix, c("_gdose", "_ca"), ".txt.gz"))
	} else{
		stopifnot(! is.null(snp.coords))
	}
	geno.ids <- rownames(X)
	snp.ids <- colnames(X)

  stopIfNotValidGenosDose(X, check.noNA=FALSE)
  stopifnot(! is.null(genos.file),
            is.data.frame(snp.coords),
            all(snp.ids %in% rownames(snp.coords)),
            ncol(snp.coords) >= 2,
            is.numeric(snp.coords[,1]),
            ! is.null(snp.info.file))
  if(! is.null(chips))
    stopifnot(is.vector(chips),
              length(chips) == length(geno.ids),
              ! is.null(names(chips)),
              all(sort(names(chips)) == sort(geno.ids)))
  if(! is.null(ped))
    stopifnot(is.data.frame(ped),
              ncol(ped) == 4,
              all(geno.ids %in% ped[,1]),
              ! is.null(ped.file))

  if(verbose > 0){
    msg <- paste0("nb of genotypes: ", nrow(X),
                  "\nnb of SNPs: ", ncol(X),
                  "\nmissing SNP genotypes: ",
                  format(100 * sum(is.na(X)) / length(X), digits=3),
                  "%")
    write(msg, stdout())
  }

  cat(sprintf(fmt="%3i", iters), file=opt.file, sep="", append=FALSE)
  cat(sprintf(fmt="%8i", Xchrom), file=opt.file, sep="", append=TRUE)
  cat(sprintf(fmt="%8i", maxlen), file=opt.file, sep="", append=TRUE)
  cat(sprintf(fmt="%8i", minlen), file=opt.file, sep="", append=TRUE)
  cat(sprintf(fmt="%8i", steps), file=opt.file, sep="", append=TRUE)
  cat(sprintf(fmt="%8i", maxhap), file=opt.file, sep="", append=TRUE)
  cat(sprintf(fmt="%8i", hapout), file=opt.file, sep="", append=TRUE)
  cat(sprintf(fmt="%8i", genout), file=opt.file, sep="", append=TRUE)
  cat(sprintf(fmt="%8i", damout), file=opt.file, sep="", append=TRUE)
  cat(sprintf(fmt="%8i", listout), file=opt.file, sep="", append=TRUE)
  tmp <- "\niters  Xchrom  maxlen  minlen  steps  maxhap  hapout  genout  damout  listout\n\n"
  cat(tmp, file=opt.file, sep="", append=TRUE)
  cat(sprintf(fmt="%9.6f", errate), file=opt.file, sep="", append=TRUE)
  tmp <- "\nerrate\n"
  cat(tmp, file=opt.file, sep="", append=TRUE)

  writeGenosFindhap(X, genos.file)

  snp.coords <- snp.coords[snp.ids,]
  writeSnpInfoFindhap(snp.coords, snp.info.file)

  if(! is.null(ped)){
		## ped <- ped[ped$id %in% geno.ids,]
    writePedFileFindhap(ped, ped.file)
  }
}

##' Read genotypes from findhap
##'
##' Read imputed genotypes/haplotypes from findhap.
##' @param file path to the file containing imputed genotypes/haplotypes
##' @param snp.ids vector of SNP identifiers
##' @param input.haplos if TRUE, the input file is supposed to contain haplotypes
##' @param output.genos if TRUE, the output will correspond to genotypes
##' @return matrix
##' @author Timothee Flutre
##' @export
readGenosFindhap <- function(file, snp.ids, input.haplos=TRUE,
                             output.genos=TRUE){
  stopifnot(file.exists(file),
            is.vector(snp.ids))
  if(! input.haplos)
    stop("support for genotype input not yet implemented")
  if(! output.genos)
    stop("support for haplotype output not yet implemented")

  tmp <- utils::read.table(file=file, header=TRUE, sep="",
                           stringsAsFactors=FALSE,
                           colClasses=c("character", "numeric",
                                        "character"))
  if(length(unique(nchar(tmp[,3]))) > 1)
    stop("badly formatted file (different numbers of calls among individuals")
  if(unique(nchar(tmp[,3])) != length(snp.ids))
    stop("number of SNP identifiers != number of SNP calls")

  out <- matrix(data=do.call(rbind, lapply(strsplit(tmp[,3], ""), as.numeric)),
                nrow=nrow(tmp), ncol=length(snp.ids),
                dimnames=list(tmp[,1], snp.ids))
  idx <- which(out == 5)
  if(length(idx) > 0)
    out[idx] <- NA

  if(output.genos){
    idx <- which(out %in% c(3,4))
    if(length(idx) > 0)
      out[idx] <- 1
    idx <- which(out %in% c(6,7,8,9))
    if(length(idx) > 0)
      out[idx] <- NA
  }

  return(out)
}

##' Read outputs from findhap
##'
##' Read output files from findhap.
##' @param out.dir directory in which the output files are saved
##' @return list
##' @author Timothee Flutre
##' @export
readOutputsFindhap <- function(out.dir){
  out <- list()

  f <- paste0(out.dir, "/stat_snp.txt")
  out$stats.snps <- utils::read.table(file=f, header=TRUE, sep="", skip=1,
                                      stringsAsFactors=FALSE)

  f <- paste0(out.dir, "/stat_snp_imp.txt")
  out$stats.snps.imp <- utils::read.table(file=f, header=TRUE, sep="", skip=1,
                                          stringsAsFactors=FALSE)

  f <- paste0(out.dir, "/stat_anim.txt")
  out$stats.inds <- utils::read.table(file=f, header=TRUE, sep="", skip=1,
                                      stringsAsFactors=FALSE)

  f <- paste0(out.dir, "/stat_anim_imp.txt")
  out$stats.inds.imp <- utils::read.table(file=f, header=TRUE, sep="", skip=1,
                                          stringsAsFactors=FALSE)

  f <- paste0(out.dir, "/genotypes_imp.txt")
  out$genos.imp <- readGenosFindhap(file=f, snp.ids=out$stats.snps.imp[,1])

  f <- paste0(out.dir, "/genotypes_imp_chip0.txt")
  if(file.exists(f)){
    out$genos.imp.ungen <- readGenosFindhap(file=f,
                                            snp.ids=out$stats.snps.imp[,1])
  }

  f <- paste0(out.dir, "/parentage_test.txt")
  if(file.exists(f)){
    out$parentage.test <- readLines(f)
  }

  f <- paste0(out.dir, "/ref_pop.txt")
  if(file.exists(f)){
  }

  return(out)
}

##' Genotype imputation via findhap
##'
##' Impute SNP genotypes via findhap (\href{http://dx.doi.org/10.3168/jds.2012-5702}{VanRaden et al, 2013}).
##' @param X matrix of bi-allelic SNP genotypes encoded in number of copies of the 2nd allele, i.e. as allele doses in \{0,1,2\}, with genotypes in rows and SNPs in columns; missing values should be encoded as NA; if not NULL, will be used in priority even if \code{vcf.file} is not NULL
##' @param chips if several microarrays were used, provide a vector with chip numbers which names are genotype identifiers (same as in \code{X})
##' @param snp.coords data.frame with SNP identifiers as row names, and two compulsory columns in that order, chromosome identifiers and coordinate/position identifiers; the maximum length of SNP identifiers is 50 characters; chromosome identifiers should be numeric; other optional column(s) contain the SNP order on the microarray(s) (maximum 10); compulsory if \code{X} is specified
##' @param vcf.file path to the VCF file (if the bgzip index doesn't exist in the same directory, it will be created); used only if \code{X=NULL}
##' @param yieldSize number of records to yield each time the file is read from  (see \code{?TabixFile})
##' @param ped data frame of pedigree with four columns in that order, genotype identifiers (same as in \code{X}), parent1 identifiers (considered as father/sire), parent2 identifiers (considered as mother/dam), and sex (as M or F)
##' @param work.dir directory in which the input and output files will be saved
##' @param task.id identifier of the task (used in temporary and output file names)
##' @param params.findhap list of additional parameters to pass to findhap
##' @param nb.threads number of threads to pass to findhap
##' @param clean remove files: none, some (temporary only), all (temporary and results)
##' @param verbose verbosity level (0/1)
##' @return list
##' @author Timothee Flutre
##' @seealso \code{\link{writeInputsFindhap}}, \code{\link{readOutputsFindhap}}, \code{\link{runFastphase}}
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
##' snp.coords$chr <- as.numeric(sub("chr", "", snp.coords$chr))
##' out.fh <- runFindhap(X=X.na, snp.coords=snp.coords, clean="all")
##' X.impfh <- outfh$genos.imp
##'
##' ## assess imputation accuracy
##' genomes$haplos[[1]][1:4, 1:7]
##' genomes$genos[1:2, 1:7]
##' X.na[1:2, 1:6]
##' X.impfh[1:2, 1:6]
##' 100 * sum(X.impfh != X) / sum(is.na(X.na))
##' }
##' @export
runFindhap <- function(X=NULL, chips=NULL, snp.coords=NULL,
                       vcf.file=NULL, yieldSize=10000,
                       ped=NULL,
                       work.dir=getwd(), task.id="findhap",
                       params.findhap=NULL, nb.threads=1,
                       clean="none", verbose=1){
  exe.name <- "findhap"
  stopifnot(file.exists(Sys.which(exe.name)))
  stopifnot(! all(is.null(X), is.null(vcf.file)))
  stopifnot(is.character(task.id),
            length(task.id) == 1,
            all(! grepl(" ", task.id)),
            clean %in% c("none", "some", "all"))
  if(! is.null(params.findhap))
    stopifnot(is.list(params.findhap))

  out <- list()

  current.dir <- getwd()
  work.dir <- path.expand(work.dir)
  tmp.dir <- paste0(work.dir, "/", task.id, "_runFindhap")
  if(dir.exists(tmp.dir)){
    msg <- paste0("temporary directory ", tmp.dir, " already exists")
    stop(msg)
  } else
    dir.create(tmp.dir)
  setwd(tmp.dir)

  if(verbose > 0){
    msg <- "prepare input files..."
    write(msg, stdout())
  }
  tmp.files <- c()
  tmp.files["genos"] <- "genotypes.txt"
  tmp.files["snp.info"] <- "chromosome.data"
  tmp.files["ped"] <- "pedigree.file"
  tmp.files["opt"] <- "findhap.options"
  if(is.null(params.findhap))
    params.findhap <- list()
  if(! "iters" %in% names(params.findhap))
    params.findhap$iters <- 4
  if(! "Xchrom" %in% names(params.findhap))
    params.findhap$Xchrom <- 0
  if(! "maxlen" %in% names(params.findhap))
    params.findhap$maxlen <- 1000
  if(! "minlen" %in% names(params.findhap))
    params.findhap$minlen <- 100
  if(! "steps" %in% names(params.findhap))
    params.findhap$steps <- 3
  if(! "maxhap" %in% names(params.findhap))
    params.findhap$maxhap <- 5000
  if(! "hapout" %in% names(params.findhap))
    params.findhap$hapout <- 1
  if(! "genout" %in% names(params.findhap))
    params.findhap$genout <- 1
  if(! "damout" %in% names(params.findhap))
    params.findhap$damout <- 1
  if(! "listout" %in% names(params.findhap))
    params.findhap$listout <- -1
  if(! "errate" %in% names(params.findhap))
    params.findhap$errate <- 0.003

  ## TODO: use temporary IDs for the genotypes because of the limit of length 10

  writeInputsFindhap(opt.file=tmp.files["opt"],
                     X=X, genos.file=tmp.files["genos"], chips=chips,
                     snp.coords=snp.coords,
										 vcf.file=vcf.file, yieldSize=yieldSize,
                     snp.info.file=tmp.files["snp.info"],
                     ped=ped, ped.file=tmp.files["ped"],
                     iters=params.findhap$iters,
                     Xchrom=params.findhap$Xchrom,
                     maxlen=params.findhap$maxlen,
                     minlen=params.findhap$minlen,
                     steps=params.findhap$steps,
                     maxhap=params.findhap$maxhap,
                     hapout=params.findhap$hapout,
                     genout=params.findhap$genout,
                     damout=params.findhap$damout,
                     listout=params.findhap$listout,
                     errate=params.findhap$errate,
                     verbose=verbose)

  if(verbose > 0){
    msg <- paste0("run ", exe.name, "...")
    write(msg, stdout())
  }
  args <- nb.threads
  if(verbose > 0){
    msg <- paste(exe.name, args)
    write(msg, stdout())
  }
  tmp.files["stdouterr"] <- paste0(task.id, "_stdouterr.txt")
  ret <- system2(command=exe.name, args=args,
                 stdout=tmp.files["stdouterr"], stderr=tmp.files["stdouterr"],
                 wait=TRUE)

  ## if(ret != 0){
  ##   msg <- paste0("command '", exe.name, "' returned '", ret, "'")
  ##   warning(msg)
  ## } else{
  ##   if(verbose > 0){
  ##     msg <- "load files..."
  ##     write(msg, stdout())
  ##   }
  ##   out <- readOutputsFindhap(out.dir)
  ##   if(verbose > 0){
  ##     msg <- paste0("nb of SNPs with no Mendelian errors after imputation: ",
  ##                   sum(out$stats.snps.imp$MendelianErr == 0))
  ##     write(msg, stdout())
  ##   }
	## }

	## if(file.exists(tmp.files["stdouterr"]))
	## 	out$stdouterr <- readLines(tmp.files["stdouterr"])
	## if(clean != "none"){
	## 	for(f in tmp.files)
	## 		if(file.exists(f))
	## 			file.remove(f)
  ##   if(clean == "all")
  ##     unlink(tmp.dir, recursive=TRUE)
  ## }

  setwd(current.dir)

  return(out)
}
