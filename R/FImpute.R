## Contains functions useful to wrap FImpute.

##' Write pedigree for FImpute
##'
##' Write a pedigree into a file formatted for FImpute.
##' @param ped data frame of pedigree with four columns in that order, genotype identifiers (same as in \code{X}), parent1 identifiers (considered as father/sire), parent2 identifiers (considered as mother/dam), and sex (as M or F); for the parents, use \code{0} for founders
##' @param ped.file path to the file in which the pedigree will be saved
##' @return nothing
##' @author Timothee Flutre
##' @export
writePedFileFimpute <- function(ped, ped.file){
  stopifnot(is.data.frame(ped),
            ncol(ped) == 4,
            is.character(ped.file))

  colnames(ped) <- c("ID", "Sire", "Dam", "Sex")

  utils::write.table(x=ped, file=ped.file, quote=FALSE, sep="\t",
                     row.names=FALSE, col.names=TRUE)
}

##' Write SNP information for FImpute
##'
##' Write SNP information into a file formatted for FImpute.
##' @param snp.coords data.frame with SNP identifiers as row names, and two compulsory columns in that order, chromosome identifiers and coordinate/position identifiers; the maximum length of SNP identifiers is 50 characters; chromosome identifiers should be numeric; other optional column(s) contain the SNP order on the microarray(s) (maximum 10)
##' @param snp.info.file path to the file in which SNP information will be saved
##' @return nothing
##' @author Timothee Flutre
##' @export
writeSnpInfoFimpute <- function(snp.coords, snp.info.file){
  stopifnot(is.data.frame(snp.coords),
            ! is.null(rownames(snp.coords)),
            all(nchar(rownames(snp.coords)) <= 50),
            ncol(snp.coords) >= 2,
            is.numeric(snp.coords[,1]),
            ! is.null(snp.info.file))

  snp.coords <- cbind(rownames(snp.coords), snp.coords)
  rownames(snp.coords) <- NULL
  colnames(snp.coords) <- c("SNP_ID", "Chr", "Pos")
  if(ncol(snp.coords) == 3){
    snp.coords$Chip1 <- 1:nrow(snp.coords)
  }

  utils::write.table(x=snp.coords, file=snp.info.file, quote=FALSE, sep="\t",
                     row.names=FALSE, col.names=TRUE)
}

##' Write SNP genotypes for FImpute
##'
##' Write SNP genotypes into a file formatted for FImpute.
##' @param X matrix of bi-allelic SNP genotypes encoded in number of copies of the 2nd allele, i.e. as allele doses in {0,1,2}, with genotypes in rows and SNPs in columns; missing values should be encoded as NA; the maximum length of genotypes identifiers is 30 characters
##' @param genos.file path to the file in which genotypes will be saved
##' @param chips if several microarrays were used, provide a vector with chip numbers which names are genotype identifiers (same as in \code{X})
##' @return nothing
##' @author Timothee Flutre
##' @export
writeGenosFimpute <- function(X, genos.file, chips=NULL){
  stopIfNotValidGenosDose(X, check.noNA=FALSE)
  stopifnot(all(nchar(rownames(X)) <= 30),
            ! is.null(genos.file))
  if(! is.null(chips))
    stopifnot(is.vector(chips),
              length(chips) == nrow(X),
              ! is.null(names(chips)),
              all(sort(names(chips)) == sort(rownames(X))))

  if(is.null(chips))
    chips <- stats::setNames(rep(1, nrow(X)), rownames(X))

  X[is.na(X)] <- 5

  tmp <- cbind(rownames(X), chips,
               apply(X, 1, paste, collapse=""))
  rownames(tmp) <- NULL
  colnames(tmp) <- c("ID", "Chip", "Call...")

  utils::write.table(x=tmp, file=genos.file, quote=FALSE, sep="\t",
                     row.names=FALSE, col.names=TRUE)
}

##' Write inputs for FImpute
##'
##' Write inputs into files formatted for FImpute.
##' @param ctl.file path to the control file in which FImpute's configuration will be saved
##' @param X matrix of bi-allelic SNP genotypes encoded in number of copies of the 2nd allele, i.e. as allele doses in {0,1,2}, with genotypes in rows and SNPs in columns; missing values should be encoded as NA; the maximum length of genotypes identifiers is 30 characters
##' @param genos.file path to the file in which genotypes will be saved
##' @param chips if several microarrays were used, provide a vector with chip numbers which names are genotype identifiers (same as in \code{X})
##' @param snp.coords data.frame with SNP identifiers as row names, and two compulsory columns in that order, chromosome identifiers and coordinate/position identifiers; the maximum length of SNP identifiers is 50 characters; chromosome identifiers should be numeric; other optional column(s) contain the SNP order on the microarray(s) (maximum 10)
##' @param snp.info.file path to the file in which SNP information will be saved
##' @param ped data frame of pedigree with four columns in that order, genotype identifiers (same as in \code{X}), parent1 identifiers (considered as father/sire), parent2 identifiers (considered as mother/dam), and sex (as M or F)
##' @param ped.file path to the file in which the pedigree will be saved
##' @param title character identifying the analysis
##' @param out.dir directory in which the output files will be saved
##' @param hap.lib.file path to a file containing a haplotype library
##' @param add.ungen to add ungenotyped genotypes in the imputation process, provide a list such as \code{list(min.fsize=4, output.min.fsize=4, output.min.call.rate=0.9)}
##' @param parentage.test check for parentage errors
##' @param ped.depth to set a maximum number of generations to be traced for family imputation, provide a number
##' @param turn.off.fam turn off family imputation
##' @param turn.off.pop turn off population imputation
##' @param save.partial save partial calls
##' @param save.genos save genotypes instead of haplotypes (encoded as {0,1,2,5/NA})
##' @param save.hap.lib save the haplotype library built from reference individuals
##' @param random.fill ramdom filling (imputation) based on allele frequency; useful to assess minimum accuracy
##' @param nb.jobs number of jobs
##' @param verbose verbosity level (0/1)
##' @return nothing
##' @author Timothee Flutre
##' @export
writeInputsFimpute <- function(ctl.file,
                                  X, genos.file, chips=NULL,
                                  snp.coords, snp.info.file,
                                  ped=NULL, ped.file=NULL,
                                  title="FImpute",
                                  out.dir,
                                  hap.lib.file=NULL, add.ungen=NULL,
                                  parentage.test=NULL, ped.depth=NULL,
                                  turn.off.fam=FALSE, turn.off.pop=FALSE,
                                  save.partial=FALSE, save.genos=FALSE,
                                  save.hap.lib=FALSE, random.fill=FALSE,
                                  nb.jobs=1, verbose=1){
  stopIfNotValidGenosDose(X, check.noNA=FALSE)
  stopifnot(all(nchar(rownames(X)) <= 30),
            ! is.null(genos.file),
            is.data.frame(snp.coords),
            all(colnames(X) %in% rownames(snp.coords)),
            all(nchar(rownames(snp.coords)) <= 50),
            ncol(snp.coords) >= 2,
            is.numeric(snp.coords[,1]),
            ! is.null(snp.info.file),
            is.character(title),
            length(title) == 1,
            is.character(out.dir),
            length(out.dir) == 1,
            is.logical(turn.off.fam),
            is.logical(turn.off.pop),
            is.logical(save.partial),
            is.logical(save.genos),
            is.logical(save.hap.lib),
            is.logical(random.fill),
            is.numeric(nb.jobs),
            nb.jobs > 0)
  if(! is.null(chips))
    stopifnot(is.vector(chips),
              length(chips) == nrow(X),
              ! is.null(names(chips)),
              all(sort(names(chips)) == sort(rownames(X))))
  if(! is.null(ped))
    stopifnot(is.data.frame(ped),
              ncol(ped) == 4,
              all(rownames(X) %in% rownames(ped)),
              ! is.null(ped.file))
  if(! is.null(hap.lib.file))
    stop("support for 'hap.lib.file' not yet implemented")
  if(! is.null(add.ungen))
    stopifnot(is.list(add.ungen),
              ! is.null(names(add.ungen)),
              all(c("min.fsize", "output.min.fsize", "output.min.call.rate")
                  %in% names(add.ungen)))
  if(! is.null(parentage.test))
    stop("support for 'parentage.test' not yet implemented")
  if(! is.null(ped.depth))
    stopifnot(is.numeric(ped.depth))

  txt <- paste0("title=\"", title, "\";\n")
  cat(txt, file=ctl.file)

  writeGenosFimpute(X, genos.file)
  txt <- paste0("genotype_file=\"", genos.file, "\";\n")
  cat(txt, file=ctl.file, append=TRUE)

  snp.coords <- snp.coords[colnames(X),]
  writeSnpInfoFimpute(snp.coords, snp.info.file)
  txt <- paste0("snp_info_file=\"", snp.info.file, "\";\n")
  cat(txt, file=ctl.file, append=TRUE)

  if(! is.null(ped)){
    writePedFileFimpute(ped, ped.file)
    txt <- paste0("ped_file=\"", ped.file, "\";\n")
    cat(txt, file=ctl.file, append=TRUE)
  }

  ## txt <- paste0("hap_lib_file=\"", , "\";\n")
  ## cat(txt, file=ctl.file, append=TRUE)

  if(dir.exists(out.dir))
    unlink(out.dir)
  txt <- paste0("output_folder=\"", out.dir, "\";\n")
  cat(txt, file=ctl.file, append=TRUE)

  if(! is.null(add.ungen)){
    txt <- paste0("add_ungen",
                  " /min_fsize=", add.ungen[["min.fsize"]],
                  " /output_min_fsize=", add.ungen[["output.min.fsize"]],
                  " /output_min_call_rate=", add.ungen[["output.min.call.rate"]],
                  ";\n")
    cat(txt, file=ctl.file, append=TRUE)
  }

  ## txt <- paste0("parentage_test", " /off", ";\n")
  ## cat(txt, file=ctl.file, append=TRUE)

  ## txt <- paste0("exclude_snp=\"", , "\";\n")
  ## cat(txt, file=ctl.file, append=TRUE)

  ## txt <- paste0("exclude_chr=\"", , "\";\n")
  ## cat(txt, file=ctl.file, append=TRUE)

  ## txt <- paste0("exclude_chip=\"", , "\";\n")
  ## cat(txt, file=ctl.file, append=TRUE)

  txt <- paste0("njob=", nb.jobs, ";\n")
  cat(txt, file=ctl.file, append=TRUE)

  ## txt <- paste0("chmod=\"", , "\";\n")
  ## cat(txt, file=ctl.file, append=TRUE)

  if(! is.null(ped.depth)){
    txt <- paste0("ped_depth=", ped.depth, ";\n")
    cat(txt, file=ctl.file, append=TRUE)
  }

  ## txt <- paste0("min_nprg_imp=\"", , "\";\n")
  ## cat(txt, file=ctl.file, append=TRUE)

  ## txt <- paste0("min_nsib_imp=\"", , "\";\n")
  ## cat(txt, file=ctl.file, append=TRUE)

  ## txt <- paste0("min_segm_len_fam=\"", , "\";\n")
  ## cat(txt, file=ctl.file, append=TRUE)

  ## txt <- paste0("trim_segm_fam=\"", , "\";\n")
  ## cat(txt, file=ctl.file, append=TRUE)

  ## txt <- paste0("ref=\"", , "\";\n")
  ## cat(txt, file=ctl.file, append=TRUE)

  ## txt <- paste0("target=\"", , "\";\n")
  ## cat(txt, file=ctl.file, append=TRUE)

  ## txt <- paste0("sw_shrink_factor=\"", , "\";\n")
  ## cat(txt, file=ctl.file, append=TRUE)

  ## txt <- paste0("sw_overlap=\"", , "\";\n")
  ## cat(txt, file=ctl.file, append=TRUE)

  ## txt <- paste0("sw_min_size=\"", , "\";\n")
  ## cat(txt, file=ctl.file, append=TRUE)

  ## txt <- paste0("sw_max_size=\"", , "\";\n")
  ## cat(txt, file=ctl.file, append=TRUE)

  ## txt <- paste0("trim_segm_pop=\"", , "\";\n")
  ## cat(txt, file=ctl.file, append=TRUE)

  if(turn.off.fam){
    txt <- paste0("turnoff_fam;\n")
    cat(txt, file=ctl.file, append=TRUE)
  }

  if(turn.off.pop){
    txt <- paste0("turnoff_pop;\n")
    cat(txt, file=ctl.file, append=TRUE)
  }

  if(save.partial){
    txt <- paste0("save_partial;\n")
    cat(txt, file=ctl.file, append=TRUE)
  }

  if(save.genos){
    txt <- paste0("save_genotype;\n")
    cat(txt, file=ctl.file, append=TRUE)
  }

  if(save.hap.lib){
    txt <- paste0("save_hap_lib", " /diplotype", ";\n")
    cat(txt, file=ctl.file, append=TRUE)
  }

  if(random.fill){
    txt <- paste0("random_fill;\n")
    cat(txt, file=ctl.file, append=TRUE)
  }

  ## TODO: gzip the output?
  ## txt <- paste0("system=\"", , "\";\n")
  ## cat(txt, file=ctl.file, append=TRUE)
}

##' Read genotypes from FImpute
##'
##' Read imputed genotypes/haplotypes from FImpute.
##' @param file path to the file containing imputed genotypes/haplotypes
##' @param snp.ids vector of SNP identifiers
##' @param input.haplos if TRUE, the input file is supposed to contain haplotypes
##' @param output.genos if TRUE, the output will correspond to genotypes
##' @return matrix
##' @author Timothee Flutre
##' @export
readGenosFimpute <- function(file, snp.ids, input.haplos=TRUE,
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

##' Read outputs from FImpute
##'
##' @param out.dir directory in which the output files are saved
##' Read output files from FImpute.
##' @return nothing
##' @author Timothee Flutre
##' @export
readOutputsFimpute <- function(out.dir){
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
  out$genos.imp <- readGenosFimpute(file=f, snp.ids=out$stats.snps.imp[,1])

  f <- paste0(out.dir, "/genotypes_imp_chip0.txt")
  if(file.exists(f)){
    out$genos.imp.ungen <- readGenosFimpute(file=f,
                                            snp.ids=out$stats.snps.imp[,1])
  }

  f <- paste0(out.dir, "/org_vs_imp.txt")
  if(file.exists(f)){
  }

  f <- paste0(out.dir, "/ref_pop.txt")
  if(file.exists(f)){
  }

  return(out)
}

##' Genotype imputation via FImpute
##'
##' Impute SNP genotypes via FImpute (\href{http://dx.doi.org/10.1186/1471-2164-15-478}{Sargolzaei et al, 2014}).
##' @param X matrix of bi-allelic SNP genotypes encoded in number of copies of the 2nd allele, i.e. as allele doses in {0,1,2}, with genotypes in rows and SNPs in columns
##' @param chips if several microarrays were used, provide a vector with chip numbers which names are genotype identifiers (same as in \code{X})
##' @param snp.coords data.frame with SNP identifiers as row names, and two compulsory columns in that order, chromosome identifiers and coordinate/position identifiers; the maximum length of SNP identifiers is 50 characters; chromosome identifiers should be numeric; other optional column(s) contain the SNP order on the microarray(s) (maximum 10)
##' @param ped data frame of pedigree with four columns in that order, genotype identifiers (same as in \code{X}), parent1 identifiers (considered as father/sire), parent2 identifiers (considered as mother/dam), and sex (as M or F)
##' @param work.dir directory in which the input and output files will be saved
##' @param task.id identifier of the task (used in temporary and output file names)
##' @param params.fimpute list of additional parameters to pass to FImpute
##' @param clean remove files
##' @param verbose verbosity level (0/1)
##' @return list
##' @author Timothee Flutre
##' @seealso \code{\link{writeInputsFimpute}}, \code{\link{readOutputsFimpute}}, \code{\link{runFastphase}}
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
##' out.FI <- runFimpute(X=X.na, snp.coords=snp.coords, clean=TRUE)
##' X.imp.FI <- out.FI$genos.imp
##'
##' ## assess imputation accuracy
##' genomes$haplos[[1]][1:4, 1:7]
##' genomes$genos[1:2, 1:7]
##' X.na[1:2, 1:6]
##' X.imp.FI[1:2, 1:6]
##' 100 * sum(X.imp.FI != X) / sum(is.na(X.na))
##' }
##' @export
runFimpute <- function(X, chips=NULL, snp.coords, ped=NULL, work.dir=getwd(),
                       task.id="fimpute", params.fimpute=NULL, clean=FALSE,
                       verbose=1){
  exe.name <- "FImpute"
  stopifnot(file.exists(Sys.which(exe.name)))
  stopIfNotValidGenosDose(X, check.noNA=FALSE)
  stopifnot(is.character(task.id),
            length(task.id) == 1,
            all(! grepl(" ", task.id)))
  if(! is.null(params.fimpute))
    stopifnot(is.list(params.fimpute))

  out <- list()

  work.dir <- path.expand(work.dir)
  task.id <- paste0(work.dir, "/", task.id)
  out.dir <- paste0(task.id, "_outputs")

  if(verbose > 0){
    msg <- "prepare input files..."
    write(msg, stdout())
  }
  tmp.files <- c()
  tmp.files["ctl"] <- paste0(task.id, "_input-control.txt")
  tmp.files["genos"] <- paste0(task.id, "_input-genos.txt")
  tmp.files["snp.info"] <- paste0(task.id, "_input-snp-info.txt")
  tmp.files["ped"] <- paste0(task.id, "_input-ped.txt")
  if(is.null(params.fimpute))
    params.fimpute <- list(hap.lib.file=NULL,
                           turn.off.fam=FALSE,
                           turn.off.pop=FALSE,
                           save.partial=FALSE,
                           save.genos=FALSE,
                           save.hap.lib=FALSE,
                           random.fill=FALSE,
                           nb.jobs=1)
  writeInputsFimpute(ctl.file=tmp.files["ctl"],
                     X=X, genos.file=tmp.files["genos"], chips=chips,
                     snp.coords=snp.coords,
                     snp.info.file=tmp.files["snp.info"],
                     ped=ped, ped.file=tmp.files["ped"],
                     out.dir=out.dir,
                     hap.lib.file=params.fimpute$hap.lib.file,
                     turn.off.fam=params.fimpute$turn.off.fam,
                     turn.off.pop=params.fimpute$turn.off.pop,
                     save.partial=params.fimpute$save.partial,
                     save.genos=params.fimpute$save.genos,
                     save.hap.lib=params.fimpute$save.hap.lib,
                     random.fill=params.fimpute$random.fill,
                     nb.jobs=params.fimpute$nb.jobs,
                     verbose=verbose-1)

  if(verbose > 0){
    msg <- paste0("run ", exe.name, "...")
    write(msg, stdout())
  }
  args <- paste(tmp.files["ctl"], "-o")
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
    out <- readOutputsFimpute(out.dir)
    if(clean)
      for(f in tmp.files)
        if(file.exists(f))
          file.remove(f)
  }

  return(out)
}
