## Contains functions useful for quantitative genetics.

##' Convert genotypes
##'
##' Convert SNP genotype data from alleles (say, "AA" and "AT") to minor allele doses (here, 0 and 1 if "T" is the minor allele). Not particularly efficient, but at least it exists.
##' @param x data.frame with SNPs in rows and individuals in columns, the SNP identifiers being in the first column
##' @param na.string a character to be interpreted as NA values
##' @param verbose verbosity level (0/1)
##' @return list of a matrix (allele doses, SNPs in columns and individuals in rows) and a vector (minor alleles)
##' @author Timothee Flutre
##' @export
alleles2dose <- function(x, na.string="--", verbose=1){
  stopifnot(is.data.frame(x),
            ! is.null(colnames(x)))

  snp.names <- x[,1]
  P <- length(snp.names)
  ind.names <- colnames(x)[-1]
  N <- length(ind.names)
  if(verbose > 0){
    txt <- paste0(P, " SNPs and ", N, " individuals")
    write(txt, stdout())
  }

  geno.doses <- matrix(data=NA, nrow=N, ncol=P,
                       dimnames=list(ind=ind.names, snp=snp.names))
  alleles <- matrix(data=NA, nrow=P, ncol=2,
                    dimnames=list(snp.names, c("minor", "major")))

  for(p in 1:P){ # for each SNP
    raw.genos <- as.character(unlist(x[p, -1]))
    raw.genos[raw.genos == na.string] <- NA
    if(all(is.na(raw.genos))){
      next
    }
    tmp <- do.call(c, strsplit(raw.genos[! is.na(raw.genos)], ""))
    distinct.alleles <- sort(unique(tmp))
    allele.counts <- sort(table(tmp))
    if(length(distinct.alleles) > 2){ # SNP with more than 2 alleles
      stop("SNP ", paste0(x[p,1], " has more than 2 alleles"))
    } else if(length(distinct.alleles) == 2){ # SNP with exactly 2 alleles
      alleles[p, "minor"] <- names(allele.counts)[1]
      alleles[p, "major"] <- names(allele.counts)[2]
      raw.genos <- gsub(pattern=paste0(alleles[p, "minor"], alleles[p, "minor"]),
                        replacement="2", x=raw.genos)
      raw.genos <- gsub(pattern=paste(paste0(alleles[p, "minor"], alleles[p, "major"]),
                            paste0(alleles[p, "major"], alleles[p, "minor"]), sep="|"),
                        replacement="1", x=raw.genos)
      raw.genos <- gsub(pattern=paste0(alleles[p, "major"], alleles[p, "major"]),
                        replacement="0", x=raw.genos)
    } else if(length(distinct.alleles) == 1){ # SNP with only 1 allele
      alleles[p, "major"] <- names(allele.counts)[1]
      raw.genos <- gsub(pattern=paste0(alleles[p, "major"], alleles[p, "major"]),
                        replacement="0", x=raw.genos)
    }
    raw.genos <- as.numeric(raw.genos)
    if(! all(names(table(raw.genos, useNA="no")) %in% c("0", "1", "2")))
      stop("check SNP ", paste0(x[p,1], " (row ", p, ")"))
    geno.doses[,p] <- raw.genos
  }

  return(list(geno.doses=geno.doses,
              alleles=alleles))
}

.isValidGenosDose <- function(X, check.coln=TRUE, check.rown=TRUE,
                              check.na=TRUE){
  all(! is.null(X),
      is.matrix(X),
      sum(X < 0, na.rm=TRUE) == 0,
      sum(X > 2, na.rm=TRUE) == 0,
      ifelse(check.coln,
             ! is.null(colnames(X)) & ! anyDuplicated(colnames(X)),
             TRUE),
      ifelse(check.rown,
             ! is.null(rownames(X)),
             TRUE),
      ifelse(check.na, sum(is.na(X)) == 0, TRUE))
}

##' Convert genotypes
##'
##' Convert SNP genotypes from minor allele doses (say, 0 and 1) into alleles (say, "AA" and "AT" with "T" being the minor allele). Not particularly efficient, but at least it exists.
##' @param X matrix of SNP genotypes encoded in number of copies of the 2nd allele, i.e. as allele doses in {0,1,2}, with individuals in rows and SNPs in columns
##' @param tX matrix with SNPs in rows and individuals in columns
##' @param alleles data.frame with SNPs in rows (names as row names) and alleles in columns (first is "minor", second is "major")
##' @param na.string a character used to replace NA values (for instance "??")
##' @return data.frame with SNPs in rows and individuals in columns
##' @author Timothee Flutre
##' @export
dose2alleles <- function(X=NULL, tX=NULL, alleles, na.string){
  stopifnot(xor(is.null(X), is.null(tX)),
            is.data.frame(alleles),
            ncol(alleles) == 2,
            colnames(alleles) == c("minor", "major"),
            ! is.null(row.names(alleles)))
  if(! is.null(X)){
    stopifnot(.isValidGenosDose(X, check.na=FALSE))
    tX <- t(X)
  }

  out <- as.data.frame(tX, row.names=rownames(tX), col.names=colnames(tX))
  for(i in 1:nrow(tX)){ # for each SNP
    idx <- which(tX[i,] == 0)
    if(length(idx) > 0)
      out[i, idx] <- paste0(alleles[i, "major"], alleles[i, "major"],
                            collapse="")
    idx <- which(tX[i,] == 1)
    if(length(idx) > 0)
      out[i, idx] <- paste0(alleles[i, "major"], alleles[i, "minor"],
                            collapse="")
    idx <- which(tX[i,] == 2)
    if(length(idx) > 0)
      out[i, idx] <- paste0(alleles[i, "minor"], alleles[i, "minor"],
                            collapse="")
    idx <- which(is.na(tX[i,]))
    if(length(idx) > 0)
      out[i, idx] <- na.string
  }

  return(out)
}

##' Missing genotypes
##'
##' Calculate the frequency of missing genotypes for each SNP.
##' @param X matrix of SNP genotypes encoded in number of copies of the 2nd allele, i.e. as allele doses in {0,1,2}, with individuals in rows and SNPs in columns
##' @return vector
##' @author Timothee Flutre
##' @export
calcFreqMissGenos <- function(X){
  stopifnot(.isValidGenosDose(X, check.coln=FALSE, check.rown=FALSE,
                              check.na=FALSE))
  N <- nrow(X)
  out <- apply(X, 2, function(x){sum(is.na(x)) / N})
  return(out)
}

##' Missing genotypes
##'
##' Plot missing SNP genotypes as a grid.
##' Data will be represented in black if missing, white otherwise.
##' @param X matrix of SNP genotypes encoded in number of copies of the 2nd allele, i.e. as allele doses in {0,1,2}, with individuals in rows and SNPs in columns; NA if missing
##' @param main an overall title for the plot
##' @param xlab a title for the x axis
##' @param ylab a title for the y axis
##' @return nothing
##' @author Timothee Flutre
##' @export
plotGridMissGenos <- function(X, main="Missing genotypes", xlab="Individuals",
                              ylab="SNPs"){
  stopifnot(.isValidGenosDose(X, check.coln=FALSE, check.rown=FALSE,
                              check.na=FALSE))

  graphics::image(1:nrow(X), 1:ncol(X), is.na(X), col=c("white","black"),
        main=main, xlab=xlab, ylab=ylab)
}

##' Allele frequencies
##'
##' Estimate the frequency of the second allele for each SNP.
##' Note that the "second" allele may not be the "minor" allele (the least frequent).
##' @param X matrix of SNP genotypes encoded in number of copies of the 2nd allele, i.e. as allele doses in {0,1,2}, with individuals in rows and SNPs in columns; missing values should be encoded as NA in order to be ignored
##' @return vector
##' @seealso \code{\link{estimMaf}}
##' @author Timothee Flutre
##' @export
estimAf <- function(X){
  stopifnot(.isValidGenosDose(X, check.coln=FALSE, check.rown=FALSE,
                              check.na=FALSE))

  N <- nrow(X)
  P <- ncol(X)

  afs <- colMeans(X, na.rm=TRUE) / 2

  return(afs)
}

##' Minor allele frequencies
##'
##' Estimate the frequency of the minor allele for each SNP.
##' @param X matrix of SNP genotypes encoded in number of copies of the 2nd allele, i.e. as allele doses in {0,1,2}, with individuals in rows and SNPs in columns; missing values should be encoded as NA in order to be ignored; if X is NULL, afs should be specified
##' @param afs vector with 2nd allele frequencies; if NULL, X should be specified, and the frequencies will be estimated with \code{\link{estimAf}}
##' @return vector
##' @seealso \code{\link{estimAf}}
##' @author Timothee Flutre
##' @export
estimMaf <- function(X=NULL, afs=NULL){
  stopifnot(xor(is.null(X), is.null(afs)))

  if(is.null(afs))
    afs <- estimAf(X)

  mafs <- apply(rbind(afs, 1 - afs), 2, min)

  return(mafs)
}

##' Minor allele frequencies
##'
##' Plot the histogram of the minor allele frequency per SNP.
##' @param X matrix of SNP genotypes encoded in number of copies of the 2nd allele, i.e. as allele doses in {0,1,2}, with individuals in rows and SNPs in columns (optional if maf is not null); if X is not NULL, minor allele frequencies will be estimated via \code{\link{estimMaf}}
##' @param maf vector of minor allele frequencies (optional if X is not null)
##' @param main string for the main title
##' @param xlim limits of the x-axis
##' @param col color for the bars
##' @param border color for the border of the bars
##' @param las see ?par
##' @param breaks see ?hist
##' @param add.ml.beta if TRUE, add a curve corresponding to a Beta distribution with parameters inferred by maximum likelihood
##' @param verbose verbosity level (0/1)
##' @param ... arguments to be passed to hist()
##' @return nothing
##' @author Timothee Flutre
##' @export
plotHistMinAllelFreq <- function(X=NULL, maf=NULL, main=NULL, xlim=c(0,0.5),
                                 col="grey", border="white", las=1,
                                 breaks="FD", add.ml.beta=FALSE,
                                 verbose=1, ...){
  stopifnot(! is.null(X) || ! is.null(maf),
            is.logical(add.ml.beta))

  if(! is.null(X) & is.null(maf)){
    stopifnot(.isValidGenosDose(X, check.coln=FALSE, check.rown=FALSE,
                                check.na=FALSE))
    N <- nrow(X)
    P <- ncol(X)
    if(verbose > 0){
      txt <- paste0(P, " SNPs and ", N, " individuals")
      write(txt, stdout())
    }
    maf <- estimMaf(X)
  }

  if(is.null(main))
    main <- paste0("MAFs of ", length(maf), " SNPs")

  tmp <- graphics::hist(x=maf, xlab="Minor allele frequency", ylab="Number of SNPs",
              main=main, xlim=xlim, col=col, border=border, las=las,
              breaks=breaks, ...)

  if(add.ml.beta){
    fn.beta <- function(theta, x){
      sum(-stats::dbeta(x, theta[1], theta[2], log=TRUE))
    }
    fit <- stats::optim(par=c(0.3,0.3), fn=fn.beta, x=maf, method="L-BFGS-B",
                 lower=c(0.0001,0.0001), upper=c(1,1))
    ## http://stackoverflow.com/a/20078645/597069
    x <- seq(min(maf), max(maf), length=100)
    y <- stats::dbeta(x, fit$par[1], fit$par[2])
    y <- y * diff(tmp$mids[1:2]) * length(maf)
    graphics::lines(x, y, col="red")
    graphics::legend("topright", legend=paste0("Beta(", format(fit$par[1], digits=2),
                                     ", ", format(fit$par[2], digits=2),
                                     ") fitted via\nmaximum likelihood"),
           col="red", lty=1, bty="n")
  }
}

##' Minor allele frequencies
##'
##' Discard the SNPs with a minor allele frequency below the given threshold.
##' @param X matrix of SNP genotypes encoded in number of copies of the 2nd allele, i.e. as allele doses in {0,1,2}, with individuals in rows and SNPs in columns
##' @param mafs vector of minor allele frequencies; if NULL, will be estimated with \code{\link{estimMaf}}
##' @param thresh threshold on minor allele frequencies below which SNPs are ignored
##' @param verbose verbosity level (0/1)
##' @return matrix similar to X but possibly with less columns
##' @author Timothee Flutre
##' @export
discardSnpsLowMaf <- function(X, mafs=NULL, thresh=0.01, verbose=1){
  stopifnot(.isValidGenosDose(X, check.coln=FALSE, check.rown=FALSE,
                              check.na=FALSE))

  if(is.null(mafs))
    mafs <- estimMaf(X=X)

  snps.low <- mafs < thresh
  if(any(snps.low)){
    if(verbose > 0){
      txt <- paste0("skip ", sum(snps.low), " SNPs with MAF below ", thresh)
      write(txt, stdout())
    }
    idx.rm <- which(snps.low)
    X <- X[, -idx.rm, drop=FALSE]
  }

  return(X)
}

##' Convert genotypes
##'
##' Convert SNP genotypes to the file formats used by BimBam, as specified in its \href{http://www.haplotype.org/download/bimbam-manual.pdf}{manual}.
##' @param X matrix of SNP genotypes encoded in number of copies of the 2nd allele, i.e. as allele doses in {0,1,2}, with individuals in rows and SNPs in columns
##' @param tX matrix with SNPs in rows and individuals in columns
##' @param alleles data.frame with SNPs in rows (names as row names) and alleles in columns (first is "minor", second is "major")
##' @param file write the genotype data to this file if non NULL (for instance 'genotypes_bimbam.txt' or 'genotypes_bimbam.txt.gz', but don't use \code{\link[base]{gzfile}})
##' @param format BimBam's format in which the data should be saved (mean/basic)
##' @return invisible data.frame
##' @author Timothee Flutre
##' @export
dose2bimbam <- function(X=NULL, tX=NULL, alleles, file=NULL, format="mean"){
  stopifnot(xor(is.null(X), is.null(tX)),
            is.data.frame(alleles),
            ncol(alleles) == 2,
            colnames(alleles) == c("minor", "major"),
            ! is.null(row.names(alleles)),
            format %in% c("mean", "basic"))
  if(! is.null(X)){
    stopifnot(.isValidGenosDose(X, check.coln=FALSE, check.rown=FALSE))
    tX <- t(X)
  }

  if(format == "mean"){
    out <- cbind(alleles, tX)
    if(! is.null(file))
      utils::write.table(x=out, file=file, quote=FALSE, sep="\t",
                         row.names=TRUE, col.names=FALSE)
  } else if(format == "basic"){
    out <- dose2alleles(tX=tX, alleles=alleles, na.string="??")
    if(! is.null(file)){
      file <- removeFileExtension(file, "gz")
      cat(paste0(ncol(tX), "\n", nrow(tX), "\n"), file=file, sep="")
      sep <- ", "
      cat(paste0("IND", sep, paste0(colnames(out), collapse=sep), "\n"),
          file=file, sep="", append=TRUE)
      tmp <- do.call(rbind, lapply(1:nrow(out), function(i){
        paste0(rownames(out)[i], sep, paste0(out[i,], collapse=sep), "\n")
      }))
      cat(tmp, file=file, sep="", append=TRUE)
      system2(command="gzip", args=file)
    }
  }

  invisible(out)
}

##' Genetic relatedness
##'
##' Reformat the output of \code{related::coancestry} (http://frasierlab.wordpress.com/software/) into a matrix.
##' By default, off-diagonal elements  correspond to coancestry coefficients between two individuals, and each diagonal element corresponds to (1 + f) / 2 where f corresponds to the inbreeding coefficient of the given individual.
##' To learn more about genetic relatedness, see Weir et al (2006), Astle & Balding (2009) and Legarra (2016).
##' @param x list returned by \code{related::coancestry}
##' @param estim.coancestry name of the coancestry estimator (e.g. "dyadml") as used in \code{related::coancestry}
##' @param estim.inbreeding name of the inbreeding estimator (e.g. "LR") as used in \code{related::coancestry}
##' @param rel.type type of relatedness, "coancestries" or "relationships" (i.e. 2 x coancestries)
##' @param debug boolean (TRUE to check the output matrix is indeed symmetric)
##' @return matrix
##' @author Timothee Flutre
##' @export
coancestry2relmat <- function(x, estim.coancestry, estim.inbreeding=NULL,
                              rel.type="coancestries", debug=FALSE){
  stopifnot(estim.coancestry %in% colnames(x$relatedness),
            rel.type %in% c("coancestries", "relationships"))
  if(! is.null(estim.inbreeding)){
    stopifnot("inbreeding" %in% names(x))
    stopifnot(estim.inbreeding %in% colnames(x$inbreeding))
  }

  ind.ids <- unique(c(x$relatedness$ind1.id,
                      x$relatedness$ind2.id))
  nb.inds <- length(ind.ids)
  mat <- matrix(NA, nrow=nb.inds, ncol=nb.inds,
              dimnames=list(ind.ids, ind.ids))

  diag(mat) <- 1 / 2
  if(! is.null(estim.inbreeding)){
    stopifnot(nrow(x$inbreeding) == nb.inds)
    for(i in 1:nrow(x$inbreeding))
      mat[x$inbreeding$ind.id[i], x$inbreeding$ind.id[i]] <-
        (1 + x$inbreeding[[estim.inbreeding]][i]) / 2
  }

  for(i in 1:nrow(x$relatedness))
    mat[x$relatedness$ind1.id[i], x$relatedness$ind2.id[i]] <-
      x$relatedness[[estim.coancestry]][i]
  idx.na.upper <- which(upper.tri(mat) & is.na(mat))
  idx.equiv.lower <- which(lower.tri(t(mat)) & is.na(t(mat)))
  mat[idx.na.upper] <- mat[idx.equiv.lower]
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]

  if(debug){ # check that the matrix is symmetric
    for(i in 1:(nrow(mat)-1))
      for(j in (i+1):ncol(mat))
        if(mat[i,j] != mat[j,i])
          stop(paste0("matrix not symmetric at row=", i, ", col=", j))
  }

  if(rel.type == "relationships")
    mat <- 2 * mat

  return(mat)
}

##' Haplotype list to matrix
##'
##' Convert a list of haplotypes into a matrix.
##' @param haplos list of matrices (one per chromosome, with individuals in rows and markers in columns)
##' @return matrix
##' @author Timothee Flutre
##' @export
haplosList2Matrix <- function(haplos){
  stopifnot(is.list(haplos),
            all(sapply(haplos, class) == "matrix"))

  nb.chroms <- length(haplos)
  nb.haplos <- nrow(haplos[[1]]) # 2 x nb of individuals
  nb.snps <- sapply(haplos, ncol)
  P <- sum(nb.snps) # nb of SNPs

  H <- matrix(data=NA, nrow=nb.haplos, ncol=P,
              dimnames=list(rownames(haplos[[1]]),
                            do.call(c, lapply(haplos, colnames))))
  names(colnames(H)) <- NULL # remove "chr" names

  H[, 1:nb.snps[1]] <- haplos[[1]]
  for(chr in 2:nb.chroms)
    H[, (cumsum(nb.snps)[chr-1]+1):(cumsum(nb.snps)[chr])] <- haplos[[chr]]

  return(H)
}

##' Site frequency spectrum
##'
##' Convert the SFS of independent replicates into a matrix of allele dosage.
##' @param seg.sites list returned by scrm()
##' @param ind.ids vector with the identifiers of the individuals
##' @param snp.ids vector with the identifiers of the SNPs (if NULL, the SNP identifiers from seg.sites will be used if they aren't NULL, too)
##' @return matrix with diploid individuals in rows and SNPs in columns
##' @author Timothee Flutre
##' @export
segSites2allDoses <- function(seg.sites, ind.ids=NULL, snp.ids=NULL){
  stopifnot(is.list(seg.sites),
            length(unique(sapply(seg.sites, nrow))) == 1)
  if(! is.null(ind.ids))
    stopifnot(is.vector(ind.ids),
              is.character(ind.ids),
              length(ind.ids) == nrow(seg.sites[[1]]) / 2)
  if(! is.null(snp.ids))
    stopifnot(is.vector(snp.ids),
              is.character(snp.ids),
              length(snp.ids) == sum(sapply(seg.sites, ncol)))

  nb.inds <- nrow(seg.sites[[1]]) / 2 # nb of diploid individuals
  nb.snps <- sum(sapply(seg.sites, ncol)) # nb of SNPs
  X <- matrix(data=NA, nrow=nb.inds, ncol=nb.snps)
  if(! is.null(ind.ids))
    rownames(X) <- ind.ids
  if(! is.null(snp.ids)){
    colnames(X) <- snp.ids
  } else{
    tmp <- do.call(c, lapply(seg.sites, colnames))
    if(all(! is.null(tmp), length(tmp) == ncol(X)))
      colnames(X) <- tmp
  }

  j <- 1
  for(x in seq_along(seg.sites)){
    X[,j:(j+ncol(seg.sites[[x]])-1)] <-
      do.call(rbind, lapply(seq(1, 2*nb.inds, by=2), function(i){
        colSums(seg.sites[[x]][c(i,i+1),])
      }))
    j <- j + ncol(seg.sites[[x]])
  }

  return(X)
}

##' Site frequency spectrum
##'
##' Make a data.frame of SNP coordinates from the SFS of independent replicates.
##' @param seg.sites list returned by scrm()
##' @param snp.ids vector of identifiers (one per SNP)
##' @param chrom.len chromosome length (same for all)
##' @param prefix character string
##' @return data.frame with SNPs in rows and 2 columns (chr, pos)
##' @author Timothee Flutre
##' @export
segSites2snpCoords <- function(seg.sites, snp.ids, chrom.len, prefix="chr"){
  stopifnot(is.list(seg.sites),
            length(unique(sapply(seg.sites, nrow))) == 1)

  nb.chrs <- length(seg.sites)
  nb.snps.per.chr <- sapply(seg.sites, ncol)
  nb.snps <- sum(nb.snps.per.chr)

  ## fill up the output data.frame
  snp.coords <- data.frame(chr=rep(NA, nb.snps),
                           pos=-1,
                           row.names=snp.ids)
  snp.coords$chr <- rep(paste0(prefix, 1:nb.chrs), nb.snps.per.chr)
  snp.coords$pos <- as.numeric(do.call(c, lapply(seg.sites, colnames)))

  ## convert genomic positions from float to integer
  snp.coords$pos <- floor(snp.coords$pos)
  for(chr in levels(snp.coords$chr)){
    tmp <- snp.coords$pos[snp.coords$chr == chr]
    stopifnot(anyDuplicated(tmp) == 0,
              min(tmp) > 0,
              max(tmp) <= chrom.len)
  }

  return(snp.coords)
}

##' Genotypes
##'
##' Simulate SNP genotypes as allele dose additively encoded, i.e. 0,1,2.
##' @param nb.genos number of genotypes (i.e. individuals)
##' @param nb.snps number of SNPs
##' @param geno.ids vector of genotype identifiers (if NULL, will be "geno001", etc)
##' @param snp.ids vector of SNP identifiers (if NULL, will be "snp001", etc)
##' @param mafs vector of minor allele frequencies; by default, they are uniformly distributed between 0.05 and 0.5
##' @return matrix with individuals in rows and SNPs in columns
##' @author Peter Carbonetto [aut], Timothee Flutre [ctb]
##' @export
simulGenosDose <- function(nb.genos, nb.snps, geno.ids=NULL, snp.ids=NULL, mafs=NULL){
  if(! is.null(geno.ids))
    stopifnot(length(geno.ids) == nb.genos)
  if(! is.null(snp.ids))
    stopifnot(length(snp.ids) == nb.snps)
  if(! is.null(mafs))
    stopifnot(length(mafs) == nb.snps)

  if(is.null(geno.ids))
    geno.ids <- sprintf(fmt=paste0("geno%0", floor(log10(nb.genos))+1, "i"),
                        1:nb.genos)

  snp.ids <- sprintf(fmt=paste0("snp%0", floor(log10(nb.snps))+1, "i"),
                     1:nb.snps)
  if(is.null(mafs))
    mafs <- stats::setNames(0.05 + 0.45 * stats::runif(nb.snps), snp.ids)

  X <- matrix(data=(stats::runif(nb.genos * nb.snps) < mafs) +
                (stats::runif(nb.genos * nb.snps) < mafs),
              nrow=nb.genos, ncol=nb.snps, byrow=TRUE,
              dimnames=list(geno.ids, snp.ids))

  return(X)
}

##' Coalescent with recombination
##'
##' Simulate haplotypes according to an approximation to the coalescent with recombination named the Sequential Coalescent with Recombination Model. Requires the scrm package (Staab et al, 2014).
##' @param nb.inds diploids (thus nb of haplotypes is 2 * nb.inds)
##' @param ind.ids vector of identifiers (one per individual)
##' @param nb.reps number of independent loci that will be produced (could be seen as distinct chromosomes)
##' @param pop.mut.rate theta = 4 N0 mu
##' @param pop.recomb.rate rho = 4 N0 r
##' @param chrom.len in bp
##' @param other character vector of length 1 with other parameters to the simulator (e.g. time-specific parameters such as "-G 6.93 -eG 0.2 0.0 -eN 0.3 0.5")
##' @param nb.pops number of populations (\code{\link{kmeans}} will then be used to pair haplotypes into diploid individuals)
##' @param mig.rate migration rate = 4 N0 m (will be symmetric)
##' @param get.trees get gene genealogies in the Newick format
##' @param get.tmrca get time to most recent common ancestor and local tree lengths
##' @param verbose verbosity level (0/1/2)
##' @return list with haplotypes (list), genotypes as allele doses (matrix) and SNP coordinates (data.frame)
##' @author Timothee Flutre
##' @examples
##' \dontrun{## simulate haplotypes and genotypes in a single population
##' nb.genos <- 200
##' Ne <- 10^4
##' chrom.len <- 10^5
##' mu <- 10^(-8)
##' c <- 10^(-8)
##' genomes <- simulCoalescent(nb.inds=nb.genos,
##'                            pop.mut.rate=4 * Ne * mu * chrom.len,
##'                            pop.recomb.rate=4 * Ne * c * chrom.len,
##'                            chrom.len=chrom.len)
##' }
##' @export
simulCoalescent <- function(nb.inds=100,
                            ind.ids=NULL,
                            nb.reps=10,
                            pop.mut.rate=40,
                            pop.recomb.rate=40,
                            chrom.len=10^4,
                            other=NULL,
                            nb.pops=1,
                            mig.rate=5,
                            get.trees=FALSE,
                            get.tmrca=FALSE,
                            verbose=1){
  requireNamespaces("scrm")
  stopifnot(nb.inds > nb.pops)
  if(! is.null(other))
    stopifnot(is.character(other),
              length(other) == 1)

  out <- list()

  if(is.null(ind.ids))
    ind.ids <- sprintf(fmt=paste0("ind%0", floor(log10(nb.inds))+1, "i"),
                       1:nb.inds)

  ## simulate according to the SCRM
  nb.samples <- nb.inds * 2 # e.g. 2 chr1 in ind1, 2 chr1 in ind2, etc
  cmd <- paste0(nb.samples, " ", nb.reps)
  cmd <- paste0(cmd, " -t ", pop.mut.rate)
  cmd <- paste0(cmd, " -r ", pop.recomb.rate, " ", chrom.len)
  if(! is.null(other))
    cmd <- paste0(cmd, " ", other)
  if(get.trees)
    cmd <- paste0(cmd, " -T")
  if(get.tmrca)
    cmd <- paste0(cmd, " -L")
  cmd <- paste0(cmd, " -SC abs") # absolute seq positions in bp
  cmd <- paste0(cmd, " -oSFS") # print site freq spectrum, requires -t
  if(nb.pops > 1){
    cmd <- paste0(cmd, " -I ", nb.pops)
    nb.inds.per.pop <- rep(0, nb.pops)
    for(p in 1:(nb.pops-1)){
      nb.inds.per.pop[p] <- floor(nb.inds / nb.pops)
      cmd <- paste0(cmd, " ", 2 * nb.inds.per.pop[p])
    }
    nb.inds.per.pop[nb.pops] <- nb.inds - sum(nb.inds.per.pop)
    cmd <- paste0(cmd, " ", 2 * nb.inds.per.pop[nb.pops])
    cmd <- paste0(cmd, " ", mig.rate)
  }
  if(verbose > 1)
    message(cmd)
  sum.stats <- scrm::scrm(cmd)
  if(verbose > 1)
    print(utils::str(sum.stats))

  prefix <- "chr"
  names(sum.stats$seg_sites) <- paste0(prefix, 1:nb.reps)

  ## make a data.frame with SNP coordinates
  nb.snps.per.chr <- sapply(sum.stats$seg_sites, ncol)
  nb.snps <- sum(nb.snps.per.chr)
  snp.ids <- sprintf(fmt=paste0("snp%0", floor(log10(nb.snps))+1, "i"),
                     1:nb.snps)
  snp.coords <- segSites2snpCoords(sum.stats$seg_sites, snp.ids, chrom.len,
                                   prefix)
  out[["snp.coords"]] <- snp.coords

  ## randomize haplotypes to make diploid individuals, and set names
  out[["haplos"]] <- list()
  idx <- sample.int(nb.samples)
  if(nb.pops > 1){
    H <- haplosList2Matrix(sum.stats$seg_sites)
    kmH <- stats::kmeans(H, nb.pops)
  }
  for(chr in 1:nb.reps){
    if(nb.pops > 1)
      idx <- do.call(c, lapply(1:nb.pops, function(pop.idx){
        sample(x=which(kmH$cluster == pop.idx),
               sum(kmH$cluster == pop.idx))
      }))
    out$haplos[[chr]] <- sum.stats$seg_sites[[chr]][idx,]
    rownames(out$haplos[[chr]]) <- paste0(rep(ind.ids, each=2), "_h", 1:2)
    colnames(out$haplos[[chr]]) <-
      snp.ids[(ifelse(chr == 1, 1, 1 + cumsum(nb.snps.per.chr)[chr-1])):
      (cumsum(nb.snps.per.chr)[chr])]
  }
  names(out$haplos) <- names(sum.stats$seg_sites)

  ## make a matrix with genotypes encoded as allele dose
  X <- segSites2allDoses(out$haplos, ind.ids, snp.ids)
  if(verbose > 0){
    txt <- paste0("nb of SNPs: ", nb.snps)
    write(txt, stdout())
    print(sapply(sum.stats$seg_sites, ncol))
  }
  out[["genos"]] <- X

  if(get.trees)
    out[["trees"]] <- sum.stats$trees
  if(get.tmrca)
    out[["tmrca"]] <- sum.stats$tmrca

  return(out)
}

##' Genomes
##'
##' From a set of individual genomes simulated together, split them into a training (for estimation) and testing (for prediction) sets.
##' @param genomes list, e.g. returned by \code{\link{simulCoalescent}}
##' @param nb.inds.pred number of individuals for whom prediction will be performed
##' @return list composed of two lists, genomes and genomes.pred
##' @author Timothee Flutre
##' @export
splitGenomesTrainTest <- function(genomes, nb.inds.pred){
  stopifnot(is.list(genomes),
            all(c("haplos", "genos") %in% names(genomes)),
            is.list(genomes$haplos),
            is.matrix(genomes$genos),
            nb.inds.pred < nrow(genomes$haplos[[1]]))

  out <- list(genomes=NULL,
              genomes.pred=NULL)

  if(nb.inds.pred == 0){
    out$genomes <- genomes
  } else{
    ind.names <- getIndNamesFromHaplos(genomes$haplos)
    ind.names.train <- ind.names[1:(length(ind.names) - nb.inds.pred)]
    ind.names.test <- ind.names[(length(ind.names) - nb.inds.pred + 1):
                                length(ind.names)]
    out$genomes <- list(haplos=getHaplosInds(genomes$haplos, ind.names.train),
                        genos=subset(genomes$genos,
                                     ind.names %in% ind.names.train))
    out$genomes.pred <- list(haplos=getHaplosInds(genomes$haplos,
                                                  ind.names.test),
                             genos=subset(genomes$genos,
                                          ind.names %in% ind.names.test))
  }

  return(out)
}

##' Pi
##'
##' Calculate the average number of differences between pairs of sequences, named pi in the population genetics literature (Tajima, 1983; Wakeley, 2008).
##' @param haplos.chr matrix of 0's and 1's with haplotypes in rows and sites in columns
##' @return numeric(1)
##' @author Timothee Flutre
##' @export
calcAvgPwDiffBtwHaplos <- function(haplos.chr){
  stopifnot(is.matrix(haplos.chr),
            all(haplos.chr %in% c(0,1)))

  pi <- 0

  n <- nrow(haplos.chr)
  k <- as.matrix(stats::dist(haplos.chr, "manhattan"))
  k[lower.tri(k)] <- 0

  ## naive implementation:
  ## for(i in 1:(n-1))
  ##   for(j in (i+1):n)
  ##     pi <- pi + k[i,j]

  ## faster implementation:
  pi <- sum(rowSums(k))

  pi <- pi / choose(n, 2)

  return(pi)
}

##' Haplotypes
##'
##' Plot haplotypes as an image.
##' @param haplos matrix of haplotypes with haplotypes (2 per individual) in rows and sites in columns
##' @param main main title
##' @return nothing
##' @author Timothee Flutre
##' @export
plotHaplosMatrix <- function(haplos, main="Haplotypes"){
  stopifnot(is.matrix(haplos),
            ! is.null(dimnames(haplos)))

  opar <- graphics::par(mar=c(1,6,5,2))

  graphics::image(t(haplos)[,nrow(haplos):1], axes=FALSE, col=c("white","black"))

  graphics::title(main=main, line=3)
  graphics::axis(side=3, at=seq(1, ncol(haplos), length.out=7) / ncol(haplos),
       labels=colnames(haplos)[seq(1, ncol(haplos), length.out=7)])
  graphics::axis(side=2, at=rev(seq(1, nrow(haplos), length.out=10) / nrow(haplos)),
       labels=rownames(haplos)[seq(1, nrow(haplos), length.out=10)],
       las=1, padj=0)

  on.exit(graphics::par(opar))
}

##' Individual names
##'
##' Return the identifiers of all individuals
##' @param haplos list of matrices (one per chromosome, with individuals in rows)
##' @return vector
##' @author Timothee Flutre
##' @export
getIndNamesFromHaplos <- function(haplos){
  stopifnot(is.list(haplos))

  ind.names <- unique(do.call(c, lapply(haplos, rownames)))
  ind.names <- lapply(strsplit(ind.names, "_"), function(ind.name){
    paste(ind.name[1:(length(ind.name)-1)], collapse="_")
  })
  ind.names <- unique(do.call(c, ind.names))

  return(ind.names)
}

##' Haplotypes of an individual
##'
##' Retrieve both haplotypes for all chromosomes of a given individual.
##' @param haplos list with one matrix per chromosome, with haplotype in rows and SNPs in columns, as returned by \code{\link{simulCoalescent}}
##' @param ind.name identifier of the individual to retrieve
##' @return list similar to "haplos"
##' @author Timothee Flutre
##' @export
getHaplosInd <- function(haplos, ind.name){
  stopifnot(is.list(haplos),
            is.character(ind.name))

  ## for each chromosome, retrieve both haplotypes of the given individual
  haplos.ind <- lapply(haplos, function(haplos.chr){
    idx <- grep(paste0("^", ind.name, "_h[12]"), rownames(haplos.chr))
    if(length(idx) == 2){
      haplos.chr[idx, , drop=FALSE]
    } else{
      msg <- paste0("can't find haplotypes for '", ind.name, "'")
      stop(msg)
    }
  })

  return(haplos.ind)
}

##' Haplotypes of several individuals
##'
##' Retrieve both haplotypes for all chromosomes of certain individuals.
##' @param haplos list with one matrix per chromosome, with haplotype in rows and SNPs in columns, as returned by \code{\link{simulCoalescent}}
##' @param ind.names identifier of the individuals to retrieve
##' @return list similar to "haplos"
##' @author Timothee Flutre
##' @export
getHaplosInds <- function(haplos, ind.names){
  stopifnot(is.list(haplos),
            length(unique(sapply(haplos, class))) == 1,
            unique(sapply(haplos, class)) == "matrix",
            length(unique(sapply(haplos, nrow))) == 1, # same nb of individuals
            all(sapply(haplos, function(x){! is.null(rownames(x))})), # has row names
            is.character(ind.names))

  ## for each chromosome, retrieve both haplotypes of the given individuals
  haplos.inds <- lapply(haplos, function(haplos.chr){
    idx <- do.call(c, lapply(ind.names, grep, rownames(haplos.chr)))
    if(length(idx) == 0){
      msg <- paste0("can't find haplotypes for the given individuals")
      stop(msg)
    } else{
      haplos.chr[idx, , drop=FALSE]
    }
  })

  return(haplos.inds)
}

##' Gamete
##'
##' Make a gamete for a given chromosome by recombining two parental haplotypes (always starting with the first haplotype).
##' @param haplos.par.chr matrix containing both haplotypes of a parent for a given chromosome (must have dimnames, such as "ind37_h1" in rows and "snp00265" in columns)
##' @param loc.crossovers positions of the crossing overs (the coordinate of the first nucleotide is assumed to be 1, and crossing-overs are assumed to occur right after the given localization)
##' @return vector
##' @author Timothee Flutre
##' @export
makeGameteSingleIndSingleChrom <- function(haplos.par.chr, loc.crossovers){
  stopifnot(is.matrix(haplos.par.chr),
            ! is.null(dimnames(haplos.par.chr)),
            is.vector(loc.crossovers))

  nb.snps <- ncol(haplos.par.chr)
  loc.crossovers <- sort(loc.crossovers)
  stopifnot(max(loc.crossovers) < ncol(haplos.par.chr))

  ind.name <- strsplit(rownames(haplos.par.chr)[1],
                       "_")[[1]]
  ind.name <- paste(ind.name[1:(length(ind.name)-1)], collapse="_")
  gam.chr <- matrix(rep(NA, nb.snps), nrow=1,
                    dimnames=list(ind.name, colnames(haplos.par.chr)))

  ## always start with the first haplotype
  if(loc.crossovers[1] != 0)
    loc.crossovers <- c(0, loc.crossovers)
  for(l in 2:length(loc.crossovers))
    gam.chr[(loc.crossovers[l-1]+1):loc.crossovers[l]] <-
      haplos.par.chr[ifelse(l %% 2 == 0, 1, 2),
      (loc.crossovers[l-1]+1):loc.crossovers[l]]

  ## if necessary, fill the gamete until the end
  if(loc.crossovers[length(loc.crossovers)] < nb.snps)
    gam.chr[(loc.crossovers[length(loc.crossovers)]+1):nb.snps] <-
    haplos.par.chr[ifelse((length(loc.crossovers)+1) %% 2 == 0, 1, 2),
    (loc.crossovers[length(loc.crossovers)]+1):nb.snps]

  stopifnot(sum(is.na(gam.chr)) == 0)

  return(gam.chr)
}

##' Gamete
##'
##' Make a gamete for all chromosomes by recombining two parental haplotypes.
##' @param haplos.par list of matrices containing both haplotypes of a parent for each chromosome
##' @param loc.crossovers list of vectors with positions of the crossing overs (the coordinate of the first nucleotide is assumed to be 1, and crossing-overs are assumed to occur right after the given localization)
##' @return list of vectors (one per chromosome)
##' @author Timothee Flutre
##' @export
makeGameteSingleInd <- function(haplos.par, loc.crossovers){
  stopifnot(is.list(haplos.par),
            is.list(loc.crossovers))

  gam <- Map(makeGameteSingleIndSingleChrom, haplos.par, loc.crossovers)

  return(gam)
}

##' Fecundation
##'
##' Perform fecundation by uniting the two gametes.
##' @param gam1 list of 1-row matrices (one per chromosome)
##' @param gam2 list of 1-row matrices (one per chromosome)
##' @param child.name identifier of the child
##' @return list of 2-row matrices (one per chromosome)
##' @author Timothee Flutre
##' @export
fecundation <- function(gam1, gam2, child.name){
  stopifnot(is.list(gam1),
            is.list(gam2),
            length(gam1) == length(gam2))

  nb.chroms <- length(gam1)

  child <- lapply(1:nb.chroms, function(c){
    tmp <- rbind(gam1[[c]], gam2[[c]])
    rownames(tmp) <- c(paste0(child.name, "_h1"), paste0(child.name, "_h2"))
    tmp
  })
  names(child) <- names(gam1)

  return(child)
}

##' Cross
##'
##' Make a cross. If two different individuals are given, then a fecundation is made with one gamete from each individual. If the same individual is given twice, an autofecondation is made with two different gametes from this individual. If a single individual is given, a haplodiploidization is made with a single gamete from this individual.
##' @param haplos.par1 list of matrices (one per chromosome) for the first parent
##' @param loc.crossovers.par1 list of vectors (one per chromosome) for the first parent
##' @param haplos.par2 list of matrices (one per chromosome) for the second parent
##' @param loc.crossovers.par2 list of vectors (one per chromosome) for the second parent
##' @param child.name identifier of the child (if NULL, <parent1>"-x-"<parent2> or <parent1>"-hd")
##' @param verbose verbosity level (0/1)
##' @return list of matrices (one per chromosome) for the child
##' @author Timothee Flutre
##' @export
makeCross <- function(haplos.par1,
                      loc.crossovers.par1,
                      haplos.par2=NULL,
                      loc.crossovers.par2=NULL,
                      child.name=NULL,
                      verbose=1){
  stopifnot(is.list(haplos.par1),
            is.list(loc.crossovers.par1),
            all(names(haplos.par1) == names(loc.crossovers.par1)))
  if(! is.null(haplos.par2))
    stopifnot(is.list(haplos.par2),
              is.list(loc.crossovers.par2),
              all(names(haplos.par2) == names(loc.crossovers.par2)))

  haplos.child <- NULL

  ## find which type of cross is wanted
  cross.type <- NULL
  name.par1 <- getIndNamesFromHaplos(haplos.par1)
  stopifnot(length(name.par1) == 1)
  if(! is.null(haplos.par2)){
    name.par2 <- getIndNamesFromHaplos(haplos.par2)
    stopifnot(length(name.par2) == 1)
    if(name.par2 == name.par1){
      cross.type <- "autofecondation"
      if(verbose > 0){
        msg <- paste0("make an autofecondation of '", name.par1, "'")
        write(msg, stdout())
      }
    } else{
      cross.type <- "fecondation"
      if(verbose > 0){
        msg <- paste0("make a fecondation of '", name.par1, "' and '",
                      name.par2, "'")
        write(msg, stdout())
      }
    }
  } else{
    cross.type <- "haplodiploidization"
    if(verbose > 0){
      msg <- paste0("make a haplodiploidization of '", name.par1, "'")
      write(msg, stdout())
      }
  }

  ## make the gamete(s) and then the cross
  gam.par1 <- makeGameteSingleInd(haplos.par1, loc.crossovers.par1)
  if(cross.type == "haplodiploidization"){
    if(is.null(child.name))
      child.name <- paste0(name.par1, "-hd")
    haplos.child <- fecundation(gam.par1, gam.par1, child.name)
  } else{
    gam.par2 <- makeGameteSingleInd(haplos.par2, loc.crossovers.par2)
    if(is.null(child.name))
      child.name <- paste0(name.par1, "-x-", name.par2)
    haplos.child <- fecundation(gam.par1, gam.par2, child.name)
  }

  return(haplos.child)
}

##' Crossing-overs
##'
##' Draw the number and location of crossing-overs per gamete.
##' @param crosses data.frame with three columns, parent1, parent2, child; if parent 1 and 2 are the same, it will be an autofecondation; if parent2 is NA, it will be a haplodiploidization
##' @param nb.snps vector with the nb of SNPs per chromosome, which names are chromosome names
##' @param lambda mean number of crossing-overs (parameter of a Poisson)
##' @return list of lists (one per cross, then one per parent, then one per chromosome) whose names are crosses$child, in the same order
##' @author Timothee Flutre
##' @export
drawLocCrossovers <- function(crosses, nb.snps, lambda=2){
  stopifnot(is.data.frame(crosses),
            ncol(crosses) >= 3,
            all(c("parent1", "parent2", "child") %in% colnames(crosses)),
            sum(is.na(crosses$parent1)) == 0,
            sum(is.na(crosses$child)) == 0,
            is.vector(nb.snps),
            ! is.null(names(nb.snps)))

  loc.crossovers <- list()

  ## draw the number of crossing-overs for each gamete
  parent.names <- c(crosses$parent1, crosses$parent2)
  parent.names <- parent.names[! is.na(parent.names)]
  nb.gametes <- length(parent.names)
  nb.chroms <- length(nb.snps)
  nb.crossovers <- stats::rpois(nb.gametes * nb.chroms, lambda)
  nb.crossovers[nb.crossovers == 0] <- 1 # at least 1 per chromosome

  ## draw the location of each crossing-over
  nb.crosses <- nrow(crosses)
  chrom.names <- names(nb.snps)
  cross.idx <- 0
  gam.idx <- 0
  while(cross.idx < nb.crosses){
    cross.idx <- cross.idx + 1
    child <- list()
    gam.idx <- gam.idx + 1
    gam.par1 <- lapply(1:nb.chroms, function(c){
      sort(sample.int(n=nb.snps[c] - 1,
                      size=nb.crossovers[(gam.idx-1)*nb.chroms + c]))
    })
    names(gam.par1) <- chrom.names
    child[[crosses$parent1[cross.idx]]] <- gam.par1
    if(! is.na(crosses$parent2[cross.idx])){ # (auto)fecondation
      gam.idx <- gam.idx + 1
      gam.par2 <- lapply(1:nb.chroms, function(c){
        sort(sample.int(n=nb.snps[c] - 1,
                        size=nb.crossovers[(gam.idx-1)*nb.chroms + c]))
      })
      names(gam.par2) <- chrom.names
      child[[crosses$parent2[cross.idx]]] <- gam.par2
    }
    loc.crossovers[[crosses$child[cross.idx]]] <- child
  }

  return(loc.crossovers)
}

##' Crosses
##'
##' Make crosses (fecundation, autofecondation, haplodiploidization).
##' @param haplos list of matrices (one per chromosome)
##' @param crosses data.frame with three columns, parent1, parent2, child (no duplicate); if parent 1 and 2 are the same, it will be an autofecondation; if parent2 is NA, it will be a haplodiploidization
##' @param loc.crossovers list of lists (one per cross, then one per parent, then one per chromosome) whose names are crosses$child, in the same order; if NULL, draw many crossing-overs localizations at once (as Poisson with parameter 2, assuming all chromosomes have the same length)
##' @param nb.cores the number of cores to use, i.e. at most how many child processes will be run simultaneously (not on Windows)
##' @param verbose verbosity level (0/1)
##' @return list of matrices (one per chromosome) with child haplotypes in rows and SNPs in columns
##' @author Timothee Flutre
##' @export
makeCrosses <- function(haplos, crosses, loc.crossovers=NULL,
                        nb.cores=1, verbose=1){
  stopifnot(is.list(haplos),
            length(unique(sapply(haplos, class))) == 1,
            unique(sapply(haplos, class)) == "matrix",
            length(unique(sapply(haplos, nrow))) == 1, # same nb of individuals
            is.data.frame(crosses),
            ncol(crosses) >= 3,
            all(c("parent1", "parent2", "child") %in% colnames(crosses)),
            sum(is.na(crosses$parent1)) == 0,
            sum(is.na(crosses$child)) == 0,
            anyDuplicated(crosses$child) == 0)
  haplos.ind.names <- getIndNamesFromHaplos(haplos)
  parent.names <- c(crosses$parent1, crosses$parent2)
  parent.names <- parent.names[! is.na(parent.names)]
  stopifnot(all(parent.names %in% haplos.ind.names))
  if(! is.null(loc.crossovers))
    stopifnot(is.list(loc.crossovers),
              all(names(loc.crossovers) == crosses$child))

  if(Sys.info()["sysname"] == "Windows")
    nb.cores <- 1

  nb.crosses <- nrow(crosses)
  nb.chroms <- length(haplos)
  chrom.names <- names(haplos)

  if(is.null(loc.crossovers)){
    if(verbose > 0){
      msg <- "draw locations of crossing-overs ..."
      write(msg, stdout())
    }
    nb.snps <- sapply(haplos, ncol) # per chromosome
    loc.crossovers <- drawLocCrossovers(crosses, nb.snps)
  }

  if(verbose > 0){
    msg <- "make all crosses ..."
    write(msg, stdout())
  }
  tmp <- parallel::mclapply(1:nb.crosses, function(i){
    if(is.na(crosses$parent2[i])){ # haplodiploidization
      makeCross(haplos.par1=getHaplosInd(haplos, crosses$parent1[i]),
                loc.crossovers.par1=loc.crossovers[[i]][[crosses$parent1[i]]],
                child.name=crosses$child[i],
                verbose=verbose-1)
    } else #(auto)fecondation
      makeCross(haplos.par1=getHaplosInd(haplos, crosses$parent1[i]),
                loc.crossovers.par1=loc.crossovers[[i]][[crosses$parent1[i]]],
                haplos.par2=getHaplosInd(haplos, crosses$parent2[i]),
                loc.crossovers.par2=loc.crossovers[[i]][[crosses$parent2[i]]],
                child.name=crosses$child[i],
                verbose=verbose-1)
  }, mc.cores=nb.cores)

  if(verbose > 0){
    msg <- "gather all children's haplotypes per chromosome ..."
    write(msg, stdout())
  }
  haplos.children <- lapply(1:nb.chroms, function(c){
    do.call(rbind, parallel::mclapply(tmp, `[[`, c, mc.cores=nb.cores))
  })
  names(haplos.children) <- chrom.names

  return(haplos.children)
}

.isValidSnpCoords <- function(snp.coords){
  all(is.data.frame(snp.coords),
      ! is.null(rownames(snp.coords)),
      ncol(snp.coords) >= 2,
      "chr" %in% colnames(snp.coords),
      "coord" %in% colnames(snp.coords) | "pos" %in% colnames(snp.coords))
}

.df2gr <- function(snp.coords.df){
  requireNamespaces(c("GenomicRanges", "S4Vectors", "IRanges"))
  stopifnot(.isValidSnpCoords(snp.coords.df))
  snp.coords.gr <-
    GenomicRanges::GRanges(seqnames=S4Vectors::Rle(snp.coords.df$chr),
                           ranges=IRanges::IRanges(start=snp.coords.df$coord,
                                                   end=snp.coords.df$coord))
  names(snp.coords.gr) <- rownames(snp.coords.df)
  return(snp.coords.gr)
}

##' Distance between SNP pairs
##'
##' For each SNP pair, return the number of "blocks" (i.e. nucleotides) between both SNPs via the \code{\link[GenomicRanges]{distance}} function.
##' @param snp.pairs data.frame with two columns "loc1" and "loc2"
##' @param snp.coords data.frame with SNP identifiers as row names, and two columns, "chr" and "coord" or "pos"
##' @param nb.cores the number of cores to use
##' @param verbose verbosity level (0/1)
##' @return vector
##' @author Timothee Flutre
##' @export
distSnpPairs <- function(snp.pairs, snp.coords, nb.cores=1, verbose=1){
  requireNamespaces(c("GenomicRanges", "S4Vectors", "IRanges"))
  stopifnot(is.data.frame(snp.pairs),
            ncol(snp.pairs) >= 2,
            all(c("loc1", "loc2") %in% colnames(snp.pairs)),
            .isValidSnpCoords(snp.coords))
  snp.pairs$loc1 <- as.character(snp.pairs$loc1)
  snp.pairs$loc2 <- as.character(snp.pairs$loc2)
  stopifnot(all(unique(unlist(snp.pairs[, c("loc1", "loc2")])) %in%
                rownames(snp.coords)))
  if(! "coord" %in% colnames(snp.coords))
    colnames(snp.coords)[colnames(snp.coords) == "pos"] <- "coord"

  if(verbose > 0)
    message("make GRanges ...")
  snp.granges <- .df2gr(snp.coords)

  if(verbose > 0)
    message("calculate pairwise distances ...")
  dist.loc <- GenomicRanges::distance(x=snp.granges[snp.pairs$loc1],
                                      y=snp.granges[snp.pairs$loc2])

  return(dist.loc)
}

##' Genomic relatedness
##'
##' Estimate genetic relationships between individuals from their SNP genotypes.
##' Note that "relationships" are estimated, and not "coancestries" which are equal to 2 times "relationhips".
##' @param X matrix of SNP genotypes encoded in number of copies of the 2nd allele, i.e. as allele doses in {0,1,2}, with individuals in rows and SNPs in columns; missing values should be encoded as NA in order to be ignored
##' @param afs vector with 2nd allele frequencies (calculated with \code{\link{estimAf}} if NULL)
##' @param thresh threshold on minor allele frequencies below which SNPs are ignored (e.g. 0.01; NULL to skip this step)
##' @param relationships relationship to estimate (additive/dominant/gaussian) where "gaussian" corresponds to the Gaussian kernel from Endelman (2011)
##' @param method if additive relationships, can be "vanraden1" (first method in VanRaden, 2008), "habier" (similar to 'vanraden1' without giving more importance to rare alleles; from Habier et al, 2007), "astle-balding" (two times equation 2.2 in Astle & Balding, 2009), "yang" (similar to 'astle-balding' but without ignoring sampling error per SNP; from Yang et al, 2010), "zhou" (centering the genotypes and not assuming that rare variants have larger effects; from Zhou et al, 2013) or "center-std"; if dominant relationships, can be "vitezica" (classical/statistical parametrization from Vitezica et al, 2013) or "su" (from Su et al, PLoS One, 2012)
##' @param theta smoothing parameter for "gauss"
##' @param verbose verbosity level (0/1)
##' @return matrix
##' @author Timothee Flutre
##' @export
estimGenRel <- function(X, afs=NULL, thresh=NULL, relationships="additive",
                        method="vanraden1", theta=0.5, verbose=1){
  stopifnot(.isValidGenosDose(X),
            relationships %in% c("additive", "dominant", "gaussian"))
  if(relationships == "additive")
    stopifnot(method %in% c("vanraden1", "habier", "astle-balding", "yang",
                            "zhou", "center-std"))
  if(relationships == "dominant")
    stopifnot(method %in% c("vitezica", "su"))
  if(! is.null(thresh))
    stopifnot(thresh >= 0, thresh <= 0.5)
  if(relationships == "gauss"){
    stopifnot(! is.null(theta),
              is.numeric(theta),
              length(theta) == 1,
              theta > 0,
              theta <= 1)
  }
  if(! is.null(afs)){
    stopifnot(! is.null(names(afs)),
              all(colnames(X) %in% names(afs)),
              all(names(afs) %in% colnames(X)))
    afs <- afs[colnames(X)] # re-order if necessary
  }

  gen.rel <- NULL # to be filled and returned

  N <- nrow(X) # nb of individuals
  P <- ncol(X) # nb of SNPs

  idx.rm <- c()

  ## discard SNPs with missing data
  snps.na <- apply(X, 2, function(x){
    any(is.na(x))
  })
  if(any(snps.na)){
    if(verbose > 0){
      txt <- paste0("skip ", sum(snps.na), " SNPs with missing data")
      write(txt, stdout())
    }
    idx.rm <- which(snps.na)
    X <- X[, -idx.rm, drop=FALSE]
    P <- ncol(X)
    if(! is.null(afs))
      afs <- afs[-idx.rm]
  }

  ## estimate AFs
  if(is.null(afs))
    afs <- estimAf(X=X)

  ## discard SNPs with low MAFs
  if(! is.null(thresh)){
    mafs <- estimMaf(afs=afs)
    X <- discardSnpsLowMaf(X=X, mafs=mafs, thresh=thresh, verbose=verbose)
    P <- ncol(X)
    afs <- afs[colnames(X)]
  }

  ## estimate relationships (not coancestries)
  if(verbose > 0){
    txt <- paste0("use ", ncol(X), " SNPs")
    write(txt, stdout())
  }
  if(relationships == "additive"){
    if(method == "vanraden1"){
      ## implementation as in VanRaden (2008)
      ## M <- X - 1 # recode genotypes as {-1,0,1}
      ## Pmat <- matrix(rep(1, N)) %*% (2 * (afs - 0.5))
      ## Z <- M - Pmat

      ## implementation as in Vitezica et al (2013)
      tmp <- matrix(rep(1, N)) %*% (2 * afs)
      Z <- X - tmp

      gen.rel <- tcrossprod(Z, Z) / (2 * sum(afs * (1 - afs)))
    } else if(method == "habier"){
      gen.rel <- tcrossprod(X, X) / (2 * sum(afs * (1 - afs)))
    } else if(method == "astle-balding"){
      tmp1 <- sweep(x=X, MARGIN=2, STATS=2*afs, FUN="-")
      tmp2 <- sweep(x=tmp1, MARGIN=2, STATS=2*afs*(1-afs), FUN="/")
      gen.rel <- (1/P) * tcrossprod(tmp1, tmp2)
    } else if(method == "yang"){
      tmp1 <- sweep(x=X, MARGIN=2, STATS=2*afs, FUN="-")
      tmp2 <- sweep(x=tmp1, MARGIN=2, STATS=2*afs*(1-afs), FUN="/")
      gen.rel <- (1/P) * tcrossprod(tmp1, tmp2)
      tmp3 <- X^2 - sweep(x=X, MARGIN=2, STATS=1+2*afs, FUN="*")
      tmp4 <- sweep(x=tmp3, MARGIN=2, STATS=2*afs^2, FUN="+")
      tmp5 <- sweep(x=tmp4, MARGIN=2, STATS=2*afs*(1-afs), FUN="/")
      diag(gen.rel) <- 1 + (1/P) * rowSums(tmp5)
    } else if(method == "zhou"){
      tmp <- scale(x=X, center=TRUE, scale=FALSE)
      gen.rel <- tcrossprod(tmp, tmp) / P
    } else if(method == "center-std"){
      tmp <- scale(x=X, center=TRUE, scale=TRUE)
      gen.rel <- tcrossprod(tmp, tmp) / P
    }
  } else if(relationships == "dominant"){
    if(method == "vitezica"){
      ## caution, compared to Vitezica et al (2013), the X matrix encodes
      ## genotypes in terms of nb of copies of the A2 allele, not A1
      is.0 <- (X == 0) # homozygotes for the first allele
      is.1 <- (X == 1) # heterozygotes
      is.2 <- (X == 2) # homozygotes for the second allele
      W <- matrix(NA, nrow=nrow(X), ncol=ncol(X), dimnames=dimnames(X))
      for(i in 1:N){
        W[i, is.0[i,]] <- - 2 * afs[is.0[i,]]^2
        W[i, is.1[i,]] <- 2 * afs[is.1[i,]] * (1 - afs[is.1[i,]])
        W[i, is.2[i,]] <- - 2 * (1 - afs[is.2[i,]])^2
      }
      gen.rel <- tcrossprod(W, W) / sum((2 * afs * (1 - afs))^2)
    } else if(method == "su"){
      H <- X
      H[H != 1] <- 0 # recode genotypes as {0,1,0}
      H <- sweep(x=H, MARGIN=2, STATS=2*afs*(1-afs), FUN="-")
      gen.rel <- tcrossprod(H, H) /
        (2 * sum(afs * (1 - afs) * (1 - 2 * afs * (1 - afs))))
    }
  } else if(relationships == "gaussian"){
    M <- X - 1 # recode genotypes as {-1,0,1}
    gen.dist <- as.matrix(stats::dist(x=M, method="euclidean")) / (2 * sqrt(P))
    gen.rel <- exp(-(gen.dist / theta)^2)
  }

  return(gen.rel)
}

##' Pairwise linkage disequilibrium
##'
##' Estimates linkage disequilibrium between pairs of SNPs when the observations are the genotypes of individuals, not their gametes (i.e. the gametic phases are unknown).
##' When ignoring kinship and population structure, the estimator of Rogers and Huff (Genetics, 2009) can be used.
##' When kinship and/or population structure are controlled for, the estimator of Mangin et al (Heredity, 2012) is used via their LDcorSV package.
##' @param X matrix of SNP genotypes encoded in number of copies of the 2nd allele, i.e. as allele doses in {0,1,2}, with individuals in rows and SNPs in columns
##' @param K matrix of kinship
##' @param pops vector of characters indicating the population of each individual
##' @param snp.coords data.frame with SNP identifiers as row names, and two columns, "chr" and "pos"
##' @param only.chr identifier of a given chromosome
##' @param only.pop identifier of a given population
##' @param use.ldcorsv required if K and/or pops are not NULL; otherwise use the square of \code{\link{cor}}
##' @param verbose verbosity level (0/1)
##' @return data frame
##' @author Timothee Flutre
##' @export
estimLd <- function(X, K=NULL, pops=NULL, snp.coords,
                    only.chr=NULL, only.pop=NULL,
                    use.ldcorsv=FALSE, verbose=1){
  if(use.ldcorsv & ! requireNamespace("LDcorSV", quietly=TRUE))
    stop("Pkg 'LDcorSV' needed for this function to work.",
         call.=FALSE)
  stopifnot(.isValidGenosDose(X),
            .isValidSnpCoords(snp.coords))
  if(! is.null(K))
    stopifnot(use.ldcorsv,
              is.matrix(K),
              nrow(K) == ncol(K),
              nrow(K) == nrow(X),
              ! is.null(dimnames(K)),
              all(rownames(K) == colnames(K)),
              all(rownames(K) == rownames(X)))
  W.s <- NA
  if(! is.null(pops)){
    stopifnot(use.ldcorsv,
              length(pops) == nrow(X),
              ! is.null(names(pops)),
              names(pops) == rownames(X))
    W.s <- stats::model.matrix(~ as.factor(pops))[, -1]
    rownames(W.s) <- names(pops)
  }
  if(! is.null(only.chr))
    if(! only.chr %in% snp.coords$chr)
      stop(paste0("chr '", only.chr, "' absent from snp.coords"))
  if(! is.null(only.pop))
    if(! only.pop %in% pops)
      stop(paste0("pop '", only.pop, "' absent from pops"))

  ld <- NULL

  subset.snps <- 1:ncol(X)
  subset.inds <- 1:nrow(X)
  if(! is.null(only.chr) | ! is.null(only.pop)){
    if(verbose > 0)
      write("extract relevant individuals and SNPs ...", stdout())
    if(! is.null(only.chr))
      subset.snps <- which(snp.coords$chr == only.chr)
    if(! is.null(only.pop))
      subset.inds <- which(pops == only.pop)
  }
  X <- X[subset.inds, subset.snps]
  if(! is.null(K))
    K <- K[subset.inds, subset.inds]

  if(verbose > 0)
    write("estimate pairwise LD ...", stdout())
  if(is.null(K)){
    if(is.null(pops)){
      if(use.ldcorsv){
        ld <- LDcorSV::LD.Measures(donnees=X,
                                   V=NA,
                                   S=NA,
                                   data="G", supinfo=FALSE, na.presence=FALSE)
      } else{
        tmp <- stats::cor(X)^2
        tmp[upper.tri(tmp)] <- NA
        diag(tmp) <- NA
        ld <- data.frame(t(utils::combn(colnames(X), 2)),
                         tmp[! is.na(tmp)],
                         stringsAsFactors=TRUE)
        colnames(ld) <- c("loc1", "loc2", "cor2")
      }
    } else{ # if(is.null(K) & ! is.null(pops))
      if(! is.null(only.pop)){
        ld <- LDcorSV::LD.Measures(donnees=X,
                                   V=NA,
                                   S=NA,
                                   data="G", supinfo=FALSE, na.presence=FALSE)
      } else{
        ld <- LDcorSV::LD.Measures(donnees=X,
                                   V=NA,
                                   S=W.s,
                                   data="G", supinfo=FALSE, na.presence=FALSE)
      }
    }
  } else{ # if(! is.null(K))
    if(is.null(pops)){
      ld <- LDcorSV::LD.Measures(donnees=X,
                                 V=K,
                                 S=NA,
                                 data="G", supinfo=FALSE, na.presence=FALSE)
    } else{ # if(! is.null(K) & ! is.null(pops))
      if(! is.null(only.pop)){
        ld <- LDcorSV::LD.Measures(donnees=X,
                                   V=K,
                                   S=NA,
                                   data="G", supinfo=FALSE, na.presence=FALSE)
      } else{
        ld <- LDcorSV::LD.Measures(donnees=X,
                                   V=K,
                                   S=W.s,
                                   data="G", supinfo=FALSE, na.presence=FALSE)
      }
    }
  }

  return(ld)
}

##' Pairwise linkage disequilibrium
##'
##' Plots the linkage disequilibrium between pairs of SNPs, as a blue density or black points, with a red loess.
##' Possibility to add two analytical approximations of E[r^2] at equilibrium (see McVean, Handbook of Stat Gen, 2007): 1 / (1 + 4 Ne c x) by Sved (1971) and (10 + 4 Ne c x) / (22 + 13 * 4 Ne c x + (4 Ne c x)^2) by Ohta and Kimura (1971).
##' @param x vector of distances between SNPs (see \code{\link{distSnpPairs}})
##' @param y vector of LD estimates (see \code{\link{estimLd}})
##' @param estim estimator of pairwise LD corresponding to the values in y (r2/r)
##' @param main main title
##' @param use.density if TRUE, uses smoothScatter; otherwise, use scatter.smooth
##' @param xlab label for the x axis
##' @param ylab label for the y axis
##' @param span the parameter alpha which controls the degree of smoothing (see \code{\link[stats]{loess}})
##' @param degree the degree of the polynomials to be used (see \code{\link[stats]{loess}})
##' @param evaluation number of points at which to evaluate the smooth curve (see \code{\link[stats]{loess.smooth}})
##' @param sample.size nb of sampled haplotypes, n, used to estimate the pairwise LD; if not NULL, n is used to plot the horizontal line at the r value above which the null hypothesis "D=0" is rejected, where r = D / sqrt(f_A f_a f_B f_b) and X^2 = n r^2 is the test statistic asymptotically following a Chi2(df=1) (see McVean, Handbook of Statistical Genetics, 2007, and Pritchard and Przewoski, AJHG, 2001)
##' @param add.ohta.kimura add the analytical approximation by Ohta and Kimura (1971); requires Ne and c
##' @param add.sved add the analytical approximation by Sved (1971); requires Ne and c
##' @param Ne effective population size
##' @param c recomb rate in events per base per generation
##' @return nothing
##' @author Timothee Flutre
##' @export
plotLd <- function(x, y, estim="r2", main,
                   use.density=TRUE,
                   xlab="Physical distance (bp)",
                   ylab=paste0("Linkage disequilibrium (", estim, ")"),
                   span=1/10, degree=1, evaluation=50,
                   sample.size=NULL,
                   add.ohta.kimura=TRUE, add.sved=TRUE, Ne=NULL, c=NULL){
  stopifnot(is.vector(x),
            is.vector(y),
            estim %in% c("r2","r"),
            is.character(main))
  if(! is.null(sample.size))
    stopifnot(is.numeric(sample.size))
  if(any(c(add.ohta.kimura, add.sved)))
    stopifnot(! is.null(Ne),
              ! is.null(c))

  ## plot the pairwise estimates of LD
  lpars <- list(col="red", cex=2)
  if(use.density){
    graphics::smoothScatter(x, y,
                  main=main,
                  xlab=xlab,
                  ylab=ylab,
                  las=1)
    pred <- stats::loess.smooth(x, y, span=span, degree=degree, evaluation=evaluation)
    do.call(graphics::lines, c(list(pred), lpars))
  } else{
    stats::scatter.smooth(x, y, lpars=lpars,
                   main=main, xlab=xlab, ylab=ylab, las=1,
                   span=span, degree=degree, evaluation=evaluation)
  }

  ## add the "significance" horizontal line
  ## reject H0:"D=0" at 5% if X2 = n x hat(r^2) >= Chi2(1)
  if(! is.null(sample.size)){
    X2 <- stats::qchisq(p=0.05, df=1, lower.tail=FALSE)
    tmp <- X2 / sample.size
    if(estim == "r")
      tmp <- sqrt(tmp)
    graphics::abline(h=tmp,
           col=ifelse(use.density, "black", "blue"),
           lty=2, lwd=2)
  }

  ## add analytical approximations
  if(any(c(add.ohta.kimura, add.sved)))
    scaled.dist <- 4 * Ne * c * x
  if(add.ohta.kimura){
    ok <- (10 + scaled.dist) / (22 + 13 * scaled.dist + scaled.dist^2)
    if(estim == "r")
      ok <- sqrt(ok)
    graphics::points(x, ok, pch=".", col="purple", cex=1.2)
  }
  if(add.sved){
    sved <- 1 / (1 + scaled.dist)
    if(estim == "r")
      sved <- sqrt(sved)
    graphics::points(x, sved, pch=".", col="green", cex=1.2)
  }

  ## add the legend
  legs <- "loess"
  cols <- "red"
  ltys <- 1
  lwds <- c(2)
  if(! is.null(sample.size)){
    legs <- c(legs,
              as.expression(bquote(paste("r | D=0, n=", .(sample.size),
                                         ", ", alpha, "=5%"))))
    cols <- c(cols, ifelse(use.density, "black", "blue"))
    ltys <- c(ltys, 2)
    lwds <- c(lwds, 2)
  }
  if(add.ohta.kimura){
    legs <- c(legs, "Ohta & Kimura (1971)")
    cols <- c(cols, "purple")
    ltys <- c(ltys, 1)
    lwds <- c(lwds, 2)
  }
  if(add.sved){
    legs <- c(legs, "Sved (1971)")
    cols <- c(cols, "green")
    ltys <- c(ltys, 1)
    lwds <- c(lwds, 2)
  }
  graphics::legend("topright", legend=legs, col=cols, lty=ltys, lwd=lwds, bty="n")
}

##' Distance between consecutive SNPs
##'
##' For each pair of consecutive SNPs, return the number of "blocks" (i.e. nucleotides) between both SNPs (to be coherent with \code{\link{distSnpPairs}}).
##' @param snp.coords data.frame with SNP identifiers as row names, and with two columns "chr" and "coord" or "pos"
##' @param only.chr identifier of a given chromosome
##' @param nb.cores the number of cores to use
##' @return list with one component per chromosome
##' @author Timothee Flutre
##' @export
distConsecutiveSnps <- function(snp.coords, only.chr=NULL, nb.cores=1){
  requireNamespaces("parallel")
  stopifnot(.isValidSnpCoords(snp.coords))
  if(! "coord" %in% colnames(snp.coords))
    colnames(snp.coords)[colnames(snp.coords) == "pos"] <- "coord"

  snp.coords$chr <- as.character(snp.coords$chr)
  chr.ids <- unique(snp.coords$chr)
  if(! is.null(only.chr)){
    stopifnot(only.chr %in% chr.ids)
    chr.ids <- only.chr
  }

  snp.dists <- parallel::mclapply(chr.ids, function(chr.id){
    coords <- snp.coords$coord[snp.coords$chr == chr.id]
    names(coords) <- rownames(snp.coords)[snp.coords$chr == chr.id]
    dis <- coords[2:length(coords)] - coords[1:(length(coords)-1)] - 1
    names(dis) <- paste(names(dis), names(coords)[-length(coords)], sep="-")
    dis
  }, mc.cores=nb.cores)
  names(snp.dists) <- chr.ids

  return(snp.dists)
}

##' Thin SNPs
##'
##' Thin SNPs according to various methods: based on their index or based on their genomic coordinate.
##' @param method index or coord
##' @param threshold keep every "threshold" SNPs (if method="index"), keep SNPs with more than "threshold" base pairs between them (if method="coord")
##' @param snp.coords data.frame with SNP identifiers as row names, and two columns, "chr" and "coord" or "pos"; SNPs will be sorted according to their coordinates per chromosome (use \code{gtools::mixedsort} if you want to also sort chromosomes)
##' @param only.chr identifier of a given chromosome
##' @return vector of SNP identifiers
##' @author Timothee Flutre
##' @export
thinSnps <- function(method, threshold, snp.coords, only.chr=NULL){
  requireNamespaces("GenomicRanges")
  stopifnot(method %in% c("index", "coord"))
  stopifnot(! is.null(snp.coords),
            .isValidSnpCoords(snp.coords),
            threshold == floor(threshold))
  if(! "coord" %in% colnames(snp.coords))
    colnames(snp.coords)[colnames(snp.coords) == "pos"] <- "coord"

  snp.coords$chr <- as.character(snp.coords$chr)
  chr.ids <- unique(snp.coords$chr)
  if(! is.null(only.chr)){
    stopifnot(only.chr %in% chr.ids)
    chr.ids <- only.chr
  }

  out.snp.ids <- do.call(c, lapply(chr.ids, function(chr.id){
    tmp <- snp.coords[snp.coords$chr == chr.id,]
    tmp <- tmp[order(tmp$coord),]
    if(method == "index"){
      idx <- seq(1, nrow(tmp), threshold)
      rownames(tmp)[idx]
    } else if(method == "coord"){
      tiles <- GenomicRanges::tileGenome(seqlengths=stats::setNames(max(tmp$coord),
                                                             chr.id),
                                         tilewidth=threshold)
      tmp.gr <- .df2gr(tmp)
      ovl <- GenomicRanges::findOverlaps(tiles, tmp.gr)
      idx <- sapply(as.list(ovl), `[`, 1)
      names(tmp.gr[idx[! is.na(idx)]])
    }
  }))

  return(out.snp.ids)
}

##' Convert genotypes
##'
##'
##' @param X matrix of SNP genotypes encoded in number of copies of the 2nd allele, i.e. as allele doses in {0,1,2}, with individuals in rows and SNPs in columns
##' @param snp.coords data.frame with SNP identifiers as row names, and two columns, "chr" and "coord" or "pos"
##' @param alleles data.frame with SNPs in rows (names as row names) and alleles in columns (first is "minor", second is "major")
##' @param file write the genotype data to this file (for instance 'genotypes_fastphase.txt' or 'genotypes_fastphase.txt.gz', but don't use \code{\link[base]{gzfile}})
##' @return nothing
##' @author Timothee Flutre
##' @export
writeInputsForFastPhase <- function(X, snp.coords, alleles, file){
  stopifnot(file.exists(Sys.which("fcgene")),
            .isValidGenosDose(X, check.na=FALSE),
            .isValidSnpCoords(snp.coords))

  ## save SNP genotypes as "r-formatted" file for fcGENE
  tmp.geno.file <- tempfile(pattern="genotypes_", fileext=".txt")
  sep <- "\t"
  txt <- paste0("SAMPLE_ID", sep, paste0(colnames(X), collapse=sep), "\n")
  cat(txt, file=tmp.geno.file, sep="")
  tmp <- do.call(rbind, lapply(1:nrow(X), function(i){
    paste0(rownames(X)[i], sep, paste0(X[i,], collapse=sep), "\n")
  }))
  cat(tmp, file=tmp.geno.file, sep="", append=TRUE)

  ## save "snpinfo"
  tmp.snpinfo.file <- tempfile(pattern="snpinfo_", fileext=".txt")
  snp.info <- data.frame(snpid=rownames(snp.coords),
                         rsid=rownames(snp.coords),
                         position=snp.coords$pos,
                         allele1=alleles$minor,
                         allele2=alleles$major,
                         stringsAsFactors=FALSE)
  utils::write.table(x=snp.info, file=tmp.snpinfo.file, quote=FALSE, sep="\t",
                     row.names=FALSE)

  ## convert with fcGENE
  tmp.out.prefix <- tempfile(pattern="output")
  args <- paste0("--rgeno ", tmp.geno.file,
                 " --snpinfo ", tmp.snpinfo.file,
                 " --oformat fastphase",
                 " --out ", tmp.out.prefix)
  retVal <- system2(command="fcgene", args=args, stdout=TRUE, stderr=TRUE)

  ## rename (and compress if necessary)
  use.gzip <- grepl(".gz", file)
  if(use.gzip)
    file <- removeFileExtension(file, "gz")
  file.rename(from=paste0(tmp.out.prefix, ".inp"), to=file)
  if(use.gzip)
    system2(command="gzip", args=file)

  ## clean
  tmp <- file.remove(c(tmp.geno.file, tmp.snpinfo.file,
                       paste0(tmp.out.prefix, "_fcgene.log"),
                       paste0(tmp.out.prefix, "_pedinfo.txt"),
                       paste0(tmp.out.prefix, "_snpinfo.txt")))
}

##' Imputation
##'
##' Extract the haplotypes from the output file of fastPHASE.
##' @param file character
##' @param snp.ids vector of SNP identifiers
##' @return matrix with genotypes in rows (2 per individual) and SNPs in columns
##' @author Timothee Flutre
##' @export
readGenosFastPhase <- function(file, snp.ids){
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

##' Imputation
##'
##' Extract the haplotypes from the output file of fastPHASE.
##' @param file character
##' @param snp.ids vector of SNP identifiers
##' @return matrix with haplotypes in rows (2 per individual) and SNPs in columns
##' @author Timothee Flutre
##' @export
readHaplosFastPhase <- function(file, snp.ids){
  stopifnot(file.exists(file),
            is.vector(snp.ids))

  lines <- readLines(file, warn=FALSE)
  idx <- (which(grepl("BEGIN GENOTYPES", lines))+1):(which(grepl("END GENOTYPES", lines))-1)
  lines <- lines[idx]
  lines <- gsub(" ", "", lines)
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

##' Imputation
##'
##' Impute SNP genotypes via fastPHASE (Scheet and Stephens, 2006).
##' @param X NxP matrix of SNP genotypes encoded in number of copies of the 2nd allele, i.e. as allele doses in {0,1,2}, with individuals in rows and SNPs in columns
##' @param snp.coords data.frame with SNP identifiers as row names, and two columns, "chr" and "coord" or "pos"
##' @param alleles data.frame with SNPs in rows (names as row names) and alleles in columns (first is "minor", second is "major")
##' @param out.dir directory in which the output files will be saved
##' @param task.id identifier of the task (used in temporary and output file names)
##' @param nb.starts number of random starts of the EM algorithm
##' @param nb.iters number of iterations of the EM algorithm
##' @param nb.samp.haplos number of haplotypes sampled from the posterior distribution obtained from a particular random start of the EM algorithm
##' @param estim.haplos estimate haplotypes by minimizing individual error
##' @param nb.clusters number of clusters
##' @param seed seed for random number generation
##' @param clean remove files
##' @param verbose verbosity level (0/1)
##' @return list with matrix of genotypes before imputation, and matrix of haplotypes (2 per individual, in rows) after imputation
##' @author Timothee Flutre
##' @export
fastphase <- function(X, snp.coords, alleles, out.dir=getwd(),
                      task.id="fastphase", nb.starts=20, nb.iters=25,
                      nb.samp.haplos=50, estim.haplos=FALSE,
                      nb.clusters=10, seed=1859, clean=FALSE,
                      verbose=1){
  stopifnot(file.exists(Sys.which("fastPHASE")),
            .isValidGenosDose(X, check.na=FALSE))

  if(verbose > 0){
    txt <- "prepare input files..."
    write(txt, stdout())
  }
  tmp.files <- c()
  tmp.files <- c(tmp.files,
                 genos=paste0(out.dir, "/genos_", task.id, ".txt"))
  writeInputsForFastPhase(X, snp.coords, alleles, tmp.files["genos"])

  if(verbose > 0){
    txt <- "run fastPHASE..."
    write(txt, stdout())
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
  if(verbose > 0)
    message(paste("fastPHASE", args))
  system2(command="fastPHASE", args=args, stdout=TRUE, stderr=TRUE)

  if(verbose > 0){
    txt <- "load files..."
    write(txt, stdout())
  }
  tmp.files["haplos"] <- paste0(out.dir, "/", task.id, "_hapguess_switch.out")
  out <- list(genos=readGenosFastPhase(tmp.files["genos"],
                                       rownames(snp.coords)),
              haplos=readHaplosFastPhase(tmp.files["haplos"],
                                         rownames(snp.coords)))

  if(clean)
    for(f in tmp.files)
      file.remove(f)

  return(out)
}

##' Animal model
##'
##' Given T traits, I genotypes, Q covariates and N=I*Q phenotypes per trait, simulate phenotypes via the following "animal model": \eqn{Y = W C + Z G_A + Z G_D + E}, where Y is N x T; W is N x Q; Z is N x I; G_A ~ Normal_IxT(0, A, V_{G_A}) with A being IxI and V_{G_A} being TxT; G_D ~ Normal_IxT(0, D, V_{G_D}) with D being IxI and V_{G_D} being TxT; E ~ Normal_NxT(0, Id_N, V_E) with Id_N being NxN and V_E being TxT. Note that using the matrix Normal (MN) is equivalent to using the multivariate Normal (MVN) via the vec operator and Kronecker product: \eqn{Y ~ MN(M, U, V) <=> vec(Y) ~ MVN(vec(M), V \otimes U)}.
##' @param T number of traits
##' @param Q number of covariates (as "fixed effects", e.g. replicates)
##' @param mu T-vector of overall means (one per trait), i.e. C[1,1:T]
##' @param mean.C mean of the univariate Normal prior on C[2:Q,1:T] (ignored if Q=1)
##' @param sd.C std dev of the univariate Normal prior on C[2:Q,1:T] (ignored if Q=1)
##' @param A IxI matrix of additive genetic relationships between genotypes (see \code{\link{estimGenRel}} with VanRaden's estimator)
##' @param V.G.A TxT matrix of additive genetic variance-covariance between traits (e.g. 15 when T=1)
##' @param scale.hC.G.A scale of the half-Cauchy prior for sqrt{V_{G_A}} (e.g. 5; used if V.G.A=NULL and T=1)
##' @param nu.G.A degrees of freedom of the Wishart prior for V_{G_A} (used if V.G.A=NULL and T>1)
##' @param D IxI matrix of dominant genetic relationships between genotypes (see \code{\link{estimGenRel}} with Vitezica's estimator)
##' @param V.G.D TxT matrix of dominant genetic variance-covariance between traits (e.g. 3 when T=1; used if D!=NULL)
##' @param scale.hC.G.D scale of the half-Cauchy prior for sqrt{V_{G_D}} (e.g. 5; used if D!=NULL, V.G.D=NULL and T=1)
##' @param nu.G.D degrees of freedom of the Wishart prior for V_{G_D} (used if D!=NULL, V.G.D=NULL and T>1)
##' @param V.E TxT matrix of error variance-covariance between traits (used if T=1 and err.df=Inf)
##' @param scale.hC.E scale of the half-Cauchy prior for sqrt{V_E} (e.g. 5; used if V.E=NULL and T=1 and err.df=Inf)
##' @param nu.E degrees of freedom of the Wishart prior for V_E (used if V.E=NULL and T>1)
##' @param err.df degrees of freedom of the Student's t-distribution of the errors (e.g. 3; Inf means Normal distribution; will be Inf if T>1)
##' @param perc.NA percentage of missing phenotypes, at random
##' @param seed seed for the pseudo-random number generator
##' @return list
##' @author Timothee Flutre
##' @examples
##' ## simulate genotypes
##' set.seed(1859)
##' X <- simulGenosDose(nb.genos=200, nb.snps=2000)
##'
##' ## estimate the additive genetic relationships
##' A <- estimGenRel(X, relationships="additive", method="vanraden1", verbose=0)
##' A <- as.matrix(Matrix::nearPD(A)$mat) # not always necessary
##'
##' # one trait with heritability h2=0.75, no covariate, Normal errors, no NA
##' model <- simulAnimalModel(T=1, Q=1, A=A, V.G.A=15, V.E=5)
##'
##' # one trait with heritability h2=0.75, three covariates, Normal errors, no NA
##' model <- simulAnimalModel(T=1, Q=3, A=A, V.G.A=15, V.E=5)
##'
##' # one trait with heritability h2=0.75, no covariate, Student errors, no NA
##' model <- simulAnimalModel(T=1, Q=1, A=A, V.G.A=15, err.df=3)
##'
##' # one trait with heritability drawn at random, no covariate, Normal errors, no NA
##' model <- simulAnimalModel(T=1, Q=1, A=A, scale.hC.G.A=10, V.E=5)
##'
##' # two traits with heritabilities drawn at random, no covariate, Normal errors, no NA
##' model <- simulAnimalModel(T=2, Q=1, A=A, nu.G.A=5, nu.E=3)
##' @export
simulAnimalModel <- function(T=1,
                             Q=3, mu=rep(50,T), mean.C=5, sd.C=2,
                             A, V.G.A=NULL, scale.hC.G.A=NULL, nu.G.A=T,
                             D=NULL, V.G.D=NULL, scale.hC.G.D=NULL, nu.G.D=T,
                             V.E=NULL, scale.hC.E=NULL, nu.E=T,
                             err.df=Inf, perc.NA=0, seed=NULL){
  requireNamespaces("MASS")
  stopifnot(length(mu) == T,
            is.matrix(A),
            nrow(A) == ncol(A),
            ! is.null(rownames(A)),
            ! is.null(colnames(A)),
            rownames(A) == colnames(A))
  if(! is.null(V.G.A)){
    stopifnot(! is.matrix(V.G.A),
              nrow(V.G.A) == ncol(V.G.A))
  } else{
    if(T == 1){
      stopifnot(! is.null(scale.hC.G.A))
    } else{
      stopifnot(! is.null(nu.G.A))
    }
  }
  if(! is.null(D)){
    stopifnot(is.matrix(D),
              nrow(D) == ncol(D),
              ! is.null(rownames(D)),
              ! is.null(colnames(D)),
              rownames(D) == colnames(D),
              rownames(D) == rownames(A))
    if(! is.null(V.G.D)){
      stopifnot(! is.matrix(V.G.D),
                nrow(V.G.D) == ncol(V.G.D))
    } else{
      if(T == 1){
        stopifnot(! is.null(scale.hC.G.D))
      } else{
        stopifnot(! is.null(nu.G.D))
      }
    }
  }
  if(T == 1 & is.infinite(err.df) & is.null(V.E))
    stopifnot(! is.null(scale.hC.E))
  if(T > 1 & is.null(V.E))
    stopifnot(! is.null(nu.E))
  stopifnot(perc.NA >= 0, perc.NA <= 100)
  if(! is.null(seed))
    set.seed(seed)

  I <- nrow(A)
  N <- Q * I

  ## incidence matrix of covariates
  levels.years <- as.character(seq(from=2010, to=2010+Q-1))
  if(N %% Q == 0){
    years <- rep(levels.years, each=N / Q)
  } else
    years <- sort(sample(x=levels.years, size=N, replace=TRUE))
  years <- as.factor(years)
  if(Q == 1){
    W <- matrix(1, nrow=N, ncol=Q)
  } else
    W <- stats::model.matrix(~ years)
  dat <- data.frame(year=years)

  ## "fixed" effects
  if(Q == 1){
    C <- matrix(data=mu, nrow=Q, ncol=T)
  } else
    C <- matrix(data=c(mu, stats::rnorm(n=(Q-1)*T, mean=mean.C, sd=sd.C)),
                byrow=TRUE, nrow=Q, ncol=T)

  ## incidence matrix of genetic "random" effects
  levels.genos <- rownames(A)
  genos <- rep(NA, N)
  for(year in levels.years)
    genos[years == year] <- levels.genos[1:sum(years == year)]
  genos <- as.factor(genos)
  Z <- stats::model.matrix(~ genos - 1)
  dat$geno <- genos

  ## additive genetic component
  G.A <- matrix(0, I, T)
  if(T == 1){
    if(is.null(V.G.A)){
      sqrt.V.G.A <- abs(stats::rcauchy(n=1, location=0, scale=scale.hC.G.A))
      V.G.A <- sqrt.V.G.A^2
    }
    G.A <- matrix(MASS::mvrnorm(n=1, mu=rep(0, I), Sigma=V.G.A * A))
  } else{ # T > 1
    if(is.null(V.G.A))
      V.G.A <- stats::rWishart(n=1, df=nu.G.A, Sigma=diag(T))[,,1]
    G.A <- rmatnorm(n=1, M=matrix(data=0, nrow=I, ncol=T),
                    U=A, V=V.G.A)[,,1]
  }
  rownames(G.A) <- rownames(A)

  ## dominant genetic component
  G.D <- matrix(0, I, T)
  if(! is.null(D)){
    if(T == 1){
      if(is.null(V.G.D)){
        sqrt.V.G.D <- abs(stats::rcauchy(n=1, location=0, scale=scale.hC.G.D))
        V.G.D <- sqrt.V.G.D^2
      }
      G.D <- matrix(MASS::mvrnorm(n=1, mu=rep(0, I), Sigma=V.G.D * D))
    } else{ # T > 1
      if(is.null(V.G.D))
        V.G.D <- stats::rWishart(n=1, df=nu.G.D, Sigma=diag(T))[,,1]
      G.D <- rmatnorm(n=1, M=matrix(data=0, nrow=I, ncol=T),
                      U=D, V=V.G.D)[,,1]
    }
    rownames(G.D) <- rownames(D)
  }

  ## errors
  E <- matrix(0, N, T)
  if(T == 1){
    if(is.infinite(err.df)){
      if(is.null(V.E)){
        sqrt.V.E <- abs(stats::rcauchy(n=1, location=0, scale=scale.hC.E))
        V.E <- sqrt.V.E^2
      }
      E <- matrix(stats::rnorm(n=N, mean=0, sd=sqrt(V.E)))
    } else
      E <- matrix(stats::rt(n=N, df=err.df, ncp=0))
  } else{ # T > 1
    if(is.null(V.E))
      V.E <- stats::rWishart(n=1, df=nu.E, Sigma=diag(T))[,,1]
    E <- rmatnorm(n=1, M=matrix(data=0, nrow=N, ncol=T),
                  U=diag(N), V=V.E)[,,1]
  }

  ## phenotypes
  Y <- W %*% C + Z %*% G.A + Z %*% G.D + E
  if(perc.NA > 0){
    idx <- sample.int(n=N*T, size=floor(perc.NA/100 * N*T))
    Y[idx] <- NA
  }
  for(t in 1:T)
    dat[[paste0("response", t)]] <- Y[,t]

  return(list(Y=Y,
              W=W, C=C,
              Z=Z, V.G.A=V.G.A, G.A=G.A,
              V.G.D=V.G.D, G.D=G.D,
              V.E=V.E,
              dat=dat))
}

##' Animal model
##'
##' Given T=1 trait, I genotypes, Q covariates and N=I*Q phenotypes per trait, compute the BLUEs and BLUPs of the following "animal model" via Henderson's mixed model equations: y = W c + Z g_A + e, where y is N x 1; W is N x Q; Z is N x I; g_A ~ Normal_Ix1(0, sigma_A^2 A); e ~ Normal_Nx1(0, V.E Id_N).
##' @param y vector of phenotypes
##' @param W incidence matrix of fixed effects
##' @param Z incidence matrix of random effects, the breeding values
##' @param sigma.A2 variance component of the breeding values
##' @param Ainv inverse of A, the matrix of additive genetic relationships
##' @param V.E variance component of the errors
##' @return vector of length QI containing the BLUEs of c and the BLUPs of g_A
##' @author Timothee Flutre
##' @examples
##' ## simulate genotypes
##' set.seed(1859)
##' X <- simulGenosDose(nb.genos=200, nb.snps=2000)
##'
##' ## estimate the additive genetic relationships
##' A <- estimGenRel(X, relationships="additive", method="vanraden1", verbose=0)
##' A <- as.matrix(Matrix::nearPD(A)$mat) # not always necessary
##'
##' ## simulate phenotypes
##' model <- simulAnimalModel(T=1, Q=3, A=A, V.G.A=15, V.E=5)
##'
##' ## calculate BLUEs and BLUPs
##' fit <- mme(y=model$Y[,1,drop=FALSE], W=model$W, Z=model$Z,
##'            sigma.A2=model$V.G.A, Ainv=solve(A), V.E=model$V.E)
##' cbind(model$C, fit[1:3])
##' cor(model$G.A, fit[4:length(fit)])
##' @export
mme <- function(y, W, Z, sigma.A2, Ainv, V.E){
  stopifnot(is.matrix(y),
            is.matrix(W),
            is.matrix(Z),
            is.matrix(Ainv),
            ncol(y) == 1,
            nrow(y) == nrow(W),
            nrow(W) == nrow(Z),
            nrow(Ainv) == ncol(Z))

  lambda <- V.E / sigma.A2

  lhs <- matrix(data=NA, nrow=ncol(W)+ncol(Z), ncol=ncol(W)+ncol(Z))
  lhs[1:ncol(W), 1:ncol(W)] <- crossprod(W, W) # faster than t(W) %*% W
  lhs[(ncol(W)+1):nrow(lhs), 1:ncol(W)] <- crossprod(Z, W)
  lhs[1:ncol(W), (ncol(W)+1):ncol(lhs)] <- crossprod(W, Z)
  lhs[(ncol(W)+1):nrow(lhs), (ncol(W)+1):ncol(lhs)] <-
    crossprod(Z, Z) + lambda * Ainv

  rhs <- matrix(data=NA, nrow=nrow(lhs), ncol=1)
  rhs[1:ncol(W), 1] <- crossprod(W, y)
  rhs[(ncol(W)+1):nrow(rhs), 1] <- crossprod(Z, y)

  theta.hat <- solve(lhs, rhs) # faster than solve(lhs) %*% rhs

  return(c(theta.hat))
}

##' Animal model
##'
##' Given I genotypes, Q covariates and N=I*Q phenotypes for the trait, fit an "animal model" with the lme4 package via the following likelihood: y = W c + Z g_A + Z g_D + epsilon, where y is Nx1; W is NxQ; Z is NxI; g_A ~ Normal_I(0, sigma_A^2 A) with A the known matrix of additive genetic relationships; g_D ~ Normal_I(0, sigma_D^2 D) with D the known matrix of dominant genetic relationships; epsilon ~ Normal_N(0, sigma^2 Id_N); Cov(g_A,g_D)=0; Cov(g_A,e)=0; Cov(g_D,e)=0.
##' @param formula formula (see \code{\link[lme4]{lmer}})
##' @param dat data.frame containing the data corresponding to formula and relmat (see \code{\link[lme4]{lmer}})
##' @param relmat list containing the matrices of genetic relationships (A is compulsory but D is optional); should use the same names as the colnames in data; can be in the "matrix" class (base) or the "dsCMatrix" class (Matrix package); see \code{\link{estimGenRel}}
##' @param REML default is TRUE (use FALSE to compare models with different fixed effects)
##' @param ci.meth method to compute confidence intervals (profile/boot); if not NULL, use \code{\link[lme4]{confint.merMod}}
##' @param verbose verbosity level (0/1)
##' @return list which first component is a \code{\link[lme4]{merMod}} object and second component a data.frame with confidence intervals (if ci.meth is not NULL)
##' @author Timothee Flutre (inspired by Ben Bolker at http://stackoverflow.com/q/19327088/597069)
##' @note If A is not positive definite, an error will be raised (via \code{\link[base]{chol}}); in such cases, using the nearPD function from the Matrix package can be useful.
##' @seealso \code{\link{inlaAM}}, \code{\link{jagsAM}}, \code{\link{stanAM}}
##' @examples
##' ## simulate genotypes
##' set.seed(1859)
##' X <- simulGenosDose(nb.genos=200, nb.snps=2000)
##'
##' ## simulate phenotypes with only additive part of genotypic values
##' A <- estimGenRel(X, relationships="additive", method="vanraden1", verbose=0)
##' A <- as.matrix(Matrix::nearPD(A)$mat) # not always necessary
##' modelA <- simulAnimalModel(T=1, Q=3, A=A, V.G.A=15, V.E=5)
##'
##' ## infer with lme4
##' fitA <- lmerAM(formula=response1 ~ year + (1|geno), dat=modelA$dat,
##'                relmat=list(geno=A), verbose=0)
##' summary(fitA$merMod)
##' REMLcrit(fitA$merMod)
##' extractAIC(fitA$merMod)
##' summary(residuals(fitA$merMod)) # "deviance residuals"
##' summary(residuals(fitA$merMod) / sigma(fitA$merMod)) # "scaled/Pearson residuals"
##' c(modelA$C); modelA$V.G.A; modelA$V.E
##' fixef(fitA$merMod)
##' coefficients(summary(fitA$merMod))[, "Std. Error"]
##' vc <- as.data.frame(VarCorr(fitA$merMod))
##' c(vc[vc$grp == "geno", "vcov"], vc[vc$grp == "Residual", "vcov"])
##' blups.geno <- ranef(fitA$merMod, condVar=TRUE, drop=TRUE)$geno
##' var.blups.geno <- setNames(attr(blups.geno, "postVar"), names(blups.geno))
##'
##' ## simulate phenotypes with additive and dominant parts of genotypic values
##' D <- estimGenRel(X, relationships="dominant", method="vitezica", verbose=0)
##' D <- as.matrix(Matrix::nearPD(D)$mat) # not always necessary
##' modelAD <- simulAnimalModel(T=1, Q=3, A=A, V.G.A=15, V.E=5,
##'                             D=D, V.G.D=3)
##'
##' ## infer with lme4
##' modelAD$dat$geno.add <- modelAD$dat$geno
##' modelAD$dat$geno.dom <- modelAD$dat$geno; modelAD$dat$geno <- NULL
##' fitAD <- lmerAM(formula=response1 ~ year + (1|geno.add) + (1|geno.dom),
##'                 dat=modelAD$dat, relmat=list(geno.add=A, geno.dom=D),
##'                 verbose=0)
##' summary(fitAD$merMod)
##' c(modelAD$C); modelAD$V.G.A; modelAD$V.E; modelAD$V.G.D
##' fixef(fitAD$merMod)
##' vc <- as.data.frame(VarCorr(fitAD$merMod))
##' c(vc[vc$grp == "geno.add", "vcov"], vc[vc$grp == "Residual", "vcov"],
##'   vc[vc$grp == "geno.dom", "vcov"])
##' @export
lmerAM <- function(formula, dat, relmat, REML=TRUE, ci.meth=NULL, verbose=1){
  requireNamespaces(c("lme4", "Matrix"))
  stopifnot(is.data.frame(dat),
            all(! duplicated(colnames(dat))),
            is.list(relmat),
            all(! duplicated(names(names(relmat)))),
            all(names(relmat) %in% colnames(dat)),
            is.logical(REML))
  for(i in seq_along(relmat))
    stopifnot(is.matrix(relmat[[i]]) || class(relmat[[i]]) == "dsCMatrix")
  if(! is.null(ci.meth))
    stopifnot(ci.meth %in% c("profile", "boot"))

  if(verbose > 0)
    write("parse the formula ...", stdout())
  parsedFormula <- lme4::lFormula(formula=formula,
                                  data=dat,
                                  control=lme4::lmerControl(
                                      check.nobs.vs.nlev="ignore",
                                      check.nobs.vs.nRE="ignore"),
                                  REML=REML)

  if(verbose > 0)
    write(paste0("structure the design and covariance matrices",
                 " of the random effects ..."),
          stdout())
  relfac <- relmat
  flist <- parsedFormula$reTrms[["flist"]] # list of grouping factors
  asgn <- attr(flist, "assign")
  Ztlist <- parsedFormula$reTrms[["Ztlist"]] # list of transpose of the sparse model matrices
  for(i in seq_along(relmat)) {
    tn <- which(match(names(relmat)[i], names(flist)) == asgn)
    zn <- rownames(Ztlist[[i]])
    relmat[[i]] <- Matrix::Matrix(relmat[[i]][zn,zn], sparse=TRUE)
    relfac[[i]] <- chol(relmat[[i]])
    Ztlist[[i]] <- relfac[[i]] %*% Ztlist[[i]]
  }
  parsedFormula$reTrms[["Ztlist"]] <- Ztlist
  parsedFormula$reTrms[["Zt"]] <- do.call(Matrix::rBind, Ztlist)

  if(verbose > 0)
    write("make the deviance function ..", stdout())
  devianceFunction <- do.call(lme4::mkLmerDevfun, parsedFormula)

  if(verbose > 0)
    write("optimize the deviance function ...", stdout())
  optimizerOutput <- lme4::optimizeLmer(devianceFunction)

  if(verbose > 0)
    write("make the output ...", stdout())
  fit <- lme4::mkMerMod(rho=environment(devianceFunction),
                        opt=optimizerOutput,
                        reTrms=parsedFormula$reTrms,
                        fr=parsedFormula$fr)

  ci <- NULL
  if(! is.null(ci.meth)){
    if(verbose > 0)
      write("compute confidence intervals ...", stdout())
    suppressMessages(ci <- lme4::confint.merMod(fit, method=ci.meth, oldNames=FALSE))
  }

  return(list(merMod=fit, ci=ci))
}

##' Animal model
##'
##' Given I genotypes, Q covariates and N=I*Q phenotypes for the trait, fit an "animal model" with the INLA package via the following likelihood: y = W c + Z g_A + Z g_D + epsilon, where y is Nx1; W is NxQ; Z is NxI; g_A ~ Normal_I(0, sigma_A^2 A) with A the known matrix of additive genetic relationships; g_D ~ Normal_I(0, sigma_D^2 D) with D the known matrix of dominant genetic relationships; epsilon ~ Normal_N(0, sigma^2 Id_N); Cov(g_A,g_D)=0; Cov(g_A,e)=0; Cov(g_D,e)=0.
##' @param dat data.frame containing the data corresponding to relmat; should have a column grep-able for "response" as well as a column "geno.add" used with matrix A; if a column "geno.dom" exists, it will be used with matrix D; any other column will be interpreted as corresponding to "fixed effects"
##' @param relmat list containing the matrices of genetic relationships (see \code{\link{estimGenRel}}); additive relationships (matrix A) are compulsory, with name "geno.add"; dominant relationships (matrix D) are optional, with name "geno.dom"; can be in the "matrix" class (base) or the "dsCMatrix" class (Matrix package)
##' @param family a string indicating the likelihood family (see \code{\link[INLA]{inla}})
##' @param nb.threads maximum number of threads the inla-program will use (see \code{\link[INLA]{inla}})
##' @param verbose verbosity level (0/1)
##' @param silent if equal to TRUE, then the inla-program would be silent (see \code{\link[INLA]{inla}})
##' @return \code{\link[INLA]{inla}} object
##' @author Timothee Flutre
##' @seealso \code{\link{lmerAM}}, \code{\link{jagsAM}}, \code{\link{stanAM}}
##' @examples
##' if(require(INLA)){
##'   ## simulate genotypes
##'   set.seed(1859)
##'   X <- simulGenosDose(nb.genos=200, nb.snps=2000)
##'
##'   ## simulate phenotypes with only additive part of genotypic values
##'   A <- estimGenRel(X, relationships="additive", method="vanraden1", verbose=0)
##'   A <- as.matrix(Matrix::nearPD(A)$mat) # not always necessary
##'   modelA <- simulAnimalModel(T=1, Q=3, A=A, V.G.A=15, V.E=5)
##'
##'   ## infer with INLA
##'   library(INLA)
##'   modelA$dat$geno.add <- modelA$dat$geno; modelA$dat$geno <- NULL
##'   fitA <- inlaAM(dat=modelA$dat, relmat=list(geno.add=A))
##'   summary(fitA)
##'   c(modelA$C); 1/modelA$V.G.A; 1/modelA$V.E
##'
##'   ## simulate phenotypes with additive and dominant parts of genotypic values
##'   D <- estimGenRel(X, relationships="dominant", method="vitezica", verbose=0)
##'   D <- as.matrix(Matrix::nearPD(D)$mat) # not always necessary
##'   modelAD <- simulAnimalModel(T=1, Q=3, A=A, V.G.A=15, V.E=5,
##'                               D=D, V.G.D=3)
##'
##'   ## infer with INLA
##'   modelAD$dat$geno.add <- modelAD$dat$geno
##'   modelAD$dat$geno.dom <- modelAD$dat$geno; modelAD$dat$geno <- NULL
##'   fitAD <- inlaAM(dat=modelAD$dat, relmat=list(geno.add=A, geno.dom=D))
##'   summary(fitAD)
##'   c(modelAD$C); 1/modelAD$V.G.A; 1/modelAD$V.E; 1/modelAD$V.G.D
##' }
##' @export
inlaAM <- function(dat, relmat, family="gaussian",
                   nb.threads=1, verbose=0, silent=TRUE){
  requireNamespaces("INLA")
  stopifnot(is.data.frame(dat),
            sum(grepl("response", colnames(dat))) == 1,
            "geno.add" %in% colnames(dat),
            is.list(relmat),
            ! is.null(names(relmat)),
            "geno.add" %in% names(relmat))

  ## make formula with response, intercept and additive genotypic value
  formula <- paste0(colnames(dat)[grepl("response", colnames(dat))],
                    " ~ 1",
                    " + f(geno.add, model=\"generic0\", constr=TRUE",
                    ", hyper=list(theta=list(param=c(0.5, 0.5), fixed=FALSE))",
                    ", Cmatrix=solve(relmat[[\"geno.add\"]])",
                    ", values=levels(dat$geno.add))")

  ## add dominant genotypic value, if present
  if("geno.dom" %in% names(relmat))
    formula <- paste0(formula,
                      " + f(geno.dom, model=\"generic0\", constr=TRUE",
                      ", hyper=list(theta=list(param=c(0.5, 0.5), fixed=FALSE))",
                      ", Cmatrix=solve(relmat[[\"geno.dom\"]])",
                      ", values=levels(dat$geno.dom))")

  ## add "fixed effects", if any
  for(cn in colnames(dat))
    if(! grepl("response", cn) & cn != "geno.add" & cn != "geno.dom")
      formula <- paste0(formula, " + ", cn)

  ## finalize formula
  formula <- stats::as.formula(formula)

  ## fit the model with INLA
  fit <- INLA::inla(formula=formula,
                    family=family,
                    data=dat,
                    control.compute=list(dic=TRUE),
                    num.threads=nb.threads,
                    verbose=verbose,
                    silent=silent)

  return(fit)
}

##' Animal model
##'
##' Given I genotypes, Q covariates and N=I*Q phenotypes for the trait, fit an "animal model" with the rjags package via the following likelihood: y = W c + Z g_A + Z g_D + epsilon, where y is Nx1; W is NxQ; Z is NxI; g_A ~ Normal_I(0, sigma_A^2 A) with A the known matrix of additive genetic relationships; g_D ~ Normal_I(0, sigma_D^2 D) with D the known matrix of dominant genetic relationships; epsilon ~ Normal_N(0, sigma^2 Id_N); Cov(g_A,g_D)=0; Cov(g_A,e)=0; Cov(g_D,e)=0.
##' @param dat data.frame containing the data corresponding to relmat; should have a column grep-able for "response" as well as a column "geno.add" used with matrix A; if a column "geno.dom" exists, it will be used with matrix D; any other column will be interpreted as corresponding to "fixed effects"
##' @param relmat list containing the matrices of genetic relationships (see \code{\link{estimGenRel}}); additive relationships (matrix A) are compulsory, with name "geno.add"; dominant relationships (matrix D) are optional, with name "geno.dom"; can be in the "matrix" class (base) or the "dsCMatrix" class (Matrix package)
##' @param inits list of initial values (possible to use 1 sub-list per chain, see \code{\link[rjags]{jags.model}}); if NULL, JAGS will choose typical values (usually mean, median, or mode of the prior)
##' @param priors list of specifying priors; each component is itself a list for which "dist" specifies the distribution and "par" the parameter; a component named "fix" is used for fixed effects,  so that ("c",5) means c[q] ~ Cauchy(0,5) and ("n",10^6) means c[q] ~ Normal(0,10^6), but note that the global intercept always is Normal(mean(y),10^6); a component named "vc" is used for variance components, so that ("hc",5) means stdev ~ half-Cauchy(0,5) and ("ig",0.001) means var ~ InvGamma(shape=0.001,rate=0.001)
##' @param nb.chains number of independent chains to run (see \code{\link[rjags]{jags.model}})
##' @param nb.adapt number of iterations for adaptation (see \code{\link[rjags]{jags.model}})
##' @param burnin number of initial iterations to discard (see the update function of the rjags package)
##' @param nb.iters number of iterations to monitor (see \code{\link[rjags]{coda.samples}})
##' @param thin thinning interval for monitored iterations (see \code{\link[rjags]{coda.samples}})
##' @param progress.bar type of progress bar (text/gui/none or NULL)
##' @param rm.jags.file remove the file specifying the model written in the JAGS-dialect of the BUGS language
##' @param verbose verbosity level (0/1)
##' @return \code{\link[coda]{mcmc.list}} object
##' @author Timothee Flutre
##' @seealso \code{\link{lmerAM}}, \code{\link{inlaAM}}, \code{\link{stanAM}}
##' @examples
##' \dontrun{## simulate genotypes
##' set.seed(1859)
##' X <- simulGenosDose(nb.genos=200, nb.snps=2000)
##'
##' ## simulate phenotypes with only additive part of genotypic values
##' A <- estimGenRel(X, relationships="additive", method="vanraden1", verbose=0)
##' A <- as.matrix(Matrix::nearPD(A)$mat) # not always necessary
##' modelA <- simulAnimalModel(T=1, Q=3, A=A, V.G.A=15, V.E=5)
##'
##' ## infer with rjags
##' modelA$dat$geno.add <- modelA$dat$geno; modelA$dat$geno <- NULL
##' fitA <- jagsAM(dat=modelA$dat, relmat=list(geno.add=A))
##' plotMcmcChain(fitA[[1]], "V.G.A", 1:4, modelA$V.G.A)
##' cbind(truth=c(c(modelA$C), modelA$V.G.A, modelA$V.E),
##'       summaryMcmcChain(fitA[[1]], c("c[1]", "c[2]", "c[3]", "sigma.A2", "V.E")))
##' }
##' @export
jagsAM <- function(dat, relmat, inits=NULL,
                   priors=list(fix=list(dist="c", par=5),
                               vc=list(dist="hc", par=5)),
                   nb.chains=1, nb.adapt=10^3, burnin=10^2,
                   nb.iters=10^3, thin=10,
                   progress.bar=NULL, rm.jags.file=TRUE, verbose=0){
  requireNamespaces("rjags")
  stopifnot(is.data.frame(dat),
            sum(grepl("response", colnames(dat))) == 1,
            "geno.add" %in% colnames(dat),
            is.list(relmat),
            ! is.null(names(relmat)),
            "geno.add" %in% names(relmat),
            is.list(priors),
            all(c("fix", "vc") %in% names(priors)),
            all(c("dist", "par") %in% names(priors$fix)),
            all(c("dist", "par") %in% names(priors$vc)))
  if(! is.null(inits))
    stopifnot(is.list(inits))

  ## define model in JAGS dialect of BUGS language in separate file
  st <- Sys.time()
  jags.file <- tempfile(paste0("jagsAM_", format(st, "%Y-%m-%d-%H-%M-%S"), "_"),
                        getwd(), ".jags")
  model.txt <- paste0(
      "# Goal: fit an \"animal model\" with rjags",
      "\n# Author: Timothee Flutre (INRA)",
      "\n# Source: rutilstimflutre ", utils::packageVersion("rutilstimflutre"),
      "\n# Date: ", format(st, "%Y-%m-%d %H:%M:%S"))
  model.txt <- paste0(model.txt, "

model {
  # likelihood
  for(n in 1:N){
    y[n] ~ dnorm(W[n,] %*% c + Z[n,] %*% g.A")
  if("geno.dom" %in% names(relmat))
    model.txt <- paste0(model.txt, " + Z[n,] %*% g.D")
  model.txt <- paste0(model.txt, ", 1/V.E)
  }")
  model.txt <- paste0(model.txt, "

  # priors
  c[1] ~ dnorm(mean.mu, 10^(-6))
  for(q in 2:Q){")
  if(priors$fix$dist == "c"){
    model.txt <- paste0(model.txt, "
    c[q] ~ dt(0, par.c, 1)")
  } else if(priors$fix$dist == "n")
    model.txt <- paste0(model.txt, "
    c[q] ~ dnorm(0, 1/par.c)")
  model.txt <- paste0(model.txt, "
  }")
  if(priors$vc$dist == "hc"){
    model.txt <- paste0(model.txt, "
  sigma.A ~ dt(0, par.vc, 1)T(0,)
  sigma.A2 <- sigma.A^2")
  } else if(priors$vc$dist == "ig")
    model.txt <- paste0(model.txt, "
  tau.A ~ dgamma(par.vc, par.vc)
  sigma.A2 <- 1 / tau.A")
  model.txt <- paste0(model.txt, "
  Sigma.g.A <- sigma.A2 * A
  g.A ~ dmnorm(mean.g.A, inverse(Sigma.g.A))")
  if("geno.dom" %in% names(relmat)){
    if(priors$vc$dist == "hc"){
    model.txt <- paste0(model.txt, "
  sigma.D ~ dt(0, par.vc, 1)T(0,)
  sigma.D2 <- sigma.D^2")
    } else if(priors$vc$dist == "ig")
      model.txt <- paste0(model.txt, "
  tau.D ~ dgamma(par.vc, par.vc)
  sigma.D2 <- 1 / tau.D")
    model.txt <- paste0(model.txt, "
  Sigma.g.D <- sigma.D2 * D
  g.D ~ dmnorm(mean.g.D, inverse(Sigma.g.D))")
  }
  if(priors$vc$dist == "hc"){
    model.txt <- paste0(model.txt, "
  sigma.E ~ dt(0, par.vc, 1)T(0,)
  V.E <- sigma.E^2")
  } else if(priors$vc$dist == "ig")
    model.txt <- paste0(model.txt, "
  tau.E ~ dgamma(par.vc, par.vc)
  V.E <- 1 / tau")
  model.txt <- paste0(model.txt, "
}")
  cat(model.txt, file=jags.file)

  ## read model file and create object
  y <- dat[, grepl("response", colnames(dat))]
  W <- matrix(1, nrow=nrow(dat), ncol=1)
  colnames(W) <- "mu"
  for(j in 1:ncol(dat))
    if(! grepl("response", colnames(dat)[j]) & colnames(dat)[j] != "geno.add" &
       colnames(dat)[j] != "geno.dom"){
      W <- cbind(W, stats::model.matrix(~ dat[,j])[,-1])
    }
  colnames(W) <- gsub("dat\\[, j\\]", "", colnames(W))
  data.list <- list(N=nrow(dat), Q=ncol(W),
                    y=y, W=W, Z=stats::model.matrix(~ dat[,"geno.add"] - 1),
                    mean.mu=mean(y), par.c=priors$fix$par,
                    par.vc=priors$vc$par, A=relmat[["geno.add"]],
                    mean.g.A=rep(0, nlevels(dat[,"geno.add"])))
  if("geno.dom" %in% names(relmat)){
    data.list$D <- relmat[["geno.dom"]]
    data.list$mean.g.D <- rep(0, nlevels(dat[,"geno.add"]))
  }
  jags <- rjags::jags.model(file=jags.file, data=data.list, inits=inits,
                            n.chains=nb.chains,
                            n.adapt=nb.adapt,
                            quiet=ifelse(verbose > 0, FALSE, TRUE))
  if(rm.jags.file)
    file.remove(jags.file)

  ## update model for burn-in period
  stats::update(jags, n.iter=burnin, progress.bar=progress.bar)

  ## extract samples from model
  vn <- c("c", "sigma.A2")
  if("geno.dom" %in% names(relmat))
    vn <- c(vn, "sigma.D2")
  vn <- c(vn, "V.E")
  vn <- c(vn, "g.A")
  if("geno.dom" %in% names(relmat))
    vn <- c(vn, "g.D")
  fit <- rjags::coda.samples(model=jags,
                             variable.names=vn,
                             n.iter=nb.iters,
                             thin=thin,
                             progress.bar=progress.bar)

  return(fit)
}

##' Animal model
##'
##'
##' @param stan.file path to a file; if NULL, a temporary one will be created
##' @param errors.Student use a Student's t distribution for the errors (useful to handle outliers)
##' @param missing.phenos indicating if there is any missing phenotypes to impute jointly with the other unknown variables
##' @return invisible path to stan.file
##' @author Timothee Flutre
##' @seealso \code{\link{stanAM}}
##' @export
stanAMwriteModel <- function(stan.file=NULL,
                             errors.Student=FALSE,
                             missing.phenos=FALSE){
  st <- Sys.time()
  if(is.null(stan.file))
    stan.file <- tempfile(paste0("stanAM_", format(st, "%Y-%m-%d-%H-%M-%S"), "_"),
                          getwd(), ".stan")

  ## meta-data
  model.code <- "# model: stanAM
# copyright: Inra
# license: AGPL-3+
# persons: Timothée Flutre [cre,aut]"
  model.code <- paste0(model.code, "
# date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "
")

  ## block: data
  model.code <- paste0(model.code, "
data {
  int<lower=0> Q;         // num years
  int<lower=0> I;         // num genotypes
  real loc_mu;            // prior mean of the intercept
  real<lower=0> scale_mu; // prior sd of the intercept
  matrix[I,I] A;          // additive genetic relationships matrix")
  if(missing.phenos){
    model.code <- paste0(model.code, "
  int<lower=0> N_obs;     // num observed phenotypes
  int<lower=0> N_mis;     // num missing phenotypes
  matrix[N_obs,Q] W_obs;  // incidence matrix
  matrix[N_mis,Q] W_mis;  // incidence matrix
  matrix[N_obs,I] Z_obs;  // incidence matrix
  matrix[N_mis,I] Z_mis;  // incidence matrix
  vector[N_obs] y_obs;    // observed phenotypes")
  } else{
    model.code <- paste0(model.code, "
  int<lower=0> N;         // num phenotypic measurements
  matrix[N,Q] W;          // incidence matrix
  matrix[N,I] Z;          // incidence matrix
  vector[N] y;            // phenotypes")
  }
  if(errors.Student){
    model.code <- paste0(model.code, "
  real<lower=0> nu;       // degrees of freedom of the Student's t
")
  } else{
    model.code <- paste0(model.code, "
")
  }
  model.code <- paste0(model.code, "}\n")

  ## block: transformed data
  model.code <- paste0(model.code, "
transformed data {
  matrix[I,I] CA;
  CA <- cholesky_decompose(A);")
  model.code <- paste0(model.code, "
}
")

  ## block: parameters
  model.code <- paste0(model.code, "
parameters {
  vector[Q] c;
  real<lower=0> sigma_E;
  vector[I] g_A_z;          // primitive of g_A
  real<lower=0> sigma_g_A;")
  if(missing.phenos){
    model.code <- paste0(model.code, "
  vector[N_mis] y_mis;")
  }
  model.code <- paste0(model.code, "
}
")

  ## block: model
  model.code <- paste0(model.code, "
model {
  vector[I] g_A;
  c[1] ~ cauchy(loc_mu, scale_mu); // intercept
  for(q in 2:Q)
    c[q] ~ cauchy(0, 5);
  sigma_E ~ cauchy(0, 5);  // implicit half-Cauchy
  sigma_g_A ~ cauchy(0, 5);
  g_A_z ~ normal(0, 1);
  g_A <- sigma_g_A * (CA * g_A_z);// implies g_A ~ multi_normal(0, sigma_g_A^2 * A)")
  if(errors.Student){
    if(missing.phenos){
      model.code <- paste0(model.code, "
  y_obs ~ student_t(nu, W_obs * c + Z_obs * g_A, sigma_E);
  y_mis ~ student_t(nu, W_mis * c + Z_mis * g_A, sigma_E);
")
    } else{
      model.code <- paste0(model.code, "
  y ~ student_t(nu, W * c + Z * g_A, sigma_E);
")
    }
  } else{
    if(missing.phenos){
      model.code <- paste0(model.code, "
  y_obs ~ normal(W_obs * c + Z_obs * g_A, sigma_E);
  y_mis ~ normal(W_mis * c + Z_mis * g_A, sigma_E);
")
    } else{
      model.code <- paste0(model.code, "
  y ~ normal(W * c + Z * g_A, sigma_E);
")
    }
  }
  model.code <- paste0(model.code, "}")

  model.code <- paste0(model.code, "

generated quantities {
  vector[I] g_A;
  g_A <- sigma_g_A * (CA * g_A_z);
")
  model.code <- paste0(model.code, "}")

  cat(model.code, file=stan.file)

  invisible(stan.file)
}

##' Animal model
##'
##' Given I genotypes, Q covariates and N=I*Q phenotypes for the trait, fit an "animal model" with the rstan package via the following likelihood: y = W c + Z g_A + Z g_D + epsilon, where y is Nx1; W is NxQ; Z is NxI; g_A ~ Normal_I(0, sigma_A^2 A) with A the known matrix of additive genetic relationships; g_D ~ Normal_I(0, sigma_D^2 D) with D the known matrix of dominant genetic relationships; epsilon ~ Normal_N(0, sigma^2 Id_N); Cov(g_A,g_D)=0; Cov(g_A,e)=0; Cov(g_D,e)=0. Missing phenotypes are jointly imputed with the other unknown variables, and the errors can follow a Student's t distribution to handle outliers.
##' @param dat data.frame containing the data corresponding to relmat; should have a column grep-able for "response" as well as a column "geno.add" used with matrix A; if a column "geno.dom" exists, it will be used with matrix D; any other column will be interpreted as corresponding to "fixed effects"
##' @param relmat list containing the matrices of genetic relationships (see \code{\link{estimGenRel}}); additive relationships (matrix A) are compulsory, with name "geno.add"; dominant relationships (matrix D) are optional, with name "geno.dom"; can be in the "matrix" class (base) or the "dsCMatrix" class (Matrix package)
##' @param errors.Student use a Student's t distribution for the errors (useful to handle outliers)
##' @param nb.chains number of independent chains to run
##' @param nb.iters number of iterations to monitor
##' @param burnin number of initial iterations to discard
##' @param thin thinning interval for monitored iterations
##' @param compile.only only compile the model, don't run it
##' @param rm.stan.file remove the file specifying the model written in the STAN language
##' @param rm.sm.file remove the file corresponding to the compiled model
##' @param verbose verbosity level (0/1)
##' @return path to compiled file if \code{compile.only=TRUE}, object of class \code{\link[rstan]{stanfit-class}} otherwise
##' @author Timothee Flutre
##' @seealso \code{\link{lmerAM}}, \code{\link{inlaAM}}, \code{\link{jagsAM}}
##' @examples
##' \dontrun{## simulate genotypes
##' set.seed(1859)
##' X <- simulGenosDose(nb.genos=200, nb.snps=2000)
##'
##' ## simulate phenotypes with only additive part of genotypic values
##' A <- estimGenRel(X, relationships="additive", method="vanraden1", verbose=0)
##' A <- as.matrix(Matrix::nearPD(A)$mat) # not always necessary
##' modelA <- simulAnimalModel(T=1, Q=3, A=A, V.G.A=15, V.E=5)
##'
##' ## infer with rstan
##' modelA$dat$geno.add <- modelA$dat$geno; modelA$dat$geno <- NULL
##' fitA <- stanAM(dat=modelA$dat, relmat=list(geno.add=A))
##' fitA <- rstan::As.mcmc.list(fitA)
##' plotMcmcChain(fitA[[1]], "sigma_g_A", 1:4, sqrt(modelA$V.G.A))
##' cbind(truth=c(c(modelA$C), sqrt(modelA$V.G.A), sqrt(modelA$V.E)),
##'       summaryMcmcChain(fitA[[1]], c("c[1]", "c[2]", "c[3]", "sigma_g_A", "sigma_E")))
##' }
##' @export
stanAM <- function(dat, relmat, errors.Student=FALSE,
                   nb.chains=1, nb.iters=10^3, burnin=10^2, thin=10,
                   compile.only=FALSE, rm.stan.file=TRUE, rm.sm.file=FALSE,
                   verbose=0){
  requireNamespaces("rstan")
  stopifnot(is.data.frame(dat),
            sum(grepl("response", colnames(dat))) == 1,
            "geno.add" %in% colnames(dat),
            is.list(relmat),
            ! is.null(names(relmat)),
            "geno.add" %in% names(relmat))

  ## make the input matrices from the input data.frame
  y <- dat[, grepl("response", colnames(dat))]
  if(all(is.matrix(y), ncol(y) == 1))
    y <- y[,1]
  is_y_obs <- (! is.na(y))
  W <- matrix(1, nrow=nrow(dat), ncol=1)
  colnames(W) <- "mu"
  for(j in 1:ncol(dat))
    if(! grepl("response", colnames(dat)[j]) & colnames(dat)[j] != "geno.add" &
       colnames(dat)[j] != "geno.dom"){
      W <- cbind(W, stats::model.matrix(~ dat[,j])[,-1])
    }
  colnames(W) <- gsub("dat\\[, j\\]", "", colnames(W))
  Z <- stats::model.matrix(~ dat[,"geno.add"] - 1)
  N <- nrow(W)
  Q <- ncol(W)
  I <- ncol(Z)

  ## define model in STAN language in separate file
  st <- Sys.time()
  stan.file <- tempfile(paste0("stanAM_", format(st, "%Y-%m-%d-%H-%M-%S"), "_"),
                        getwd(), ".stan")
  stanAMwriteModel(stan.file, errors.Student, sum(is_y_obs) != N)

  ## compile, or make the input list and run
  sm.file <- tempfile(paste0("stanAM_", format(st, "%Y-%m-%d-%H-%M-%S"), "_"),
                      getwd(), "_sm.RData")
  if(compile.only){
    if(verbose > 0)
      write(paste0("compile..."), stdout())
    rt <- rstan::stanc(file=stan.file, model_name="stanAM")
    sm <- rstan::stan_model(stanc_ret=rt,
                            verbose=ifelse(verbose > 0, TRUE, FALSE))
    save(sm, file=sm.file)
    return(sm.file)
  } else{
    ldat <- list(Q=Q,
                 I=I,
                 loc_mu=mean(y, na.rm=TRUE),
                 scale_mu=5,
                 A=relmat[["geno.add"]])
    if(sum(is_y_obs) != N){
      ldat$N_obs <- sum(is_y_obs)
      ldat$N_mis <- sum(! is_y_obs)
      ldat$W_obs<- W[is_y_obs,]
      ldat$W_mis<- W[! is_y_obs,]
      ldat$Z_obs <- Z[is_y_obs,]
      ldat$Z_mis <- Z[! is_y_obs,]
      ldat$y_obs <- y[is_y_obs]
    } else{
      ldat$N <- N
      ldat$W<- W
      ldat$Z <- Z
      ldat$y <- y
    }
    if(errors.Student){
      ldat$nu <- 3
    }
    if(file.exists(sm.file)){
      if(verbose > 0)
        write(paste0("run..."), stdout())
      load(sm.file)
      fit <- rstan::sampling(sm,
                             data=ldat,
                             iter=nb.iters + burnin,
                             warmup=burnin,
                             thin=thin,
                             chains=nb.chains)
    } else
      if(verbose > 0)
        write(paste0("compile and run..."), stdout())
      fit <- rstan::stan(file=stan.file,
                         data=ldat,
                         iter=nb.iters + burnin,
                         warmup=burnin,
                         thin=thin,
                         chains=nb.chains)
  }
  if(rm.stan.file)
    file.remove(stan.file)
  if(rm.sm.file)
    file.remove(sm.file)

  return(fit)
}

##' BVSR
##'
##' Simulate phenotypes according to the following model: Y = W c + Z X_A a + epsilon where Y is N x 1; W is N x Q; c is Q x 1; Z is N x I; X_A is I x P and epsilon is N x 1 with epsilon ~ Normal_N(0, sigma^2 Id) and c ~ Normal(mean_a, sd_a) so that sd_a is large ("fixed effect"). For SNP p, a_p ~ Prob(gamma_p=1) Normal_1(0, sigma_a^2) + Prob(gamma_p=0) delta_0, where Prog(gamma_p=1) is named pi, and delta_0 is Dirac's delta function at 0. For the case where pi is small, see Guan & Stephens (2011), Carbonetto & Stephens (2012), Peltola et al (2012), Verzelen (2012).
##' @param Q number of fixed effects, i.e. the intercept plus the number of years during which individuals are phenotyped (starting in 2010)
##' @param mu overall mean
##' @param mean.c mean of the prior on c[2:Q]
##' @param sd.c std dev of the prior on c[2:Q]
##' @param X matrix of SNP genotypes encoded in number of copies of the 2nd allele, i.e. as allele doses in {0,1,2}, with individuals in rows and SNPs in columns (SNPs with missing values or low MAF should be discarded beforehand); will be used in the simulations as X_A which is the column-centered version of X when encoded in {-1,0,1}
##' @param pi proportion of marker effects (a) that are non-zero; setting pi at 1 means simulating from the additive infinitesimal model (equivalent to ridge regression)
##' @param pve.A proportion of phenotypic variance explained by SNPs with non-zero effect ("heritability"); PVE_A = V[X_A a] / V[y]; used along with option sigma.a2 to choose a value for sigma^2
##' @param sigma.a2 prior variance of the non-zero additive effects
##' @param perc.NA percentage of missing phenotypes, at random
##' @param err.df degrees of freedom of errors' Student's t-distribution
##' @param seed seed for the pseudo-random number generator
##' @return list
##' @author Timothee Flutre
##' @examples
##' ## simulate genotypes
##' set.seed(1859)
##' I <- 200
##' P <- 2000
##' X <- simulGenosDose(nb.genos=I, nb.snps=P)
##'
##' ## additive sparse genetic architecture
##' ## choose pi so that sum(gamma * (1 + log(P / sum(gamma)))) < I
##' Q <- 3
##' model <- simulBvsr(Q=Q, X=X, pi=0.01, pve.A=0.7, sigma.a2=1)
##'
##' \dontrun{if(all(require(lme4), require(varbvs))){
##'   dat <- data.frame(response=model$Y[,1],
##'                     year=factor(rep(2010:(2010+Q-1), each=I)),
##'                     geno=factor(rep(rownames(X), Q)))
##'   fit1 <- lmer(formula=response ~ year + (1|geno), data=dat)
##'   cbind(model$c, blues <- fixef(fit1))
##'   blups <- ranef(fit1, drop=TRUE)$geno[rownames(X)]
##'   cor(model$g.A, blups)
##'
##'   fit2 <- varbvs(X=model$X.A, Z=NULL, y=blups, verbose=FALSE)
##'   print(fit2s <- summary(fit2))
##'
##'   (pi.hat <- 10^(fit2s$logodds$x0) / (1 + 10^(fit2s$logodds$x0)))
##'   (pi.hat.low <- 10^(fit2s$logodds$a) / (1 + 10^(fit2s$logodds$a)))
##'   (pi.hat.high <- 10^(fit2s$logodds$b) / (1 + 10^(fit2s$logodds$b)))
##'
##'   w <- c(normalizelogweights(fit2$logw))
##'   pips <- c(fit2$alpha %*% w)
##'   cols <- rep("black", P)
##'   cols[model$gamma != 0] <- "red"
##'   plot(x=1:P, y=pips, col=cols)
##'
##'   y.pred <- predict(fit2, X=model$X.A, Z=NULL)
##'   cor(blups, y.pred)
##' }}
##' @export
simulBvsr <- function(Q=3, mu=50, mean.c=5, sd.c=2,
                      X, pi=1, pve.A=0.7, sigma.a2=1,
                      perc.NA=0, err.df=Inf, seed=NULL){
  stopifnot(.isValidGenosDose(X),
            sd.c >= 0,
            pi >= 0,
            pi <= 1,
            pve.A >= 0,
            pve.A <= 1,
            sigma.a2 >= 0,
            perc.NA >= 0,
            perc.NA <= 100)
  if(! is.null(seed))
    set.seed(seed)

  I <- nrow(X)
  P <- ncol(X)
  if(Q > 1){
    N <- Q * I
  } else
    N <- I

  ## incidence matrix of the non-genetic predictors having "fixed effects"
  if(Q > 1){
    levels.years <- as.character(seq(from=2010, to=2010+Q-1))
    if(N %% Q == 0){
      years <- rep(levels.years, each=N / Q)
    } else
      years <- sort(sample(x=levels.years, size=N, replace=TRUE))
    years <- as.factor(years)
    W <- stats::model.matrix(~ years)
  } else
    W <- matrix(data=1, nrow=N, ncol=1)

  ## "fixed effects"
  if(Q > 1){
    c <- matrix(data=c(mu, stats::rnorm(n=Q-1, mean=mean.c, sd=sd.c)),
                nrow=Q, ncol=1)
  } else if(Q == 1){
    c <- matrix(data=mu, nrow=Q, ncol=1)
  } else
    c <- matrix(data=0, nrow=1, ncol=1)

  ## incidence matrices of the genetic predictors
  levels.inds <- rownames(X)
  inds <- rep(NA, N)
  if(Q > 1){
    for(year in levels.years)
      inds[years == year] <- levels.inds[1:sum(years == year)]
  } else
    inds <- levels.inds
  inds <- as.factor(inds)
  Z <- stats::model.matrix(~ inds - 1)
  X.A <- scale(x=X - 1, center=TRUE, scale=FALSE)

  ## incidence vector of the causal genetic predictors
  gamma <- stats::setNames(object=stats::rbinom(n=P, size=1, prob=pi),
                           nm=colnames(X))

  ## "random effects"
  a <- stats::setNames(object=stats::rnorm(n=P, mean=0, sd=sqrt(sigma.a2)),
                       nm=colnames(X))
  a[gamma == 0] <- 0
  g.A <- X.A %*% a

  ## errors
  sigma2 <- ((1 - pve.A) / pve.A) * stats::var(g.A)
  if(is.infinite(err.df)){
    epsilon <- matrix(stats::rnorm(n=N, mean=0, sd=sqrt(sigma2)))
  } else
    epsilon <- matrix(stats::rt(n=N, df=err.df, ncp=0))

  ## phenotypes
  Y <- W %*% c + Z %*% X.A %*% a + epsilon
  if(perc.NA > 0){
    idx <- sample.int(n=N, size=floor(perc.NA/100 * N))
    Y[idx, 1] <- NA
  }

  return(list(Y=Y,
              W=W, c=c,
              Z=Z, X.A=X.A, gamma=gamma,
              pi=pi, pve.A=pve.A, sigma.a2=sigma.a2, sigma2=sigma2,
              a=a, g.A=g.A))
}

##' BSLMM
##'
##' Simulate phenotypes according to the Bayesian sparse linear mixed model (Zhou, Carbonetto & Stephens, 2013)): y = W alpha + Z X_c beta-tilde + Z u + epsilon where y is N x 1; W is N x Q; alpha is Q x 1; Z is N x I; X_c is I x P; u is I x 1 and epsilon is N x 1. For all p, beta-tilde_p ~ pi Norm_1(0, sigma_betat^2/sigma^2) + (1 - pi) delta_0; u ~ Norm_I(0, (sigma_u^2/sigma^2) K) and epsilon ~ Norm_N(0, sigma^2 I).
##' @param Q number of fixed effects, i.e. the intercept plus the number of years during which individuals are phenotyped (starting in 2010)
##' @param mu overall mean
##' @param mean.a mean of the prior on alpha[2:Q]
##' @param sd.a std dev of the prior on alpha[2:Q]
##' @param X matrix of SNP genotypes encoded in number of copies of the 2nd allele, i.e. as allele doses in {0,1,2}, with individuals in rows and SNPs in columns (SNPs with missing values or low MAF should be discarded beforehand).
##' @param pi proportion of beta-tilde values that are non-zero; if NULL, log(pi) is drawn from Unif(log(1/P),log(1))
##' @param h approximation to E[PVE]; if NULL, drawn from Unif(0,1)
##' @param rho approximation to E[PGE]; if NULL and pi=0, then rho=0, else rho is drawn from Unif(0,1)
##' @param tau precision of the errors (1 / sigma^2)
##' @param enforce.zhou if TRUE, X will be transformed into X_c by centering columns, and K will be calculated as X_c X_c^T / P; otherwise, X will be transformed into (-1,0,1) and K will be calculated as \code{\link{estimGenRel}} with method="vanraden1"
##' @param perc.NA percentage of missing phenotypes, at random
##' @param err.df degrees of freedom of errors' Student's t-distribution
##' @param seed seed for the pseudo-random number generator
##' @return list
##' @author Timothee Flutre
##' @examples
##' set.seed(1859)
##' X <- simulGenosDose(nb.genos=200, nb.snps=2000)
##' afs <- estimAf(X)
##'
##' ## particular case: LMM (only u contributes)
##' mod.lmm <- simulBslmm(X=X, pi=0, h=0.7, seed=3591)
##'
##' ## particular case: BVSR (only beta-tilde contributes)
##' mod.bvsr <- simulBslmm(X=X, pi=0.1, h=0.7, rho=1, seed=3591)
##'
##' ## general case: BSLMM (both u and beta-tilde contribute)
##' mod.bslmm <- simulBslmm(X=X, pi=0.1, h=0.7, rho=0.7, seed=3591)
##'
##' \dontrun{library(rrBLUP)
##' mod.lmmC <- simulBslmm(X=X, pi=0, h=0.7, enforce.zhou=FALSE, seed=3591)
##' fit.lmmC1 <- mixed.solve(y=mod.lmmC$y, Z=mod.lmmC$Z,
##'                          K=mod.lmmC$X.c %*% t(mod.lmmC$X.c),
##'                          X=mod.lmmC$W)
##' cbind(truth=c(mod.lmmC$alpha), estim=fit.lmmC1$beta)
##' cbind(truth=c(mod.lmmC$sigma2, mod.lmmC$sigma.u2 / sum(2 * afs * (1-afs))),
##'       estim=c(fit.lmmC1$Ve, fit.lmmC1$Vu))
##' fit.lmmC2 <- mixed.solve(y=mod.lmmC$y, Z=mod.lmmC$Z %*% mod.lmmC$X.c,
##'                          X=mod.lmmC$W)
##' cbind(truth=c(mod.lmmC$alpha), estim=fit.lmmC2$beta)
##' cbind(truth=c(mod.lmmC$sigma2, mod.lmmC$sigma.u2 / sum(2 * afs * (1-afs))),
##'       estim=c(fit.lmmC2$Ve, fit.lmmC2$Vu))
##' }
##' @export
simulBslmm <- function(Q=3, mu=50, mean.a=5, sd.a=2,
                       X, pi=NULL, h=NULL, rho=NULL, tau=1,
                       enforce.zhou=TRUE, perc.NA=0, err.df=Inf,
                       seed=NULL){
  requireNamespaces("MASS")
  stopifnot(.isValidGenosDose(X))
  if(! is.null(seed))
    set.seed(seed)

  I <- nrow(X)
  P <- ncol(X)
  if(Q > 1){
    N <- Q * I
  } else
    N <- I

  ## incidence matrix of the non-genetic predictors having "fixed effects"
  if(Q > 1){
    levels.years <- as.character(seq(from=2010, to=2010+Q-1))
    if(N %% Q == 0){
      years <- rep(levels.years, each=N / Q)
    } else
      years <- sort(sample(x=levels.years, size=N, replace=TRUE))
    years <- as.factor(years)
    W <- stats::model.matrix(~ years)
  } else
    W <- matrix(data=1, nrow=N, ncol=1)

  ## "fixed" effects
  if(Q > 1){
    alpha <- matrix(data=c(mu, stats::rnorm(n=Q-1, mean=mean.a, sd=sd.a)),
                    nrow=Q, ncol=1)
  } else if(Q == 1){
    alpha <- matrix(data=mu, nrow=Q, ncol=1)
  } else
    alpha <- matrix(data=0, nrow=1, ncol=1)

  ## incidence matrices of the genetic predictors
  levels.inds <- rownames(X)
  inds <- rep(NA, N)
  if(Q > 1){
    for(year in levels.years)
      inds[years == year] <- levels.inds[1:sum(years == year)]
  } else
    inds <- levels.inds
  inds <- as.factor(inds)
  Z <- stats::model.matrix(~ inds - 1)
  if(enforce.zhou){
    X.c <- scale(x=X, center=TRUE, scale=FALSE)
    K <- tcrossprod(X.c, X.c) / P
  } else{
    X.c <- X - 1
    K <- estimGenRel(X=X, relationships="additive", method="vanraden1",
                     verbose=0)
  }

  ## hyper-parameters
  s.a <- (1 / (N*P)) * sum(colSums(X.c^2))
  s.b <- (1/N)  * sum(diag(K))
  if(is.null(pi))
    pi <- exp(stats::runif(n=1, min=log(1/P), max=log(1)))
  if(is.null(h))
    h <- stats::runif(n=1, min=0, max=1)
  if(is.null(rho)){
    if(pi == 0){
      rho <- 0
    } else
      rho <- stats::runif(n=1, min=0, max=1)
  }
  sigma.betat2 <- (h * rho) / ((1 - h) * P * pi * s.a) # -> sigma_a^2 in paper
  if(is.nan(sigma.betat2))
    sigma.betat2 <- 0
  sigma.u2 <- (h * (1 - rho)) / ((1 - h) * s.b) # -> sigma_b^2 in paper
  if(is.null(tau))
    ## tau <- rgamma(n=1, shape=1, rate=0.5)
    tau <- 1

  ## sparse genetic effects
  betat <- stats::setNames(object=rep(0, P), nm=colnames(X))
  gamma <- stats::setNames(object=stats::rbinom(n=P, size=1, prob=pi), nm=colnames(X))
  betat[gamma == 1] <- stats::rnorm(n=sum(gamma == 1), mean=0,
                             sd=sqrt(sigma.betat2 * tau^(-1)))

  ## polygenic effects
  u <- stats::setNames(object=MASS::mvrnorm(n=1, mu=rep(0, I),
                                     Sigma=sigma.u2 * tau^(-1) * K),
                nm=rownames(X))

  ## errors
  if(is.infinite(err.df)){
    epsilon <- matrix(stats::rnorm(n=N, mean=0, sd=sqrt(tau^(-1))))
  } else
    epsilon <- matrix(stats::rt(n=N, df=err.df, ncp=0))

  ## phenotypes
  y <- W %*% alpha + Z %*% X.c %*% betat + Z %*% u + epsilon
  if(perc.NA > 0){
    idx <- sample.int(n=N, size=floor(perc.NA/100 * N))
    y[idx] <- NA
  }

  return(list(y=y,
              W=W, alpha=alpha,
              Z=Z, X.c=X.c, s.a=s.a, s.b=s.b, K=K,
              pi=pi, h=h, rho=rho,
              sigma.betat2=sigma.betat2, sigma.u2=sigma.u2, sigma2=1/tau,
              betat=betat, u=u))
}

##' Launch GEMMA
##'
##' See Zhou & Stephens (Nature Genetics, 2012) and Zhou et al (PLoS Genetics, 2013).
##' @param model name of the model to fit (ulmm/bslmm)
##' @param y N-vector of phenotypes
##' @param X NxP matrix of SNP genotypes encoded in number of copies of the 2nd allele, i.e. as allele doses in {0,1,2}, with individuals in rows and SNPs in columns
##' @param snp.coords data.frame with 3 columns (snp, coord, chr)
##' @param alleles data.frame with SNPs in rows (names as row names) and alleles in columns (first is "minor", second is "major")
##' @param maf SNPs with minor allele frequency strictly below this threshold will be discarded
##' @param K.c NxN kinship matrix; if NULL, will be estimated using X via \code{\link{estimGenRel}} with relationships="additive" and method="zhou"
##' @param W NxQ matrix of covariates with individuals in rows (names as row names), a first column of 1 and a second column of covariates values
##' @param out.dir directory in which the output files will be saved
##' @param task.id identifier of the task (used in temporary and output file names)
##' @param verbose verbosity level (0/1)
##' @param clean remove files: none, some (temporary only), all (temporary and results)
##' @param seed seed for the generator of pseudo-random numbers
##' @param burnin number of iterations to discard as burn-in (if model="bslmm")
##' @param nb.iters number of iterations (if model="bslmm")
##' @param thin thining (if model="bslmm")
##' @return invisible list
##' @author Timothee Flutre [aut,cre], Dalel Ahmed [ctb]
##' @examples
##' \dontrun{## simulate genotypes
##' set.seed(1859)
##' I <- 200
##' P <- 2000
##' X <- simulGenosDose(nb.genos=I, nb.snps=P)
##'
##' ## simulate phenotypes
##' model <- simulBvsr(Q=1, X=X, pi=0.01, pve.A=0.7, sigma.a2=1)
##'
##' ## test SNPs one by one with the univariate LMM
##' snp.coords <- data.frame(snp=colnames(X), coord=1:ncol(X),
##'                          chr="chr1", stringsAsFactors=FALSE)
##' alleles <- data.frame(minor=rep("a", ncol(X)),
##'                       major="A", stringsAsFactors=FALSE)
##' rownames(alleles) <- colnames(X)
##' fit.u <- gemma(model="ulmm", model$Y[,1], X, snp.coords, alleles,
##'                W=model$W, out.dir=tempdir(), clean="all")
##' cor(model$a[model$gamma == 1], fit.u$tests$beta[model$gamma == 1])
##' cols <- rep("black",ncol(X)); cols[model$gamma==1] <- "red"
##' qqplotPval(fit.u$tests$p_wald, col=cols)
##'
##' ## fit all SNPs jointly with the BSLMM
##' burnin <- 10^3
##' nb.iters <- 10^4
##' thin <- 10^2
##' fit.bs <- gemma(model="bslmm", model$Y[,1], X, snp.coords, alleles,
##'                 W=model$W, out.dir=tempdir(), clean="all",
##'                 burnin=burnin, nb.iters=nb.iters, thin=thin)
##' posterior.samples <- coda::mcmc(data=fit.bs$hyperparams, start=burnin + 1,
##'                                 end=burnin + nb.iters, thin=thin)
##' summary(posterior.samples)
##' }
##' @export
gemma <- function(model="ulmm", y, X, snp.coords, alleles, maf=0.01, K.c=NULL,
                  W, out.dir=getwd(), task.id="gemma", verbose=1, clean="none",
                  seed=1859, burnin=1000, nb.iters=7000, thin=10){
  stopifnot(file.exists(Sys.which("gemma")),
            model %in% c("ulmm", "bslmm"))
  if(is.matrix(y)){
    stopifnot(ncol(y) == 1)
    y <- as.vector(y)
  }
  stopifnot(is.vector(y),
            .isValidGenosDose(X, check.rown=FALSE),
            is.data.frame(snp.coords),
            ncol(snp.coords) == 3,
            ! is.null(row.names(alleles)),
            colnames(alleles) == c("minor","major"),
            all(colnames(snp.coords) == c("snp", "coord", "chr")),
            all(snp.coords$snp == colnames(X)),
            all(snp.coords$snp == rownames(alleles)),
            is.numeric(maf),
            length(maf) == 1,
            maf >= 0,
            maf <= 1,
            is.matrix(W),
            all(W[,1] == 1),
            file.exists(out.dir),
            is.character(task.id),
            length(task.id) == 1,
            clean %in% c("none", "some", "all"))

  output <- list()

  ## prepare input files
  tmp.files <- c()
  tmp.files <- c(tmp.files,
                 genos=paste0(out.dir, "/genos_", task.id, ".txt.gz"))
  dose2bimbam(X=X, alleles=alleles,
              file=gzfile(tmp.files["genos"]))
  tmp.files <- c(tmp.files,
                 snp.coords=paste0(out.dir, "/snp-coords_", task.id, ".txt"))
  utils::write.table(x=snp.coords,
              file=tmp.files["snp.coords"],
              quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  if(is.null(K.c))
    K.c <- estimGenRel(X=X, thresh=maf, relationships="additive",
                       method="zhou", verbose=verbose)
  tmp.files <- c(tmp.files,
                 kinship.center=paste0(out.dir, "/kinship-center_", task.id,
                                       ".txt"))
  utils::write.table(x=K.c,
              file=tmp.files["kinship.center"],
              quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  tmp.files <- c(tmp.files,
                 covars=paste0(out.dir, "/covars_", task.id, ".txt"))
  utils::write.table(x=W,
              file=tmp.files["covars"],
              quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  tmp.files <- c(tmp.files,
                 phenos=paste0(out.dir, "/phenos_", task.id, ".txt.gz"))
  utils::write.table(x=y,
              file=gzfile(tmp.files["phenos"]),
              quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

  ## prepare cmd-line and execute it
  cmd <- paste0("cd ", out.dir,
                "; gemma",
                " -g ", tmp.files["genos"],
                " -maf ", maf,
                " -p ", tmp.files["phenos"],
                " -a ", tmp.files["snp.coords"],
                " -outdir ./",
                " -o results_", task.id)
  if(model == "ulmm"){
    cmd <- paste0(cmd, " -k ", tmp.files["kinship.center"],
                  " -c ", tmp.files["covars"],
                  " -lmm 4")
  } else if(model == "bslmm"){
    cmd <- paste0(cmd, " -bslmm 1",
                  " -w ", format(x=burnin, scientific=FALSE),
                  " -s ", format(x=nb.iters, scientific=FALSE),
                  " -rpace ", thin,
                  " -seed ", seed)
  }
  if(verbose <= 0){
    tmp.files <- c(tmp.files,
                   stdoe=paste0(out.dir, "/stdouterr_gemma_", task.id, ".txt"))
    cmd <- paste0(cmd, " > ", tmp.files["stdoe"], " 2>&1")
  }
  if(verbose > 0)
    write(cmd, stdout())
  system(cmd)
  output[["cmd"]] <- cmd

  ## load output files
  f <- paste0(out.dir, "/results_", task.id, ".log.txt")
  output[["log"]] <- readLines(f)
  if(clean == "all")
    file.remove(f)
  if(model == "ulmm"){
    f <- paste0(out.dir, "/results_", task.id, ".assoc.txt")
    tmp <- utils::read.table(file=f, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    rownames(tmp) <- tmp$rs
    output[["tests"]] <- tmp
    if(clean == "all")
      file.remove(f)
  } else if(model == "bslmm"){
    f <- paste0(out.dir, "/results_", task.id, ".hyp.txt")
    output[["hyperparams"]] <- utils::read.table(file=f,
                                          sep="\t",
                                          skip=1, stringsAsFactors=FALSE,
                                          header=FALSE,
                                          col.names=c("h", "pve", "rho", "pge",
                                                      "pi", "n_gamma", ""))
    output[["hyperparams"]][7] <- NULL
    if(clean == "all")
      file.remove(f)
    f <- paste0(out.dir, "/results_", task.id, ".param.txt")
    output[["params"]] <- utils::read.table(file=f,
                                     sep="\t", skip=1,
                                     stringsAsFactors=FALSE, header=FALSE,
                                     col.names=c("chr", "rs", "ps", "n_miss",
                                                 "alpha", "beta", "gamma"))
    if(clean == "all")
      file.remove(f)
  }

  ## clean temporary files
  if(clean != "none")
    for(f in tmp.files)
      file.remove(f)

  invisible(output)
}

##' GEMMA uLMM per chromosome
##'
##' See Zhou & Stephens (Nature Genetics, 2012).
##' @param y vector of phenotypes
##' @param X matrix of SNP genotypes encoded in number of copies of the 2nd allele, i.e. as allele doses in {0,1,2}, with individuals in rows and SNPs in columns
##' @param snp.coords data.frame with SNP identifiers as row names, and two columns, "chr" and "coord" or "pos"
##' @param alleles data.frame with SNPs in rows (names as row names) and alleles in columns (first is "minor", second is "major"); will be made of A's and B's if missing
##' @param chr.ids set of chromosome identifiers to analyze (optional, all by default)
##' @param W matrix of covariates with individuals in rows (names as row names), a first column of 1 and a second column of covariates values
##' @param out.dir directory in which the data files will be found
##' @param task.id identifier of the task (used in output file names)
##' @param clean remove files: none, some (temporary only), all (temporary and results)
##' @param verbose verbosity level (0/1)
##' @return a data.frame with GEMMA's output for all chromosomes
##' @author Timothee Flutre [aut,cre], Dalel Ahmed [ctb]
##' @export
gemmaUlmmPerChr <- function(y, X, snp.coords, alleles=NULL, chr.ids=NULL, W,
                            out.dir, task.id="", clean="none", verbose=1){
  stopifnot(is.vector(y),
            .isValidGenosDose(X, check.rown=FALSE),
            .isValidSnpCoords(snp.coords),
            all(rownames(snp.coords) %in% colnames(X)),
            all(colnames(X) %in% rownames(snp.coords)),
            is.matrix(W),
            all(W[,1] == 1),
            file.exists(out.dir),
            is.character(task.id),
            length(task.id) == 1,
            clean %in% c("none", "some", "all"))
  if(! is.null(alleles))
    stopifnot(is.data.frame(alleles),
              ! is.null(rownames(alleles)),
              colnames(alleles) == c("minor","major"),
              all(colnames(X) %in% rownames(alleles)),
              all(rownames(alleles) %in% colnames(X)))

  out <- list()

  snp.coords <- snp.coords[colnames(X),]
  if(! "coord" %in% colnames(snp.coords))
    colnames(snp.coords)[colnames(snp.coords) == "pos"] <- "coord"
  if(is.null(alleles)){
    alleles <- data.frame(minor=rep("A", nrow(snp.coords)),
                          major=rep("B", nrow(snp.coords)),
                          stringsAsFactors=FALSE)
    rownames(alleles) <- rownames(snp.coords)
  } else
    alleles <- alleles[colnames(X),]
  sc.gemma <- data.frame(snp=rownames(snp.coords),
                         coord=snp.coords$coord,
                         chr=snp.coords$chr,
                         stringsAsFactors=FALSE)

  if(is.null(chr.ids))
    chr.ids <- sort(unique(as.character(snp.coords$chr)))

  for(chr.id in chr.ids){
    subset.snp.ids <- (snp.coords$chr == chr.id)
    if(verbose > 0)
      message(paste0(chr.id, ": ", sum(subset.snp.ids), " SNPs"))
    K.c <- estimGenRel(X=X[, ! subset.snp.ids], method="zhou")
    out[[chr.id]] <- gemma(model="ulmm",
                           y=y,
                           X=X[, subset.snp.ids],
                           snp.coords=sc.gemma[subset.snp.ids,],
                           alleles=alleles[subset.snp.ids,],
                           K.c=K.c, W=W, out.dir=out.dir,
                           task.id=paste0(task.id, "-", chr.id),
                           verbose=verbose-1, clean=clean)
  }
  out <- do.call(rbind, out)
  rownames(out) <- out$rs

  return(out)
}

##' QTLRel per chromosome
##'
##' Given I individuals, P SNPs, Q covariates (e.g. replicates) and N=I*Q phenotypes, the whole likelihood is: y = W alpha + Z X beta + epsilon, where y is Nx1, W is NxQ, Z is NxI, X is IxP and u = X beta.
##' The QTLRel package (Cheng et al, BMC Genetics, 2011) decomposes the inference into an \emph{ad hoc} procedure of two steps.
##' First, the variance components and fixed effects are estimated: y = W alpha + Z u + epsilon.
##' Second the allele effects are tested SNP per SNP: for all p in {1,..,P}, y = W alpha + Z x_p beta_p + Z u + epsilon, using the fixed effects and variance components estimated in the first step.
##' As the SNPs can be in linkage disequilibrium, it is advised (Yang et al, Nature Genetics, 2014) to perform the procedure chromosome per chromosome, which is the goal of this function; but I advise to also run it with chr.ids=NULL.
##' @param y vector or one-column matrix of phenotypes (the order is important and should be in agreement with the other arguments X, W and Z)
##' @param X matrix of SNP genotypes encoded in number of copies of the 2nd allele, i.e. as allele doses in {0,1,2}, with individuals in rows and SNPs in columns
##' @param snp.coords data.frame with SNP identifiers as row names and two columns named "chr" and "coord" (or "pos")
##' @param thresh threshold on minor allele frequencies below which SNPs are ignored via \code{\link{discardSnpsLowMaf}} (default=0.01; NULL to skip this step; SNPs are ignored for testing, but all are still used to calculate the matrix of additive genetic relationships as rare SNPs are informative in this purpose)
##' @param chr.ids vector of chromosome identifiers to analyze (if NULL, the regular QTLRel procedure is launched, i.e. all chromosomes are used to estimate the variance components)
##' @param W incidence matrix of covariates (should not contain the column of 1's for the intercept; use \code{\link[stats]{model.matrix}} if necessary)
##' @param Z incidence matrix relating phenotypes to individuals (if nrow(y) and nrow(X) are different, diagonal otherwise; use \code{\link[stats]{model.matrix}} if necessary)
##' @param method.A method to estimate the additive relationships (see \code{\link{estimGenRel}})
##' @param verbose verbosity level (0/1)
##' @return a list with three data.frames as components, variance.components, fixed.effects and scan
##' @author Timothee Flutre
##' @examples
##' if(require(QTLRel)){
##'   set.seed(1859)
##'   I <- 200
##'   P <- 2000
##'   Q <- 3
##'   N <- Q * I
##'   dat <- data.frame(ind=as.factor(rep(paste0("ind", 1:I), times=Q)),
##'                     year=as.factor(rep(paste0(2003:(2003+Q-1)), each=I)))
##'   W <- stats::model.matrix(~ year, dat)
##'   alpha <- rnorm(n=Q, mean=50, sd=30)
##'   X <- simulGenosDose(nb.genos=I, nb.snps=P)
##'   beta <- rnorm(n=P, mean=0, sd=2)
##'   Z <- stats::model.matrix(~ ind - 1, dat)
##'   dat$response <- as.vector(W %*% alpha + Z %*% X %*% beta + rnorm(N))
##'   snp.coords <- data.frame(chr=sample(paste0("chr",1:10), P, TRUE),
##'                            coord=sample.int(10^6, P))
##'   rownames(snp.coords) <- colnames(X)
##'   res <- qtlrelPerChr(dat$response, X, snp.coords, 0.01, "chr1", W=W[,-1], Z=Z)
##'   if(interactive())
##'      res <- qtlrelPerChr(dat$response, X, snp.coords, 0.01, NULL, W=W[,-1], Z=Z)
##'  }
##' @export
qtlrelPerChr <- function(y, X, snp.coords, thresh=0.01, chr.ids=NULL, W=NULL, Z=NULL,
                         method.A="vanraden1", verbose=1){
  requireNamespaces("QTLRel")
  if(is.matrix(y))
    stopifnot(ncol(y) == 1)
  y <- as.vector(y)
  stopifnot(.isValidGenosDose(X, check.rown=FALSE),
            .isValidSnpCoords(snp.coords),
            all(rownames(snp.coords) %in% colnames(X)),
            all(colnames(X) %in% rownames(snp.coords)))
  if(! is.null(W)){
    stopifnot(is.matrix(W),
              nrow(W) == length(y),
              ! all(W[,1] == 1))
  }
  if(! is.null(Z)){
    stopifnot(is.matrix(Z),
              nrow(Z) == length(y),
              ncol(Z) == nrow(X))
  } else{
    stopifnot(length(y) == nrow(X))
    Z <- diag(length(y))
  }

  out <- list()

  if(! "coord" %in% colnames(snp.coords))
    colnames(snp.coords)[colnames(snp.coords) == "pos"] <- "coord"
  snp.coords$chr <- as.character(snp.coords$chr)
  X <- X[, rownames(snp.coords)]

  ## discard SNPs with low MAFs
  if(! is.null(thresh)){
    X <- discardSnpsLowMaf(X=X, thresh=thresh, verbose=verbose)
    snp.coords <- snp.coords[colnames(X),]
  }

  ## infer with all chromosomes together
  if(is.null(chr.ids)){
    if(verbose > 0){
      msg <- paste0("estimate additive genetic relationships (", ncol(X),
                    " SNPs, method=", method.A, ") ...")
      write(msg, stdout())
    }
    A.mark <- estimGenRel(X=X, thresh=0, relationships="additive",
                          method=method.A, verbose=verbose - 1)
    if(is.null(W)){
      if(verbose > 0)
        write("estimate variance components ...", stdout())
      vc <- QTLRel::estVC(y=y,
                          v=list(AA=Z %*% A.mark %*% t(Z), DD=NULL, HH=NULL,
                                 AD=NULL, MH=NULL, EE=diag(length(y))))
      if(verbose > 0){
        msg <- paste0("test alleles' effect (", ncol(X), " SNPs) ...")
        write(msg, stdout())
      }
      scan <- QTLRel::scanOne(y=y,
                              gdat=Z %*% (X - 1),
                              vc=vc,
                              test="F",
                              numGeno=TRUE)
    } else{
      if(verbose > 0)
        write("estimate variance components ...", stdout())
      vc <- QTLRel::estVC(y=y,
                          x=W,
                          v=list(AA=Z %*% A.mark %*% t(Z), DD=NULL, HH=NULL,
                                 AD=NULL, MH=NULL, EE=diag(length(y))))
      if(verbose > 0){
        msg <- paste0("test alleles' effect (", ncol(X), " SNPs) ...")
        write(msg, stdout())
      }
      scan <- QTLRel::scanOne(y=y,
                              x=W,
                              gdat=Z %*% (X - 1),
                              vc=vc,
                              test="F",
                              numGeno=TRUE)
    }
    out$variance.components <- vc$par[c("AA", "EE")]
    out$fixed.effects <- vc$par[1:(length(vc$par)-2)]
    out$scan <- data.frame(chr=snp.coords$chr,
                           coord=snp.coords$coord,
                           pvalue=scan$p,
                           pve=scan$v / 100,
                           stringsAsFactors=FALSE)
    out$scan <- cbind(out$scan,
                      do.call(rbind, scan$parameters))
    rownames(out$scan) <- names(scan$p)
  } else{ ## infer chromosome per chromosome
    stopifnot(all(chr.ids %in% unique(snp.coords$chr)))
    for(chr.id in chr.ids){
      subset.snp.ids <- (snp.coords$chr == chr.id)
      if(verbose > 0)
        message(paste0(chr.id, ": ", sum(subset.snp.ids), " SNPs"))
      if(verbose > 0){
        msg <- paste0("estimate additive genetic relationships (",
                      ncol(X[, ! subset.snp.ids]), " SNPs, method=",
                      method.A, ") ...")
        write(msg, stdout())
      }
      A.mark <- estimGenRel(X=X[, ! subset.snp.ids], relationships="additive",
                            method=method.A, verbose=verbose - 1)
      if(is.null(W)){
        if(verbose > 0)
          write("estimate variance components ...", stdout())
        vc <- QTLRel::estVC(y=y,
                            v=list(AA=Z %*% A.mark %*% t(Z), DD=NULL, HH=NULL,
                                   AD=NULL, MH=NULL, EE=diag(length(y))))
        if(verbose > 0){
          msg <- paste0("test alleles' effect (", ncol(X[, subset.snp.ids]),
                        " SNPs) ...")
          write(msg, stdout())
        }
        scan <- QTLRel::scanOne(y=y,
                                gdat=Z %*% (X[, subset.snp.ids] - 1),
                                vc=vc,
                                test="F",
                                numGeno=TRUE)
      } else{
        if(verbose > 0)
          write("estimate variance components ...", stdout())
        vc <- QTLRel::estVC(y=y,
                            x=W,
                            v=list(AA=Z %*% A.mark %*% t(Z), DD=NULL, HH=NULL,
                                   AD=NULL, MH=NULL, EE=diag(length(y))))
        if(verbose > 0){
          msg <- paste0("test alleles' effect (", ncol(X[, subset.snp.ids]),
                        " SNPs) ...")
          write(msg, stdout())
        }
        scan <- QTLRel::scanOne(y=y,
                                x=W,
                                gdat=Z %*% (X[, subset.snp.ids] - 1),
                                vc=vc,
                                test="F",
                                numGeno=TRUE)
      }
      out$variance.components[[chr.id]] <- vc$par[c("AA", "EE")]
      out$fixed.effects[[chr.id]] <- vc$par[1:(length(vc$par)-2)]
      out$scan[[chr.id]] <- data.frame(snp=names(scan$p),
                                       chr=chr.id,
                                       coord=snp.coords$coord[subset.snp.ids],
                                       pvalue=scan$p,
                                       pve=scan$v / 100,
                                       stringsAsFactors=FALSE)
      out$scan[[chr.id]] <- cbind(out$scan[[chr.id]],
                                  do.call(rbind, scan$parameters))
      colnames(out$scan[[chr.id]])[(ncol(out$scan[[chr.id]]) -
                                    length(scan$parameters[[1]]) + 1):
                                   ncol(out$scan[[chr.id]])] <-
        paste0("coef", 1:length(scan$parameters[[1]]))
    } # end of for loop over chromosomes
    out$variance.components <- do.call(rbind, out$variance.components)
    rownames(out$variance.components) <- chr.ids
    out$fixed.effects <- do.call(rbind, out$fixed.effects)
    rownames(out$fixed.effects) <- chr.ids
    out$scan <- do.call(rbind, out$scan)
    rownames(out$scan) <- out$scan$snp
  } # end of chromosome-by-chromosome inference

  return(out)
}

##' Asymptotic Bayes factor
##'
##' Calculate the asymptotic Bayes factor proposed by Wakefield in Genetic Epidemiology 33:79-86 (2009, \url{http://dx.doi.org/10.1002/gepi.20359}).
##' @param theta.hat MLE of the additive genetic effect
##' @param V variance of theta.hat
##' @param W variance of the prior on theta
##' @param log10 return the log10 of the ABF
##' @return numeric
##' @author Timothee Flutre
##' @export
calcAsymptoticBayesFactorWakefield <- function(theta.hat, V, W, log10=TRUE){
  z2 <- theta.hat^2 / V # Wald statistic

  log10.ABF <- 0.5 * log10(V) - 0.5 * log10(V + W) +
    (0.5 * z2 * W / (V + W)) / log(10)

  if(log10)
    return(log10.ABF)
  else
    return(10^log10.ABF)
}

##' Exact Bayes factor
##'
##' Calculate the exact Bayes factor proposed by Servin and Stephens in PLoS Genetics 3,7 (2007, \url{http://dx.doi.org/10.1371/journal.pgen.0030114}).
##' @param G vector of genotypes
##' @param Y vector of phenotypes
##' @param sigma.a variance of the prior on the additive genetic effect
##' @param sigma.d variance of the prior on the dominance genetic effect
##' @param log10 return the log10 of the ABF
##' @return numeric
##' @author Bertrand Servin [aut], Timothee Flutre [ctb,cre]
##' @export
calcExactBayesFactorServinStephens <- function(G, Y, sigma.a, sigma.d,
                                               log10=TRUE){
  stopifnot(is.vector(G), is.vector(Y))

  subset <- stats::complete.cases(Y) & stats::complete.cases(G)
  Y <- Y[subset]
  G <- G[subset]
  stopifnot(length(Y) == length(G))

  N <- length(G)
  X <- cbind(rep(1,N), G, G == 1)
  inv.Sigma.B <- diag(c(0, 1/sigma.a^2, 1/sigma.d^2))
  inv.Omega <- inv.Sigma.B + t(X) %*% X
  inv.Omega0 <- N
  tY.Y <- t(Y) %*% Y
  log10.BF <- as.numeric(0.5 * log10(inv.Omega0) -
                           0.5 * log10(det(inv.Omega)) -
                           log10(sigma.a) - log10(sigma.d) -
                           (N/2) * (log10(tY.Y - t(Y) %*% X %*% solve(inv.Omega)
                                          %*% t(X) %*% cbind(Y)) -
                                    log10(tY.Y - N*mean(Y)^2)))

  if(log10)
    return(log10.BF)
  else
    return(10^log10.BF)
}

##' Approximate Bayes factor
##'
##' Calculate the log10(ABF) of Wen & Stephens in Annals of Applied Statistics (2014, \url{http://dx.doi.org/10.1214/13-AOAS695}) according to the "exchangeable standardized effects" model.
##' @param sstats matrix of summary statistics with one row per subgroup and three columns, "bhat", "sebhat" and "t"
##' @param phi2 prior variance of the \eqn{b_s} given \eqn{\bar{b}}; controls the prior expected degree of heterogeneity among subgroups
##' @param oma2 prior variance of \eqn{\bar{b}}; controls the prior expected size of the average effect across subgroups
##' @return numeric
##' @author Xiaoquan Wen [aut], Timothee Flutre [ctb,cre]
##' @export
calcL10ApproximateBayesFactorWenStephens <- function(sstats, phi2, oma2){
  stopifnot(is.matrix(sstats),
            all(colnames(sstats) == c("bhat","sebhat","t")),
            is.numeric(phi2), length(phi2) == 1,
            is.numeric(oma2), length(oma2) == 1)

  l10abf <- NA

  bbarhat.num <- 0
  bbarhat.denom <- 0
  varbbarhat <- 0
  l10abfs.single <- c()
  for(i in 1:nrow(sstats)){ # for each subgroup
    if(sum(is.na(sstats[i,])) == length(sstats[i,]))
      next
    bhat <- sstats[i,"bhat"]
    varbhat <- sstats[i,"sebhat"]^2
    t <- sstats[i,"t"]
    if(abs(t) > 10^(-8)){
      bbarhat.num <- bbarhat.num + bhat / (varbhat + phi2)
      bbarhat.denom <- bbarhat.denom + 1 / (varbhat + phi2)
      varbbarhat <- varbbarhat + 1 / (varbhat + phi2)
      tmp <- 0.5 * log10(varbhat) -
        0.5 * log10(varbhat + phi2) +
          (0.5 * t^2 * phi2 / (varbhat + phi2)) / log(10)
    } else
      tmp <- 0
    l10abfs.single <- c(l10abfs.single, tmp)
  }

  if(length(l10abfs.single) != 0){
    bbarhat <- ifelse(bbarhat.denom != 0, bbarhat.num / bbarhat.denom, 0)
    varbbarhat <- ifelse(varbbarhat != 0, 1 / varbbarhat, Inf)
    if(bbarhat != 0 & ! is.infinite(varbbarhat)){
      T2 <- bbarhat^2 / varbbarhat
      l10abf.bbar <- ifelse(T2 != 0,
                            0.5 * log10(varbbarhat) -
                              0.5 * log10(varbbarhat + oma2) +
                              (0.5 * T2 * oma2 / (varbbarhat + oma2)) /
                              log(10),
                            0)
      l10abf <- l10abf.bbar
      for(i in 1:length(l10abfs.single))
        l10abf <- l10abf + l10abfs.single[i]
    } else
      l10abf <- 0
  }

  return(as.numeric(l10abf))
}

##' Returns the genetic map contained in a BioMercator TXT file.
##'
##' http://moulon.inra.fr/index.php/en/tranverse-team/atelier-de-bioinformatique/projects/projets/135
##' @param file the name of the file which the data are to be read from
##' @return list
##' @author Timothee Flutre
##' @export
readBiomercator <- function(file){
    stopifnot(file.exists(file))

    gmap <- list()

    lines <- readLines(file)

    ## load meta-data
    i <- 1
    while(! grepl("chr=", lines[i])){
        tokens <- strsplit(lines[i], "=")[[1]]
        key <- gsub(" ", "", tokens[1])
        if(key %in% names(gmap))
            key <- paste(key, sum(key == names(gmap)) + 1, sep=".")
        val <- tokens[2]
        gmap[[key]] <- val
        i <- i + 1
    }

    ## load chromosomes (can have several linkage groups)
    chrs <- split(lines[-(1:(i-1))], cumsum(grepl("^chr=", lines[-(1:(i-1))])))
    gmap[["map"]] <- lapply(chrs, function(chr){
        lgs <- split(chr[-1], cumsum(grepl("^lg=", chr[-1])))
        data <- lapply(lgs, function(lg){
            tmp <- as.data.frame(do.call(rbind, strsplit(lg[-1], "\t")))
            tmp <- tmp[,-1] # discard column of identifiers
            tmp <- tmp[! duplicated(tmp[,1]),] # discard redundant marker names
            df <- as.numeric(as.character(tmp[,2]))
            names(df) <- as.character(tmp[,1])
            df[order(df)]
        })
        names(data) <- paste("lg",
                             sapply(lgs, function(lg){
                                 strsplit(lg[1], "=")[[1]][2]
                             }),
                             sep=".")
        data
    })
    names(gmap[["map"]]) <- paste("chr",
                                  sapply(chrs, function(chr){
                                      strsplit(chr[1], "=")[[1]][2]
                                  }),
                                  sep=".")

    ## add useful numbers
    gmap$nbMarkers <- sum(sapply(gmap$map, function(chr){
        sapply(chr, function(lg){length(lg)})
    }))
    gmap$nbLinkageGroups <- sum(sapply(gmap$map, function(chr){length(chr)}))
    gmap$mapSize <- sum(sapply(gmap$map, function(chr){
        sapply(chr, function(lg){lg[length(lg)]})
    }))

    txt <- paste0("map '", gmap$mapName, "':")
    txt <- paste0(txt, "\n\tnb of individuals: ", gmap$popSize)
    txt <- paste0(txt, "\n\tnb of markers: ", gmap$nbMarkers)
    txt <- paste0(txt, "\n\tnb of chromosomes: ", length(gmap$map))
    txt <- paste0(txt, "\n\tnb of linkage groups: ", gmap$nbLinkageGroups)
    txt <- paste0(txt, "\n\tmap size: ", gmap$mapSize, " ", gmap$mapUnit)
    mean.dist <- do.call(c, lapply(gmap$map, function(chr){
        sapply(chr, function(lg){
            rev(lg)[1:(length(lg)-1)] - lg[(length(lg)-1):1]
        })}))
    txt <- paste0(txt, "\n\tmean distances per linkage group:")
    message(paste0(txt))
    print(summary(mean.dist))

    return(gmap)
}

##' Plot pedigree
##'
##' Plot a pedigree using the "igraph" package.
##' This function was inspired by plot.pedigree() from the "synbreed" package (under GPL-3).
##' It add options for monoecious species and auto-fecondation.
##' @param inds identifiers of the individuals
##' @param mothers identifiers of their mother; can be NA
##' @param fathers identifiers of their father; can be NA
##' @param generations should start at 0
##' @param sexes "F" for female (circle), "M" for male (square) and "H" for hermaphrodite (triangle); can also be NA (no shape)
##' @param edge.col.mother see ?igraph.plotting
##' @param edge.col.father see ?igraph.plotting
##' @param vertex.label.color see ?igraph.plotting
##' @param vertex.color see ?igraph.plotting
##' @param vertex.size see ?igraph.plotting
##' @param vertex.shape see ?igraph.plotting
##' @param vertex.label.family see ?igraph.plotting
##' @param mult.edge.curve see ?igraph.plotting
##' @param edge.arrow.width see ?igraph.plotting
##' @param edge.arrow.size see ?igraph.plotting
##' @param ... other plotting options; see ?plot.igraph and ?igraph.plotting
##' @return invisible pedigree as an "igraph" object
##' @author Timothee Flutre
##' @export
plotPedigree <- function(inds, mothers, fathers, generations, sexes=NULL,
                         edge.col.mother="black", edge.col.father="darkgrey",
                         vertex.label.color="darkblue", vertex.color="white",
                         vertex.size=20, vertex.shape="none",
                         vertex.label.family="Helvetica", mult.edge.curve=0.25,
                         edge.arrow.width=0.75, edge.arrow.size=0.75,
                         ...){
  stopifnot(is.vector(inds),
            is.vector(mothers),
            is.vector(fathers),
            is.vector(generations))
  if(! is.null(sexes))
    stopifnot(is.vector(sexes))

  ## check inds
  inds <- as.character(inds)
  if(length(unique(inds)) != length(inds)){
    msg <- paste0(length(inds), " individuals, but only ",
                  length(unique(inds)), " unique")
    stop(msg)
  }

  ## check mothers
  mothers <- as.character(mothers)
  stopifnot(length(mothers) == length(inds),
            all(mothers[! is.na(mothers)] %in% inds))

  ## check fathers
  fathers <- as.character(fathers)
  stopifnot(length(fathers) == length(inds),
            all(fathers[! is.na(fathers)] %in% inds))

  ## check generations
  generations <- as.numeric(generations)
  stopifnot(length(generations) == length(inds))
  generations <- generations - min(generations)

  ## check sexes
  if(! is.null(sexes)){
    sexes <- as.character(sexes)
    stopifnot(length(sexes) == length(inds))
    sexes <- toupper(sexes)
    uniq.sex <- unique(sexes)
    if(any(! is.na(uniq.sex)))
      stopifnot(all(uniq.sex[! is.na(uniq.sex)] %in% c("F","M","H")))
  }

  ## sort according to generations
  idx <- order(generations)
  inds <- inds[idx]
  mothers <- mothers[idx]
  fathers <- fathers[idx]
  generations <- generations[idx]
  if(! is.null(sexes))
    sexes <- sexes[idx]

  ## add "triangle" as shape
  mytriangle <- function(coords, v=NULL, params) {
    vertex.color <- params("vertex", "color")
    if(length(vertex.color) != 1 && !is.null(v))
      vertex.color <- vertex.color[v]
    vertex.size <- 1/200 * params("vertex", "size")
    if(length(vertex.size) != 1 && !is.null(v))
      vertex.size <- vertex.size[v]
    graphics::symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
            stars=cbind(vertex.size, vertex.size, vertex.size),
            add=TRUE, inches=FALSE)
  }
  igraph::add_shape("triangle", clip=igraph::shapes("circle")$clip,
                    plot=mytriangle)

  ## make "igraph" object
  relations <- rbind(cbind(mothers[! is.na(mothers)], inds[! is.na(mothers)]),
                     cbind(fathers[! is.na(fathers)], inds[! is.na(fathers)]))
  colnames(relations) <- c("parent","child")
  v.df <- data.frame(ind=inds,
                     mother=mothers,
                     father=fathers,
                     generation=generations,
                     stringsAsFactors=FALSE)
  if(is.null(sexes)){
    v.df[["label"]] <- inds
  } else{
    v.df[["sex"]] <- sexes
    v.df[["label"]] <- paste0(inds, "\n(", sexes, ")")
  }
  v.df[["shape"]] <- rep(vertex.shape, length(inds))
  v.df[["size"]] <- rep(vertex.size, length(inds))
  v.df[["color"]] <- vertex.color
  v.df[["label.color"]] <- vertex.label.color
  v.df[["label.family"]] <- vertex.label.family
  ped.graph <- igraph::graph_from_data_frame(d=relations,
                                             directed=TRUE,
                                             vertices=v.df)

  ## check multiplicity corresponding to auto-fecondation
  has.autof <- FALSE
  if(igraph::has.multiple(ped.graph)){
    has.autof <- TRUE
    stopifnot(all(igraph::count_multiple(ped.graph) %in% c(1,2)))
  }

  ## set plot coordinates for vertices
  coords <- matrix(data=NA, nrow=length(inds), ncol=2,
                   dimnames=list(inds, c("x", "y")))
  coords[,2] <- max(generations) - generations
  ## coords[,1] <- order(generations, partial=order(inds, decreasing=TRUE)) -
  ##   cumsum(c(0, table(generations)))[generations + 1]
  coords[,1] <- order(generations) -
    cumsum(c(0, table(generations)))[generations + 1]
  coords[nrow(coords):1,1] <- unlist(tapply(coords[,1], coords[,2], function(x){
    if(length(x) == 1){
      x <- 0
    } else
      x <- rev(scale(x))
    return(x)
  }))

  ## set edge color depending on parental sex
  nb.rel.mother <- sum(! is.na(mothers))
  nb.rel.father <- sum(! is.na(fathers))
  edge.cols <- c(rep(edge.col.mother, nb.rel.mother),
                 rep(edge.col.father, nb.rel.father))

  ## tune curvature in case of auto-fecondation
  edge.curvatures <- rep(0, igraph::ecount(ped.graph))
  if(has.autof){
    idx.mult <- which(igraph::count_multiple(ped.graph) > 1)
    stopifnot(length(idx.mult) %% 2 == 0)
    edge.curvatures[idx.mult[idx.mult <= nb.rel.mother]] <- mult.edge.curve
    edge.curvatures[idx.mult[idx.mult > nb.rel.mother]] <- - mult.edge.curve
  }

  ## plot, finally
  igraph::plot.igraph(x=ped.graph,
                      layout=coords,
                      edge.color=edge.cols,
                      edge.curved=edge.curvatures,
                      edge.arrow.width=edge.arrow.width,
                      edge.arrow.size=edge.arrow.size,
                      ...)

  invisible(ped.graph)
}
