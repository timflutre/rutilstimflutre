## Contains functions useful for quantitative genetics.

##' Convert SNP genotype data from alleles (say, "AA" and "AT") to minor allele
##' doses (here, 0 and 1 if "T" is the minor allele).
##'
##' Not particularly efficient, but at least it exists.
##' @param x data.frame with SNPs in rows and individuals in columns, the SNP
##' identifiers being in the first column
##' @param na.string a character to be interpreted as NA values
##' @param verbose verbosity level (0/default=1)
##' @return list of a matrix (allele doses, SNPs in columns and individuals in
##' rows) and a vector (minor alleles)
##' @author Timothee Flutre
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
    raw.genos <- unlist(x[p, -1])
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

##' Plot missing SNP genotypes as a grid.
##'
##' Data will be represented in black if missing, white otherwise.
##' @param x matrix with SNP genotypes as allele doses (NA if missing) with
##' SNPs in columns and individuals in rows
##' @param main an overall title for the plot (default="Missing genotypes")
##' @param xlab a title for the x axis (default="Individuals")
##' @param ylab a title for the y axis (default="SNPs")
##' @return nothing
##' @author Timothee Flutre
plotGridMissGenos <- function(x, main="Missing genotypes", xlab="Individuals",
                              ylab="SNPs"){
  if(ncol(x) < nrow(x))
    warning("did you put SNPs in columns and individuals in rows?")
  image(1:nrow(x), 1:ncol(x), is.na(x), col=c("white","black"),
        main=main, xlab=xlab, ylab=ylab)
}

##' Estimate minor allele frequencies of SNPs.
##'
##' Missing values should be encoded as NA.
##' @param X matrix of SNP genotypes encoded as allele doses, with SNPs in
##' columns and individuals in rows
##' @return vector
##' @author Timothee Flutre
estimMaf <- function(X){
  stopifnot(is.matrix(X),
            sum(X < 0, na.rm=TRUE) == 0,
            sum(X > 2, na.rm=TRUE) == 0)
  N <- nrow(X)
  P <- ncol(X)
  if(P < N)
    warning("input matrix doesn't seem to have SNPs in columns and individuals in rows")

  maf <- apply(X, 2, function(x){
    x <- x[complete.cases(x)]
    tmp <- sum(x) / (2 * length(x))
    ifelse(tmp <= 0.5, tmp, 1 - tmp)
  })

  return(maf)
}

##' Plot the histogram of the minor allele frequency per SNP
##'
##' Missing values (encoded as NA) are discarded.
##' @param X matrix of SNP genotypes encoded as allele doses, with SNPs in
##' columns and individuals in rows (optional if maf is not null)
##' @param maf vector of minor allele frequencies (optional if X is not null)
##' @param main string for the main title (default="")
##' @param xlim default=c(0,0.5)
##' @param col color for the bars
##' @param border color for the border of the bars
##' @param las see ?par (default=1)
##' @param breaks see ?hist (default="FD")
##' @param verbose verbosity level (0/default=1)
##' @param ... arguments to be passed to hist()
##' @return nothing
##' @author Timothee Flutre
plotHistMinAllelFreq <- function(X=NULL, maf=NULL, main="", xlim=c(0,0.5),
                                 col="grey", border="white", las=1,
                                 breaks="FD", verbose=1, ...){
  stopifnot(! is.null(X) || ! is.null(maf))

  if(! is.null(X) & is.null(maf)){
    N <- nrow(X)
    P <- ncol(X)
    if(P < N)
      warning("input matrix doesn't seem to have SNPs in columns and individuals in rows")
    if(verbose > 0){
      txt <- paste0(P, " SNPs and ", N, " individuals")
      write(txt, stdout())
    }
    maf <- estimMaf(X)
  }

  tmp <- hist(x=maf, xlab="Minor allele frequency", ylab="Number of SNPs",
              main=main, xlim=xlim, col=col, border=border, las=las,
              breaks=breaks, ...)
}

##' Convert genotype data to the "mean genotype" file format from BimBam
##'
##' The format is specified in BimBam's manual http://www.haplotype.org/download/bimbam-manual.pdf#page=6
##' @param X matrix with individuals in rows and SNPs in columns
##' @param tX matrix with SNPs in rows and individuals in columns
##' @param alleles data.frame with SNPs in rows (names as row names) and
##' alleles in columns (first is "minor", second is "major")
##' @param file prints the genotype data to this file if non NULL (for instance
##' 'genotypes_bimbam.txt' or gzfile('genotypes_bimbam.txt.gz'))
##' @return data.frame
##' @author Timothee Flutre
dose2bimbam <- function(X=NULL, tX=NULL, alleles, file=NULL){
    stopifnot(xor(is.null(X), is.null(tX)),
              ! is.null(row.names(alleles)),
              colnames(alleles) == c("minor","major"))
    if(is.null(tX))
        tX <- t(X)
    tmp <- cbind(alleles, tX)
    if(! is.null(file))
        write.table(x=tmp, file=file, quote=FALSE, sep="\t", row.names=TRUE,
                    col.names=FALSE)
    return(tmp)
}

##' Return the additive relationship matrix, A, from the output of the coancestry() function in the "related" package.
##'
##' The "related" package can be found here: http://frasierlab.wordpress.com/software/. The A matrix is also known as the numerator relationship matrix. It is calculated as explained in chapter 2 from Mrode (2005).
##' @param x list returned by coancestry()
##' @param estim.coancestry name of the coancestry estimator (e.g. "dyadml")
##' @param estim.inbreeding name of the inbreeding estimator (e.g. "LR")
##' @param debug boolean (TRUE to check the output matrix is indeed symmetric)
##' @return matrix
##' @author Timothee Flutre
coancestry2addrel <- function(x, estim.coancestry, estim.inbreeding=NULL,
                              debug=FALSE){
  stopifnot(estim.coancestry %in% colnames(x$relatedness))
  if(! is.null(estim.inbreeding)){
    stopifnot("inbreeding" %in% names(x))
    stopifnot(estim.inbreeding %in% colnames(x$inbreeding))
  }

  ind.ids <- unique(c(x$relatedness$ind1.id,
                      x$relatedness$ind2.id))
  nb.inds <- length(ind.ids)
  A <- matrix(NA, nrow=nb.inds, ncol=nb.inds,
              dimnames=list(ind.ids, ind.ids))

  diag(A) <- 1
  if(! is.null(estim.inbreeding)){
    stopifnot(nrow(x$inbreeding) == nb.inds)
    for(i in 1:nrow(x$inbreeding))
      A[x$inbreeding$ind.id[i], x$inbreeding$ind.id[i]] <-
        1 + x$inbreeding[[estim.inbreeding]][i]
  }

  for(i in 1:nrow(x$relatedness))
    A[x$relatedness$ind1.id[i], x$relatedness$ind2.id[i]] <-
      x$relatedness[[estim.coancestry]][i]
  idx.na.upper <- which(upper.tri(A) & is.na(A))
  idx.equiv.lower <- which(lower.tri(t(A)) & is.na(t(A)))
  A[idx.na.upper] <- A[idx.equiv.lower]
  A[lower.tri(A)] <- t(A)[lower.tri(A)]

  if(debug){ # check that the matrix is symmetric
    for(i in 1:(nrow(A)-1))
      for(j in (i+1):ncol(A))
        if(A[i,j] != A[j,i])
          stop(paste0("matrix not symmetric at (", i, ",", j, ")"))
  }

  return(A)
}

##' Convert the SFS of independent replicates into a matrix of allele dosage.
##'
##' SFS stands for site frequency spectrum
##' @param seg.sites list returned by scrm()
##' @return matrix with diploid individuals in rows and SNPs in columns
##' @author Timothee Flutre
segSites2allDoses <- function(seg.sites){
  stopifnot(is.list(seg.sites),
            length(unique(sapply(seg.sites, nrow))) == 1)

  nb.inds <- nrow(seg.sites[[1]]) / 2 # nb of diploid individuals
  nb.chrs <- sum(sapply(seg.sites, ncol)) # nb of SNPs
  X <- matrix(data=NA, nrow=nb.inds, ncol=nb.chrs)

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

##' Make a data.frame of SNP coordinates from the SFS of independent replicates
##'
##' SFS stands for site frequency spectrum
##' @param seg.sites list returned by scrm()
##' @param snp.ids vector of identifiers (one per SNP)
##' @param prefix character string
##' @return data.frame with SNPs in rows and 2 columns (chr, pos)
##' @author Timothee Flutre
segSites2snpCoords <- function(seg.sites, snp.ids, prefix="chr"){
  stopifnot(is.list(seg.sites),
            length(unique(sapply(seg.sites, nrow))) == 1)

  nb.chrs <- length(seg.sites) # nb of chromosomes
  nb.snps.per.chr <- sapply(seg.sites, ncol)
  nb.snps <- sum(nb.snps.per.chr) # nb of SNPs

  snp.coords <- data.frame(chr=rep(NA, nb.snps),
                           pos=-1,
                           row.names=snp.ids)
  snp.coords$chr <- rep(paste0(prefix, 1:nb.chrs), nb.snps.per.chr)
  snp.coords$pos <- as.numeric(do.call(c, lapply(seg.sites, colnames)))

  ## convert genomic positions from float to integer
  i <- 0
  while(TRUE){
    tmp <- c(by(floor(snp.coords$pos * 10^i), factor(snp.coords$chr),
                function(x){anyDuplicated(x)}))
    if(! any(tmp))
      break
    i <- i + 1
  }
  snp.coords$pos <- floor(snp.coords$pos * 10^i)

  return(snp.coords)
}

##' Simulate according to an approximation to the coalescent with recombination named the Sequential Coalescent with Recombination Model.
##'
##' Requires the scrm package (Staab et al, 2014).
##' @param nb.inds diploids (thus nb of haplotypes is 2 * nb.inds)
##' @param ind.ids vector of identifiers (one per individual)
##' @param nb.reps number of independent loci that will be produced (could be seen as distinct chromosomes)
##' @param pop.mut.rate theta = 4 N0 mu
##' @param pop.recomb.rate rho = 4 N0 r
##' @param chrom.len in bp
##' @param nb.pops number of populations
##' @param mig.rate migration rate = 4 N0 m (symmetric)
##' @param verbose verbosity level (default=0=nothing, 1=few, 2=more)
##' @return list with haplotypes (list), genotypes as allele doses (matrix) and SNP coordinates (data.frame)
##' @author Timothee Flutre
simulCoalescent <- function(nb.inds=100,
                            ind.ids=NULL,
                            nb.reps=20,
                            pop.mut.rate=50,
                            pop.recomb.rate=5,
                            chrom.len=10^3,
                            nb.pops=1,
                            mig.rate=5,
                            verbose=0){
  if(! requireNamespace("scrm", quietly=TRUE))
    stop("Pkg scrm needed for this function to work. Please install it.",
         call.=FALSE)
  stopifnot(nb.inds > nb.pops)

  if(is.null(ind.ids))
    ind.ids <- sprintf(fmt=paste0("ind%0", floor(log10(nb.inds))+1, "i"),
                       1:nb.inds)

  ## simulate according to the SCRM
  nb.samples <- nb.inds * 2 # e.g. 2 chr1 in ind1, 2 chr1 in ind2, etc
  cmd <- paste0(nb.samples, " ", nb.reps)
  cmd <- paste0(cmd, " -t ", pop.mut.rate)
  cmd <- paste0(cmd, " -r ", pop.recomb.rate, " ", chrom.len)
  cmd <- paste0(cmd, " -T") # print genealogies in newick
  cmd <- paste0(cmd, " -L") # print TMRCA and local tree lengths
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
    print(str(sum.stats))

  ## make a data.frame with SNP coordinates
  prefix <- "chr"
  names(sum.stats$seg_sites) <- paste0(prefix, 1:nb.reps)
  nb.snps.per.chr <- sapply(sum.stats$seg_sites, ncol)
  nb.snps <- sum(nb.snps.per.chr)
  snp.ids <- sprintf(fmt=paste0("snp%0", floor(log10(nb.snps))+1, "i"),
                     1:nb.snps)
  snp.coords <- segSites2snpCoords(sum.stats$seg_sites, snp.ids, prefix)
  for(c in 1:nb.reps){
    colnames(sum.stats$seg_sites[[c]]) <-
      snp.ids[(ifelse(c == 1, 1, 1 + cumsum(nb.snps.per.chr)[c-1])):
                (cumsum(nb.snps.per.chr)[c])]
    rownames(sum.stats$seg_sites[[c]]) <-
      paste0(rep(ind.ids, each=2), rep(c("_h1", "_h2"), nb.inds))
  }

  ## make a matrix with genotypes as allele doses
  X <- segSites2allDoses(sum.stats$seg_sites)
  rownames(X) <- ind.ids
  colnames(X) <- snp.ids
  if(verbose > 0){
    txt <- paste0("nb of SNPs: ", nb.snps)
    write(txt, stdout())
    print(sapply(sum.stats$seg_sites, ncol))
  }

  return(list(haplos=sum.stats$seg_sites,
              genos=X,
              snp.coords=snp.coords))
}

##' Individual names
##'
##' Return the identifiers of all individuals
##' @param haplos list of matrices (one per chromosome, with individuals in rows)
##' @return vector
##' @author Timothee Flutre
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
getHaplosInd <- function(haplos, ind.name){
  stopifnot(is.list(haplos),
            is.character(ind.name))

  ## for each chromosome, retrieve both haplotypes of the given individual
  haplos.ind <- lapply(haplos, function(haplos.chr){
    idx <- grep(ind.name, rownames(haplos.chr))
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
##' @param child.name identifier of the child (if NULL, <parent1>"_x_"<parent2>)
##' @return list of 2-row matrices (one per chromosome)
##' @author Timothee Flutre
fecundation <- function(gam1, gam2, child.name=NULL){
  stopifnot(is.list(gam1),
            is.list(gam2),
            length(gam1) == length(gam2))

  nb.chroms <- length(gam1)
  if(is.null(child.name)){
    parent1 <- unique(sapply(gam1, rownames))
    parent2 <- unique(sapply(gam2, rownames))
    child.name <- paste0(parent1, "_x_", parent2)
  }

  child <- lapply(1:nb.chroms, function(c){
    tmp <- rbind(gam1[[c]], gam2[[c]])
    rownames(tmp) <- c(paste0(child.name, "_h1"), paste0(child.name, "_h2"))
    tmp
  })
  names(child) <- names(gam1)

  return(child)
}

##' Doubled haploids
##'
##' Make a doubled haploid from all chromosomes of the individual under consideration.
##' @param haplos.ind list of matrices (one per chromosome)
##' @param loc.crossovers list of vectors (one per chromosome)
##' @param child.name identifier of the child (if NULL, <parent1>"_x_"<parent2>)
##' @return list of matrices
##' @author Timothee Flutre
makeDoubledHaploidsSingleInd <- function(haplos.ind, loc.crossovers,
                                         child.name=NULL){
  stopifnot(is.list(haplos.ind),
            is.list(loc.crossovers),
            all(names(haplos.ind) == names(loc.crossovers)))

  gam <- makeGameteSingleInd(haplos.ind, loc.crossovers)
  fecundation(gam, gam, child.name)
}

##' Doubled haploids
##'
##' Make a doubled haploid for each individual under consideration.
##' @param haplos list of matrices (one per chromosome)
##' @param ind.names identifiers of the individuals to handle (if NULL, all individuals)
##' @param loc.crossovers list of lists (one per chromosome, then one per individual); if NULL, draw many crossing-overs localizations at once (as Poisson with parameter 2, assuming all chromosomes roughly have the same length)
##' @param child.names identifiers of the children (if NULL, <ind>"_x_"<ind>)
##' @param nb.cores the number of cores to use, i.e. at most how many child processes will be run simultaneously (not on Windows)
##' @return list of matrices
##' @author Timothee Flutre
makeDoubledHaploids <- function(haplos, ind.names=NULL, loc.crossovers=NULL,
                                child.names=NULL, nb.cores=1){
  stopifnot(is.list(haplos),
            length(unique(sapply(haplos, class))) == 1,
            unique(sapply(haplos, class)) == "matrix",
            length(unique(sapply(haplos, nrow))) == 1) # same nb of individuals
  haplos.ind.names <- getIndNamesFromHaplos(haplos)
  if(! is.null(ind.names))
    stopifnot(is.vector(ind.names),
              is.character(ind.names),
              all(ind.names %in% haplos.ind.names))
  if(! is.null(loc.crossovers))
    stopifnot(is.list(loc.crossovers),
              all(names(loc.crossovers) == ind.names),
              length(loc.crossovers[[1]]) == length(haplos)) # same nb of chromosomes
  if(Sys.info()["sysname"] == "Windows")
    nb.cores <- 1

  if(is.null(ind.names)){
    ind.names <- haplos.ind.names
  } else if(length(ind.names) < nrow(haplos[[1]])){
    haplos <- getHaplosInds(haplos, ind.names)
  }

  if(is.null(child.names))
    child.names <- paste0(ind.names, "_x_", ind.names)

  nb.chroms <- length(haplos)
  nb.inds <- length(ind.names)

  ## if required, draw locations of crossing-overs
  if(is.null(loc.crossovers)){
    nb.crossovers <- rpois(n=nb.inds * nb.chroms, lambda=2) # 2 per chromosome
    nb.crossovers[nb.crossovers == 0] <- 1 # at least 1 per chromosome
    loc.crossovers <- lapply(1:nb.inds, function(i){
      tmp <- lapply(1:nb.chroms, function(c){
        sort(sample.int(n=ncol(haplos[[c]]) - 1,
                        size=nb.crossovers[(i-1)*nb.chroms + c]))
      })
      names(tmp) <- nrow(haplos)
      tmp
    })
    names(loc.crossovers) <- ind.names
  }

  ## make haplodiploidization for each individual
  tmp <- parallel::mclapply(1:nb.inds, function(i){
    ind.name <- ind.names[i]
    haplos.ind <- getHaplosInd(haplos, ind.name)
    child.name <- child.names[i]
    makeDoubledHaploidsSingleInd(haplos.ind, loc.crossovers[[i]], child.name)
  }, mc.cores=nb.cores)

  ## gather all individuals per chromosome
  dbl.hpd <- lapply(1:nb.chroms, function(c){
    do.call(rbind, parallel::mclapply(tmp, `[[`, c, mc.cores=nb.cores))
  })
  names(dbl.hpd) <- names(haplos)

  return(dbl.hpd)
}

##' Calculate the distances between SNPs, assuming they are sorted.
##'
##' Useful before estimating pairwise linkage disequilibrium.
##' @param snp.coords data.frame with SNP identifiers as row names, and with two columns "chr" and "pos"
##' @param nb.cores the number of cores to use (default=1)
##' @return list with one component per chromosome
##' @author Timothee Flutre
snpDistances <- function(snp.coords, nb.cores=1){
  if(! requireNamespace("parallel", quietly=TRUE))
    stop("Pkg parallel needed for this function to work. Please install it.",
         call.=FALSE)
  stopifnot(is.data.frame(snp.coords),
            colnames(snp.coords) == c("chr", "pos"),
            ! is.null(rownames(snp.coords)))

  chr.names <- unique(snp.coords$chr)
  snp.dists <- parallel::mclapply(chr.names, function(chr.name){
    pos <- snp.coords$pos[snp.coords$chr == chr.name]
    names(pos) <- rownames(snp.coords)[snp.coords$chr == chr.name]
    dis <- pos[2:length(pos)] - pos[1:(length(pos)-1)]
    names(dis) <- paste(names(dis), names(pos)[-length(pos)], sep="-")
    dis
  }, mc.cores=nb.cores)
  names(snp.dists) <- chr.names

  return(snp.dists)
}

##' Estimate genetic relationships between individuals from their SNP genotypes.
##'
##' SNPs with missing data are ignored.
##' @param X matrix of SNP genotypes encoded as allele doses ({0,1,2}), with SNPs in
##' columns and individuals in rows
##' @param mafs vector with minor allele frequencies (calculated with `estimMaf` if NULL)
##' @param thresh threshold on allele frequencies below which SNPs are ignored (default=0.01, NULL to skip this step)
##' @param relationships genetic relationship to estimate (default=additive/dominance)
##' @param method if additive relationships, can be "astle-balding" (default; see equation 2.2 in Astle & Balding, 2009), "vanraden1" (first method in VanRaden, 2008), "habier" (similar to 'vanraden1' without giving more importance to rare alleles; from Habier et al, 2007), "zhou" (centering the genotypes and not assuming that rare variants have larger effects; from Zhou et al, 2013) or "center-std"
##' @param verbose verbosity level (default=1)
##' @return matrix
##' @author Timothee Flutre
estimGenRel <- function(X, mafs=NULL, thresh=0.01, relationships="additive",
                        method="astle-balding", verbose=1){
  stopifnot(is.matrix(X),
            relationships %in% c("additive", "dominance"))
  if(relationships != "additive")
    stop(paste0(relationships, " relationships is not (yet) implemented"))
  stopifnot(method %in% c("astle-balding", "vanraden1", "habier", "zhou", "center-std"))
  if(! is.null(thresh))
    stopifnot(thresh >= 0, thresh <= 0.5)

  gen.rel <- NULL # to be filled and returned

  N <- nrow(X) # nb of individuals
  P <- ncol(X) # nb of SNPs
  if(P < N)
    warning("input matrix doesn't seem to have SNPs in columns and individuals in rows")

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
    X <- X[, -idx.rm]
    P <- ncol(X)
  }

  ## estimate MAFs
  if(is.null(mafs)){
    mafs <- estimMaf(X)
    if(verbose > 1){
      txt <- paste0("allele freqs: ",
                    "min=", format(min(mafs), digits=2),
                    " Q1=", format(quantile(mafs, 0.25), digits=2),
                    " med=", format(median(mafs), digits=2),
                    " mean=", format(mean(mafs), digits=2),
                    " Q3=", format(quantile(mafs, 0.75), digits=2),
                    " max=", format(max(mafs), digits=2))
      write(txt, stdout())
    }
  }

  ## discard SNPs with low MAFs
  if(! is.null(thresh)){
    snps.low <- mafs < thresh
    if(any(snps.low)){
      if(verbose > 0){
        txt <- paste0("skip ", sum(snps.low), " SNPs with freq below ", thresh)
        write(txt, stdout())
      }
      idx.rm <- which(snps.low)
      X <- X[, -idx.rm]
      P <- ncol(X)
      mafs <- mafs[-idx.rm]
    }
  }

  ## estimate genetic relationships
  if(relationships == "additive"){
    if(method == "astle-balding"){
      tmp <- sweep(x=X, MARGIN=2, STATS=2 * mafs, FUN="-")
      tmp <- sweep(x=tmp, MARGIN=2, STATS=sqrt(4 * mafs * (1 - mafs)), FUN="/")
      gen.rel <- tcrossprod(tmp, tmp) / P
    } else if(method == "vanraden1"){
      M <- X - 1
      Pmat <- matrix(rep(1,N)) %*% (2 * (mafs - 0.5))
      Z <- M - Pmat
      gen.rel <- tcrossprod(Z, Z) / (2  * sum(mafs * (1 - mafs)))
    } else if(method == "habier"){
      gen.rel <- tcrossprod(X, X) / (2  * sum(mafs * (1 - mafs)))
    } else if(method == "zhou"){
      ## tmp <- sweep(x=X, MARGIN=2, STATS=colMeans(X), FUN="-")
      tmp <- scale(x=X, center=TRUE, scale=FALSE)
      gen.rel <- tcrossprod(tmp, tmp) / P
    } else if(method == "center-std"){
      ## X.cs <- sweep(x=X, MARGIN=2, STATS=colMeans(X), FUN="-")
      ## tmp <- sweep(x=X.cs, MARGIN=2, STATS=apply(X=X.cs, MARGIN=2, sd), FUN="/")
      tmp <- scale(x=X, center=TRUE, scale=TRUE)
      gen.rel <- tcrossprod(tmp, tmp) / P
    }
  } else if(relationships == "dominance"){
    H <- matrix(data=NA, nrow=N, ncol=P) # <- TODO
    gen.rel <- tcrossprod(H, H) / (2 * sum(mafs * (1 - mafs) * (1 - 2 * mafs * (1 - mafs))))
  }

  return(gen.rel)
}

##' Return estimates of linkage disequilibrium between pairs of SNPs.
##'
##' Requires package LDcorSV
##' @param X matrix of SNP genotypes encoded as allele doses, with SNPs in columns and individuals in rows
##' @param K matrix of kinship
##' @param pops vector of characters indicating the population of each individual
##' @param snp.coords data.frame with SNP identifiers as row names, and two columns, "chr" and "pos"
##' @param only.chr identifier of a given chromosome
##' @param only.pop identifier of a given population
##' @param verbose verbosity level (default=0/1)
##' @return data frame
##' @author Timothee Flutre
estimLd <- function(X, K=NULL, pops=NULL, snp.coords,
                    only.chr=NULL, only.pop=NULL, verbose=0){
  if(! requireNamespace("LDcorSV", quietly=TRUE))
    stop("Pkg LDcorSV needed for this function to work. Please install it.",
         call.=FALSE)
  stopifnot(is.matrix(X),
            ! is.null(dimnames(X)),
            sum(is.na(X)) == 0,
            is.data.frame(snp.coords),
            colnames(snp.coords) == c("chr", "pos"))
  if(! is.null(K))
    stopifnot(is.matrix(K),
              nrow(K) == ncol(K),
              nrow(K) == nrow(X),
              ! is.null(dimnames(K)),
              all(rownames(K) == colnames(K)),
              all(rownames(K) == rownames(X)))
  W.s <- NA
  if(! is.null(pops)){
    stopifnot(length(pops) == nrow(X),
              ! is.null(names(pops)),
              names(pops) == rownames(X))
    W.s <- model.matrix(~ as.factor(pops))[, -1]
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
  if(! is.null(only.chr))
    subset.snps <- which(snp.coords$chr == only.chr)
  subset.inds <- 1:nrow(X)
  if(! is.null(only.pop))
    subset.inds <- which(pops == only.pop)

  if(verbose > 0)
    write("estimate pairwise LD ...", stdout())
  if(is.null(K)){
    if(is.null(only.pop)){
      ld <- LDcorSV::LD.Measures(donnees=X[subset.inds, subset.snps],
                                 V=NA,
                                 S=W.s,
                                 data="G", supinfo=FALSE, na.presence=FALSE)
    } else
      ld <- LDcorSV::LD.Measures(donnees=X[subset.inds, subset.snps],
                                 V=NA,
                                 S=NA,
                                 data="G", supinfo=FALSE, na.presence=FALSE)
  } else{
    if(is.null(only.pop)){
      ld <- LDcorSV::LD.Measures(donnees=X[subset.inds, subset.snps],
                                 V=K[subset.inds, subset.inds],
                                 S=W.s,
                                 data="G", supinfo=FALSE, na.presence=FALSE)
    } else
      ld <- LDcorSV::LD.Measures(donnees=X[subset.inds, subset.snps],
                                 V=K[subset.inds, subset.inds],
                                 S=NA,
                                 data="G", supinfo=FALSE, na.presence=FALSE)
  }

  return(ld)
}

##' Simulate phenotypes from a basic "animal model" (LMM).
##'
##' y = W alpha + Z u + epsilon
##' where y is N x 1; W is N x Q; Z is N x I
##' u ~ Norm_I(0, G=sigma_u^2 A); epsilon_n ~ Student(nu, 0); Cov(u,e)=0
##' @param Q number of years
##' @param mu overall mean
##' @param mean.a mean of the prior on alpha[2:Q]
##' @param sd.a std dev of the prior on alpha[2:Q]
##' @param lambda ratio of variances (sigma_u^2 /sigma^2)
##' @param sigma.u2 genetic variance (e.g. 15)
##' @param scale.halfCauchy scale of the half-Cauchy prior for sigma_u (e.g. 5)
##' @param A matrix of additive genetic relationships
##' @param perc.NA percentage of missing phenotypes, at random
##' @param err.df degrees of freedom of errors' Student's t-distribution
##' @param seed seed for the pseudo-random number generator
##' @return list
##' @author Timothee Flutre
simulAnimalModel <- function(Q=3, mu=50, mean.a=5, sd.a=2,
                             A, lambda=3, sigma.u2=NULL, scale.halfCauchy=NULL,
                             perc.NA=0, err.df=Inf,
                             seed=NULL){
  if(! requireNamespace("MASS", quietly=TRUE))
    stop("Pkg MASS needed for this function to work. Please install it.",
         call.=FALSE)
  if(! requireNamespace("Matrix", quietly=TRUE))
    stop("Pkg Matrix needed for this function to work. Please install it.",
         call.=FALSE)
  stopifnot(xor(is.null(sigma.u2), is.null(scale.halfCauchy)),
            is.matrix(A),
            nrow(A) == ncol(A),
            ! is.null(rownames(A)),
            ! is.null(colnames(A)),
            rownames(A) == colnames(A))
  if(! is.null(seed))
    set.seed(seed)

  I <- nrow(A)
  N <- Q * I

  levels.years <- as.character(seq(from=2010, to=2010+Q-1))
  if(N %% Q == 0){
    years <- rep(levels.years, each=N / Q)
  } else
    years <- sort(sample(x=levels.years, size=N, replace=TRUE))
  years <- as.factor(years)
  W <- model.matrix(~ years)
  dat <- data.frame(year=years)

  ## "fixed" effects
  alpha <- matrix(data=c(mu, rnorm(n=Q-1, mean=mean.a, sd=sd.a)),
                  nrow=Q, ncol=1)

  levels.inds <- rownames(A)
  inds <- rep(NA, N)
  for(year in levels.years)
    inds[years == year] <- levels.inds[1:sum(years == year)]
  inds <- as.factor(inds)
  ## Z <- as.matrix(Matrix::t(as(inds, Class="sparseMatrix")))
  Z <- model.matrix(~ inds - 1)
  dat$ind <- inds

  ## additive genetic component
  if(is.null(sigma.u2)){
    sigma.u <- abs(rcauchy(n=1, location=0, scale=scale.halfCauchy))
    sigma.u2 <- sigma.u^2
  }
  G <- as.matrix(Matrix::nearPD(sigma.u2 * A)$mat)
  u <- setNames(object=MASS::mvrnorm(n=1, mu=rep(0, I), Sigma=G),
                nm=rownames(A))

  ## errors
  sigma2 <- sigma.u2 / lambda
  h2 <- sigma.u2 / (sigma.u2 + sigma2)
  if(is.infinite(err.df)){
    epsilon <- rnorm(n=N, mean=0, sd=sqrt(sigma2))
  } else
    epsilon <- rt(n=N, df=err.df, ncp=0)

  ## phenotypes
  y <- W %*% alpha + Z %*% u + epsilon
  if(perc.NA > 0){
    idx <- sample.int(n=N, size=floor(perc.NA/100 * N))
    y[idx] <- NA
  }
  dat$response <- y[,1]

  return(list(y=y,
              W=W, alpha=alpha,
              Z=Z, G=G, u=u, sigma.u2=sigma.u2,
              sigma2=sigma2,
              h2=h2,
              dat=dat))
}

##' Fit a basic "animal model" with lme4
##'
##' y = W alpha + Z u + epsilon
##' where y is N x 1; W is N x Q; Z is N x I
##' u ~ Norm_I(0, G=sigma_u^2 A); epsilon_n ~ N(0, sigma^2 I); Cov(u,e)=0
##' See http://stackoverflow.com/q/19327088/597069.
##' Still problems: see https://github.com/lme4/lme4/issues/340.
##' @param formula formula
##' @param data data.frame containing the data corresponding to formula and relmat
##' @param relmat list containing the matrix of additive genetic relationships (should use the same name as the colname in data)
##' @param REML default is TRUE, but use FALSE to compare models with different fixed effects
##' @param verbose verbosity level (0/default=1)
##' @return merMod object
##' @author Timothee Flutre
lmerAM <- function(formula, data, relmat, REML=TRUE, verbose=1){
  if(! requireNamespace("lme4", quietly=TRUE))
    stop("Pkg lme4 needed for this function to work. Please install it.",
         call.=FALSE)
  if(! requireNamespace("Matrix", quietly=TRUE))
    stop("Pkg Matrix needed for this function to work. Please install it.",
         call.=FALSE)
  for(i in seq_along(relmat))
    stopifnot(is.matrix(relmat[[i]]),
              names(relmat)[i] %in% colnames(data))

  if(verbose > 0)
    write("parse the formula", stdout())
  parsedFormula <- lme4::lFormula(formula=formula,
                                  data=data,
                                  control=lme4::lmerControl(
                                      check.nobs.vs.nlev="ignore",
                                      check.nobs.vs.nRE="ignore"),
                                  REML=REML)
  if(verbose > 0)
    write("structure the design and covariance matrices of the random effects",
          stdout())
  relfac <- relmat
  flist <- parsedFormula$reTrms[["flist"]] # list of grouping factors
  Ztlist <- parsedFormula$reTrms[["Ztlist"]] # list of transpose of the sparse model matrices
  stopifnot(all(names(relmat) %in% names(flist)))
  asgn <- attr(flist, "assign")
  for(i in seq_along(relmat)) {
    tn <- which(match(names(relmat)[i], names(flist)) == asgn)
    if(length(tn) > 1)
      stop("a relationship matrix must be associated",
           " with only one random effects term", call.=FALSE)
    relmat[[i]] <- Matrix::Matrix(relmat[[i]], sparse=TRUE)
    relfac[[i]] <- chol(relmat[[i]])
    Ztlist[[i]] <- relfac[[i]] %*% Ztlist[[i]]
  }
  parsedFormula$reTrms[["Ztlist"]] <- Ztlist
  parsedFormula$reTrms[["Zt"]] <- do.call(Matrix::rBind, Ztlist)

  if(verbose > 0)
    write("make the deviance function", stdout())
  devianceFunction <- do.call(lme4::mkLmerDevfun, parsedFormula)

  if(verbose > 0)
    write("optimize the deviance function", stdout())
  optimizerOutput <- lme4::optimizeLmer(devianceFunction)

  if(verbose > 0)
    write("make the output", stdout())
  fit <- lme4::mkMerMod(rho=environment(devianceFunction),
                        opt=optimizerOutput,
                        reTrms=parsedFormula$reTrms,
                        fr=parsedFormula$fr)
  return(fit)
}

##' Simulate phenotypes according to the BSLMM model.
##'
##' y = W alpha + Z X tilde{beta} + Z u + epsilon
##' where y is N x 1; W is N x Q; Z is N x I; X is I x P
##' tilde{beta}_p ~ pi Norm_1(0, sigma_beta^2/sigma^2) + (1 - pi) delta_0
##' u ~ Norm_I(0, (sigma_u^2/sigma^2) A)
##' epsilon ~ Norm_N(0, sigma^2 I)
##' See Zhou, Carbonetto & Stephens (2013).
##' @param Q number of years
##' @param mu overall mean
##' @param mean.a mean of the prior on alpha[2:Q]
##' @param sd.a std dev of the prior on alpha[2:Q]
##' @param X matrix of SNP genotypes encoded as allele doses, with SNPs in
##' columns and individuals in rows (SNPs with missing values or low MAF
##' should be discarded beforehand)
##' @param pi proportion of beta-tilde values that are non-zero
##' @param h approximation to E[PVE] (h and rho should be NULL or not together)
##' @param rho approximation to E[PGE]
##' @param perc.NA percentage of missing phenotypes, at random
##' @param err.df degrees of freedom of errors' Student's t-distribution
##' @param seed seed for the pseudo-random number generator
##' @return list
##' @author Timothee Flutre
simulBslmm <- function(Q=3, mu=50, mean.a=5, sd.a=2,
                       X, pi=NULL, h=NULL, rho=NULL,
                       perc.NA=0, err.df=Inf,
                       seed=NULL){
  if(! requireNamespace("MASS", quietly=TRUE))
    stop("Pkg MASS needed for this function to work. Please install it.",
         call.=FALSE)
  ## if(! requireNamespace("Matrix", quietly=TRUE))
  ##   stop("Pkg Matrix needed for this function to work. Please install it.",
  ##        call.=FALSE)
  stopifnot(xor(is.null(h) & is.null(rho), ! (is.null(h) & is.null(rho))),
            sum(is.na(X)) == 0,
            ! is.null(rownames(X)),
            ! is.null(colnames(X)))
  if(! is.null(seed))
    set.seed(seed)

  I <- nrow(X)
  P <- ncol(X)
  if(P < I)
    warning("input matrix doesn't seem to have SNPs in columns and individuals in rows")
  if(Q > 0){
    N <- Q * I
  } else
    N <- I

  ## incidence matrix of the non-genetic predictors having "fixed effects"
  if(Q > 0){
    levels.years <- as.character(seq(from=2010, to=2010+Q-1))
    if(N %% Q == 0){
      years <- rep(levels.years, each=N / Q)
    } else
      years <- sort(sample(x=levels.years, size=N, replace=TRUE))
    years <- as.factor(years)
    W <- model.matrix(~ years)
  } else
    W <- matrix(data=0, nrow=N, ncol=1)

  ## "fixed" effects
  if(Q > 0){
    alpha <- matrix(data=c(mu, rnorm(n=Q-1, mean=mean.a, sd=sd.a)),
                    nrow=Q, ncol=1)
  } else
    alpha <- matrix(data=0, nrow=1, ncol=1)

  ## incidence matrices of the genetic predictors
  levels.inds <- rownames(X)
  inds <- rep(NA, N)
  if(Q > 0){
    for(year in levels.years)
      inds[years == year] <- levels.inds[1:sum(years == year)]
  } else
    inds <- levels.inds
  inds <- as.factor(inds)
  ## Z <- as.matrix(Matrix::t(as(inds, Class="sparseMatrix")))
  Z <- model.matrix(~ inds - 1)
  X.c <- scale(x=X, center=TRUE, scale=FALSE)
  A <- tcrossprod(X.c, X.c) / P

  ## hyper-parameters
  s.a <- (1 / (N*P)) * sum(colSums(X.c^2))
  s.b <- (1/N)  * sum(diag(A))
  if(is.null(pi))
    pi <- exp(runif(n=1, min=log(1/P), max=log(1)))
  if(is.null(h))
    h <- runif(n=1, min=0, max=1)
  if(is.null(rho))
    rho <- runif(n=1, min=0, max=1)
  sigma.betat2 <- (h * rho) / ((1 - h) * P * pi * s.a)
  if(is.nan(sigma.betat2))
    sigma.betat2 <- 0
  sigma.u2 <- (h * (1 - rho)) / ((1 - h) * s.b)
  tau <- 1

  ## sparse genetic effects
  betat <- setNames(object=rep(0, P), nm=colnames(X))
  gamma <- setNames(object=rbinom(n=P, size=1, prob=pi), nm=colnames(X))
  betat[gamma == 1] <- rnorm(n=sum(gamma == 1), mean=0,
            sd=sqrt(sigma.betat2 * tau^(-1)))

  ## polygenic effects
  u <- setNames(object=MASS::mvrnorm(n=1, mu=rep(0, I),
                                     Sigma=sigma.u2 * tau^(-1) * A),
                nm=rownames(X))

  ## errors
  if(is.infinite(err.df)){
    epsilon <- matrix(rnorm(n=N, mean=0, sd=sqrt(tau^(-1))))
  } else
    epsilon <- matrix(rt(n=N, df=err.df, ncp=0))

  ## phenotypes
  y <- W %*% alpha + Z %*% X.c %*% betat + Z %*% u + epsilon
  if(perc.NA > 0){
    idx <- sample.int(n=N, size=floor(perc.NA/100 * N))
    y[idx] <- NA
  }

  return(list(y=y,
              W=W, alpha=alpha,
              Z=Z, X.c=X.c, s.a=s.a, s.b=s.b, A=A,
              pi=pi, h=h, rho=rho,
              sigma.betat2=sigma.betat2, sigma.u2=sigma.u2, sigma2=1/tau,
              betat=betat, u=u))
}

##' Launch GEMMA
##'
##' @param model name of the model to fit (default=ulmm/bslmm)
##' @param y vector of phenotypes
##' @param X matrix of SNP genotypes encoded as allele doses, with SNPs in
##' columns and individuals in rows
##' @param snp.coords data.frame with 3 columns: chromosomes' identifiers,
##' coordinates of SNPs on the chromosome and SNPs names
##' @param alleles data.frame with SNPs in rows (names as row names) and
##' alleles in columns (first is "minor", second is "major")
##' @param K.c kinship centered matrix of SNPs
##' @param W matrix of covariates with individuals in rows (names as row names), a first column of 1 and a second column of covariates values
##' @param out.dir directory in which the output files will be saved
##' @param task.id identifier of the task (used in output file names)
##' @param verbose verbosity level (default=1)
##' @param clean remove temporary files
##' @param burnin number of iterations to discard as burn-in
##' @param nb.iters number of iterations
##' @param thin thining
##' @return invisible data.frame
##' @author Timothee Flutre [cre,aut], Dalel Ahmed [ctb]
gemma <- function(model="ulmm", y, X, snp.coords, alleles, K.c=NULL, W,
                  out.dir, task.id="", verbose=1, clean=FALSE,
                  burnin=1000, nb.iters=7000, thin=10){
  stopifnot(model %in% c("ulmm", "bslmm"),
            is.vector(y),
            is.matrix(X),
            ! is.null(colnames(X)),
            anyDuplicated(colnames(X)) == 0,
            ! is.null(row.names(alleles)),
            colnames(alleles) == c("minor","major"),
            "chr" %in% colnames(snp.coords),
            "snp" %in% colnames(snp.coords),
            all(snp.coords$snp == colnames(X)),
            all(snp.coords$snp == rownames(alleles)),
            file.exists(out.dir))

  ## prepare input files
  X.bimbam <- dose2bimbam(X=X, alleles=alleles,
                          file=gzfile(paste0(out.dir,
                              "/genos_bimbam_", task.id, ".txt.gz")))
  write.table(x=snp.coords,
              file=paste0(out.dir, "/snp_coordinates_", task.id, ".txt"),
              quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  if(is.null(K.c))
    K.c <- estimGenRel(X=X, method="center")
  write.table(x=K.c,
              file=paste0(out.dir, "/kinship-center_", task.id, ".txt"),
              quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  write.table(x=W,
              file=paste0(out.dir, "/covar_gemma_", task.id, ".txt"),
              quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  write.table(x=y,
              file=gzfile(paste0(out.dir, "/phenotypes_", task.id, ".txt.gz")),
              quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

  ## prepare cmd-line and execute it
  cmd <- paste0("cd ", out.dir,
                "; gemma",
                " -g genos_bimbam_", task.id, ".txt.gz",
                " -p phenotypes_", task.id, ".txt.gz",
                " -a snp_coordinates_", task.id, ".txt",
                " -o results_simul_", task.id)
  if(model == "ulmm"){
    cmd <- paste0(cmd, " -k kinship-center_", task.id, ".txt",
                  " -c covar_gemma_", task.id, ".txt",
                  " -lmm 4")
  } else if(model == "bslmm"){
    cmd <- paste0(cmd, " -bslmm 1",
                  " -w ", format(x=burnin, scientific=FALSE),
                  " -s ", format(x=nb.iters, scientific=FALSE),
                  " -rpace ", thin)
  }
  if(verbose == 0)
    cmd <- paste0(cmd, " >& stdouterr_gemma_", task.id, ".txt")
  system(cmd)

  ## load output files
  if(model == "ulmm"){
    output <- read.table(file=paste0(out.dir, "/output/results_simul_",
                             task.id, ".assoc.txt"), sep="\t",
                         header=TRUE, stringsAsFactors=FALSE)
    rownames(output) <- output$rs
  } else if(model == "bslmm"){
    output <- list()
    output[["hyperparams"]] <- read.table(file=paste0(out.dir,
                                              "/output/results_simul_",
                                              task.id, ".hyp.txt"), sep="\t",
                                          skip=1, stringsAsFactors=FALSE,
                                          header=FALSE,
                                          col.names=c("h", "pve", "rho", "pge",
                                              "pi", "n_gamma", ""))
    output[["hyperparams"]][7] <- NULL
    output[["params"]] <- read.table(file=paste0(out.dir,
                                         "/output/results_simul_",
                                         task.id, ".param.txt"), sep="\t", skip=1,
                                     stringsAsFactors=FALSE, header=FALSE,
                                     col.names=c("chr", "rs", "ps", "n_miss",
                                         "alpha", "beta", "gamma"))
  }

  invisible(output)
}

##' Launch GEMMA uLMM chromosome per chromosome
##'
##' @param y vector of phenotypes
##' @param X matrix of SNP genotypes encoded as allele doses, with SNPs in
##' columns and individuals in rows
##' @param snp.coords data.frame with 3 columns (chr, snp, coord)
##' @param alleles data.frame with SNPs in rows (names as row names) and
##' @param chr.ids set of chromosome identifiers to analyze (optional, all by default)
##' @param W matrix of covariates with individuals in rows (names as row names), a first column of 1 and a second column of covariates values
##' @param out.dir directory in which the data files will be found
##' @param task.id identifier of the task (used in output file names)
##' @param clean remove temporary files
##' @param verbose verbosity level (default=1)
##' @return a list of the GEMMA outputs for all chromosomes
##' @author Timothee Flutre [cre,aut], Dalel Ahmed [ctb]
gemmaUlmmPerChr <- function(y, X, snp.coords, alleles, chr.ids=NULL, W, out.dir,
                            task.id="", clean=FALSE, verbose=1){
  stopifnot(is.vector(y),
            is.matrix(X),
            ! is.null(colnames(X)),
            anyDuplicated(colnames(X)) == 0,
            ! is.null(row.names(alleles)),
            colnames(alleles) == c("minor","major"),
            "chr" %in% colnames(snp.coords),
            "snp" %in% colnames(snp.coords),
            all(snp.coords$snp == colnames(X)),
            all(snp.coords$snp == rownames(alleles)),
            file.exists(out.dir))

  out <- list()
  if(is.null(chr.ids))
    chr.ids <- unique(snp.coords$chr)

  for(chr.id in chr.ids){
    subset.snp.ids <- (snp.coords$chr == chr.id)
    if(verbose > 0)
      message(paste0(chr.id, ": ", sum(subset.snp.ids), " SNPs"))
      K.c <- estimGenRel(X=X[, ! subset.snp.ids], method="center")
    out[[chr.id]] <- gemma(model="ulmm",
                           y=y,
                           X=X[, subset.snp.ids],
                           snp.coords=snp.coords[subset.snp.ids,],
                           alleles=alleles[subset.snp.ids,],
                           K.c=K.c, W=W, out.dir=out.dir,
                           task.id=paste0(task.id, "-", chr.id),
                           verbose=verbose-1, clean=clean)
  }
  out <- do.call(rbind, out)
  rownames(out) <- out$rs

  return(out)
}

##' Produce a quantile-quantile plot for p values and display its confidence
##' interval
##'
##' A quantile is an order statistic, and the j-th order statistic from a
##' Uniform(0,1) sample has a Beta(j,N-j+1) distribution (Casella & Berger,
##' 2001, 2nd edition, p230).
##' Let us assume we have N independent p values, \eqn{\{p_1,\ldots,p_N\}}, for
##' instance: pvalues <- c(runif(99000,0,1), rbeta(1000,0.5,1)). Under the
##' null, they are independent and identically uniformly distributed:
##' \eqn{\forall i \; p_i \sim \mathcal{U}_{[0,1]}}.
##' Therefore, the 95% confidence interval for the j-th quantile of the set
##' of p values can be calculated with: qbeta(0.95, j, N-j+1).
##' TODO: look at this https://github.com/stephenturner/qqman/blob/v0.0.0/qqman.r
##' @param pvalues vector of raw p values
##' @param plot.conf.int show the confidence interval (default=TRUE)
##' @param xlab a title for the x axis (see default)
##' @param ylab a title for the x axis (see default)
##' @param main an overall title for the plot (default: "Q-Q plot (<length(pvalues)> p-values)")
##' @param col plotting color for the points (default is all points in black)
##' @param ... graphical parameters other than xlim, ylim, xlab, ylab, las and col
##' @author Timothee Flutre (inspired from an anonymous comment to http://gettinggeneticsdone.blogspot.fr/2009/11/qq-plots-of-p-values-in-r-using-ggplot2.html)
qqplotPval <- function(pvalues, plot.conf.int=TRUE,
                       xlab=expression(Expected~~-log[10](italic(p)~values)),
                       ylab=expression(Observed~~-log[10](italic(p)~values)),
                       main=NULL, col=NULL){
  N <- length(pvalues)
  expected <- - log10(1:N / N)
  observed <- - log10(pvalues)
  MAX <- max(c(expected, observed))

  if(plot.conf.int){
    c95 <- rep(0, N)
    c05 <- rep(0, N)
    for(j in 1:N){
      c95[j] <- qbeta(0.95, j, N-j+1)
      c05[j] <- qbeta(0.05, j, N-j+1)
    }
    c95 <- - log10(c95)
    c05 <- - log10(c05)
    plot(expected, c95, ylim=c(0,MAX), xlim=c(0,MAX), type="l",
         axes=FALSE, xlab="", ylab="")
    par(new=T)
    plot(expected, c05, ylim=c(0,MAX), xlim=c(0,MAX), type="l",
         axes=FALSE, xlab="", ylab="")
    par(new=T)
  }

  if(is.null(main))
    main <- paste0("Q-Q plot (", N, " p-values)")

  if(is.null(col)){
    col <- rep(1, N)
  } else if(length(col) != N)
    stop("param 'col' should have the same length as 'pvalues'")

  plot(x=sort(expected), y=sort(observed),
       xlim=c(0,MAX), ylim=c(0,MAX),
       las=1, col=col[order(observed)],
       xlab=xlab, ylab=ylab, main=main)
  abline(0, 1, col="red")
}

##' Asymptotic Bayes factor
##'
##' Calculate the asymptotic Bayes factor proposed by Wakefield in Genetic Epidemiology 33:79-86 (2009, \url{http://dx.doi.org/10.1002/gepi.20359}).
##' @param theta.hat MLE of the additive genetic effect
##' @param V variance of theta.hat
##' @param W variance of the prior on theta
##' @param log10 to return the log10 of the ABF (default=TRUE)
##' @return numeric
##' @author Timothee Flutre
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
##' @param log10 to return the log10 of the ABF (default=TRUE)
##' @return numeric
##' @author Bertrand Servin [cre,aut], Timothee Flutre [ctb]
calcExactBayesFactorServinStephens <- function(G, Y, sigma.a, sigma.d,
                                               log10=TRUE){
  stopifnot(is.vector(G), is.vector(Y))

  subset <- complete.cases(Y) & complete.cases(G)
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
##' @author Xiaoquan Wen [cre,aut], Timothee Flutre [ctb]
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
##' Plot a pedigree using the "igraph" package. Inspired by plot.pedigree() from the "synbreed" package (under GPL-3). Add options for monoecious species and auto-fecondation.
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
    symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
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
