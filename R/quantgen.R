## Contains functions useful for quantitative genetics.

##' Convert genotypes
##'
##' Reformat bi-allelic SNP genotypes encoded in genotypic classes to ease subsequent manipulations.
##' @param file the name of the file which the data are to be read from, with a header line, columns separated by a tabulation, and row names as the first column
##' @param x if \code{file=NULL}, data.frame of bi-allelic SNP genotypes encoded in genotypic classes, i.e. in {AA,AB or BA,BB}, with SNPs in rows and genotypes in columns; if it is a matrix, it will be silently transformed into a data frame
##' @param na.string a character to be interpreted as NA values
##' @param sep separator of alleles within each SNP genotype to be removed
##' @param verbose verbosity level (0/1)
##' @return matrix with SNPs in rows and genotypes in columns
##' @author Eric Duchene [aut], Timothee Flutre [ctb]
##' @seealso \code{\link{genoClasses2genoDoses}}, \code{\link{genoClasses2JoinMap}}
##' @export
reformatGenoClasses <- function(file=NULL, x=NULL, na.string="--", sep="",
                                verbose=1){
  stopifnot(xor(! is.null(file), ! is.null(x)))
  if(! is.null(file)){
    stopifnot(is.character(file),
              file.exists(file))
  } else if(! is.null(x)){
    if(is.matrix(x))
      x <- as.data.frame(x)
    stopifnot(is.data.frame(x))
  }

  if(! is.null(file))
    x <- utils::read.table(file=file, header=TRUE, sep="\t", row.names=1)
  x <- convertFactorColumnsToCharacter(x)
  x <- as.matrix(x)
  nb.snps <- nrow(x)
  nb.genos <- ncol(x)
  if(verbose > 0){
    msg <- paste0("reformat ", nb.snps, " SNPs for ",
                  nb.genos, " genotypes...")
    write(msg, stdout())
  }

  ## remove the allele separator
  if(sep != ""){
    x <- apply(x, 2, gsub, pattern=sep, replacement="")
  }

  ## insure that each genotypic class is sorted by alleles
  allowed.genoClasses <- c("AA", "AC", "AG", "AT",
                           "CA", "CC", "CG", "CT",
                           "GA", "GC", "GG", "GT",
                           "TA", "TC", "TG", "TT",
                           na.string, NA)
  tmp <- data.frame(bad=c("CA","GA","TA","GC","TC","TG"),
                    good=c("AC","AG","AT","CG","CT","GT"),
                    stringsAsFactors=FALSE)
  for(i in 1:nrow(tmp)){
    is.bad <- x == tmp$bad[i]
    if(any(is.bad, na.rm=TRUE))
      x[is.bad] <- tmp$good[i]
  }

  ## replace not allowed genotypic classes with NA
  distinct.genoClasses <- unique(c(x))
  if(verbose > 0){
    msg <- paste0("distinct genotypic classes: ", length(distinct.genoClasses))
    write(msg, stdout())
  }
  is.bad <- ! distinct.genoClasses %in% allowed.genoClasses
  if(any(is.bad)){
    if(verbose > 0){
      msg <- paste0(sum(is.bad), " not allowed (e.g. ",
                    distinct.genoClasses[is.bad][1],
                    "): replaced by NA")
      write(msg, stdout())
    }
    for(not.allowed.gC in distinct.genoClasses[is.bad])
      x[x == not.allowed.gC] <- NA
  }

  ## replace na.string with NA
  idx.na.string <- which(x == na.string)
  if(length(idx.na.string) > 0)
    x[idx.na.string] <- NA

  return(x)
}

##' Convert genotypes
##'
##' Convert SNP genotypes from "genotypic classes" into "allele doses" by counting the number of minor alleles.
##' The code is not particularly efficient, but at least it exists.
##' @param x matrix or data.frame with bi-allelic SNPs in rows and genotypes in columns
##' @param na.string a character to be interpreted as NA values
##' @param snpIDs1stCol indicates if the SNP identifiers are in the first column; otherwise they should be in the row names
##' @param errorIfNotBi raises an error if a SNP is not bi-allelic; otherwise, raises a warning and set the SNP genotypes to missing
##' @param verbose verbosity level (0/1)
##' @return list of a matrix of allele doses with SNPs in columns and genotypes in rows, and a matrix with major and minor alleles
##' @author Timothee Flutre
##' @export
genoClasses2genoDoses <- function(x, na.string="--", snpIDs1stCol=TRUE,
                                  errorIfNotBi=TRUE, verbose=1){
  stopifnot(is.matrix(x) || is.data.frame(x),
            ! is.null(colnames(x)))
  if(! snpIDs1stCol)
    stopifnot(! is.null(rownames(x)))

  if(snpIDs1stCol){
    rownames(x) <- x[,1]
    x[,1] <- NULL
  }
  snp.names <- rownames(x)
  P <- length(snp.names)
  ind.names <- colnames(x)
  N <- length(ind.names)
  if(verbose > 0){
    msg <- paste0("convert to doses ", P, " SNPs and ", N, " genotypes...")
    write(msg, stdout())
  }

  geno.doses <- matrix(data=NA, nrow=N, ncol=P,
                       dimnames=list(ind.names, snp.names))
  alleles <- matrix(data=NA, nrow=P, ncol=2,
                    dimnames=list(snp.names, c("major", "minor")))

  for(p in 1:P){ # for each SNP
    raw.genos <- as.character(unlist(x[p,]))
    raw.genos[raw.genos == na.string] <- NA
    if(all(is.na(raw.genos))){
      next
    }
    tmp <- do.call(c, strsplit(raw.genos[! is.na(raw.genos)], ""))
    distinct.alleles <- sort(unique(tmp))
    allele.counts <- sort(table(tmp))
    if(length(distinct.alleles) > 2){ # SNP with more than 2 alleles
      msg <- paste0("SNP ", rownames(x)[p], " has more than 2 alleles")
      if(errorIfNotBi)
        stop(msg)
      else{
        warning(msg, immediate.=TRUE)
        raw.genos <- rep(NA, length(raw.genos))
      }
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

##' JoinMap format
##'
##' The \href{https://www.kyazma.nl/index.php/JoinMap/}{JoinMap} software needs a locus genotype file ("loc-file") containing the information (genotypes, phased or not) of the loci for a single segregating population.
##' This function converts a data.frame in the JoinMap v2 format into the JoinMap v3-v4 format.
##' It is assumed that the population is of type CP: "a population resulting from a cross between two heterogeneously heterozygous and homozygous diploid parents, linkage phases originally (possibly) unknown".
##' @param x data.frame in the JoinMap v2 format; the first four lines of a JoinMap loc-file should be absent; columns \code{<locus name>}, \code{[PHASE]} and \code{[CLAS]} are ignored; missing genotypes are coded as "--"
##' @param verbose verbosity level (0/1)
##' @return data.frame
##' @author Timothee Flutre
##' @export
updateJoinMap <- function(x, verbose=1){
  stopifnot(is.data.frame(x),
            all(colnames(x)[1:4] == c("locus", "seg", "phase", "clas")))

  x <- convertFactorColumnsToCharacter(x)

  output <- x

  if(verbose > 0)
    pb <- utils::txtProgressBar(min=0, max=nrow(x), style=3)

  for(i in 1:nrow(output)){
    stopifnot(x$seg[i] %in% c("<abxcd>", "<abxac>", "<abxab>", "<abxaa>",
                              "<aaxab>"))

    if(x$seg[i] == "<abxac>"){
      if(! all(x[i, 5: ncol(x)] %in% c("aa", "ac", "ca", "ba", "ab", "bc",
                                       "cb", "--"))){
        msg <- paste0("wrong genotype(s) at row ", i)
        stop(msg)
      }
      output$seg[i] <- "<efxeg>"
      for(j in 5:ncol(x))
        if(x[i, j] == "aa"){
          output[i, j] <- "ee"
        } else if(x[i, j] == "ac"){
          output[i, j] <- "eg"
        } else if(x[i, j] == "ca"){
          output[i, j] <- "ge"
        } else if(x[i, j] == "ba"){
          output[i, j] <- "fe"
        } else if(x[i, j] == "ab"){
          output[i, j] <- "ef"
        } else if(x[i, j] == "bc"){
          output[i, j] <- "fg"
        } else if(x[i, j] == "cb"){
          output[i, j] <- "gf"
        }

    } else if(x$seg[i] == "<abxab>"){
      if(! all(x[i, 5: ncol(x)] %in% c("aa", "ab", "ba", "bb", "--"))){
        msg <- paste0("wrong genotype(s) at row ", i)
        stop(msg)
      }
      output$seg[i] <- "<hkxhk>"
      for(j in 5:ncol(x))
        if(x[i, j] == "aa"){
          output[i, j] <- "hh"
        } else if(x[i, j] == "ab"){
          output[i, j] <- "hk"
        } else if(x[i, j] == "ba"){
          output[i, j] <- "kh"
        } else if(x[i, j] == "bb"){
          output[i, j] <- "kk"
        }

    } else if(x$seg[i] == "<abxaa>"){
      if(! all(x[i, 5: ncol(x)] %in% c("aa", "ba", "ab", "--"))){
        msg <- paste0("wrong genotype(s) at row ", i)
        stop(msg)
      }
      output$seg[i] <- "<lmxll>"
      for(j in 5:ncol(x))
        if(x[i, j] == "aa"){
          output[i, j] <- "ll"
        } else if(x[i, j] == "ba"){
          output[i, j] <- "ml"
        } else if(x[i, j] == "ab"){
          output[i, j] <- "lm"
        }

    } else if(x$seg[i] == "<aaxab>"){
      if(! all(x[i, 5: ncol(x)] %in% c("aa", "ab", "ba", "--"))){
        msg <- paste0("wrong genotype(s) at row ", i)
        stop(msg)
      }
      output$seg[i] <- "<nnxnp>"
      for(j in 5:ncol(x))
        if(x[i, j] == "aa"){
          output[i, j] <- "nn"
        } else if(x[i, j] == "ab"){
          output[i, j] <- "np"
        } else if(x[i, j] == "ba"){
          output[i, j] <- "pn"
        }

    }
    if(verbose > 0)
      utils::setTxtProgressBar(pb, i)
  }
  if(verbose > 0)
    close(pb)

  return(output)
}

##' Convert genotypes
##'
##' Convert SNP genotypes of an outbred bi-parental cross from genotypic classes into the \href{https://www.kyazma.nl/index.php/JoinMap/}{JoinMap} format.
##' For the moment, missing genotypes in parents result in the SNP being ignored, but we could imagine using genotypes in offsprings to impute such cases.
##' @param x data.frame of bi-allelic SNP genotypes, with SNPs in rows and genotypes in columns; row names should contain SNP identifiers, the first column should contain the SNP genotypes of the first parent (traditionnaly the mother), the second column should contain the SNP genotypes of the second parent (traditionnaly the father), and the remaining columns should contain the SNP genotypes of the offsprings (full siblings); if it is a matrix, it will be silently transformed into a data frame
##' @param reformat.input if TRUE, the function \code{\link{reformatGenoClasses}} will be used
##' @param na.string a character to be interpreted as NA values by \code{\link{reformatGenoClasses}}
##' @param thresh.counts threshold (per SNP) on the number of offsprings having a particular genotypic class, below which counts are converted into \code{NA}
##' @param thresh.na threshold (per SNP) on the number of offsprings having \code{NA} (applied after \code{thresh.counts})
##' @param is.F2 if TRUE, it is assumed that the two outbred parents were crossed once to make a F1 which was then autofecundated multiple times to make the offsprings; note that the SNP genotypes of the F1 should not be given (and they are often unavailable)
##' @param verbose verbosity level (0/1)
##' @return data.frame
##' @author Timothee Flutre
##' @seealso \code{\link{reformatGenoClasses}}, \code{\link{writeSegregJoinMap}}
##' @examples
##' \dontrun{
##' nb.snps <- 6
##' x <- data.frame(par1=c("AA", "GC", "CG", "AT", NA, "AA"),
##'                 par2=c("AT", "GC", "GG", "AT", "AT", "AT"),
##'                 off1=c("AA", "GG", "CG", "AA", "AA", "AT"),
##'                 off2=c("AT", "GG", "CG", "AT", "AT", "AA"),
##'                 off3=c("AT", "GG", "GG", "TT", "TT", NA),
##'                 off4=c(NA, NA, NA, NA, NA, NA),
##'                 row.names=paste0("snp", 1:nb.snps),
##'                 stringsAsFactors=FALSE)
##' genoClasses2JoinMap(x=x, reformat.input=TRUE, thresh.na=2, verbose=1)
##' }
##' @export
genoClasses2JoinMap <- function(x, reformat.input=TRUE, na.string="--",
                                thresh.counts=NULL, thresh.na=NULL,
                                is.F2=FALSE, verbose=1){
  stopifnot(is.logical(reformat.input))
  if(! is.null(thresh.counts))
    stopifnot(is.numeric(thresh.counts),
              thresh.counts < ncol(x) - 2)
  if(! is.null(thresh.na))
    stopifnot(is.numeric(thresh.na),
              thresh.na < ncol(x) - 2)
  if(reformat.input){
    x <- as.data.frame(reformatGenoClasses(x=x, na.string=na.string,
                                           verbose=verbose))
  } else{
    if(is.matrix(x))
      x <- as.data.frame(x)
    stopifnot(is.data.frame(x),
              ! is.null(rownames(x)))
  }
  x <- convertFactorColumnsToCharacter(x)
  if(verbose > 0){
    msg <- paste0("nb of offsprings: ", ncol(x) - 2)
    write(msg, stdout())
  }

  output <- data.frame(p1=x[,1],
                       p2=x[,2],
                       p1.A=NA,
                       p1.B=NA,
                       p2.C=NA,
                       p2.D=NA,
                       seg.pars=NA,
                       seg.offs=NA,
                       seg=NA,
                       stringsAsFactors=FALSE)
  colnames(output)[1:2] <- colnames(x)[1:2]
  output <- cbind(output, x[,-c(1,2)])

  ## -------------------------------------------------------
  ## vectorized version:

  ## determine parental alleles
  output[, c("p1.A","p1.B")] <- do.call(rbind, strsplit(x[,1], ""))
  output[, c("p2.C","p2.D")] <- do.call(rbind, strsplit(x[,2], ""))

  ## assess conditions before identifying parental segregations
  par1.hom <- output[,"p1.A"] == output[,"p1.B"]
  par2.hom <- output[,"p2.C"] == output[,"p2.D"]
  is.na.pars <- apply(x[, 1:2], 1, function(in.row){any(is.na(in.row))})
  if(nrow(x) == 1){
    tmp <- apply(output[, 3:6], 1, unique)
    par.alleles <- list()
    par.alleles[[colnames(tmp)]] <- tmp[,1]
  } else
    par.alleles <- apply(output[, 3:6], 1, unique)
  par.alleles <- lapply(par.alleles, sort, na.last=NA)
  two.par.alleles <- sapply(par.alleles, length) == 2

  ## assess conditions before identifying offspring segregations
  if(nrow(x) == 1){
    tmp <- apply(x[, 3:ncol(x)], 1, table, useNA="always")
    nb.gclasses.offs <- list()
    nb.gclasses.offs[[colnames(tmp)]] <- tmp[,1]
  } else
    nb.gclasses.offs <- apply(x[, 3:ncol(x)], 1, table, useNA="always")
  if(! is.null(thresh.counts)){
    below.thresh <- sapply(nb.gclasses.offs, function(nb.gc){
      any(nb.gc < thresh.counts)
    })
    if(any(below.thresh))
      nb.gclasses.offs <- lapply(nb.gclasses.offs, function(nb.gc){
        to.rmv <- ! is.na(names(nb.gc)) & nb.gc < thresh.counts
        nb.gc[is.na(names(nb.gc))] <- nb.gc[is.na(names(nb.gc))] +
          sum(nb.gc[to.rmv])
        nb.gc[! to.rmv]
      })
  }
  pass.na <- rep(TRUE, length(nb.gclasses.offs))
  if(! is.null(thresh.na))
    pass.na <- sapply(nb.gclasses.offs, function(nb.gc){
      nb.gc[is.na(names(nb.gc))] < thresh.na
    })

  if(is.F2){
    ## identify proper segregation per SNP based on parents
    proper.seg.pars <- ! is.na.pars & two.par.alleles

    ## identify proper segregation per SNP based on offsprings
    proper.seg.offs <- sapply(nb.gclasses.offs, length) >= 2 &
      sapply(nb.gclasses.offs, length) <= 4 & pass.na

    ## identify segregations per SNP coherent for parents and offsprings
    proper.seg <- proper.seg.offs & proper.seg.pars
    output$seg[proper.seg] <- "F2"

    ## code offspring genotypes when unproper segregation
    output[! proper.seg, 10:ncol(output)] <- NA

    ## code offspring genotypes when heterozygotes
    offs.het <- apply(output[proper.seg, 10:ncol(output)], 2, function(in.col){
      sapply(strsplit(in.col, split=""), function(alleles){
        sum(! is.na(unique(alleles))) == 2
      })
    })
    output[proper.seg, 10:ncol(output)][offs.het] <- "H"

    ## possible.F1 <- lapply(1:sum(proper.seg), function(i){
    ##   par.all <- output[proper.seg, 3:6][i,]
    ##   unique(sort(c(paste(sort(par.all[c("p1.A","p2.C")]), collapse=""),
    ##                 paste(sort(par.all[c("p1.A","p2.D")]), collapse=""),
    ##                 paste(sort(par.all[c("p1.B","p2.C")]), collapse=""),
    ##                 paste(sort(par.all[c("p1.B","p2.D")]), collapse=""))))
    ## })
    ## names(possible.F1) <- rownames(output[proper.seg,])

    ## code offspring genotypes when homozygotes from parent 1
    offs.hom1 <- matrix(data=FALSE, nrow=sum(proper.seg),
                        ncol=length(10:ncol(output)))
    rownames(offs.hom1) <- names(proper.seg[proper.seg])
    colnames(offs.hom1) <- colnames(output)[10:ncol(output)]
    if(verbose > 0){
      pb <- utils::txtProgressBar(min=0, max=nrow(offs.hom1), style=3)
      pb.tmp <- stats::setNames(object=cut(x=1:nrow(offs.hom1), breaks=10,
                                           labels=FALSE),
                                nm=1:nrow(offs.hom1)) # update pb <= 10 times
    }
    for(i in 1:nrow(offs.hom1)){ # for each SNP
      offs.genos <- unlist(output[proper.seg, 10:ncol(output)][i,])
      idx.hom <- which(! is.na(offs.genos) & ! offs.het[i,])
      if(length(idx.hom) > 0){
        ## identify from which parent(s) each allele comes from
        offs.genos.hom <- offs.genos[idx.hom]
        uniq.offs.homs <- unique(sort(unlist(offs.genos.hom)))
        uniq.offs.homs.alleles <- sapply(strsplit(uniq.offs.homs, ""), `[`, 1)
        origin.uniq.offs.homs.alleles <-
          lapply(uniq.offs.homs.alleles, function(a){
            tmp <- c()
            if(a %in% output[proper.seg, c("p1.A","p1.B")][i,])
              tmp <- c(tmp, "par1")
            if(a %in% output[proper.seg, c("p2.C","p2.D")][i,])
              tmp <- c(tmp, "par2")
            tmp
          })
        names(origin.uniq.offs.homs.alleles) <- uniq.offs.homs.alleles
        ## assign one allele to parent 1
        has.par1 <- sapply(origin.uniq.offs.homs.alleles, function(pars){
          "par1" %in% pars
        })
        allele.par1 <- names(has.par1[has.par1])[1] # arbitrary choice
        output$seg.pars[proper.seg][i] <- paste0("A=", allele.par1)
        ## fill offs.hom1
        geno.hom.par1 <- paste(rep(allele.par1, 2), collapse="")
        is.hom.par1 <- offs.genos.hom == geno.hom.par1
        if(any(is.hom.par1, na.rm=TRUE))
          offs.hom1[i, names(is.hom.par1[is.hom.par1])] <- TRUE
      }
      if(verbose > 0){
        ## utils::setTxtProgressBar(pb, i)
        if(i %in% as.numeric(names(pb.tmp)[cumsum(table(pb.tmp))]))
          utils::setTxtProgressBar(pb, which(1:nrow(offs.hom1) == i))
      }
    }
    if(verbose > 0)
      close(pb)
    output[proper.seg, 10:ncol(output)][offs.hom1] <- "A"

    ## code offspring genotypes when homozygotes from parent 2
    offs.na <- is.na(output[proper.seg, 10:ncol(output)])
    output[proper.seg, 10:ncol(output)][! offs.na &
                                        ! offs.het &
                                        ! offs.hom1] <- "B"

  } else{ # is.F2 == FALSE
    ## identify proper segregation per SNP based on parents
    proper.seg.pars <- ! (par1.hom & par2.hom) & ! is.na.pars & two.par.alleles

    ## identify proper segregation per SNP based on offsprings
    proper.seg.offs <- sapply(nb.gclasses.offs, length) >= 3 &
      sapply(nb.gclasses.offs, length) <= 4 & pass.na

    ## among proper segregations, determine type per SNP based on parents
    pars.hkhk <- proper.seg.pars & x[,1] == x[,2] & (! par1.hom)
    pars.lmll <- proper.seg.pars & x[,1] != x[,2] & (! par1.hom) & par2.hom
    pars.nnnp <- proper.seg.pars & x[,1] != x[,2] & par1.hom & (! par2.hom)
    output[pars.hkhk, "seg.pars"] <- "<hkxhk>"
    output[pars.lmll, "seg.pars"] <- "<lmxll>"
    output[pars.nnnp, "seg.pars"] <- "<nnxnp>"

    ## among proper segregations, determine type per SNP based on offsprings
    offs.hkhk <- proper.seg.offs & sapply(nb.gclasses.offs, length) == 4
    offs.lmll.or.nnnp <- proper.seg.offs & sapply(nb.gclasses.offs, length) == 3
    output[offs.hkhk, "seg.offs"] <- "<hkxhk>"
    output[offs.lmll.or.nnnp, "seg.offs"] <- "<lmxll>_<nnxnp>"

    ## identify segregations per SNP coherent for parents and offsprings
    proper.seg <- proper.seg.offs & proper.seg.pars

    ## determine coherent segregations among parents and offsprings
    is.hkhk <- (proper.seg.pars & pars.hkhk) &
      (proper.seg.offs & offs.hkhk)
    is.hkhk[is.na(is.hkhk)] <- FALSE
    is.lmll <- (proper.seg.pars & pars.lmll) &
      (proper.seg.offs & offs.lmll.or.nnnp)
    is.lmll[is.na(is.lmll)] <- FALSE
    is.nnnp <- (proper.seg.pars & pars.nnnp) &
      (proper.seg.offs & offs.lmll.or.nnnp)
    is.nnnp[is.na(is.nnnp)] <- FALSE
    if(verbose > 0){
      msg <- paste0("coherent segregations: ", sum(is.hkhk) + sum(is.lmll) +
                                               sum(is.nnnp),
                    "\n<hkxhk>=", sum(is.hkhk),
                    " <lmxll>=", sum(is.lmll),
                    " <nnxnp>=", sum(is.nnnp))
      write(msg, stdout())
    }

    ## set proper segregation classes
    output[is.hkhk, 1:2] <- "hk"
    output[is.lmll, 1] <- "lm"
    output[is.lmll, 2] <- "ll"
    output[is.nnnp, 1] <- "nn"
    output[is.nnnp, 2] <- "np"

    ## set segregation type per SNP
    idx <- which(output[,1] %in% c("hk", "lm", "nn"))
    output$seg[idx] <- paste0("<", output[idx,1], "x", output[idx,2], ">")

    ## set segregation classes of offsprings when seg type is <hkxhk>
    gclasses.hh <- paste0(output[is.hkhk, "p1.A"], output[is.hkhk, "p2.C"])
    tmp1 <- t(apply(cbind(gclasses.hh, output[is.hkhk, 10:ncol(output)]), 1,
                   function(in.row){in.row[-1] == in.row[1]}))
    tmp1[is.na(tmp1)] <- FALSE
    if(any(tmp1))
      output[is.hkhk, 10:ncol(output)][tmp1] <- "hh"
    gclasses.kk <- paste0(output[is.hkhk, "p1.B"], output[is.hkhk, "p2.D"])
    tmp2 <- t(apply(cbind(gclasses.kk, output[is.hkhk, 10:ncol(output)]), 1,
                   function(in.row){in.row[-1] == in.row[1]}))
    tmp2[is.na(tmp2)] <- FALSE
    if(any(tmp2))
      output[is.hkhk, 10:ncol(output)][tmp2] <- "kk"
    gclasses.hk <- paste0(output[is.hkhk, "p1.A"], output[is.hkhk, "p2.D"])
    tmp3 <- t(apply(cbind(gclasses.hk, output[is.hkhk, 10:ncol(output)]), 1,
                   function(in.row){in.row[-1] == in.row[1]}))
    tmp3[is.na(tmp3)] <- FALSE
    if(any(tmp3))
      output[is.hkhk, 10:ncol(output)][tmp3] <- "hk"
    ## handle bad offspring segregations and thresh.counts != NULL
    tmp4 <- ! (tmp1 | tmp2 | tmp3)
    if(any(tmp4))
      output[is.hkhk, 10:ncol(output)][tmp4] <- NA

    ## set segregation classes of offsprings when seg type is <lmxll> or <nnxnp>
    lmll.or.nnnp <- is.lmll | is.nnnp
    tmp1 <- t(apply(cbind(x[,1], output[, 10:ncol(output)])[lmll.or.nnnp,], 1,
                    function(in.row){in.row[-1] == in.row[1]}))
    tmp1[is.na(tmp1)] <- FALSE
    if(any(tmp1))
      output[lmll.or.nnnp, 10:ncol(output)][tmp1] <-
        rep(output[lmll.or.nnnp, 1], ncol(tmp1))[tmp1]
    tmp2 <- t(apply(cbind(x[,2], output[, 10:ncol(output)])[lmll.or.nnnp,], 1,
                    function(in.row){in.row[-1] == in.row[1]}))
    tmp2[is.na(tmp2)] <- FALSE
    if(any(tmp2))
      output[lmll.or.nnnp, 10:ncol(output)][tmp2] <-
        rep(output[lmll.or.nnnp, 2], ncol(tmp2))[tmp2]
    ## handle bad offspring segregations and thresh.counts != NULL
    tmp3 <- ! xor(tmp1, tmp2)
    if(any(tmp3))
      output[lmll.or.nnnp, 10:ncol(output)][tmp3] <- NA
  }

  ## -------------------------------------------------------
  ## un-vectorized version:

  ## if(verbose > 0)
  ##   pb <- utils::txtProgressBar(min=0, max=nrow(x), style=3)

  ## for(i in 1:nrow(x)){
  ##   if(verbose > 0)
  ##     utils::setTxtProgressBar(pb, i)

  ##   ## count nb genotypic classes in offsprings
  ##   counts <- table(as.character(x[i, 3:ncol(x)]), useNA="always")

  ##   ## discard if no segregation
  ##   if(length(counts) <= 2)
  ##     next
  ##   ## discard if more than 2 parental alleles
  ##   if(length(unique(as.character(output[i,3:6]))) > 2)
  ##     next

  ##   ## if no missing genotype in parents, determine the segregation type
  ##   ## by looking only at the parental genotypes
  ##   if(! any(is.na(x[i, 1:2]))){

  ##     ## hkxhk
  ##     if(all(x[i,1] == x[i,2],
  ##            output[i,"p1.A"] != output[i,"p1.B"])){
  ##       output[i,1] <- "hk"
  ##       output[i,2] <- "hk"

  ##       tmp <- paste0(output[i,"p1.A"], output[i,"p2.C"])
  ##       idx <- which(output[i, 8:ncol(output)] == tmp)
  ##       if(length(idx) > 0)
  ##         output[i, 7+idx] <- rep("hh", length(idx))
  ##       tmp <- paste0(output[i,"p1.B"], output[i,"p2.D"])
  ##       idx <- which(output[i, 8:ncol(output)] == tmp)
  ##       if(length(idx) > 0)
  ##         output[i, 7+idx] <- rep("kk", length(idx))

  ##     } else{ ## lmxll or nnxnp
  ##       if(all(output[i,"p1.A"] != output[i,"p1.B"],
  ##              output[i,"p2.C"] == output[i,"p2.D"])){
  ##         output[i,1] <- "lm"
  ##         output[i,2] <- "ll"

  ##       } else{
  ##         output[i,1] <- "nn"
  ##         output[i,2] <- "np"
  ##       }
  ##     }

  ##     ## fill offspring columns
  ##     idx <- which(output[i, 8:ncol(output)] == x[i,1])
  ##     if(length(idx) > 0)
  ##       output[i, 7+idx] <- rep(output[i,1], length(idx))
  ##     idx <- which(output[i, 8:ncol(output)] == x[i,2])
  ##     if(length(idx) > 0)
  ##       output[i, 7+idx] <- rep(output[i,2], length(idx))
  ##   }
  ## }
  ## if(verbose > 0)
  ##   close(pb)

  ## ## fill the "seg" column
  ## idx <- which(output[,1] %in% c("nn", "lm", "hk"))
  ## output$seg[idx] <- paste0("<", output[idx,1], "x", output[idx,2], ">")
  ## if(verbose > 0){
  ##   msg <- paste0("proper segregations: ", sum(! is.na(output$seg)),
  ##                 " / ", nrow(output))
  ##   write(msg, stdout())
  ## }

  return(output)
}

##' Filter for segregation
##'
##' Filter genotypes based on the chi2 statistic to test for segregation distortion.
##' @param x data.frame similar to the output from \code{\link{genoClasses2JoinMap}}, with locus in rows and where the first columns corresponding to parents are discarded, i.e. the first column should be named "seg" and the following should correspond to offsprings
##' @param thresh.pval significance threshold at which to control the FWER (Bonferroni-adjusted p values) and the FDR (Benjamini-Hochberg-adjusted p values); the number of non-rejected null hypotheses will be shown only if \code{verbose} is non-zero
##' @param return.counts if TRUE, counts used to calculate the chi2 statistic are returned
##' @param is.F2 if TRUE, it is assumed that the two outbred parents were crossed once to make a F1 which was then autofecundated multiple times to make the offsprings; note that the SNP genotypes of the F1 should not be given (and they are often unavailable)
##' @param nb.cores the number of cores to use, i.e. at most how many child processes will be run simultaneously (not on Windows)
##' @param verbose verbosity level (0/1)
##' @return data.frame with one row per locus and several columns (chi2, p value, Bonferroni-adjusted p value, Benjamini-Hochberg-adjusted p value), as well as counts (if \code{return.counts=TRUE})
##' @author Timothee Flutre
##' @seealso \code{\link{genoClasses2JoinMap}}, \code{\link{updateJoinMap}}, \code{\link{plotHistPval}}, \code{\link{qqplotPval}}
##' @export
filterSegreg <- function(x, thresh.pval=0.05, return.counts=FALSE,
                         is.F2=FALSE, nb.cores=1, verbose=1){
  stopifnot(is.data.frame(x),
            colnames(x)[1] == "seg",
            is.numeric(thresh.pval),
            length(thresh.pval) == 1,
            all(thresh.pval >= 0, thresh.pval <= 1))
  if(is.F2){
    stopifnot(all(x$seg %in% c("F2", "NA", NA)),
              all(unlist(x[,-1]) %in% c("A", "H", "B", "NA", NA)))
  } else
    stopifnot(all(x$seg %in% c("<nnxnp>", "<lmxll>", "<hkxhk>", "<efxeg>",
                               "<abxcd>", "NA", NA)))

  x <- convertFactorColumnsToCharacter(x)

  idx.rows <- (1:nrow(x))[! is.na(x$seg)]
  if(verbose > 0){
    msg <- paste0("nb of segregations to test: ", length(idx.rows))
    write(msg, stdout())
  }

  info.seg <- list("<nnxnp>"=list(classes=sort(c("nn", "np")),
                                  exp=rep(0.5, 2)),
                   "<lmxll>"=list(classes=sort(c("ll", "lm")),
                                  exp=rep(0.5, 2)),
                   "<hkxhk>"=list(classes=sort(c("hh", "hk", "kk")),
                                  exp=c(0.25, 0.5, 0.25)),
                   "<efxeg>"=list(classes=sort(c("ee", "eg", "ef", "fg")),
                                  exp=rep(0.25, 4)),
                   "<abxcd>"=list(classes=sort(c("ac", "ad", "bc", "bd")),
                                  exp=rep(0.25, 4)),
                   "F2"=list(classes=c("A", "H", "B"),
                             exp=c(0.25, 0.5, 0.25)))
  if(is.F2){
    info.seg <- info.seg[names(info.seg) %in% "F2"]
  } else
    info.seg <- info.seg[names(info.seg) %in% x$seg]

  output <- data.frame(seg=x$seg,
                       nb.classes=NA,
                       row.names=rownames(x),
                       stringsAsFactors=FALSE)
  for(cn in c(paste0("class", 1:4), paste0("obs", 1:4), paste0("exp", 1:4),
              "chi2", "pvalue"))
    output[[cn]] <- NA

  ## -------------------------------------------------------
  ## parallelized and vectorized version:

  tmp.out <- do.call(rbind, parallel::mclapply(names(info.seg), function(seg){
    idx.seg <- which(x$seg == seg)
    nb.classes <- length(info.seg[[seg]]$classes)
    output$nb.classes[idx.seg] <- nb.classes
    output[idx.seg, 2+(1:nb.classes)] <- rep(info.seg[[seg]]$classes,
                                             each=length(idx.seg))
    obs.counts <- t(apply(x[idx.seg, -1, drop=FALSE], 1, function(obs.classes){
      table(factor(obs.classes, levels=info.seg[[seg]]$classes), useNA="no")
    }))
    output[idx.seg, 2+4+(1:nb.classes)] <- obs.counts
    exp.counts <- rep(info.seg[[seg]]$exp, each=length(idx.seg)) *
      rowSums(obs.counts)
    output[idx.seg, 2+4+4+(1:nb.classes)] <- exp.counts
    output[idx.seg,]
  }, mc.cores=nb.cores))
  output <- rbind(tmp.out, output[is.na(x$seg),])
  output <- output[rownames(x),]

  ## -------------------------------------------------------
  ## vectorized version:

  ## zz <- lapply(names(info.seg), function(seg){
  ##   idx.seg <- which(x$seg == seg)
  ##   nb.classes <- length(info.seg[[seg]]$classes)
  ##   output$nb.classes[idx.seg] <<- nb.classes
  ##   output[idx.seg, 2+(1:nb.classes)] <<- rep(info.seg[[seg]]$classes,
  ##                                             each=length(idx.seg))
  ##   obs.counts <- t(apply(x[idx.seg, -1, drop=FALSE], 1, function(obs.classes){
  ##     table(factor(obs.classes, levels=info.seg[[seg]]$classes), useNA="no")
  ##   }))
  ##   output[idx.seg, 2+4+(1:nb.classes)] <<- obs.counts
  ##   exp.counts <- rep(info.seg[[seg]]$exp, each=length(idx.seg)) *
  ##     rowSums(obs.counts)
  ##   output[idx.seg, 2+4+4+(1:nb.classes)] <<- exp.counts
  ## })

  ## -------------------------------------------------------
  ## un-vectorized version:

  ## if(verbose > 0){
  ##   pb <- utils::txtProgressBar(min=0, max=length(idx.rows), style=3)
  ##   tmp <- stats::setNames(object=cut(x=idx.rows, breaks=10, labels=FALSE),
  ##                          nm=idx.rows) # to update pb no more than 10 times
  ## }
  ## for(i in idx.rows){
  ##   nb.classes <- length(seg2classes[[x$seg[i]]])
  ##   output[i, "nb.classes"] <- nb.classes
  ##   output[i, 2+(1:nb.classes)] <- seg2classes[[x$seg[i]]]
  ##   obs.classes <- factor(x=x[i,-1], levels=seg2classes[[x$seg[i]]])
  ##   counts <- table(obs.classes, useNA="no")
  ##   output[i, 2+4+(1:nb.classes)] <- counts

  ##   ## expected counts assuming no segregation distortion
  ##   if(x$seg[i] %in% c("<nnxnp>", "<lmxll>")){
  ##     output[i, 2+4+4+(1:2)] <- rep(0.5 * sum(counts), 2)
  ##   } else if(x$seg[i] == "<hkxhk>"){
  ##     output[i, 2+4+4+(1:3)] <- c(0.25 * sum(counts),
  ##                                 0.5 * sum(counts),
  ##                                 0.25 * sum(counts))
  ##   } else if(x$seg[i] %in% c("<efxeg>", "<abxcd>"))
  ##     output[i, 2+4+4+(1:4)] <- c(0.25 * sum(counts),
  ##                                 0.25 * sum(counts),
  ##                                 0.25 * sum(counts),
  ##                                 0.25 * sum(counts))

  ##   if(verbose > 0){
  ##     ## utils::setTxtProgressBar(pb, i) # update pb at each loop iteration
  ##     if(i %in% as.numeric(names(tmp)[cumsum(table(tmp))]))
  ##       utils::setTxtProgressBar(pb, which(idx.rows == i))
  ##   }
  ## }
  ## if(verbose > 0)
  ##   close(pb)

  ## calculate chi2 and p value
  output$chi2 <- rowSums((output[,2+4+(1:4)] - output[,2+4+4+(1:4)])^2 /
                         output[,2+4+4+(1:4)], na.rm=TRUE)
  output$chi2[is.na(output$nb.classes)] <- NA
  output$pvalue <- stats::pchisq(q=output$chi2, df=output$nb.classes - 1,
                                 lower.tail=FALSE)

  ## adjust p values
  output$pvalue.bonf <- stats::p.adjust(p=output$pvalue, method="bonferroni")
  output$pvalue.bh <- stats::p.adjust(p=output$pvalue, method="BH")

  if(verbose > 0){
    idx <- ! is.na(x$seg)
    msg <- paste0("nb of non-rejected H0 at ",
                  format(100*thresh.pval, digits=2), "%",
                  " (Bonferroni): ", sum(output$pvalue.bonf[idx] >
                                         thresh.pval))
    msg <- paste0(msg, "\nnb of non-rejected H0 at ",
                  format(100*thresh.pval, digits=2), "%",
                  " (BH): ", sum(output$pvalue.bh[idx] >
                                 thresh.pval))
    write(msg, stdout())
  }

  idx <- 1:ncol(output)
  if(! return.counts)
    idx <- -(1:(grep("chi2", colnames(output)) - 1))
  return(output[, idx])
}

##' Test for linkage
##'
##' Computes the LOD score(s) given recombination fraction(s) as well as number of non-recombinants and recombinants, testing for "ref.frac = 0".
##' @param rec.frac vector of recombination fraction(s); if above 0.5, will use 1 - rec.frac; will be recycled to match nb.nonrec and nb.rec
##' @param nb.nonrec vector of number(s) of non-recombinant offsprings; will be recycled to match rec.frac
##' @param nb.rec vector of number(s) of recombinant offsprings; will be recycled to match rec.frac
##' @return vector of LOD score(s)
##' @author Timothee Flutre
##' @examples
##' \dontrun{## from "Genetic map construction with R/qtl" by Karl Broman (2012)
##' ## page 10: http://www.rqtl.org/tutorials/geneticmaps.pdf
##' n <- 300
##' nb.rec <- 27 + 27
##' nb.nonrec <- 300 - nb.rec
##' (est.rec.frac <- nb.rec / n) # estimator for a backcross; 0.18
##' lodLinkage(est.rec.frac, nb.nonrec, nb.rec) # ~28.9
##'
##' ## show how input rec.frac is recycled
##' lodLinkage(seq(0, 0.5, 0.05), nb.nonrec, nb.rec)
##' plot(lodLinkage(seq(0, 0.5, 0.05), nb.nonrec, nb.rec), type="b",
##'      xlab="recombination fraction", ylab="LOD score", las=1)
##'
##' ## show how inputs nb.nonrec and nb.rec are recycled
##' lodLinkage(est.rec.frac, n - 60, 60)
##' lodLinkage(est.rec.frac, c(n - 54, n - 60), c(54, 60))
##' }
##' @export
lodLinkage <- function(rec.frac, nb.nonrec, nb.rec){
  stopifnot(is.numeric(rec.frac),
            is.numeric(nb.nonrec),
            is.numeric(nb.rec),
            length(nb.nonrec) == length(nb.rec))

  lod <- rep(NA, max(c(length(rec.frac), length(nb.nonrec), length(nb.rec))))

  ## handle rec.frac
  is.not.NA <- ! is.na(rec.frac)
  stopifnot(rec.frac[is.not.NA] >= 0,
            rec.frac[is.not.NA] <= 1)

  ## handle nb.nonrec and nb.rec
  if(any(! is.not.NA)){
    if(length(nb.nonrec) == 1){
      nb.nonrec <- rep(nb.nonrec, sum(is.not.NA))
      nb.rec <- rep(nb.rec, sum(is.not.NA))
    } else if(length(nb.nonrec) == length(rec.frac)){
      nb.nonrec <- nb.nonrec[is.not.NA]
      nb.rec <- nb.rec[is.not.NA]
    } else{
      msg <- "can't handle ref.frac with NA's and length(nb.nonrec) != length(rec.frac)"
      stop(msg)
    }
  }

  lod[is.not.NA] <- nb.nonrec * log10(1 - rec.frac[is.not.NA]) +
    nb.rec * log10(rec.frac[is.not.NA]) +
    (nb.nonrec + nb.rec) * log10(2)

  return(lod)
}

##' JoinMap/MapQTL to R/qtl
##'
##' Return the correspondence in terms of genotype coding between the format used by \href{https://www.kyazma.nl/index.php/JoinMap/}{JoinMap} and the one used by \href{https://cran.r-project.org/package=qtl}{qtl}.
##' @return data frame
##' @author Timothee Flutre
##' @export
correspondenceJoinMap2qtl <- function(){
  out <- data.frame(
      segreg=c(rep("lmxll", 2),
               rep("nnxnp", 2),
               rep("abxcd", 4),
               rep("efxeg", 4),
               rep("hkxhk", 4)),
      phase=c("0-", "1-",
              "-0", "-1",
              "00", "01", "10", "11",
              "00", "01", "10", "11",
              "00", "01", "10", "11"),
      mother.AB=c("lm", "ml",
                  "nn", "nn",
                  "ab", "ab", "ba", "ba",
                  "ef", "ef", "fe", "fe",
                  "hk", "hk", "kh", "kh"),
      father.CD=c("ll", "ll",
                  "np", "pn",
                  "cd", "dc", "cd", "dc",
                  "eg", "ge", "eg", "ge",
                  "hk", "kh", "hk", "kh"),
      F1.AC.joinmap=c("ll", "ml",
                      "nn", "np",
                      "ac", "ad", "bc", "bd",
                      "ee", "eg", "fe", "fg",
                      "hh", "hk", "kh", "kk"),
      F1.AC.qtl=c(5, 5,
                  7, 7,
                  1, 1, 1, 1,
                  1, 1, 1, 1,
                  1, 9, 9, 1),
      F1.AD.joinmap=c("ll", "ml",
                      "np", "nn",
                      "ad", "ac", "bd", "bc",
                      "eg", "ee", "fg", "fe",
                      "hk", "hh", "kk", "kh"),
      F1.AD.qtl=c(5, 5,
                  8, 8,
                  3, 3, 3, 3,
                  3, 3, 3, 3,
                  10, 3, 3, 10),
      F1.BC.joinmap=c("ml", "ll",
                      "nn", "np",
                      "bc", "bd", "ac", "ad",
                      "fe", "fg", "ee", "eg",
                      "kh", "kk", "hh", "hk"),
      F1.BC.qtl=c(6, 6,
                  7, 7,
                  2, 2, 2, 2,
                  2, 2, 2, 2,
                  10, 2, 2, 10),
      F1.BD.joinmap=c("ml", "ll",
                      "np", "nn",
                      "bd", "bc", "ad", "ac",
                      "fg", "fe", "eg", "ee",
                      "kk", "kh", "hk", "hh"),
      F1.BD.qtl=c(6, 6,
                  8, 8,
                  4, 4, 4, 4,
                  4, 4, 4, 4,
                  4, 9, 9, 4),
      stringsAsFactors=FALSE)

  return(out)
}

##' Genotype coding
##'
##' Convert genotypes encoded in the \href{https://www.kyazma.nl/index.php/JoinMap/}{JoinMap} format for cross type "CP" with known parental phases into the format used by the \href{https://cran.r-project.org/package=qtl}{qtl} package.
##' @param x data frame in the JoinMap format, for instance from \code{\link{genoClasses2JoinMap}}; any column named "clas" will be discarded from the start; columns named "seg" and "phase" should be present; a column named "locus" should be present unless row names are given; all other columns are supposed to correspond to genotypes
##' @param verbose verbosity level (0/1/2)
##' @return matrix in the "qtl" format
##' @author Timothee Flutre
##' @seealso \code{\link{correspondenceJoinMap2qtl}}, \code{\link{updateJoinMap}}
##' @export
phasedJoinMapCP2qtl <- function(x, verbose=1){
  seg.types <- c("<abxcd>", "<efxeg>", "<hkxhk>", "<lmxll>", "<nnxnp>")
  phase.types <- c("{00}", "{01}", "{10}", "{11}",
                   "{0-}", "{1-}", "{-0}", "{-1}")
  geno.types <- c("ac", "ad", "bc", "bd", "ee", "eg", "ef", "fg",
                  "hh", "hk", "kk", "h-", "k-", "lm", "ll", "nn", "np")
  stopifnot(is.data.frame(x),
            ! is.null(colnames(x)),
            all(c("seg", "phase") %in% colnames(x)),
            all(x$seg %in% c(seg.types, NA)),
            all(x$phase %in% c(phase.types, NA)))

  if(verbose > 0){
    msg <- paste0("convert genotypes encoded in the JoinMap format",
                  " for cross type \"CP\" with known parental phases",
                  " into the format used by the 'qtl' package")
    write(msg, stdout())
  }

  ## reformat the input data frame
  x <- convertFactorColumnsToCharacter(x)
  if("clas" %in% colnames(x))
    x <- x[, - which(colnames(x) == "clas")]
  if("locus" %in% colnames(x)){
    rownames(x) <- x$locus
    x <- x[, - which(colnames(x) == "locus")]
  }
  tmp <- x[! is.na(x$seg), - which(colnames(x) %in% c("seg", "phase"))]
  stopifnot(all(unlist(tmp) %in% c(geno.types, NA)))
  is.seg.NA <- is.na(x$seg)
  if(any(is.seg.NA))
    x$phase[which(is.seg.NA)] <- NA
  is.phase.NA <- is.na(x$phase)
  if(any(is.phase.NA))
    x$seg[which(is.phase.NA)] <- NA

  ## make the empty output matrix
  locus.names <- rownames(x)
  nb.locus <- length(locus.names)
  geno.names <- colnames(x)[! colnames(x) %in% c("seg", "phase")]
  nb.genos <- length(geno.names)
  out <- matrix(data=NA, nrow=nb.locus, ncol=nb.genos,
                dimnames=list(locus.names, geno.names))

  ## fill the output matrix per seg-phase
  ## parent 1 (mother): AB
  ## parent 2 (father): CD
  ## offspring (F1): AC or AD or BC or BD
  ## https://github.com/kbroman/qtl/blob/master/R/read.cross.mq.R#L301
  pairs.seg.phase <- paste(x$seg, x$phase, sep="_")
  for(seph in unique(pairs.seg.phase)){
    if(seph == "NA_NA")
      next
    idx.rows <- which(pairs.seg.phase == seph)
    seg <- x$seg[idx.rows[1]]
    phase <- x$phase[idx.rows[1]]

    if(phase == "{0-}"){
      if(seg == "<lmxll>"){ # AB=lm ; CD=ll
        idx <- which(x[idx.rows, geno.names] == "ll")
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 5
        idx <- which(x[idx.rows, geno.names] == "lm")
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 6
      } else{
        msg <- paste("unrecognized segregation", seg, "with phase", phase)
        stop(msg)
      }
    } else if(phase == "{1-}"){
      if(seg == "<lmxll>"){ # AB=ml ; CD=ll
        idx <- which(x[idx.rows, geno.names] == "ll")
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 6
        idx <- which(x[idx.rows, geno.names] == "lm")
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 5
      } else{
        msg <- paste("unrecognized segregation", seg, "with phase", phase)
        stop(msg)
      }
    } else if(phase == "{-0}"){
      if(seg == "<nnxnp>"){ # AB=nn ; CD=np
        idx <- which(x[idx.rows, geno.names] == "nn")
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 7
        idx <- which(x[idx.rows, geno.names] == "np")
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 8
      } else{
        msg <- paste("unrecognized segregation", seg, "with phase", phase)
        stop(msg)
      }
    } else if(phase == "{-1}"){
      if(seg == "<nnxnp>"){ # AB=nn ; CD=pn
        idx <- which(x[idx.rows, geno.names] == "nn")
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 8
        idx <- which(x[idx.rows, geno.names] == "np")
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 7
      } else{
        msg <- paste("unrecognized segregation", seg, "with phase", phase)
        stop(msg)
      }
    } else if(phase == "{00}"){
      if(seg == "<abxcd>"){
        idx <- which(x[idx.rows, geno.names] == "ac") # AC
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 1
        idx <- which(x[idx.rows, geno.names] == "ad") # AD
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 3
        idx <- which(x[idx.rows, geno.names] == "bc") # BC
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 2
        idx <- which(x[idx.rows, geno.names] == "bd") # BD
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 4
      } else if(seg == "<efxeg>"){
        idx <- which(x[idx.rows, geno.names] == "ee") # AC
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 1
        idx <- which(x[idx.rows, geno.names] == "eg") # AD
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 3
        idx <- which(x[idx.rows, geno.names] == "ef") # BC
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 2
        idx <- which(x[idx.rows, geno.names] == "fg") # BD
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 4
      } else if(seg == "<hkxhk>"){
        idx <- which(x[idx.rows, geno.names] == "hh") # AC
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 1
        idx <- which(x[idx.rows, geno.names] == "hk") # AD or BC
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 10
        idx <- which(x[idx.rows, geno.names] == "kk") # BD
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 4
        idx <- which(x[idx.rows, geno.names] == "h-") # not BD
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 14
        idx <- which(x[idx.rows, geno.names] == "k-") # not AC
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 11
      } else{
        msg <- paste("unrecognized segregation", seg, "with phase", phase)
        stop(msg)
      }
    } else if(phase == "{01}"){
      if(seg == "<abxcd>"){
        idx <- which(x[idx.rows, geno.names] == "ad") # AC
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 1
        idx <- which(x[idx.rows, geno.names] == "ac") # AD
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 3
        idx <- which(x[idx.rows, geno.names] == "bd") # BC
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 2
        idx <- which(x[idx.rows, geno.names] == "bc") # BD
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 4
      } else if(seg == "<efxeg>"){
        idx <- which(x[idx.rows, geno.names] == "eg") # AC
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 1
        idx <- which(x[idx.rows, geno.names] == "ee") # AD
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 3
        idx <- which(x[idx.rows, geno.names] == "fg") # BC
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 2
        idx <- which(x[idx.rows, geno.names] == "ef") # BD
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 4
      } else if(seg == "<hkxhk>"){
        idx <- which(x[idx.rows, geno.names] == "hk") # AC or BD
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 9
        idx <- which(x[idx.rows, geno.names] == "hh") # AD
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 3
        idx <- which(x[idx.rows, geno.names] == "kk") # BC
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 2
        idx <- which(x[idx.rows, geno.names] == "h-") # not BC
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 12
        idx <- which(x[idx.rows, geno.names] == "k-") # not AD
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 13
      } else{
        msg <- paste("unrecognized segregation", seg, "with phase", phase)
        stop(msg)
      }
    } else if(phase == "{10}"){
      if(seg == "<abxcd>"){
        idx <- which(x[idx.rows, geno.names] == "bc") # AC
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 1
        idx <- which(x[idx.rows, geno.names] == "bd") # AD
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 3
        idx <- which(x[idx.rows, geno.names] == "ac") # BC
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 2
        idx <- which(x[idx.rows, geno.names] == "ad") # BD
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 4
      } else if(seg == "<efxeg>"){
        idx <- which(x[idx.rows, geno.names] == "ef") # AC
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 1
        idx <- which(x[idx.rows, geno.names] == "fg") # AD
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 3
        idx <- which(x[idx.rows, geno.names] == "ee") # BC
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 2
        idx <- which(x[idx.rows, geno.names] == "eg") # BD
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 4
      } else if(seg == "<hkxhk>"){
        idx <- which(x[idx.rows, geno.names] == "hk") # AC or BD
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 9
        idx <- which(x[idx.rows, geno.names] == "kk") # AD
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 3
        idx <- which(x[idx.rows, geno.names] == "hh") # BC
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 2
        idx <- which(x[idx.rows, geno.names] == "h-") # not AD
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 13
        idx <- which(x[idx.rows, geno.names] == "k-") # not BC
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 12
      } else{
        msg <- paste("unrecognized segregation", seg, "with phase", phase)
        stop(msg)
      }
    } else if(phase == "{11}"){
      if(seg == "<abxcd>"){
        idx <- which(x[idx.rows, geno.names] == "bd") # AC
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 1
        idx <- which(x[idx.rows, geno.names] == "bc") # AD
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 3
        idx <- which(x[idx.rows, geno.names] == "ad") # BC
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 2
        idx <- which(x[idx.rows, geno.names] == "ac") # BD
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 4
      } else if(seg == "<efxeg>"){
        idx <- which(x[idx.rows, geno.names] == "fg") # AC
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 1
        idx <- which(x[idx.rows, geno.names] == "ef") # AD
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 3
        idx <- which(x[idx.rows, geno.names] == "eg") # BC
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 2
        idx <- which(x[idx.rows, geno.names] == "ee") # BD
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 4
      } else if(seg == "<hkxhk>"){
        idx <- which(x[idx.rows, geno.names] == "kk") # AC
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 1
        idx <- which(x[idx.rows, geno.names] == "hk") # AD or BC
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 10
        idx <- which(x[idx.rows, geno.names] == "hh") # BD
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 4
        idx <- which(x[idx.rows, geno.names] == "h-") # not AC
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 11
        idx <- which(x[idx.rows, geno.names] == "k-") # not BD
        if(length(idx) > 0)
          out[idx.rows,][idx] <- 14
      } else{
        msg <- paste("unrecognized segregation", seg, "with phase", phase)
        stop(msg)
      }
    }
  }

  return(out)
}

##' Read genotypes for JoinMap/MapQTL
##'
##' Read genotype data from a file in the \href{https://www.kyazma.nl/index.php/JoinMap/}{JoinMap} format ("loc").
##' @param file name of the file from which the data will be read (can be compressed with gzip)
##' @param na.string a character to be interpreted as NA values
##' @param verbose verbosity level (0/1)
##' @return data frame in the JoinMap format
##' @author Timothee Flutre
##' @seealso \code{\link{writeSegregJoinMap}}
##' @note the code is heavily inspired from the "read.cross.mq.loc" function in the "qtl" package, which I also wrote
##' @export
readSegregJoinMap <- function(file, na.string="--", verbose=1){
  stopifnot(file.exists(file))

  if(verbose > 0){
    msg <- paste0("read genotypes from JoinMap file '",
                  file, "'...")
    write(msg, stdout())
  }

  lines <- readLines(con=file, warn=FALSE)

  ## drop comments
  lines <- vapply(strsplit(lines, ";"), "[", "", 1)

  ## drop empty lines
  lines <- lines[!is.na(lines)]
  blank <- grep("^\\s*$", lines)
  if(length(blank) > 0)
    lines <- lines[-blank]

  grab.param <- function(lines, param, longname, filetype){
    ## stuff for error message
    if(missing(filetype)) filetype <- ""
    else filetype <- paste(" in", filetype, "file")
    if(missing(longname)) longname <- param

    g <- grep(param, lines)
    if(length(g) == 0)
      stop("Cannot find ", longname, " in ", filetype)

    ## remove white space
    line <- gsub("\\s+", "", lines[g])

    result <- strsplit(line, "=")[[1]][2]

    return(list(result, g)) # g is the line number
  }

  ## extract the population name
  res <- grab.param(lines, "^name", "population name", "loc")
  pop.name <- res[[1]]
  todrop <- res[[2]]
  if(verbose > 0){
    msg <- paste0("population name: ", pop.name)
    write(msg, stdout())
  }

  ## extract the population type
  res <- grab.param(lines, "^popt", "population type", "loc")
  pop.type <- res[[1]]
  todrop <- c(todrop, res[[2]])
  pop.types <- c("BC1","F2","RIx","DH","DH1","DH2","HAP","HAP1","CP","BCpxFy",
                 "IMxFy")
  if(! pop.type %in% pop.types){
    msg <- paste("unknown population type", pop.type)
    stop(msg, call.=FALSE)
  }
  if(verbose > 0){
    msg <- paste0("population type: ", pop.type)
    write(msg, stdout())
  }

  ## extract the number of loci
  res <- grab.param(lines, "^nloc", "Number of loci", "loc")
  nb.loci <- as.numeric(res[[1]])
  todrop <- c(todrop, res[[2]])
  if(verbose > 0){
    msg <- paste0("nb of loci: ", nb.loci)
    write(msg, stdout())
  }
  seg <- rep(NA, nb.loci)
  phase <- rep(NA, nb.loci)
  classif <- rep(NA, nb.loci)

  ## extract the number of individuals
  res <- grab.param(lines, "^nind", "Number of individuals", "loc")
  nb.inds <- as.numeric(res[[1]])
  todrop <- c(todrop, res[[2]])
  if(verbose > 0){
    msg <- paste0("nb of individuals: ", nb.inds)
    write(msg, stdout())
  }

  genotypes <- matrix(NA, nrow=nb.loci, ncol=nb.inds)

  ## extract the individual names (if any)
  idx <- grep("individual names:", lines)
  if(length(idx) > 0){
    range.idx <- (idx+1):(idx+nb.inds)
    if(any(range.idx > length(lines))){
      msg <- "can't retrieve the individual names"
      stop(msg)
    }
    ind.names <- lines[range.idx]
    stopifnot(! anyDuplicated(ind.names))
    colnames(genotypes) <- ind.names
    todrop <- c(todrop, c(idx, range.idx))
  } else if(length(idx) > 1){
    msg <- "individual names defined multiple times"
    stop(msg)
  }

  lines <- lines[-todrop]
  if(length(lines) != nb.loci){
    msg <- paste0("the number of loci indicated in the header (", nb.loci,
                  " seems to be wrong")
    stop(msg)
  }
  spl <- strsplit(lines, "\\s+")
  spl.lengths <- vapply(spl, length, 1)
  if(any(spl.lengths > nb.inds+4))
    stop("lines should have no more than ", nb.inds+4, " columns\n",
         "Problems in lines", seq(along=spl.lengths)[spl.lengths > nb.inds+4])

  rownames(genotypes) <- vapply(spl, "[", "", 1)

  for(line.id in 1:length(lines)){
    tokens <- spl[[line.id]]

    if(length(tokens) > nb.inds + 1){
      for(i in 2:(length(tokens)-nb.inds)){
        if(length(grep(pattern="\\{", x=tokens[i])) > 0){
          phase[line.id] <- tokens[i]
        } else if(length(grep(pattern="\\(", x=tokens[i])) > 0){
          classif[line.id] <- tokens[i]
        } else if(length(grep(pattern="<", x=tokens[i])) > 0){
          seg[line.id] <- tokens[i]
        }
      }
    }

    genotypes[line.id,] <- tokens[(length(tokens)-nb.inds+1):length(tokens)]
  }

  ## convert all missing data to NA
  is.miss <- genotypes == na.string
  if(any(is.miss, na.rm=TRUE))
    genotypes[which(is.miss)] <- NA

  out <- data.frame(dummy=rep(NA, nrow(genotypes)))
  if(! all(is.na(seg)))
    out <- cbind(out, seg=seg)
  if(! all(is.na(phase)))
    out <- cbind(out, phase=phase)
  if(! all(is.na(classif)))
    out <- cbind(out, classif=classif)
  out <- cbind(out, as.data.frame(genotypes))
  out <- out[,-1]
  out <- convertFactorColumnsToCharacter(out)

  return(out)
}

##' JoinMap segregation types
##'
##' Get segregation types for each marker from a data frame of marker genotypes in the JoinMap format.
##' @param genos data frame containing genotypes encoded in the JoinMap format,  similar to the output from \code{\link{genoClasses2JoinMap}}
##' @param na.string character indicating a missing genotype
##' @return vector that can be used in \code{\link{writeSegregJoinMap}}
##' @author Charlotte Brault, Timothee Flutre
##' @export
getJoinMapSegregs <- function(genos, na.string="--"){
  stopifnot(is.data.frame(genos))

  genos <- convertFactorColumnsToCharacter(genos)

  segregs <- apply(genos, 1, function(genos.i){
    seg <- NA
    if(all(genos.i %in% c("ac", "ad", "bc", "bd", na.string))){
      seg <- "<abxcd>"
    } else if(all(genos.i %in% c("ef", "eg", "fg", "ee", na.string))){
      seg <- "<efxeg>"
    } else if(all(genos.i %in% c("hh", "kk", "hk", na.string))){
      seg <- "<hkxhk>"
    } else if(all(genos.i %in% c("lm", "ll", na.string))){
      seg <- "<lmxll>"
    } else if(all(genos.i %in% c("nn", "np", na.string))){
      seg <- "<nnxnp>"
    } else
      seg <- NA
    seg
  })
  if(! is.null(rownames(genos)))
    names(segregs) <- rownames(genos)

  miss.seg <- is.na(segregs)
  if(any(miss.seg)){
    msg <- paste0(sum(miss.seg), "marker",
                  ifelse(sum(miss.seg) > 0, "s have", " has"),
                  " missing segregation")
    warning(msg)
  }

  return(segregs)
}

##' Write genotypes for JoinMap/MapQTL
##'
##' Write genotype data in the \href{https://www.kyazma.nl/index.php/JoinMap/}{JoinMap} format into a "loc" file used by this software.
##' @param pop.name name of the population
##' @param pop.type type of the population
##' @param locus vector of locus names; should be shorter or equal to 20 characters, otherwise a warning will be issued; should be in the same order as the rows of the "genos" argument
##' @param segregs vector of segregation types; should be in the same order as the rows of the "genos" argument
##' @param genos data frame containing genotypes encoded in the JoinMap format,  similar to the output from \code{\link{genoClasses2JoinMap}}
##' @param phases vector of phase types (optional; should be in the same order as the rows of the "genos" argument)
##' @param classifs vector of classification types (optional; should be in the same order as the rows of the "genos" argument)
##' @param file name of the file in which the data will be saved (will be compressed if ends with ".gz" and "gzip" is available in the PATH)
##' @param save.ind.names if TRUE, individual names will be saved at the end of the file
##' @param na.string character to replace NA's in "genos"
##' @param verbose verbosity level (0/1)
##' @return nothing
##' @author Timothee Flutre
##' @seealso \code{\link{genoClasses2JoinMap}}, \code{\link{readSegregJoinMap}}, \code{\link{writeGenMapJoinMap}}, \code{\link{writePhenoJoinMap}}
##' @examples
##' \dontrun{## make fake data
##' nb.snps <- 6
##' x <- data.frame(par1=c("AA", "GC", "CG", "AT", NA, "AA"),
##'                 par2=c("AT", "GC", "GG", "AT", "AT", "AT"),
##'                 off1=c("AA", "GG", "CG", "AA", "AA", "AT"),
##'                 off2=c("AT", "GG", "CG", "AT", "AT", "AA"),
##'                 off3=c("AT", "GG", "GG", "TT", "TT", NA),
##'                 off4=c(NA, NA, NA, NA, NA, NA),
##'                 row.names=paste0("snp", 1:nb.snps),
##'                 stringsAsFactors=FALSE)
##' jm <- genoClasses2JoinMap(x=x, reformat.input=TRUE, thresh.na=2, verbose=1)
##' idx <- which(! is.na(jm$seg))
##' writeSegregJoinMap(pop.name="test", pop.type="CP", locus=rownames(jm[idx,]),
##'                    segregs=jm[idx,"seg"], genos=jm[idx,-c(1:9)], file="test.txt")
##' }
##' @export
writeSegregJoinMap <- function(pop.name, pop.type="CP",
                               locus, segregs, genos,
                               phases=NULL, classifs=NULL,
                               file, save.ind.names=TRUE, na.string="--",
                               verbose=1){
  requireNamespace("tools", quietly=TRUE)
  stopifnot(is.character(pop.name),
            is.character(pop.type),
            pop.type %in% c("BC1", "F2", "RIx", "DH", "DH1", "DH2", "HAP",
                            "HAP1", "CP", "BCpxFy", "IMxFy"),
            is.character(locus),
            is.character(segregs),
            is.data.frame(genos),
            all(! c("locus", "seg", "phase", "clas") %in% colnames(genos)),
            nrow(genos) == length(locus),
            nrow(genos) == length(segregs),
            is.character(file),
            is.logical(save.ind.names),
            is.character(na.string))
  if(! is.null(phases))
    stopifnot(is.character(phases),
              nrow(genos) == length(phases))
  if(! is.null(classifs))
    stopifnot(is.character(classifs),
              nrow(genos) == length(classifs))
  locus.ids.too.long <- nchar(locus) > 20
  if(any(locus.ids.too.long)){
    msg <- paste0(sum(locus.ids.too.long), " locus identifier",
                  ifelse(sum(locus.ids.too.long) > 1, "s", ""),
                  " longer than 20 characters")
    warning(msg)
  }
  if(verbose > 0){
    msg <- paste0("write locus genotypes in the JoinMap/MapQTL format...",
                  "\nnb of locus: ", nrow(genos),
                  "\nnb of individuals: ", ncol(genos))
    write(msg, stdout())
  }

  file.ext <- tools::file_ext(file)
  if(file.ext == "gz")
    file <- sub("\\.gz$", "", file)
  con <- file(description=file, open="w")

  txt <- c("; written by the `writeSegregJoinMap` function",
           paste0("; from the `rutilstimflutre` package version ",
                  utils::packageVersion("rutilstimflutre")),
           "",
           paste0("name = ", pop.name),
           paste0("popt = ", pop.type),
           paste0("nloc = ", nrow(genos)),
           paste0("nind = ", ncol(genos)),
           "")
  writeLines(text=txt, con=con)

  tmp <- cbind(locus, segregs)
  if(! is.null(phases))
    tmp <- cbind(tmp, phases)
  if(! is.null(classifs))
    tmp <- cbind(tmp, classifs)
  rownames(tmp) <- NULL
  genos[is.na(genos)] <- na.string
  tmp <- cbind(tmp, genos)
  suppressWarnings(utils::write.table(x=tmp, file=con, append=TRUE,
                                      quote=FALSE, sep="\t",
                                      row.names=FALSE, col.names=FALSE))

  if(save.ind.names){
    stopifnot(all(nchar(colnames(genos)) <= 20))
    txt <- c("",
             "individual names: ",
             colnames(genos))
    writeLines(text=txt, con=con)
  }

  close(con)

  if(file.ext == "gz"){
    if(file.exists(Sys.which("gzip"))){ # should work on Linux and Mac OS
      if(file.exists(paste0(file, ".gz")))
        file.remove(paste0(file, ".gz"))
        system(command=paste("gzip", file))
    }
  }
}

##' Write genetic map for JoinMap/MapQTL
##'
##' Writes the given genetic map in the JoinMap/MapQTL format.
##' @param genmap data frame with at least three columns named "linkage.group", "locus" and "genetic.distance", and one row per locus
##' @param file name of the file in which the data will be saved
##' @param verbose verbosity level (0/1)
##' @return nothing
##' @author Timothee Flutre
##' @seealso \code{\link{writePhenoJoinMap}}, \code{\link{writeSegregJoinMap}}
##' @export
writeGenMapJoinMap <- function(genmap, file, verbose=1){
  stopifnot(is.data.frame(genmap),
            all(c("linkage.group", "locus", "genetic.distance") %in% colnames(genmap)),
            all(nchar(genmap$linkage.group) <= 20),
            is.numeric(genmap$genetic.distance),
            is.character(file))

  genmap <- convertFactorColumnsToCharacter(genmap)
  linkgroups <- unique(genmap$linkage.group)
  nb.linkgroups <- length(linkgroups)
  if(verbose > 0){
    msg <- paste0("write genetic map in the JoinMap/MapQTL format...",
                  "\nnb of linkage groups: ", nb.linkgroups,
                  "\nnb of locus: ", nrow(genmap))
    write(msg, stdout())
  }

  con <- file(description=file, open="w")

  txt <- c("; written by the `writeGenMapJoinMap` function",
           paste0("; from the `rutilstimflutre` package version ",
                  utils::packageVersion("rutilstimflutre")))
  writeLines(text=txt, con=con)

  for(lg in linkgroups){
    txt <- c("",
             paste("group", lg, sep="\t"))
    idx <- which(genmap$linkage.group == lg)
    tmp <- genmap[idx,]
    tmp <- tmp[order(tmp$genetic.distance),]
    txt <- c(txt,
             paste(tmp$locus, tmp$genetic.distance, sep="\t"))
    writeLines(text=txt, con=con)
  }

  close(con)
}

##' Write phenotypes for JoinMap/MapQTL
##'
##' Writes the given phenotypes in the JoinMap/MapQTL format.
##' @param phenos data frame with one genotype per row and traits in columns (one colun can contain the genotype identifiers)
##' @param file name of the file in which the data will be saved
##' @param alias.miss alias used to replace missing values (which should be encoded as NA in "genmap")
##' @param verbose verbosity level (0/1)
##' @return nothing
##' @author Timothee Flutre
##' @seealso \code{\link{writeGenMapJoinMap}}, \code{\link{writeSegregJoinMap}}
##' @export
writePhenoJoinMap <- function(phenos, file, alias.miss=".", verbose=1){
  stopifnot(is.data.frame(phenos),
            is.character(file))

  phenos <- convertFactorColumnsToCharacter(phenos)
  if(verbose > 0){
    msg <- paste0("write phenotypes in the JoinMap/MapQTL format...",
                  "\nnb of traits: ", ncol(phenos),
                  "\nnb of individuals: ", nrow(phenos))
    write(msg, stdout())
  }

  con <- file(description=file, open="w")

  txt <- c("; written by the `writePhenoJoinMap` function",
           paste0("; from the `rutilstimflutre` package version ",
                  utils::packageVersion("rutilstimflutre")),
           "",
           paste0("ntrt = ", ncol(phenos)),
           paste0("nind = ", nrow(phenos)),
           paste0("miss = ", alias.miss),
           "")
  writeLines(text=txt, con=con)

  tmp <- phenos
  if(any(is.na(tmp)))
    tmp[is.na(tmp)] <- alias.miss
  suppressWarnings(utils::write.table(x=tmp, file=con, append=TRUE,
                                      quote=FALSE, sep="\t",
                                      row.names=FALSE, col.names=TRUE))

  close(con)
}

##' Info about a given genetic map
##'
##' Returns some information on a given genetic map, notably the relations between linkage groups and chromosomes.
##' Beforehand, the genetic map can be filtered.
##' @param genmap data frame with at least three columns named "linkage.group", "locus" and "chr", and one row per locus
##' @param min.mrks.per.lg to filter the genetic map, any linkage group with strictly less than this threshold will be discarded
##' @param min.mrks.per.chr.in.lg to filter the genetic map, for linkage groups with markers belonging to several chromosomes, any marker belonging to chromosomes with stricly less than this threshold will be discarded
##' @param chrs.todrop optional vector of chromosome names to drop from the analysis
##' @param verbose verbosity level (0/1/2)
##' @return list with several summaries, as well as a clean version of the genetic map with one linkage group per chromosome (among the most represented)
##' @author Timothee Flutre
##' @export
infoGeneticMap <- function(genmap, min.mrks.per.lg=0,
                           min.mrks.per.chr.in.lg=0, chrs.todrop=NULL,
                           verbose=1){
  stopifnot(is.data.frame(genmap),
            all(c("linkage.group", "locus", "chr") %in%
                colnames(genmap)),
            ! any(is.na(genmap$linkage.group)),
            ! anyDuplicated(genmap$locus),
            is.numeric(genmap$genetic.distance))
  if(! is.null(chrs.todrop))
    stopifnot(is.character(chrs.todrop))

  if(! is.null(chrs.todrop)){
    if(verbose > 0){
      msg <- paste0("discard ", length(chrs.todrop), " chromosome",
                    ifelse(length(chrs.todrop) == 1, "", "s"), "...")
      write(msg, stdout())
    }
    genmap <- genmap[! genmap$chr %in% chrs.todrop,]
  }

  if(min.mrks.per.lg > 0){
    if(verbose > 0){
      msg <- paste0("discard linkage groups with strictly less than ",
                    min.mrks.per.lg, " markers...")
      write(msg, stdout())
    }
    tmp <- tapply(1:nrow(genmap), factor(genmap$linkage.group), length)
    has.enough.mrks <- tmp >= min.mrks.per.lg
    if(! all(has.enough.mrks)){
      genmap <- genmap[genmap$linkage.group %in% names(tmp)[has.enough.mrks],]
    }
  }

  if(min.mrks.per.chr.in.lg > 0){
    if(verbose > 0){
      msg <- paste0("discard chromosomes in linkage groups if strictly less than ",
                    min.mrks.per.chr.in.lg, " markers...\n")
      write(msg, stdout())
    }
    tmp <- tapply(genmap$chr, factor(genmap$linkage.group), function(x){
      tmp2 <- table(x)
      names(tmp2[tmp2 >= min.mrks.per.chr.in.lg])
    })
    tmp[sapply(tmp, length) == 0] <- NULL
    if(length(tmp) == 0)
      stop("all data were filtered out!")
    genmap <- do.call(rbind, lapply(names(tmp), function(lg){
      genmap[genmap$linkage.group == lg & genmap$chr %in% tmp[[lg]],]
    }))
  }

  if(verbose > 0){
    msg <- paste0("total nb of markers: ", nrow(genmap))
    write(msg, stdout())
  }

  ## get chr and lg names, and sort them in natural order if possible
  chr.names <- unique(genmap$chr)
  if(requireNamespace("gtools", quietly=TRUE))
    chr.names <- gtools::mixedsort(chr.names)
  nb.chrs <- length(chr.names)
  lg.names <- unique(genmap$linkage.group)
  if(requireNamespace("gtools", quietly=TRUE))
    lg.names <- gtools::mixedsort(lg.names)
  nb.lgs <- length(lg.names)

  tab.chr <- table(genmap$chr, useNA="always")
  tab.chr <- tab.chr[chr.names]
  if(verbose > 0){
    msg <- paste0("\nnb of chromosomes: ", nb.chrs,
                  "\nnb of markers per chromosome:")
    write(msg, stdout())
    print(tab.chr)
  }

  tab.lg <- table(genmap$linkage.group, useNA="always")
  tab.lg <- tab.lg[order(tab.lg, decreasing=TRUE)]
  if(verbose > 0){
    msg <- paste0("\nnb of linkage groups: ", nb.lgs,
                  "\nnb of markers per linkage group:")
    write(msg, stdout())
    print(tab.lg)
  }

  ## stats per lg
  lg2chr <- tapply(1:nrow(genmap), factor(genmap$linkage.group),
                   function(idx){
                     stats::setNames(genmap$chr[idx], genmap$locus[idx])
                   })
  tab.lg2chr <- lapply(lg2chr, function(x){
    table(factor(x, levels=chr.names))
  })
  tab.lg2chr <- tab.lg2chr[names(tab.lg)[order(tab.lg, decreasing=TRUE)]
                           [-length(tab.lg)]]
  mat.lg2chr <- t(do.call(rbind, tab.lg2chr))
  tab.lg2chr <- lapply(tab.lg2chr, function(x){
    x <- x[order(x, decreasing=TRUE)]
    x[x > 0]
  })
  chrs.per.lg <- sapply(tab.lg2chr, length)

  ## stats per chr
  chr2lg <- list()
  for(chr in chr.names){
    idx <- which(genmap$chr == chr)
    chr2lg[[chr]] <- stats::setNames(genmap$linkage.group[idx],
                                     genmap$locus[idx])
  }
  tab.chr2lg <- lapply(chr2lg, function(x){
    table(factor(x, levels=lg.names))
  })
  tab.chr2lg <- tab.chr2lg[chr.names]
  mat.chr2lg <- t(do.call(rbind, tab.chr2lg))
  tab.chr2lg <- lapply(tab.chr2lg, function(x){
    x <- x[order(x, decreasing=TRUE)]
    x[x > 0]
  })
  lgs.per.chr <- sapply(tab.chr2lg, length)

  if(verbose > 0){
    msg <- "\nnb of chromosomes per linkage group:"
    write(msg, stdout())
    print(chrs.per.lg)
    msg <- "\nnb of linkage groups per chromosome:"
    write(msg, stdout())
    print(lgs.per.chr)
    if(verbose > 1){
      msg <- "\nnb of markers per chromosome in each linkage group:"
      write(msg, stdout())
      print(tab.lg2chr)
      msg <- "\nnb of markers per linkage group for each chromosome:"
      write(msg, stdout())
      print(tab.chr2lg)
    }
  }

  genmap.clean <- genmap[genmap$linkage.group %in%
                         sapply(tab.chr2lg, function(x){
                           names(x)[1]
                         }),]
  rownames(genmap.clean) <- NULL

  return(list(clean=genmap.clean,
              tab.chr=tab.chr, tab.lg=tab.lg,
              lg2chr=lg2chr, tab.lg2chr=tab.lg2chr, mat.lg2chr=mat.lg2chr,
              chr2lg=chr2lg, tab.chr2lg=tab.chr2lg, mat.chr2lg=mat.chr2lg))
}

##' Stacked barplot of markers
##'
##' Make a barplot of markers per linkage group (or chromosome), stacked per chromosome (or linkage group).
##' @param counts matrix with the number of markers, with chromosomes in rows and linkage groups in columns, or vice-versa; should have row and column names
##' @param las see \code{\link[graphics]{par}}
##' @param border see \code{\link[graphics]{barplot}}
##' @param col see \code{\link[graphics]{barplot}}
##' @param leg.txt see \code{\link[graphics]{barplot}}; if NULL, the stack labels will appear in the middle of each bar
##' @param args.leg see \code{\link[graphics]{barplot}}
##' @param ... arguments to be passed to \code{\link[graphics]{barplot}}
##' @return see \code{\link[graphics]{barplot}}
##' @author Timothee Flutre
##' @export
barplotGeneticMap <- function(counts,
                              las=1, border=NA,
                              col=grDevices::rainbow(nrow(counts)),
                              leg.txt=rownames(counts),
                              args.leg=list(x="topright", bty="n", border=NA,
                                            fill=rev(grDevices::rainbow(nrow(counts)))),
                              ...){
  stopifnot(is.matrix(counts),
            ! is.null(rownames(counts)),
            ! is.null(colnames(counts)))

  bp.x <- graphics::barplot(counts, las=las, border=border,
                            col=col,
                            legend.text=leg.txt,
                            args.legend=args.leg,
                            xaxt="n",
                            ...)
  graphics::text(x=bp.x, y=graphics::par("usr")[3] - 1, srt=45, adj=1,
                 labels=colnames(counts), xpd=TRUE)

  ## add stack labels
  if(is.null(leg.txt)){
    nb.nonzeros <- apply(counts, 2, function(x){
      sum(x != 0)
    })
    labs <- do.call(c, lapply(1:ncol(counts), function(j){
      x <- counts[,j]
      rownames(counts)[x != 0]
    }))
    txt.x <- rep(bp.x, nb.nonzeros)
    txt.y <- do.call(c, lapply(1:ncol(counts), function(j){
      x <- counts[,j]
      is.nonzero <- x != 0
      tmp <- rep(NA, sum(is.nonzero))
      tmp[1] <- x[is.nonzero][1] / 2
      i <- 2
      while(i <= sum(is.nonzero)){
        tmp[i] <- sum(x[is.nonzero][1:(i-1)]) + x[is.nonzero][i] / 2
        i <- i + 1
      }
      tmp
    }))
    graphics::text(x=txt.x, y=txt.y, labels=labs, srt=45)
  }

  invisible(bp.x)
}

##' Genotype coding
##'
##' Write marker genotypes formatted for CarthaGene in a file.
##' @param x matrix of characters containing marker genotypes in the CarthaGene format, for instance from \code{\link{joinMap2backcross}}; row names should correspond to marker names
##' @param file path to the file to which the data will be written
##' @param type type of the data ("f2 backcross")
##' @param aliases aliases
##' @return nothing
##' @author Timothee Flutre
##' @seealso \code{\link{joinMap2backcross}}
##' @export
writeCartagene <- function(x, file, type="f2 backcross", aliases=NULL){
  stopifnot(is.matrix(x),
            ! is.null(rownames(x)),
            type %in% c("f2 backcross"))
  if(! is.null(aliases))
    stopifnot(is.vector(aliases),
              is.character(aliases))
  if(file.exists(file))
    file.remove(file)

  con <- file(file, open="a")

  ## first line: type of the data
  txt <- paste0("data type ", type, "\n")
  cat(txt, file=con, append=TRUE)

  ## second line: size of the data
  txt <- paste0(ncol(x), " ", nrow(x), " 0 0")
  if(! is.null(aliases))
    for(a in aliases)
      txt <- paste0(" ", a)
  txt <- paste0(txt, "\n")
  cat(txt, file=con, append=TRUE)

  ## last section: marker data
  txt <- apply(cbind(paste0("*", rownames(x)),
                     apply(x, 1, paste, collapse="")),
               1, paste, collapse="\t")
  rownames(txt) <- NULL
  cat(txt, file=con, append=TRUE, sep="\n")

  close(con)
}

##' Genotype coding
##'
##' Convert unphased genotypes encoded in the \href{https://www.kyazma.nl/index.php/JoinMap/}{JoinMap} format into the backcross configuration according to the \href{http://www7.inra.fr/mia/T/CarthaGene/}{CarthaGene} or R/qtl formats.
##' Genotypes are assumed to come from a bi-parental cross with outbred parents (population type "CP" in JoinMap) and the goal is to build two parental genetic maps via a pseudo-test-cross strategy (Grattapaglia and Sederoff, 1994).
##' This function hence encodes genotypes as two backcrosses, one per parent.
##' All locus will be duplicated and the suffix "_m" (for "mirror") added to their identifiers.
##' @param x data frame in the JoinMap format, for instance from \code{\link{genoClasses2JoinMap}}; any column named "phase" or "clas" will be discarded from the start; a column named "seg" should be present (markers with missing segregation will be discarded); a column named "locus" should be present unless row names are given; all other columns are supposed to correspond to genotypes
##' @param alias.hom alias to specify homozygotes ("A" for CarthaGene, 1 for R/qtl); should be different from alias.dup
##' @param alias.het alias to specify heterozygotes ("H" for CarthaGene, 2 for R/qtl); should be different from alias.dup
##' @param alias.dup alias used when duplicating the locus; should be different from both alias.hom and alias.het, but should be of the same mode (i.e. character or numeric)
##' @param alias.miss alias to specify missing values
##' @param parent.names vector containing the names of the parents
##' @param verbose verbosity level (0/1/2)
##' @return list of two matrices, one per parent
##' @author Timothee Flutre [aut], Agnes Doligez [ctb]
##' @seealso \code{\link{genoClasses2JoinMap}}, \code{\link{openCarthagene}}, \code{\link{setJoinMapPhasesFromParentalLinkGroups}}
##' @export
joinMap2backcross <- function(x, alias.hom="A", alias.het="H", alias.dup="Z",
                              alias.miss="-", parent.names=c("parent1", "parent2"),
                              verbose=1){
  seg.types <- c("<abxcd>", "<efxeg>", "<hkxhk>", "<lmxll>", "<nnxnp>")
  stopifnot(is.data.frame(x),
            ncol(x) >= 5,
            ! is.null(colnames(x)),
            "seg" %in% colnames(x),
            all(x$seg %in% c(seg.types, NA)),
            is.vector(alias.hom),
            mode(alias.hom) %in% c("character", "numeric"),
            length(alias.hom) == 1,
            is.vector(alias.het),
            mode(alias.het) %in% c("character", "numeric"),
            length(alias.het) == 1,
            is.vector(alias.dup),
            mode(alias.dup) %in% c("character", "numeric"),
            length(alias.dup) == 1,
            alias.hom != alias.dup,
            alias.het != alias.dup,
            all(c(mode(alias.hom), mode(alias.dup)) == mode(alias.dup)),
            is.vector(alias.miss),
            length(alias.miss) == 1,
            is.character(parent.names),
            length(parent.names) == 2)
  if(! "locus" %in% colnames(x))
    stopifnot(! is.null(rownames(x)))

  if(verbose > 0){
    msg <- "convert unphased genotypes from JoinMap into the backcross\n(i.e. \"parental\") configuration..."
    write(msg, stdout())
  }

  ## reformat the input
  x <- convertFactorColumnsToCharacter(x)
  idx <- which(colnames(x) %in% c("phase", "clas"))
  if(length(idx) > 0)
    x <- x[, -idx]
  idx.loc <- which(colnames(x) == "locus")
  idx.seg <- which(colnames(x) == "seg")
  if(length(idx.loc) == 0){
    x <- cbind(x, locus=rownames(x))
    idx.loc <- ncol(x)
  }
  x <- x[, c(idx.loc, idx.seg, setdiff(1:ncol(x), c(idx.loc, idx.seg)))]
  rownames(x) <- NULL

  ## clean the input
  miss.segreg <- is.na(x$seg)
  if(any(miss.segreg)){
    if(verbose > 0){
      msg <- paste0("discard ", sum(miss.segreg),
                    " locus with missing segregation...")
      write(msg, stdout())
    }
    x <- x[! miss.segreg,]
  }
  if(verbose > 0){
    msg <- paste0("convert ", nrow(x), " locus...")
    write(msg, stdout())
  }

  ## make one matrix per parent
  idx.seg <- lapply(seg.types, function(seg){
    which(x$seg == seg)
  })
  names(idx.seg) <- seg.types
  idx.loc.p1 <- c(idx.seg[["<abxcd>"]], idx.seg[["<efxeg>"]],
                  idx.seg[["<hkxhk>"]], idx.seg[["<lmxll>"]])
  idx.loc.p2 <- c(idx.seg[["<abxcd>"]], idx.seg[["<efxeg>"]],
                  idx.seg[["<hkxhk>"]], idx.seg[["<nnxnp>"]])
  p1 <- as.matrix(x[idx.loc.p1, -1])
  rownames(p1) <- x[idx.loc.p1, "locus"]
  p2 <- as.matrix(x[idx.loc.p2, -1])
  rownames(p2) <- x[idx.loc.p2, "locus"]

  ## convert genotypes' format
  seg <- "<abxcd>"
  idx.p1 <- which(p1[,"seg"] == seg)
  if(length(idx.p1) > 0)
    p1[idx.p1, -1] <-
      t(matrix(apply(p1[idx.p1, -1, drop=FALSE], 2, function(genos){
        genos <- gsub("ac|ad", alias.hom, genos)
        gsub("bc|bd", alias.het, genos)
      }))) # assume parent 2 is homozygous "aa"
  idx.p2 <- which(p2[,"seg"] == seg)
  if(length(idx.p2) > 0)
    p2[idx.p2, -1] <-
      t(matrix(apply(p2[idx.p2, -1, drop=FALSE], 2, function(genos){
        genos <- gsub("ac|bc", alias.hom, genos)
        gsub("ad|bd", alias.het, genos)
      }))) # assume parent 1 is homozygous "cc"

  seg <- "<efxeg>"
  idx.p1 <- which(p1[,"seg"] == seg)
  if(length(idx.p1) > 0)
    p1[idx.p1, -1] <-
      t(matrix(apply(p1[idx.p1, -1, drop=FALSE], 2, function(genos){
        genos <- gsub("ee|eg", alias.hom, genos)
        gsub("ef|fg", alias.het, genos)
      }))) # assume parent 2 is homozygous "ee"
  idx.p2 <- which(p2[,"seg"] == seg)
  if(length(idx.p2) > 0)
    p2[idx.p2, -1] <-
      t(matrix(apply(p2[idx.p2, -1, drop=FALSE], 2, function(genos){
        genos <- gsub("ee|ef", alias.hom, genos)
        gsub("eg|fg", alias.het, genos)
      }))) # assume parent 1 is homozygous "ee"

  seg <- "<hkxhk>"
  idx.p1 <- which(p1[,"seg"] == seg)
  if(length(idx.p1) > 0)
    p1[idx.p1, -1] <-
      t(matrix(apply(p1[idx.p1, -1, drop=FALSE], 2, function(genos){
        genos <- gsub("hh", alias.hom, genos)
        genos <- gsub("hk", alias.miss, genos)
        gsub("kk", alias.het, genos)
      }))) # "hk" as missing because phase is unknown
  idx.p2 <- which(p2[,"seg"] == seg)
  if(length(idx.p2) > 0)
    p2[idx.p2, -1] <-
      t(matrix(apply(p2[idx.p2, -1, drop=FALSE], 2, function(genos){
        genos <- gsub("hh", alias.hom, genos)
        genos <- gsub("hk", alias.miss, genos)
        gsub("kk", alias.het, genos)
      }))) # "hk" as missing because phase is unknown

  seg <- "<lmxll>"
  idx.p1 <- which(p1[,"seg"] == seg)
  if(length(idx.p1) > 0)
    p1[idx.p1, -1] <-
      t(matrix(apply(p1[idx.p1, -1, drop=FALSE], 2, function(genos){
        genos <- gsub("ll", alias.hom, genos)
        gsub("lm", alias.het, genos)
      })))

  seg <- "<nnxnp>"
  idx.p2 <- which(p2[,"seg"] == seg)
  if(length(idx.p2) > 0)
    p2[idx.p2, -1] <-
      t(matrix(apply(p2[idx.p2, -1, drop=FALSE], 2, function(genos){
        genos <- gsub("nn", alias.hom, genos)
        gsub("np", alias.het, genos)
      })))

  if(verbose > 0){
    msg <- paste0("nb of locus: ", parent.names[1], "=", nrow(p1),
                  " ", parent.names[2], "=", nrow(p2))
    write(msg, stdout())
  }

  if(verbose > 0){
    msg <- "duplicate all locus and invert the coding..."
    write(msg, stdout())
  }
  p1 <- p1[, - which(colnames(p1) == "seg")]
  p1.m <- apply(p1, 2, function(genos){
    genos <- gsub(alias.hom, alias.dup, genos)
    genos <- gsub(alias.het, alias.hom, genos)
    gsub(alias.dup, alias.het, genos)
  })
  rownames(p1.m) <- paste0(rownames(p1), "_m")

  p2 <- p2[, - which(colnames(p2) == "seg")]
  p2.m <- apply(p2, 2, function(genos){
    genos <- gsub(alias.hom, alias.dup, genos)
    genos <- gsub(alias.het, alias.hom, genos)
    gsub(alias.dup, alias.het, genos)
  })
  rownames(p2.m) <- paste0(rownames(p2), "_m")

  ## make the output
  out <- list()
  out[[1]] <- rbind(p1, p1.m)
  if(mode(out[[1]]) != mode(alias.hom))
    mode(out[[1]]) <- mode(alias.hom)
  if(! is.na(alias.miss)){
    idx.na <- which(is.na(out[[1]]))
    if(length(idx.na) > 0)
      out[[1]][idx.na] <- alias.miss
  }
  out[[2]] <- rbind(p2, p2.m)
  if(mode(out[[2]]) != mode(alias.hom))
    mode(out[[2]]) <- mode(alias.hom)
  if(! is.na(alias.miss)){
    idx.na <- which(is.na(out[[2]]))
    if(length(idx.na) > 0)
      out[[2]][idx.na] <- alias.miss
  }
  names(out) <- parent.names

  return(out)
}

##' Plot physical versus genetic distances
##'
##' Plots physical versus genetic distances
##' @param x data frame with markers in rows and at least three columns named "chr", "physical.distance", "genetic.distance"
##' @param chrs chromosome(s) to plot (all if NULL); will be sorted by mixedsort if the "gtools" package is available; if several chromosomes are to be plotted, there will be one plot per chromosome, with two plots per row as specified to \code{\link[graphics]{par}}
##' @param type type of plot to draw (see \code{\link[graphics]{plot}})
##' @param las style of axis labels (see \code{\link[graphics]{par}})
##' @param xlab title for the x axis
##' @param ylab title for the y axis
##' @param ... arguments to be passed to \code{\link[graphics]{plot}}, except "x" (will be the physical distance), "y" (will be the genetic distance), "main" (will be the chromosome identifier), "type", "las", "xlab" and "ylab"
##' @return nothing
##' @author Timothee Flutre
##' @export
plotPhyVsGenDistances <- function(x, chrs=NULL, type="b", las=1,
                                  xlab="physical distance (in bp)",
                                  ylab="genetic distance (in cM)",
                                  ...){
  stopifnot(is.data.frame(x),
            all(c("chr", "physical.distance", "genetic.distance") %in%
                colnames(x)))
  if(! is.null(chrs))
    stopifnot(all(chrs %in% unique(x$chr)))

  if(is.null(chrs))
    chrs <- unique(x$chr)
  if(requireNamespace("gtools", quietly=TRUE))
    chrs <- gtools::mixedsort(chrs)

  def.par <- graphics::par(mfrow=c(1,1))
  if(length(chrs) > 1)
    def.par <- graphics::par(mfrow=c(length(chrs) / 2,2))

  for(chr in chrs){
    idx <- which(x$chr == chr)
    tmp <- x[idx,]
    tmp <- tmp[order(tmp$physical.distance),]
    graphics::plot(tmp$physical.distance, tmp$genetic.distance,
                   main=chr, type=type, las=las, xlab=xlab, ylab=ylab, ...)
  }

  on.exit(graphics::par(def.par))
}

##' Genotype coding
##'
##' Set linkage phases in the \href{https://www.kyazma.nl/index.php/JoinMap/}{JoinMap} format from two sets of parental linkage groups.
##' This function is tested only for segregation types <hkxhk>, <lmxll> and <nnxnp>, that is, in the case of bi-allelic SNPs segregating in a cross of outbred parents.
##' @param x data frame in the JoinMap format, for instance from \code{\link{genoClasses2JoinMap}}; row names should contain locus names; the first column should be "seg"; any column(s) "phase" or "clas" already existing will be discarded; other columns should contain genotype data
##' @param lg.par1 linkage groups of the first parent (usually the mother) as a data frame with one row per marker and at least two columns named "linkage.group" and "locus"; "mirror" markers should have suffix "_m" as done by \code{\link{joinMap2backcross}}
##' @param lg.par2 linkage groups of the second parent (usually the father) as a data frame with one row per marker and at least two columns named "linkage.group" and "locus"; "mirror" markers should have suffix "_m" as done by \code{\link{joinMap2backcross}}
##' @return data frame in the JoinMap format with a "phase" column
##' @author Timothee Flutre
##' @export
setJoinMapPhasesFromParentalLinkGroups <- function(x, lg.par1, lg.par2){
  stopifnot(is.data.frame(x),
            ! is.null(rownames(x)),
            "seg" %in% colnames(x),
            is.data.frame(lg.par1),
            ncol(lg.par1) >= 2,
            all(c("linkage.group", "locus") %in% colnames(lg.par1)),
            is.data.frame(lg.par2),
            ncol(lg.par2) >= 2,
            all(c("linkage.group", "locus") %in% colnames(lg.par2)))
  if(any(! x$seg %in% c("<hkxhk>", "<lmxll>", "<nnxnp>", NA)))
    warning("some segregation types are not in <hkxhk>, <lmxll> or <nnxnp>")

  ## reformat the inputs
  lg.par1 <- convertFactorColumnsToCharacter(lg.par1)
  lg.par2 <- convertFactorColumnsToCharacter(lg.par2)
  stopifnot(all(unique(lg.par1$linkage.group) %in%
                unique(lg.par2$linkage.group)),
            all(unique(lg.par2$linkage.group) %in%
                unique(lg.par1$linkage.group)))
  x <- convertFactorColumnsToCharacter(x)
  if("phase" %in% colnames(x))
    x <- x[, - which(colnames(x) == "phase")]
  if("clas" %in% colnames(x))
    x <- x[, - which(colnames(x) == "clas")]

  ## add phases for each parent
  addPhase <- function(lg, par.num){
    lg$locus.init <- gsub("_m", "", lg$locus)
    lg[[paste0("phase.par", par.num)]] <- 0
    lg[[paste0("phase.par", par.num)]][grepl("_m", lg$locus)] <- 1
    return(lg)
  }
  lg.par1 <- addPhase(lg.par1, 1)
  lg.par2 <- addPhase(lg.par2, 2)

  ## merge parental phases
  lg <- merge(lg.par1, lg.par2, by="locus.init", all.x=TRUE,
              all.y=TRUE) # column "locus.init" will be sorted
  lg$phase.par1[is.na(lg$phase.par1)] <- "-"
  lg$phase.par2[is.na(lg$phase.par2)] <- "-"
  lg$phase <- paste0("{", lg$phase.par1, lg$phase.par2, "}")

  ## add phases to the JoinMap data frame
  x.phased <- merge(cbind(locus.init=rownames(x), x),
                    lg[, c("locus.init", "phase")],
                    by="locus.init", all.x=TRUE, all.y=TRUE)
  rownames(x.phased) <- x.phased$locus.init
  x.phased <- cbind(seg=x.phased$seg,
                    phase=x.phased$phase,
                    x.phased[, -c(1,2,ncol(x.phased))])
  x.phased <- convertFactorColumnsToCharacter(x.phased)

  return(x.phased)
}

##' Genotype coding
##'
##' Convert genotypes encoded in the \href{https://www.kyazma.nl/index.php/JoinMap/}{JoinMap} format into a design matrix.
##' Genotypes are assumed to come from a bi-parental family with heterozygous parents (population type CP in JoinMap).
##' * When phase information is not used, see equation 3 of \href{http://dx.doi.org/10.1186/1471-2156-12-82}{Wang (2011)} for the parameterization, as well as \href{https://doi.org/10.1177/1471082X16644998}{Chiquet et al (2016)} for the constraints.
##' Additive and dominance effects are handled, but epistasis is ignored.
##' * When phase information is used, each marker is coded with four columns, one per parental haplotype.
##' Only additive effects are handled.
##' @param jm data frame in the JoinMap format, for instance from \code{\link{genoClasses2JoinMap}}; the first three columns should be "locus", "seg" and "phase"; any "clas" column will be discarded; if \code{use.phase=FALSE}, "phase" will also be discarded; no missing data is allowed
##' @param use.phase if TRUE, phase information will be used
##' @param parameterization parameterization (allele/F_infinity/allele-count); ignored if \code{use.phase=TRUE}
##' @param constraints constraints (NULL/mu-zero/last-zero); ignored if \code{use.phase=TRUE}
##' @param rm.col.zeros if TRUE, the columns full of zeros will be removed
##' @param rm.dom if TRUE, the columns corresponding to dominance effects will be removed; ignored if \code{use.phase=TRUE}
##' @param verbose verbosity level (0/1)
##' @return matrix, with genotypes in rows and markers in columns
##' @author Timothee Flutre
##' @seealso \code{\link{genoClasses2JoinMap}}, \code{\link{writeSegregJoinMap}}, \code{\link{updateJoinMap}}
##' @examples
##' \dontrun{nb.snps <- 6
##' x <- data.frame(par1=c("AA", "GC", "CG", "AT", NA, "AA"),
##'                 par2=c("AT", "GC", "GG", "AT", "AT", "AT"),
##'                 off1=c("AA", "GG", "CG", "AA", "AA", "AT"),
##'                 off2=c("AT", "GG", "CG", "AT", "AT", "AA"),
##'                 off3=c("AT", "GG", "GG", "TT", "TT", NA),
##'                 off4=c(NA, NA, NA, NA, NA, NA),
##'                 row.names=paste0("snp", 1:nb.snps),
##'                 stringsAsFactors=FALSE)
##' jm <- genoClasses2JoinMap(x=x, reformat.input=TRUE, thresh.na=2, verbose=1)
##' jm <- jm[! is.na(jm$seg),
##'          c("seg", grep("^off", colnames(jm), value=TRUE))]
##' jm <- cbind(locus=rownames(jm), jm)
##' any.na <- apply(jm, 2, function(x) any(is.na(x)))
##' jm <- jm[, - which(any.na)]
##' joinMap2designMatrix(jm)
##' }
##' @export
joinMap2designMatrix <- function(jm, use.phase=FALSE,
                                 parameterization="allele",
                                 constraints=NULL, rm.col.zeros=TRUE,
                                 rm.dom=FALSE,
                                 verbose=1){
  stopifnot(is.data.frame(jm),
            all(colnames(jm)[1:2] == c("locus", "seg")),
            all(jm$seg %in% c("<abxcd>", "<efxeg>", "<hkxhk>", "<lmxll>",
                              "<nnxnp>")),
            is.logical(use.phase))
  if("clas" %in% colnames(jm))
    jm[["clas"]] <- NULL
  if(use.phase){
    stopifnot(colnames(jm)[3] == "phase",
              all(! is.na(jm[, -c(1:3)])))
  } else{
    if("phase" %in% colnames(jm))
      jm[["phase"]] <- NULL
    stopifnot(all(! is.na(jm[, -c(1:2)])))
  }
  stopifnot(parameterization %in% c("allele", "F_infinity", "allele-count"),
            is.logical(rm.col.zeros),
            is.logical(rm.dom))
  if(! is.null(constraints))
    stopifnot(constraints %in% c("mu-zero", "last-zero"))

  out <- NULL

  nb.locus <- nrow(jm)
  if(use.phase){
    nb.inds <- ncol(jm) - 3
  } else
    nb.inds <- ncol(jm) - 2
  if(verbose > 0){
    msg <- paste0("nb of locus: ", nb.locus,
                  "\nnb of inds: ", nb.inds)
    write(msg, stdout())
  }
  alleles <- lapply(strsplit(jm$seg, ""), function(x){
    sort(unique(x[! x %in% c("<","x",">")]))
  })
  nb.alleles <- sapply(alleles, length)

  if(verbose > 0){
    msg <- "for each locus, make the unconstrained design matrix..."
    write(msg, stdout())
  }
  X <- list()
  for(l in 1:nb.locus){ # slow...

    if(use.phase){

      locus <- jm$locus[l]
      X[[l]] <- matrix(data=0, nrow=nb.inds, ncol=4)
      rownames(X[[l]]) <- colnames(jm)[-(1:3)]
      colnames(X[[l]]) <- paste0(locus,
                                 rep(paste0(rep(paste0(".par", 1:2), each=2),
                                            rep(paste0(".haplo", 1:2), 2))))
      par.phases <- unique(strsplit(gsub("\\{|\\}", "", jm$phase[l]), "")[[1]])
      if(length(par.phases) == 1){ # {00} or {11}
        pat <- paste0(locus, ".par[12].haplo", as.numeric(par.phases) + 1)
      } else if(par.phases[1] == "-"){ # {-0} or {-1}
        pat <- paste0(locus, ".par2.haplo", as.numeric(par.phases[2]) + 1)
      } else if(par.phases[2] == "-"){ # {0-} or {1-}
        pat <- paste0(locus, ".par1.haplo", as.numeric(par.phases[1]) + 1)
      } else if(par.phases[1] == "0"){ # {01}
        pat <- paste(paste0(locus, c(".par1.haplo1", ".par2.haplo2")),
                     collapse="|")
      } else{ # {10}
        pat <- paste(paste0(locus, c(".par1.haplo2", ".par2.haplo1")),
                     collapse="|")
      }
      X[[l]][, grep(pat, colnames(X[[l]]))] <- 1

    } else{

      if(parameterization == "allele"){ # equation 3 of Wang (2011)
        ## set up the design matrix
        locus <- jm$locus[l]
        m <- nb.alleles[l]
        X[[l]] <- matrix(data=0, nrow=nb.inds, ncol=m + m*(m+1)/2)
        rownames(X[[l]]) <- colnames(jm)[-(1:2)]
        tmp <- c()
        for(j in 1:m)
          tmp <- c(tmp, paste0(locus, ".", alleles[[l]][j]))
        for(j in 1:m){
          for(k in j:m)
            tmp <- c(tmp, paste0(locus, ".", alleles[[l]][j], alleles[[l]][k]))
        }
        colnames(X[[l]]) <- tmp

        ## fill the design matrix
        ## caution, not all cases are present in a bi-parental family!
        if(jm$seg[l] == "<abxcd>"){ # a b c d aa ab ac ad bb bc bd cc cd dd
          idx <- which(jm[l,-(1:2)] == "ac")
          if(length(idx) > 0)
            X[[l]][idx, c(1,3,7)] <- matrix(c(1,1,1), nrow=length(idx),
                                            ncol=3, byrow=TRUE)
          idx <- which(jm[l,-(1:2)] == "ad")
          if(length(idx) > 0)
            X[[l]][idx, c(1,4,8)] <- matrix(c(1,1,1), nrow=length(idx),
                                            ncol=3, byrow=TRUE)
          idx <- which(jm[l,-(1:2)] == "bc")
          if(length(idx) > 0)
            X[[l]][idx, c(2,3,10)] <- matrix(c(1,1,1), nrow=length(idx),
                                             ncol=3, byrow=TRUE)
          idx <- which(jm[l,-(1:2)] == "bd")
          if(length(idx) > 0)
            X[[l]][idx, c(2,4,11)] <- matrix(c(1,1,1), nrow=length(idx),
                                             ncol=3, byrow=TRUE)
        } else if(jm$seg[l] == "<efxeg>"){ # e f g ee ef eg ff fg gg
          idx <- which(jm[l,-(1:2)] == "ee")
          if(length(idx) > 0)
            X[[l]][idx, c(1,4)] <- matrix(c(2,1), nrow=length(idx), ncol=2,
                                          byrow=TRUE)
          idx <- which(jm[l,-(1:2)] == "eg")
          if(length(idx) > 0)
            X[[l]][idx, c(1,3,6)] <- matrix(c(1,1,1), nrow=length(idx), ncol=3,
                                            byrow=TRUE)
          idx <- which(jm[l,-(1:2)] == "ef")
          if(length(idx) > 0)
            X[[l]][idx, c(1,2,5)] <- matrix(c(1,1,1), nrow=length(idx), ncol=3,
                                            byrow=TRUE)
          idx <- which(jm[l,-(1:2)] == "fg")
          if(length(idx) > 0)
            X[[l]][idx, c(2,3,8)] <- matrix(c(1,1,1), nrow=length(idx), ncol=3,
                                            byrow=TRUE)
        } else if(jm$seg[l] == "<hkxhk>"){ # h k hh hk kk
          idx <- which(jm[l,-(1:2)] == "hh")
          if(length(idx) > 0)
            X[[l]][idx, c(1,3)] <- matrix(c(2,1), nrow=length(idx), ncol=2,
                                          byrow=TRUE)
          idx <- which(jm[l,-(1:2)] == "hk")
          if(length(idx) > 0)
            X[[l]][idx, c(1,2,4)] <- matrix(c(1,1,1), nrow=length(idx), ncol=3,
                                            byrow=TRUE)
          idx <- which(jm[l,-(1:2)] == "kk")
          if(length(idx) > 0)
            X[[l]][idx, c(2,5)] <- matrix(c(2,1), nrow=length(idx), ncol=2,
                                          byrow=TRUE)
        } else if(jm$seg[l] == "<lmxll>"){ # l m ll lm mm
          idx <- which(jm[l,-(1:2)] == "ll")
          if(length(idx) > 0)
            X[[l]][idx, c(1,3)] <- matrix(c(2,1), nrow=length(idx), ncol=2,
                                          byrow=TRUE)
          idx <- which(jm[l,-(1:2)] == "lm")
          if(length(idx) > 0)
            X[[l]][idx, c(1,2,4)] <- matrix(c(1,1,1), nrow=length(idx), ncol=3,
                                            byrow=TRUE)
        } else if(jm$seg[l] == "<nnxnp>"){ # n p nn np pp
          idx <- which(jm[l,-(1:2)] == "nn")
          if(length(idx) > 0)
            X[[l]][idx, c(1,3)] <- matrix(c(2,1), nrow=length(idx), ncol=2,
                                          byrow=TRUE)
          idx <- which(jm[l,-(1:2)] == "np")
          if(length(idx) > 0)
            X[[l]][idx, c(1,2,4)] <- matrix(c(1,1,1), nrow=length(idx), ncol=3,
                                            byrow=TRUE)
        }
      } else{
        msg <- paste0("'", parameterization, "' parameterization",
                      " not yet implemented")
        stop(msg)
      }

    } # end of if use.phase=FALSE

  } # end for loop over locus
  names(X) <- jm$locus

  if(all(! use.phase, ! is.null(constraints))){
    if(verbose > 0){
      msg <- "apply constraints..."
      write(msg, stdout())
    }
    if(constraints == "last-zero"){
      idx <- which(jm$seg == "<abxcd>")
      if(length(idx) > 0)
        for(l in idx)
          X[[l]] <- X[[l]][, -c(4,8,11,13,14)]
      idx <- which(jm$seg == "<efxeg>")
      if(length(idx) > 0)
        for(l in idx)
          X[[l]] <- X[[l]][, -c(3,6,8,9)]
      idx <- which(jm$seg == "<hkxhk>")
      if(length(idx) > 0)
        for(l in idx)
          X[[l]] <- X[[l]][, -c(2,4,5)]
      idx <- which(jm$seg == "<lmxll>")
      if(length(idx) > 0)
        for(l in idx)
          X[[l]] <- X[[l]][, -c(2,4,5)]
      idx <- which(jm$seg == "<nnxnp>")
      if(length(idx) > 0)
        for(l in idx)
          X[[l]] <- X[[l]][, -c(2,4,5)]
    }
  }
  out <- do.call(cbind, X)
  if(any(use.phase, is.null(constraints))){
    out <- cbind(rep(1, nb.inds), out)
    colnames(out)[1] <- "intercept"
  }

  if(rm.col.zeros){
    idx <- which(apply(out, 2, function(col.j){all(col.j == 0)}))
    if(length(idx) > 0){
      if(verbose > 0){
        msg <- paste0("remove ", length(idx), " column",
                      ifelse(length(idx) > 1, "s", ""),
                      " full of zeros...")
        write(msg, stdout())
      }
      out <- out[, -idx]
    }
  }

  if(all(! use.phase, rm.dom)){
    idx <- which(sapply(strsplit(colnames(out), ""), function(x){
      x[length(x)-1] != "."
    }))
    if(colnames(out)[1] == "intercept")
      idx <- idx[-1]
    if(verbose > 0){
      msg <- paste0("remove ", length(idx), " columns encoding dominance...")
      write(msg, stdout())
    }
    out <- out[, -idx]
  }

  return(out)
}

getSegregatingLocusPerParent <- function(tX){
  stopifnot(is.matrix(tX),
            ncol(tX) == 2)

  out <- list()

  obs1 <- ! is.na(tX[,1])
  obs2 <- ! is.na(tX[,2])
  par1.het <- tX[,1] == 1
  par2.het <- tX[,2] == 1

  hkhk <- obs1 & par1.het & obs2 & par2.het
  lmll <- obs1 & par1.het & ((obs2 & (! par2.het)) | (! obs2))
  nnnp <- ((obs1 & (! par1.het)) | (! obs1)) & obs2 & par2.het

  out$parent1 <- which(hkhk | lmll)
  out$parent2 <- which(hkhk | nnnp)

  return(out)
}

##' Convert genotypes
##'
##' Convert SNP genotypes of a bi-parental cross from "allele doses" for usage in the "ASMap" package.
##' @param X matrix of bi-allelic SNP genotypes encoded, for each SNP, in number of copies of its second allele, i.e. as allele doses in {0,1,2}, with genotypes in rows and SNPs in columns; the "second" allele is arbitrary, it can correspond to the minor (least frequent) or the major (most frequent) allele
##' @param tX matrix with SNPs in rows and genotypes in columns
##' @return list of two data frames, one per parental map
##' @author Timothee Flutre
##' @export
genoDoses2ASMap <- function(X=NULL, tX=NULL){
  stopifnot(xor(is.null(X), is.null(tX)))
  if(! is.null(X)){
    stopIfNotValidGenosDose(X=X, check.noNA=FALSE, check.notImputed=TRUE)
    tX <- t(X)
  }

  out <- list()

  convertAndDuplicate <- function(tX.par){
    tmp <- apply(tX.par, 2, function(x){
      x[is.na(x)] <- "U"
      gsub("0|2", "A", gsub("1", "B", x))
    })
    tmp.m <- apply(tX.par, 2, function(x){
      x[is.na(x)] <- "U"
      gsub("0|2", "B", gsub("1", "A", x))
    })
    rownames(tmp.m) <- paste0(rownames(tmp.m), "_m")
    return(as.data.frame(rbind(tmp, tmp.m), stringsAsFactors=FALSE))
  }

  idx.loc.pars <- getSegregatingLocusPerParent(tX[, 1:2])
  out$parent1 <- convertAndDuplicate(tX[idx.loc.pars[[1]], -(1:2)])
  out$parent2 <- convertAndDuplicate(tX[idx.loc.pars[[2]], -(1:2)])

  return(out)
}

##' Set-up a R/qtl "cross" object
##'
##' Set up a "cross" object from the \href{https://cran.r-project.org/package=qtl}{qtl} package.
##' @param gendat genotype data in the qtl format as a matrix with genotypes in rows and markers in columns; if genmap is provided, only the markers also on the map will be kept
##' @param cross.type type of cross ("bc" for a backcross, "4way" for a cross between two outbred parents)
##' @param genmap genetic map as a data frame with one row per marker and at least two columns named "linkage.group" and "locus"; if a third one, "genetic.distance", is absent, it will be created and filled with fake, incremental values; if NULL, a fake genetic map will be created with all markers in the same linkage group; it is assumed that all linkage groups correspond to autosomes
##' @param phenos data frame containing the phenotypes; if not NULL, the genotype identifiers should be in a column named "Genotype" or as row names; otherwise, will be created with such a column
##' @return object of class "cross" as defined in the qtl package
##' @author Timothee Flutre
##' @seealso \code{\link{joinMap2backcross}}
##' @export
setupQtlCrossObject <- function(gendat, cross.type, genmap=NULL, phenos=NULL){
  stopifnot(is.matrix(gendat),
            ! is.null(colnames(gendat)),
            ! is.null(rownames(gendat)),
            cross.type %in% c("bc", "4way"))
  if(! is.null(genmap))
    stopifnot(is.data.frame(genmap),
              ncol(genmap) >= 2,
              all(c("linkage.group", "locus") %in% colnames(genmap)))
  if(! is.null(phenos))
    stopifnot(is.data.frame(phenos),
              any("Genotype" %in% colnames(phenos),
                  ! is.null(rownames(phenos))))

  out <- list()
  class(out) <- c(cross.type, "cross")

  ## reformat input genetic map
  if(is.null(genmap)){
    marker.names <- colnames(gendat)
    genmap <- data.frame(linkage.group=rep(1, length(marker.names)),
                         locus=marker.names)
  }
  if(! "genetic.distance" %in% colnames(genmap))
    genmap$genetic.distance <-
      unlist(tapply(genmap$locus, factor(genmap$linkage.group), function(x){
        seq(0.0, 100, length.out=length(x))
      }))
  genmap <- convertFactorColumnsToCharacter(genmap)
  if(! is.character(genmap$linkage.group))
    genmap$linkage.group <- as.character(genmap$linkage.group)
  if(requireNamespace("gtools", quietly=TRUE))
    genmap <- genmap[gtools::mixedorder(genmap$linkage.group),]
  lg.names <- unique(as.character(genmap$linkage.group))
  genmap.list <- lapply(lg.names, function(lg.name){
    tmp <- genmap[genmap$linkage.group == lg.name, "genetic.distance"]
    names(tmp) <- genmap[genmap$linkage.group == lg.name, "locus"]
    tmp
  })
  names(genmap.list) <- lg.names
  lg2loc <- lapply(genmap.list, names)

  ## reformat input genotype data
  is.marker.on.map <- colnames(gendat) %in% genmap$locus
  if(! all(is.marker.on.map)){
    if(all(! is.marker.on.map)){
      msg <- "all markers are absent from the provided map"
      stop(msg)
    }
    msg <- paste0("keep only the ", sum(is.marker.on.map),
                  " markers also present on the genetic map")
    warning(msg)
    gendat <- gendat[, is.marker.on.map]
  }
  gendat.list <- lapply(lg.names, function(lg.name){
    gendat[, lg2loc[[lg.name]]]
  })
  names(gendat.list) <- lg.names

  ## set up the "geno" component
  out$geno <- lapply(lg.names, function(lg.name){
    list(data=gendat.list[[lg.name]],
         map=genmap.list[[lg.name]])
  })
  names(out$geno) <- lg.names
  for(lg.name in lg.names)
    class(out$geno[[lg.name]]) <- "A"

  ## set up the "pheno" component
  if(is.null(phenos)){
    out$pheno <- data.frame(Genotype=rownames(gendat))
  } else{
    stopifnot(nrow(phenos) == nrow(gendat))
    geno.ids.from.phenos <- NULL
    if("Genotype" %in% colnames(phenos)){
      stopifnot(! anyDuplicated(phenos$Genotype))
      geno.ids.from.phenos <- phenos$Genotype
    } else
      geno.ids.from.phenos <- rownames(phenos)
    stopifnot(length(geno.ids.from.phenos) == nrow(gendat),
              all(sort(geno.ids.from.phenos) == sort(rownames(gendat))))
    out$pheno <- phenos
  }

  return(out)
}

##' Haplotypes
##'
##' Check that the input is a valid list of matrices of bi-allelic marker haplotypes.
##' @param haplos input list
##' @param check.hasColNames logical
##' @param check.noNA logical
##' @author Timothee Flutre
##' @export
stopIfNotValidHaplos <- function(haplos, check.hasColNames=FALSE,
                                 check.noNA=TRUE){
  stopifnot(! missing(haplos),
            ! is.null(haplos),
            is.list(haplos),
            all(sapply(haplos, is.matrix)),
            all(sapply(haplos, is.numeric)),
            length(unique(sapply(haplos, nrow))) == 1,
            ifelse(check.hasColNames,
                   ! any(sapply(haplos, function(x){is.null(colnames(x))})),
                   TRUE),
            ifelse(check.noNA,
                   ! any(sapply(haplos, function(x){any(is.na(x))})),
                   TRUE),
            all(sapply(haplos, function(x){all(x %in% c(0,1))})))
}

##' Genotypes
##'
##' Check that the input is a valid matrix of bi-allelic SNP genotypes coded in allele doses.
##' @param X input matrix
##' @param check.hasColNames logical
##' @param check.hasRowNames logical
##' @param check.noNA logical
##' @param check.isDose logical
##' @param check.notImputed logical
##' @author Timothee Flutre
##' @export
stopIfNotValidGenosDose <- function(X, check.hasColNames=TRUE,
                                    check.hasRowNames=TRUE,
                                    check.noNA=TRUE,
                                    check.isDose=TRUE,
                                    check.notImputed=FALSE){
  stopifnot(! missing(X),
            ! is.null(X),
            is.matrix(X),
            is.numeric(X),
            ifelse(check.hasColNames,
                   ! is.null(colnames(X)) & ! anyDuplicated(colnames(X)),
                   TRUE),
            ifelse(check.hasRowNames,
                   ! is.null(rownames(X)),
                   TRUE),
            ifelse(check.noNA, ! any(is.na(X)), TRUE),
            ifelse(check.isDose,
                   sum(X < 0, na.rm=TRUE) == 0,
                   TRUE),
            ifelse(check.isDose,
                   sum(X > 2, na.rm=TRUE) == 0,
                   TRUE),
            ifelse(check.notImputed,
                   sum(X == 0, na.rm=TRUE) + sum(X == 1, na.rm=TRUE) +
                   sum(X == 2, na.rm=TRUE) + sum(is.na(X)) ==
                   length(X),
                   TRUE))
}

##' Convert genotypes
##'
##' Convert SNP genotypes from "allele doses" into "genotypic classes".
##' Let us denote by "A" (or A1) the first allele, and "B" (or A2) the second allele.
##' The specific nucleotides (A, C, G or T) corresponding to "A" and "B" are specified via the \code{allele} argument to the function.
##' Moreover, whether A is the major (most frequent) allele or not, doesn't matter for the function to work.
##' At a given SNP, genotypes being homozygous for the first allele are coded as "0" in allele dose, and will be converted into "AA" in genotypic class.
##' Those being heterozygous are "1" and will be "AB".
##' Finally, those being homozygous for the second allele are "2" and will be "BB".
##' The code is not particularly efficient, but at least it exists.
##' @param X matrix of bi-allelic SNP genotypes encoded, for each SNP, in number of copies of its second allele, i.e. as allele doses in {0,1,2}, with genotypes in rows and SNPs in columns; the "second" allele is arbitrary, it corresponds to the second column of \code{alleles}, which can be the minor or the major allele
##' @param tX matrix with SNPs in rows and genotypes in columns
##' @param alleles data.frame with SNPs in rows (names as row names) and alleles in columns (exactly 2 columns are required); the second column should correspond to the allele which number of copies is counted at each SNP in \code{X}; only the rows corresponding to SNPs in X or tX will be kept
##' @param na.string character used to replace NA values
##' @param nb.cores the number of cores to use, i.e. at most how many child processes will be run simultaneously (not on Windows)
##' @param verbose verbosity level (0/1)
##' @return data.frame with SNPs in rows and genotypes in columns
##' @author Timothee Flutre
##' @seealso \code{\link{genoClasses2genoDoses}}
##' @export
genoDoses2genoClasses <- function(X=NULL, tX=NULL, alleles, na.string="--",
                                  nb.cores=1, verbose=1){
  stopifnot(xor(is.null(X), is.null(tX)),
            is.data.frame(alleles),
            ncol(alleles) == 2,
            ! is.null(row.names(alleles)),
            all(! is.na(alleles[,1])),
            all(! is.na(alleles[,2])))
  if(! is.null(X)){
    stopIfNotValidGenosDose(X=X, check.noNA=FALSE)
    tX <- t(X)
  }
  alleles <- convertFactorColumnsToCharacter(alleles)
  alleles <- alleles[rownames(tX),]

  if(verbose > 0){
    msg <- paste0("convert ", nrow(tX), " SNPs for ",
                  ncol(tX), " genotypes...")
    write(msg, stdout())
  }

  ## -------------------------------------------------------
  ## vectorized version:

  alleles$gclass <- paste0(alleles[,1], alleles[,2])
  if(verbose > 1)
    print(table(alleles$gclass))

  out <- parallel::mclapply(unique(alleles$gclass), function(gclass){
    idx.gclass <- which(alleles$gclass == gclass)
    tmp <- c(tX[idx.gclass,]) # work with a vector

    ## convert homozygous for the first allele (AA or A1A1)
    idx <- which(tX[idx.gclass,] == 0)
    if(length(idx) > 0)
      tmp[idx] <- paste0(alleles[idx.gclass[1], 1],
                         alleles[idx.gclass[1], 1], collapse="")

    ## convert heterozygous (AB or A1A2)
    idx <- which(tX[idx.gclass,] == 1)
    if(length(idx) > 0)
      tmp[idx] <- paste0(alleles[idx.gclass[1], 1],
                         alleles[idx.gclass[1], 2], collapse="")

    ## convert homozygous for the second allele (BB or A2A2)
    idx <- which(tX[idx.gclass,] == 2)
    if(length(idx) > 0)
      tmp[idx] <- paste0(alleles[idx.gclass[1], 2],
                         alleles[idx.gclass[1], 2], collapse="")

    ## convert missing (NA)
    idx <- which(is.na(tX[idx.gclass,]))
    if(length(idx) > 0)
      tmp[idx] <- na.string

    ## transform into a matrix
    matrix(data=tmp,
           nrow=nrow(tX[idx.gclass, , drop=FALSE]),
           ncol=ncol(tX[idx.gclass, , drop=FALSE]),
           dimnames=dimnames(tX[idx.gclass,, drop=FALSE]))
  }, mc.cores=nb.cores)

  ## reformat the output
  out <- do.call(rbind, out)
  out <- out[rownames(tX), colnames(tX)]

  ## -------------------------------------------------------
  ## un-vectorized version:

  ## out <- as.data.frame(tX, row.names=rownames(tX), col.names=colnames(tX))

  ## if(verbose > 0){
  ##   idx.rows <- 1:nrow(tX)
  ##   pb <- utils::txtProgressBar(min=0, max=length(idx.rows), style=3)
  ##   tmp <- stats::setNames(object=cut(x=idx.rows, breaks=10, labels=FALSE),
  ##                          nm=idx.rows) # to update pb no more than 10 times
  ## }

  ## ## for each SNP
  ## for(i in 1:nrow(tX)){

  ##   ## convert homozygous for the first allele (AA or A1A1)
  ##   idx <- which(tX[i,] == 0)
  ##   if(length(idx) > 0)
  ##     out[i, idx] <- paste0(alleles[i, 1], alleles[i, 1], collapse="")

  ##   ## convert heterozygous (AB or A1A2)
  ##   idx <- which(tX[i,] == 1)
  ##   if(length(idx) > 0)
  ##     out[i, idx] <- paste0(alleles[i, 1], alleles[i, 2], collapse="")

  ##   ## convert homozygous for the second allele (BB or A2A2)
  ##   idx <- which(tX[i,] == 2)
  ##   if(length(idx) > 0)
  ##     out[i, idx] <- paste0(alleles[i, 2], alleles[i, 2], collapse="")

  ##   ## convert missing (NA)
  ##   idx <- which(is.na(tX[i,]))
  ##   if(length(idx) > 0)
  ##     out[i, idx] <- na.string

  ##   if(verbose > 0){
  ##     ## utils::setTxtProgressBar(pb, i) # update pb at each loop iteration
  ##     if(i %in% as.numeric(names(tmp)[cumsum(table(tmp))]))
  ##       utils::setTxtProgressBar(pb, which(idx.rows == i))
  ##   }
  ## }
  ## if(verbose > 0)
  ##   close(pb)

  return(out)
}

##' Index the dose of each SNP genotype
##'
##' Return if a given SNP enotype is 0 or 1, even if it was imputed.
##' @param X matrix of bi-allelic SNP genotypes encoded, for each SNP, in number of copies of its second allele, i.e. as allele doses in {0,1,2}, with genotypes in rows and SNPs in columns; the "second" allele is arbitrary, which can be the minor or the major allele
##' @param boundaries vector with interval boundaries of allele dose
##' @return matrix of logicals with two columns, is.0 and is.1
##' @author Timothee Flutre
##' @export
indexGenoDoses <- function(X, boundaries=seq(from=0, to=2, length.out=4)){
  stopIfNotValidGenosDose(X=X, check.hasColNames=FALSE,
                          check.hasRowNames=FALSE, check.noNA=FALSE,
                          check.notImputed=FALSE)
  stopifnot(is.vector(boundaries),
            length(boundaries) == 4)

  out <- matrix(data=rep(NA, 2 * length(X)),
                ncol=2)
  colnames(out) <- c("is.0", "is.1")

  out[,"is.0"] <- (X <= boundaries[2]) # homozygotes for the first allele (ref)
  out[,"is.1"] <- (X > boundaries[2] & X <= boundaries[3]) # heterozygotes
  ## out[,"is.2"] <- (X > boundaries[3]) # homozygotes for the second allele (alt)

  return(out)
}

##' Convert genotypes
##'
##' Convert SNP genotypes from "allele doses" into the \href{http://samtools.github.io/hts-specs/}{VCF} format.
##' @param X matrix of bi-allelic SNP genotypes encoded, for each SNP, in number of copies of its second allele, i.e. as allele doses in {0,1,2}, with genotypes in rows and SNPs in columns; the "second" allele corresponds to the second column of \code{alleles} and will be interpreted as 'alt'
##' @param snp.coords data.frame with SNP identifiers as row names, and two columns, "chr" and "pos" (or "coord"); columns as factors will be converted into characters
##' @param alleles data.frame (or matrix) with SNPs in rows (names as row names) and alleles in columns (exactly 2 columns are required); the first column will be interpreted as 'ref' and the second column, which should correspond to the allele which number of copies is counted at each SNP in \code{X}, will be interpreted as 'alt'; columns as factors will be converted into characters
##' @param file.date date to indicate into the object
##' @param si output from the "seqinfo" function from the GenomeInfoDb package
##' @param verbose verbosity level (0/1)
##' @return CollapsedVCF (see pkg \href{http://bioconductor.org/packages/VariantAnnotation/}{VariantAnnotation})
##' @author Timothee Flutre
##' @export
genoDoses2Vcf <- function(X, snp.coords, alleles,
                          file.date=format(Sys.Date(), "%Y%m%d"),
                          si=NULL,
                          verbose=1){
  requireNamespaces(c("IRanges", "S4Vectors", "Biostrings",
                      "VariantAnnotation"))
  stopIfNotValidGenosDose(X=X, check.noNA=FALSE)
  stopifnot(.isValidSnpCoords(snp.coords),
            all(colnames(X) %in% rownames(snp.coords)),
            is.data.frame(alleles) | is.matrix(alleles),
            ncol(alleles) == 2,
            ! is.null(row.names(alleles)),
            all(colnames(X) %in% rownames(alleles)),
            is.character(file.date))
  if(! is.null(si)){
    requireNamespace("GenomeInfoDb", quietly=TRUE)
    stopifnot(methods::is(si, "Seqinfo"))
  }

  ind.ids <- rownames(X)
  nb.inds <- length(ind.ids)
  snp.ids <- colnames(X)
  nb.snps <- length(snp.ids)
  snp.coords <- droplevels(snp.coords[snp.ids,])
  stopifnot(all(! is.na(snp.coords[,1])),
            all(! is.na(snp.coords[,2])))
  snp.coords <- convertFactorColumnsToCharacter(snp.coords)
  if("coord" %in% colnames(snp.coords))
    colnames(snp.coords)[colnames(snp.coords) == "coord"] <- "pos"
  if(is.matrix(alleles))
    alleles <- as.data.frame(alleles)
  alleles <- droplevels(alleles[snp.ids,])
  alleles <- convertFactorColumnsToCharacter(alleles)
  stopifnot(all(! is.na(alleles[,1])),
            all(! is.na(alleles[,2])))
  if(verbose > 0){
    msg <- paste0("convert SNP genotypes of ", nb.inds, " genotypes at ",
                  nb.snps, " SNPs into VCF ...")
    write(msg, stdout())
  }

  gr <- seqIdStartEnd2GRanges(seq.id=snp.coords$chr,
                              seq.start=snp.coords$pos,
                              seq.end=snp.coords$pos,
                              subseq.name=snp.ids)
  Df <- S4Vectors::DataFrame(Samples=1:nb.inds,
                             row.names=ind.ids)
  vcf <- VariantAnnotation::VCF(rowRanges=gr, colData=Df)

  VariantAnnotation::header(vcf) <-
    VariantAnnotation::VCFHeader(samples=ind.ids)
  VariantAnnotation::meta(VariantAnnotation::header(vcf)) <-
    IRanges::DataFrameList(
                 META=S4Vectors::DataFrame(Value=c("VCFv4.2",
                                                   file.date),
                                           row.names=c("fileformat",
                                                       "fileDate")))
  ## also possible to add reference and assembly to META
  VariantAnnotation::geno(VariantAnnotation::header(vcf)) <-
    S4Vectors::DataFrame(Number="1", Type="String",
                         Description="Genotype",
                         row.names="GT")
  VariantAnnotation::ref(vcf) <- Biostrings::DNAStringSet(alleles[,1])
  VariantAnnotation::alt(vcf) <- Biostrings::DNAStringSetList(
                                                 as.list(alleles[,2]))
  VariantAnnotation::fixed(vcf)[c("REF", "ALT")]
  GenomeInfoDb::seqlevels(vcf) <-
    GenomeInfoDb::sortSeqlevels(GenomeInfoDb::seqlevels(vcf))
  if(! is.null(si)){
    GenomeInfoDb::seqlevels(si) <-
      GenomeInfoDb::sortSeqlevels(GenomeInfoDb::seqlevels(si))
    stopifnot(all(GenomeInfoDb::seqlevels(vcf) ==
                  GenomeInfoDb::seqlevels(si)))
    GenomeInfoDb::seqinfo(vcf) <- si
  }

  idx <- indexGenoDoses(X)
  X[idx[,"is.0"]] <- "0/0"
  X[idx[,"is.1"]] <- "0/1"
  X[! idx[,"is.0"] & ! idx[,"is.1"]] <- "1/1"
  is.NA <- is.na(X)
  if(any(is.NA))
    X[is.NA] <- "./."
  VariantAnnotation::geno(vcf) <- S4Vectors::SimpleList(GT=t(X))

  return(vcf)
}

##' Missing genotypes
##'
##' Calculate the frequency of missing genotypes for each marker.
##' @param X matrix of marker genotypes, with genotypes in rows and markers in columns; missing values should be encoded as NA; if not NULL, will be used in priority even if \code{vcf.file} is not NULL
##' @param vcf.file path to the VCF file (if the bgzip index doesn't exist in the same directory, it will be created); used only if \code{X=NULL}
##' @param yieldSize number of records to yield each time the VCF file is read from (see ?TabixFile)
##' @return vector
##' @author Timothee Flutre
##' @examples
##' \dontrun{## simulate fake SNP genotypes
##' set.seed(1859)
##' X <- simulGenosDose(nb.genos=200, nb.snps=10^3)
##'
##' ## randomly create missing SNP genotypes
##' idx <- sample.int(n=length(X), size=10^3)
##' X[idx] <- NA
##'
##' miss.snp <- calcFreqMissSnpGenosPerSnp(X=X)
##' hist(miss.snp, xlab="proportion of missing data at a given SNP,\nmeasured across all individuals", ylab="number of SNPs")
##' }
##' @export
calcFreqMissSnpGenosPerSnp <- function(X=NULL, vcf.file=NULL, yieldSize=10000){
  stopifnot(! all(is.null(X), is.null(vcf.file)))

  out <- c()

  if(! is.null(X)){
    stopIfNotValidGenosDose(X, check.hasColNames=FALSE,
                            check.hasRowNames=FALSE,
                            check.noNA=FALSE, check.isDose=FALSE)
    N <- nrow(X) # number of genotypes
    out <- colSums(is.na(X)) / N
  } else{
    requireNamespaces(c("Rsamtools", "VariantAnnotation"))
    stopifnot(file.exists(vcf.file),
              ! is.na(yieldSize))
    if(! file.exists(paste0(vcf.file, ".tbi")))
      Rsamtools::indexTabix(file=vcf.file, format="vcf")
    tabix.file <- Rsamtools::TabixFile(file=vcf.file,
                                       yieldSize=yieldSize)
    open(tabix.file)
    vcf <- VariantAnnotation::readVcf(file=tabix.file, genome="")
    while(nrow(vcf)){
      out <- c(out, calcFreqNaVcf(vcf=vcf, with.coords=FALSE))
      vcf <- VariantAnnotation::readVcf(file=tabix.file, genome="")
    }
    close(tabix.file)
  }

  return(out)
}

##' Missing genotypes
##'
##' Calculate the frequency of missing marker genotypes for each genotype.
##' @param X matrix of marker genotypes, with genotypes in rows and markers in columns; missing values should be encoded as NA
##' @return vector
##' @author Timothee Flutre
##' @export
calcFreqMissSnpGenosPerGeno <- function(X){
  stopIfNotValidGenosDose(X, check.hasColNames=FALSE,
                          check.hasRowNames=FALSE,
                          check.noNA=FALSE, check.isDose=FALSE)

  P <- ncol(X) # number of markers
  out <- rowSums(is.na(X)) / P

  return(out)
}

##' Missing genotypes
##'
##' Plot missing marker genotypes as a grid, represented in black if missing, white otherwise.
##' @param X matrix of marker genotypes, with genotypes in rows and markers in columns; missing values should be encoded as NA
##' @param main an overall title for the plot
##' @param xlab a title for the x axis
##' @param ylab a title for the y axis
##' @return nothing
##' @author Timothee Flutre
##' @export
plotGridMissGenos <- function(X, main="Missing marker genotypes",
                              xlab="Genotypes", ylab="markers"){
  stopIfNotValidGenosDose(X, check.hasColNames=FALSE,
                          check.hasRowNames=FALSE,
                          check.noNA=FALSE, check.isDose=FALSE)

  graphics::image(1:nrow(X), 1:ncol(X), is.na(X), col=c("white","black"),
        main=main, xlab=xlab, ylab=ylab)
}

##' Missing genotypes
##'
##' Discard all markers with at least one missing genotype.
##' @param X matrix of marker genotypes, with genotypes in rows and markers in columns; missing values should be encoded as NA
##' @param verbose verbosity level (0/1)
##' @return matrix similar to X but possibly with less columns
##' @author Timothee Flutre
##' @seealso \code{\link{imputeGenosWithMean}}
##' @export
discardMarkersMissGenos <- function(X, verbose=1){
  stopIfNotValidGenosDose(X, check.hasColNames=FALSE, check.hasRowNames=FALSE,
                          check.noNA=FALSE, check.isDose=FALSE)

  if(verbose > 0){
    msg <- "discard markers with missing data ..."
    write(msg, stdout())
  }

  tokeep <- apply(X, 2, function(x){
    sum(is.na(x)) == 0
  })
  if(verbose > 0){
    msg <- paste0("nb of markers to keep: ", sum(tokeep))
    write(msg, stdout())
  }

  X.out <- X[, tokeep]

  return(X.out)
}

##' Quantiles of allelic R2 after binning
##'
##' Calculate the quantiles of allelic R2 after binning.
##' @param snp.data numeric vector of SNP data (one value per SNP) with SNP identifiers as names
##' @param snp.bins factor as returned by \code{\link{cut}} applied to SNP data (percentage of missing data, minor allele frequencies, etc) with SNP identifiers as names
##' @param quant.probs numeric vector of probabilities with values in [0,1] to be used with \code{\link[stats]{quantile}}
##' @param verbose verbosity level (0/1)
##' @return list with the quantiles per bin, as well as the midpoint of each bin
##' @author Timothee Flutre
##' @export
quantilesBinnedSnpData <- function(snp.data, snp.bins,
                                   quant.probs=seq(0, 1, 0.25),
                                   verbose=1){
  stopifnot(is.vector(snp.data),
            is.numeric(snp.data),
            ! is.null(names(snp.data)),
            is.factor(snp.bins),
            length(snp.bins) == length(snp.data),
            ! is.null(names(snp.bins)),
            all(names(snp.bins) == names(snp.data)),
            is.vector(quant.probs),
            is.numeric(quant.probs),
            all(quant.probs >= 0),
            all(quant.probs <= 1),
            ! anyDuplicated(quant.probs))

  ## quantiles of SNP data
  quant.probs <- sort(quant.probs)
  quant.bin.snp.dat <-
    tapply(snp.data, list(snp.bins), function(subset.dat){
      stats::quantile(subset.dat, probs=quant.probs, na.rm=TRUE)
    })
  for(i in which(sapply(quant.bin.snp.dat, is.null)))
    quant.bin.snp.dat[[i]] <- rep(NA, length(quant.probs))
  quant.bin.snp.dat <- do.call(rbind, quant.bin.snp.dat)
  quant.bin.snp.dat <- quant.bin.snp.dat[! is.na(quant.bin.snp.dat[,1]),]

  ## midpoint of each bin
  bin.lowers <- regmatches(levels(snp.bins),
                           regexec("\\(([0-9\\.]{1,4}),",
                                   levels(snp.bins)))
  bin.lowers <- as.numeric(sapply(bin.lowers, function(x){x[length(x)]}))
  bin.uppers <- regmatches(levels(snp.bins),
                           regexec(",([0-9\\.]{1,4})\\]",
                                   levels(snp.bins)))
  bin.uppers <- as.numeric(sapply(bin.uppers, function(x){x[length(x)]}))
  bin.lengths <- bin.uppers - bin.lowers
  bin.mids <- bin.lowers + (bin.lengths / 2)
  names(bin.mids) <- levels(snp.bins)
  bin.mids <- bin.mids[rownames(quant.bin.snp.dat)]

  return(list(quant.bin.snp.dat=quant.bin.snp.dat,
              bin.mids=bin.mids))
}

##' Convert imputed genotypes to 0/1/2
##'
##' Converts imputed genotypes to 0/1/2.
##' @param X matrix of bi-allelic SNP genotypes encoded in allele doses in [0,2], with genotypes in rows and SNPs in columns; missing values should be encoded as NA
##' @return matrix
##' @author Timothee Flutre
##' @export
convertImputedTo012 <- function(X){
  stopIfNotValidGenosDose(X=X, check.hasColNames=FALSE,
                          check.hasRowNames=FALSE, check.noNA=FALSE,
                          check.notImputed=FALSE)

  boundaries <- seq(from=0, to=2, length.out=4)
  is.0 <- (X <= boundaries[2]) # homozygotes for the first allele (ref)
  is.1 <- (X > boundaries[2] & X <= boundaries[3]) # heterozygotes
  is.2 <- (X > boundaries[3]) # homozygotes for the second allele (alt)
  idx <- indexGenoDoses(X)
  X[idx[,"is.0"]] <- 0
  X[idx[,"is.1"]] <- 1
  X[! idx[,"is.0"] & ! idx[,"is.1"]] <- 2
  is.NA <- is.na(X)
  if(any(is.NA))
  X[is.NA] <- NA

  return(X)
}

##' Allele frequencies
##'
##' Estimate allele frequencies of bi-allelic SNPs.
##' If previous sample sizes and counts per SNP are provided, it updates them.
##' @param X matrix of bi-allelic SNP genotypes encoded in allele doses in [0,2], with genotypes in rows and SNPs in columns; missing values should be encoded as NA
##' @param allow.updating if TRUE, the sample sizes and counts per SNP will be kept to allow further updating
##' @param prev.nb.genos vector of number of genotypes (sample sizes per SNP); the length and names should correspond to the columns of the X matrix
##' @param prev.nb.refalls vector of number of reference alleles (counts per SNP); the length and names should correspond to the columns of the X matrix
##' @return vector, with sample sizes and counts per SNP as attributes if updating is allowed
##' @seealso \code{\link{estimSnpMaf}}
##' @author Timothee Flutre
##' @export
estimSnpAf <- function(X, allow.updating=FALSE,
                       prev.nb.genos=NULL, prev.nb.refalls=NULL){
  stopIfNotValidGenosDose(X, check.hasColNames=FALSE,
                          check.hasRowNames=FALSE, check.noNA=FALSE)
  stopifnot(is.logical(allow.updating),
            xor(all(is.null(prev.nb.genos),
                    is.null(prev.nb.refalls)),
                all(! is.null(prev.nb.genos),
                    ! is.null(prev.nb.refalls))))
  if(! is.null(prev.nb.genos)){
    stopifnot(is.vector(prev.nb.genos),
              length(prev.nb.genos) == ncol(X),
              ! is.null(colnames(X)),
              ! is.null(names(prev.nb.genos)),
              all(names(prev.nb.genos) == colnames(X)),
              all(! is.na(prev.nb.genos)),
              all(prev.nb.genos >= 0, na.rm=TRUE))
    stopifnot(is.vector(prev.nb.refalls),
              length(prev.nb.refalls) == ncol(X),
              ! is.null(names(prev.nb.refalls)),
              all(names(prev.nb.refalls) == colnames(X)),
              all(! is.na(prev.nb.refalls)),
              all(prev.nb.refalls >= 0, na.rm=TRUE))
  }

  ## compute inputs allowing allele frequencies to be easily updated
  tot.nb.genos <- colSums(! is.na(X))
  tot.nb.refalls <- colSums(X, na.rm=TRUE)

  ## update inputs if necessary
  if(! is.null(prev.nb.genos)){
    tot.nb.genos <- tot.nb.genos + prev.nb.genos
    tot.nb.refalls <- tot.nb.refalls + prev.nb.refalls
  }

  ## compute allele frequencies
  afs <- tot.nb.refalls / (2 * tot.nb.genos)

  ## keep inputs to allow update later on
  if(allow.updating){
    attr(x=afs, which="nb.genos") <- tot.nb.genos
    attr(x=afs, which="nb.refalls") <- tot.nb.refalls
  }

  return(afs)
}

##' Minor allele frequencies
##'
##' Estimate minor allele frequencies of bi-allelic SNPs.
##' @param X matrix of bi-allelic SNP genotypes encoded in allele doses in [0,2], with genotypes in rows and SNPs in columns; missing values should be encoded as NA
##' @param afs vector of allele frequencies; if NULL, X should be specified, and the allele frequencies will be estimated with \code{\link{estimSnpAf}}
##' @return vector
##' @seealso \code{\link{estimSnpAf}}
##' @author Timothee Flutre
##' @export
estimSnpMaf <- function(X=NULL, afs=NULL){
  stopifnot(xor(is.null(X), is.null(afs)))

  if(is.null(afs))
    afs <- estimSnpAf(X)

  mafs <- apply(rbind(afs, 1 - afs), 2, min)

  return(mafs)
}

##' Minor allele frequencies
##'
##' Recode bi-allelic SNP genotypes in terms of the number of copies of the minor allele.
##' Let be a bi-allelic SNP with two alleles, A, the first allele, and B, the second allele.
##' Here, "first" and "second" are arbitrary, meaning that "A" as "B" can be the minor (least frequent) allele.
##' The SNP genotypes are originally coded in number of copies of the second allele: that is, if a genotype is "AA", its allele dose is \code{0}; if "AB", then \code{1}; and if "BB", then \code{2}.
##' Let us now assume that, in fact, "A" is the minor (least frequent) allele, then, when recoding SNP genotypes in terms of minor alleles, "AA" will correspond to \code{2}, "AB" to \code{1} and "BB" to \code{0}.
##' @param X matrix of bi-allelic SNP genotypes encoded, for each SNP, in number of copies of its second allele, i.e. as allele doses in [0,2], with genotypes in rows and SNPs in columns; missing values should be encoded as NA; the "second" allele is arbitrary, it corresponds to the second column of \code{alleles}, which can be the minor or the major allele
##' @param alleles data.frame with SNPs in rows (names as row names) and alleles in columns (exactly 2 columns are required); the second column should correspond to the allele which number of copies is counted at each SNP in \code{X}
##' @param verbose verbosity level (0/1)
##' @return list with a matrix of SNP genotypes encoded, for each SNP, in number of copies of its minor allele, with genotypes in rows and SNPs in columns, a data.frame of alleles which columns are named "minor" and "major", and a vector of logicals indicating which SNPs have been recoded
##' @author Timothee Flutre
##' @export
recodeGenosMinorSnpAllele <- function(X, alleles, verbose=1){
  stopIfNotValidGenosDose(X, check.noNA=FALSE)
  stopifnot(is.data.frame(alleles),
            ncol(alleles) == 2,
            ! is.null(row.names(alleles)),
            all(colnames(X) %in% rownames(alleles)))

  alleles <- alleles[colnames(X),] # put in same order
  alleles <- convertFactorColumnsToCharacter(alleles)
  out <- list(X=X,
              alleles=alleles,
              recoded=stats::setNames(rep(FALSE, ncol(X)),
                                      colnames(X)))
  colnames(out$alleles) <- c("major", "minor")

  afs <- estimSnpAf(X=X)
  idx.snps <- which(afs > 0.5)

  if(length(idx.snps) > 0){
    if(verbose > 0){
      msg <- paste0(length(idx.snps), " SNP",
                    ifelse(length(idx.snps) > 1, "s", ""),
                    " to recode out of ", ncol(X))
      write(msg, stdout())
    }
    out$X[, idx.snps] <- abs(X[, idx.snps] - 2)
    out$alleles[idx.snps, "minor"] <- alleles[idx.snps, 1]
    out$alleles[idx.snps, "major"] <- alleles[idx.snps, 2]
    out$recoded[idx.snps] <- TRUE
  }

  return(out)
}

##' Genotypic classes
##'
##' Count the genotypic classes per SNP (not imputed, i.e. 0, 1, 2 and NA).
##' @param X matrix of bi-allelic SNP genotypes encoded in allele doses in {0,1,2}, with genotypes in rows and SNPs in columns; missing values should be encoded as NA
##' @return matrix with 4 columns and as many rows as there are SNPs
##' @author Timothee Flutre
##' @export
countGenotypicClasses <- function(X){
  stopIfNotValidGenosDose(X, check.hasColNames=FALSE, check.hasRowNames=FALSE,
                          check.noNA=FALSE, check.notImputed=TRUE)

  counts <- t(apply(X, 2, function(Xp){
    c(sum(Xp == 0, na.rm=TRUE), sum(Xp == 1, na.rm=TRUE),
      sum(Xp == 2, na.rm=TRUE), sum(is.na(Xp)))
  }))
  ## colnames(counts) <- c("homozygotes.first", "heterozygotes",
  ##                       "homozygotes.second", "missing")
  colnames(counts) <- c("0", "1", "2", "NA")

  return(counts)
}

##' Chi-squared for Hardy-Weinberg
##'
##' Calculate the chi2 statistic to test for Hardy-Weinberg equilibrium.
##' See \href{http://www.karger.com/?doi=10.1159/000108939}{Graffelman & Camarena (2007)} and \href{https://dx.doi.org/10.1002/gepi.20612}{Shriner (2011)}.
##' @param X matrix of bi-allelic SNP genotypes encoded as allele doses in {0,1,2}, with genotypes in rows and SNPs in columns; missing values should be encoded as NA
##' @param mafs vector of minor allele frequencies; if NULL, the frequencies will be estimated by feeding \code{X} to \code{\link{estimSnpMaf}}
##' @param c continuity correction (\code{0} means "no correction"; usually, \code{0.5} is used, from \href{http://dx.doi.org/10.2307/2983604}{Yates (1934)})
##' @param thresh.c threshold on minor allele frequencies below which the continuity correction isn't applied (used when \code{c > 0}); see \href{https://dx.doi.org/10.1016\%2Fj.ajhg.2009.11.019}{Graffelman (2010)}
##' @param calc.with.D calculate the chi2 statistic with D, the deviation from independence for the heterozygote, as in equation 1 from Graffelman & Camarena (2007), which only requires the number of heterozygotes; otherwise, use equation 4
##' @param calc.pval calculate the p values associated with the test statistics (Chi-squared distribution with one degree of freedom)
##' @return matrix
##' @author Timothee Flutre
##' @seealso \code{\link{estimSnpMaf}}, \code{\link{countGenotypicClasses}}
##' @examples
##' \dontrun{set.seed(1859)
##' library(scrm)
##' nb.genos <- 100
##' Ne <- 10^4
##' chrom.len <- 10^5
##' mu <- 10^(-8)
##' c <- 10^(-8)
##' genomes <- simulCoalescent(nb.inds=nb.genos,
##'                            pop.mut.rate=4 * Ne * mu * chrom.len,
##'                            pop.recomb.rate=4 * Ne * c * chrom.len,
##'                            chrom.len=chrom.len)
##' X <- genomes$genos
##' out <- chiSqSnpGenos(X)
##' head(out)
##' sum(p.adjust(out[,"pvalue"], "BH") <= 0.05)
##' ## library(HardyWeinberg) # available on CRAN
##' ## cts <- countGenotypicClasses(X=X)[, -4]
##' ## colnames(cts) <- c("AA","AB","BB")
##' ## out2 <- HWChisqMat(X=cts, cc=0.5, verbose=FALSE)
##' ## lapply(out2, head)
##' ## HWTernaryPlot(X=cts, region=2)
##' }
##' @export
chiSqSnpGenos <- function(X, mafs=NULL, c=0.5, thresh.c=0.01, calc.with.D=TRUE,
                          calc.pval=TRUE){
  stopIfNotValidGenosDose(X=X, check.hasColNames=TRUE, check.hasRowNames=FALSE,
                          check.noNA=FALSE, check.isDose=TRUE,
                          check.notImputed=TRUE)
  if(! is.null(mafs))
    stopifnot(is.vector(mafs),
              is.numeric(mafs),
              ! is.null(names(mafs)),
              all(mafs <= 0.5),
              all(colnames(X) %in% names(mafs)))
  stopifnot(is.numeric(c),
            length(c) == 1,
            all(c >= 0, c <= 1),
            is.numeric(thresh.c),
            length(thresh.c) == 1,
            all(thresh.c >= 0, thresh.c <= 1),
            is.logical(calc.with.D))

  out <- NULL

  N <- nrow(X)
  if(is.null(mafs)){
    mafs <- estimSnpMaf(X=X)
  } else
    mafs <- mafs[colnames(X)] # put in same order
  p.hat <- mafs
  q.hat <- 1 - p.hat
  out <- as.matrix(p.hat)
  colnames(out) <- "maf"

  ## calculate the chi2 test statistics
  chi2 <- stats::setNames(rep(NA, ncol(X)), colnames(X))
  if(calc.with.D){
    ## observed counts of heterozygotes
    obs.AB <- apply(X, 2, function(Xp){
      sum(Xp == 1, na.rm=TRUE)
    })

    ## expected counts of heterozygotes under HWE
    N.noNA <- N - apply(X, 2, function(Xp){sum(is.na(Xp))})
    exp.AB <- 2 * N.noNA * p.hat * q.hat

    D <- 0.5 * (obs.AB - exp.AB)
    numerator <- D^2
    if(c > 0){
      idx <- which(mafs >= thresh.c)
      if(length(idx) > 0)
        numerator[idx] <- numerator[idx] -
          2 * c * abs(D[idx]) * (1 - p.hat[idx] * q.hat[idx]) +
          c^2 * (1 - (3 /2) * p.hat[idx] * q.hat[idx])
    }
    chi2 <- numerator / (p.hat^2 * q.hat^2 * N.noNA)

  } else{ # calc.with.D=FALSE
    observed.counts <- countGenotypicClasses(X=X)
    obs.AA <- observed.counts[, "0"]
    obs.AB <- observed.counts[, "1"]
    obs.BB <- observed.counts[, "2"]
    N <- nrow(X)
    N.noNA <- N - observed.counts[, "NA"]
    exp.AA <- N.noNA * q.hat^2
    exp.AB <- 2 * N.noNA * p.hat * q.hat
    exp.BB <- N.noNA * p.hat^2
    chi2 <- (obs.AA - exp.AA)^2 / exp.AA +
      (obs.AB - exp.AB)^2 / exp.AB +
      (obs.BB - exp.BB)^2 / exp.BB
    if(c > 0){ # eq 5 in Graffelman & Camarena (2007)
      idx <- which(mafs >= thresh.c)
      if(length(idx) > 0)
        chi2[idx] <- (abs(obs.AA[idx] - exp.AA[idx]) - c)^2 / exp.AA[idx] +
          (abs(obs.AB[idx] - exp.AB[idx]) - c)^2 / exp.AB[idx] +
          (abs(obs.BB[idx] - exp.BB[idx]) - c)^2 / exp.BB[idx]
    }
  }
  out <- cbind(out, chi2=chi2)

  if(calc.pval)
    out <- cbind(out, pvalue=stats::pchisq(q=chi2, df=1, lower.tail=FALSE))

  return(out)
}

##' Allele frequencies
##'
##' Plot the histogram of the allele frequency per SNP.
##' @param X matrix of bi-allelic SNP genotypes encoded, for each SNP, in number of copies of its minor allele, i.e. as allele doses in [0,2], with genotypes in rows and SNPs in columns; if X is not NULL, minor allele frequencies will be estimated via \code{\link{estimSnpAf}}
##' @param afs vector of allele frequencies (not used if X is not NULL)
##' @param main string for the main title
##' @param xlim limits of the x-axis
##' @param col color for the bars
##' @param border color for the border of the bars
##' @param las see ?par
##' @param breaks see ?hist
##' @param add.ml.beta if TRUE, add a curve corresponding to a Beta distribution with parameters inferred by maximum likelihood
##' @param verbose verbosity level (0/1)
##' @param ... arguments to be passed to hist()
##' @return invisible list with the output of \code{hist} and \code{optim} (if \code{add.ml.beta=TRUE})
##' @author Timothee Flutre
##' @export
plotHistAllelFreq <- function(X=NULL, afs=NULL, main=NULL, xlim=c(0,1),
                              col="grey", border="white", las=1,
                              breaks="FD", add.ml.beta=FALSE,
                              verbose=1, ...){
  stopifnot(xor(! is.null(X), ! is.null(afs)),
            is.logical(add.ml.beta))

  if(! is.null(X)){
    stopIfNotValidGenosDose(X, check.hasColNames=FALSE,
                            check.hasRowNames=FALSE,
                            check.noNA=FALSE)
    N <- nrow(X)
    P <- ncol(X)
    if(verbose > 0){
      msg <- paste0(P, " SNPs and ", N, " genotypes")
      write(msg, stdout())
    }
    afs <- estimSnpAf(X)
  } else # is.null(X)
    stopifnot(is.numeric(afs),
              all(afs >= 0, afs <= 1))

  if(is.null(main))
    main <- paste0("AFs of ", length(afs), " markers")

  tmp <- graphics::hist(x=afs, xlab="Allele frequency",
                        ylab="Number of markers",
                        main=main, xlim=xlim, col=col, border=border,
                        las=las, breaks=breaks, ...)

  fit <- NULL
  if(add.ml.beta){
    fn.beta <- function(theta, x){
      sum(-stats::dbeta(x, theta[1], theta[2], log=TRUE))
    }
    fit <- stats::optim(par=c(0.3,0.3), fn=fn.beta, x=afs, method="L-BFGS-B",
                 lower=c(0.0001,0.0001), upper=c(1,1))
    ## http://stackoverflow.com/a/20078645/597069
    x <- seq(min(afs), max(afs), length=100)
    y <- stats::dbeta(x, fit$par[1], fit$par[2])
    y <- y * diff(tmp$mids[1:2]) * length(afs)
    graphics::lines(x, y, col="red")
    graphics::legend("topright", legend=paste0("Beta(", format(fit$par[1], digits=2),
                                     ", ", format(fit$par[2], digits=2),
                                     ") fitted via\nmaximum likelihood"),
           col="red", lty=1, bty="n")
  }

  invisible(list(out.hist=tmp, out.optim=fit))
}

##' Minor allele frequencies
##'
##' Plot the histogram of the minor allele frequency per SNP.
##' @param X matrix of bi-allelic SNP genotypes encoded, for each SNP, in number of copies of its minor allele, i.e. as allele doses in [0,2], with genotypes in rows and SNPs in columns; if X is not NULL, minor allele frequencies will be estimated via \code{\link{estimSnpAf}}
##' @param mafs vector of minor allele frequencies (not used if X is not NULL)
##' @param main string for the main title
##' @param xlim limits of the x-axis
##' @param col color for the bars
##' @param border color for the border of the bars
##' @param las see ?par
##' @param breaks see ?hist
##' @param add.ml.beta if TRUE, add a curve corresponding to a Beta distribution with parameters inferred by maximum likelihood
##' @param verbose verbosity level (0/1)
##' @param ... arguments to be passed to hist()
##' @return invisible list with the output of \code{hist} and \code{optim} (if \code{add.ml.beta=TRUE})
##' @author Timothee Flutre
##' @export
plotHistMinAllelFreq <- function(X=NULL, mafs=NULL, main=NULL, xlim=c(0,0.5),
                                 col="grey", border="white", las=1,
                                 breaks="FD", add.ml.beta=FALSE,
                                 verbose=1, ...){
  stopifnot(xor(! is.null(X), ! is.null(mafs)),
            is.logical(add.ml.beta))

  if(! is.null(X)){
    stopIfNotValidGenosDose(X, check.hasColNames=FALSE,
                            check.hasRowNames=FALSE,
                            check.noNA=FALSE)
    N <- nrow(X)
    P <- ncol(X)
    if(verbose > 0){
      msg <- paste0(P, " SNPs and ", N, " genotypes")
      write(msg, stdout())
    }
    mafs <- estimSnpAf(X)
    if(any(mafs > 0.5, na.rm=TRUE)){
      msg <- paste0(sum(mafs > 0.5), " SNPs are not encoded",
                    " in terms of their minor allele",
                    "\nuse recodeGenosMinorSnpAllele()")
      stop(msg)
    }
  } else{ # is.null(X)
    stopifnot(is.vector(mafs),
              is.numeric(mafs),
              all(mafs[! is.na(mafs)] <= 0.5))
  }

  if(is.null(main))
    main <- paste0("MAFs of ", length(mafs), " markers")

  tmp <- graphics::hist(x=mafs, xlab="Minor allele frequency",
                        ylab="Number of markers",
                        main=main, xlim=xlim, col=col, border=border,
                        las=las, breaks=breaks, ...)

  fit <- NULL
  if(add.ml.beta){
    fn.beta <- function(theta, x){
      sum(-stats::dbeta(x, theta[1], theta[2], log=TRUE))
    }
    fit <- stats::optim(par=c(0.3,0.3), fn=fn.beta, x=mafs, method="L-BFGS-B",
                 lower=c(0.0001,0.0001), upper=c(1,1))
    ## http://stackoverflow.com/a/20078645/597069
    x <- seq(min(mafs), max(mafs), length=100)
    y <- stats::dbeta(x, fit$par[1], fit$par[2])
    y <- y * diff(tmp$mids[1:2]) * length(mafs)
    graphics::lines(x, y, col="red")
    graphics::legend("topright", legend=paste0("Beta(", format(fit$par[1], digits=2),
                                     ", ", format(fit$par[2], digits=2),
                                     ") fitted via\nmaximum likelihood"),
           col="red", lty=1, bty="n")
  }

  invisible(list(out.hist=tmp, out.optim=fit))
}

##' Minor allele frequencies
##'
##' Discard the SNPs with a minor allele frequency below the given threshold, as well as monomorphic SNPs.
##' @param X matrix of bi-allelic SNP genotypes encoded in allele doses in [0,2], with genotypes in rows and SNPs in columns; missing values should be encoded as NA
##' @param mafs vector of minor allele frequencies; if NULL, will be estimated with \code{\link{estimSnpMaf}}
##' @param thresh threshold on minor allele frequencies strictly below which SNPs are ignored
##' @param rmv.mono if TRUE, monomorphic SNPs are also removed
##' @param verbose verbosity level (0/1)
##' @return matrix similar to X but possibly with less columns
##' @author Timothee Flutre
##' @export
discardSnpsLowMaf <- function(X, mafs=NULL, thresh=0.01, rmv.mono=FALSE,
                              verbose=1){
  stopIfNotValidGenosDose(X, check.hasColNames=FALSE, check.hasRowNames=FALSE,
                          check.noNA=FALSE)
  if(! is.null(mafs))
    stopifnot(is.vector(mafs),
              is.numeric(mafs),
              ! is.null(names(mafs)),
              all(mafs <= 0.5),
              all(names(mafs) %in% colnames(X)),
              all(colnames(X) %in% names(mafs)))
  stopifnot(is.numeric(thresh),
            thresh >= 0,
            thresh <= 0.5,
            is.logical(rmv.mono))

  if(is.null(mafs)){
    mafs <- estimSnpMaf(X=X)
  } else
    mafs <- mafs[colnames(X)] # put in same order

  snps.low <- mafs < thresh
  if(any(snps.low)){
    if(verbose > 0){
      msg <- paste0("discard ", sum(snps.low), " SNPs with MAF < ", thresh)
      write(msg, stdout())
    }
    idx.rm <- which(snps.low)
    X <- X[, -idx.rm, drop=FALSE]
  }

  if(rmv.mono){
    nb.alleles <- apply(X, 2, unique)
    nb.alleles <- sapply(nb.alleles, length)
    snps.mono <- nb.alleles == 1
    if(any(snps.mono)){
      if(verbose > 0){
        msg <- paste0("discard ", sum(snps.mono), " monomorphic SNPs")
        write(msg, stdout())
      }
      idx.rm <- which(snps.mono)
      X <- X[, -idx.rm, drop=FALSE]
    }
  }

  return(X)
}

##' Genetic relatedness
##'
##' Reformat the output of \code{related::coancestry} (http://frasierlab.wordpress.com/software/) into a matrix.
##' By default, off-diagonal elements  correspond to coancestry coefficients between two genotypes, and each diagonal element corresponds to (1 + f) / 2 where f corresponds to the inbreeding coefficient of the given genotype.
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

##' Convert haplotypes
##'
##' Convert a list of haplotypes into a matrix.
##' @param haplos list of matrices (one per chromosome, with genotypes in rows and SNPs in columns)
##' @return matrix
##' @author Timothee Flutre
##' @export
haplosList2Matrix <- function(haplos){
  stopifnot(is.list(haplos),
            all(sapply(haplos, methods::is, "matrix")))

  nb.chroms <- length(haplos)
  nb.haplos <- nrow(haplos[[1]]) # 2 x nb of genotypes
  nb.snps <- sapply(haplos, ncol)
  P <- sum(nb.snps) # nb of SNPs

  H <- matrix(data=NA, nrow=nb.haplos, ncol=P)
  rownames(H) <- rownames(haplos[[1]])
  tmp <- lapply(haplos, colnames)
  names(tmp) <- NULL
  colnames(H) <- do.call(c, tmp)

  H[, 1:nb.snps[1]] <- haplos[[1]]
  if(nb.chroms > 1)
    for(chr in 2:nb.chroms)
      H[, (cumsum(nb.snps)[chr-1]+1):(cumsum(nb.snps)[chr])] <- haplos[[chr]]

  return(H)
}

##' Haplotypes
##'
##' Permute alleles in haplotypes.
##' In the \href{https://cran.r-project.org/package=scrm}{scrm} package, haplotypes are made of 0's and 1's, the former indicating ancestral alleles and the latter derived alleles.
##' However, when one wants to simulate realistic data, one often converts haplotypes into genotypes encoded as allele dose.
##' At this step, one may not want to always count the number of copies of the derived allele.
##' It hence is useful to permute alleles beforehand.
##' @param haplos list of haplotypes returned by \code{scrm}, each component of which corresponds to a matrix with haplotypes in rows and SNP in columns
##' @param snps.toperm vector of SNP identifiers corresponding to the SNPs to which allele permutation will be performed (column names of \code{haplos})
##' @param verbose verbosity level (0/1)
##' @return list of haplotypes
##' @author Timothee Flutre
##' @seealso \code{\link{simulCoalescent}}, \code{\link{segSites2allDoses}}, \code{\link{haplosAlleles2num}}
##' @export
permuteAllelesInHaplosNum <- function(haplos, snps.toperm, verbose=0){
  stopIfNotValidHaplos(haplos=haplos, check.hasColNames=FALSE, check.noNA=TRUE)
  stopifnot(is.vector(snps.toperm),
            all(snps.toperm %in% do.call(c, lapply(haplos, colnames))))

  if(verbose > 0){
    msg <- "permute alleles in haplotypes ..."
    write(msg, stdout())
  }

  out <- lapply(haplos, function(mat){
    if(any(colnames(mat) %in% snps.toperm)){
      col.idx <- which(colnames(mat) %in% snps.toperm)
      row.idx.ancestral <- which(mat[, col.idx] == 0)
      row.idx.derived <- which(mat[, col.idx] == 1)
      mat[, col.idx][row.idx.ancestral] <- 1
      mat[, col.idx][row.idx.derived] <- 0
    }
    return(mat)
  })

  return(out)
}

##' Genotypes
##'
##' Permute alleles in genotypes once alleles have been permuted in the corresponding haplotypes.
##' @param X matrix of bi-allelic SNP genotypes encoded, for each SNP, in number of copies of its second allele, i.e. as allele doses in {0,1,2}, with genotypes in rows and SNPs in columns; the "second" allele is arbitrary, it can be the minor or the major allele
##' @param snps.toperm vector of SNP identifiers corresponding to the SNPs to which allele permutation will be performed (column names of \code{X})
##' @param verbose verbosity level (0/1)
##' @return matrix of genotypes
##' @author Timothee Flutre
##' @seealso \code{\link{permuteAllelesInHaplosNum}}
##' @export
permuteAllelesInGenosDose <- function(X, snps.toperm, verbose=0){
  stopIfNotValidGenosDose(X=X, check.hasColNames=TRUE, check.noNA=FALSE)
  stopifnot(is.vector(snps.toperm),
            all(snps.toperm %in% colnames(X)))

  if(verbose > 0){
    msg <- "permute alleles in genotypes ..."
    write(msg, stdout())
  }

  out <- X

  if(any(colnames(X) %in% snps.toperm)){
    col.idx <- which(colnames(X) %in% snps.toperm)
    row.idx.0 <- which(X[, col.idx] == 0)
    row.idx.2 <- which(X[, col.idx] == 2)
    out[, col.idx][row.idx.0] <- 2
    out[, col.idx][row.idx.2] <- 0
  }

  return(out)
}

##' Site frequency spectrum
##'
##' Convert the SFS of independent replicates into a matrix of allele doses.
##' @param seg.sites list of haplotypes returned by \code{scrm::scrm}, each component of which corresponds to a matrix with haplotypes in rows and SNP in columns
##' @param ind.ids vector with the identifiers of the genotypes
##' @param snp.ids vector with the identifiers of the SNPs (if NULL, the SNP identifiers from seg.sites will be used if they aren't NULL, too)
##' @param verbose verbosity level (0/1)
##' @return matrix with diploid genotypes in rows and SNPs in columns
##' @author Timothee Flutre
##' @export
segSites2allDoses <- function(seg.sites, ind.ids=NULL, snp.ids=NULL,
                              verbose=0){
  stopIfNotValidHaplos(haplos=seg.sites, check.hasColNames=FALSE,
                       check.noNA=TRUE)
  if(! is.null(ind.ids))
    stopifnot(is.vector(ind.ids),
              is.character(ind.ids),
              length(ind.ids) == nrow(seg.sites[[1]]) / 2)
  if(! is.null(snp.ids))
    stopifnot(is.vector(snp.ids),
              is.character(snp.ids),
              length(snp.ids) == sum(sapply(seg.sites, ncol)))

  if(verbose > 0){
    msg <- "convert haplotypes into genotypes encoded as allele dose ..."
    write(msg, stdout())
  }

  nb.inds <- nrow(seg.sites[[1]]) / 2 # nb of diploid genotypes
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
  for(x in seq_along(seg.sites)){ # for loop over "chromosomes"
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
##' Make a data.frame of SNP coordinates (1-based) from the SFS of independent replicates.
##' @param seg.sites list of haplotypes returned by \code{scrm}, each component of which corresponds to a matrix with haplotypes in rows and SNP in columns
##' @param snp.ids vector of identifiers (one per SNP)
##' @param chrom.len chromosome length (same for all)
##' @param prefix character string
##' @param verbose verbosity level (0/1)
##' @return data.frame with SNPs in rows and 2 columns (chr, pos)
##' @author Timothee Flutre
##' @export
segSites2snpCoords <- function(seg.sites, snp.ids, chrom.len, prefix="chr",
                               verbose=1){
  stopIfNotValidHaplos(haplos=seg.sites, check.hasColNames=TRUE,
                       check.noNA=FALSE)

  if(verbose > 0){
    msg <- "make a data.frame with SNP coordinates ..."
    write(msg, stdout())
  }

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
  snp.coords$pos <- ceiling(snp.coords$pos)
  for(chr in unique(snp.coords$chr)){
    ## message(chr) # when debugging
    idx <- which(snp.coords$chr == chr)
    pos <- snp.coords$pos[idx]
    if(any(pos == 0))
      pos <- pos + 1
    if(max(pos) > chrom.len){
      pos[which.max(pos)] <- chrom.len
    }
    while(anyDuplicated(pos)){
      ## message("dedup") # when debugging
      i.dup <- which(duplicated(pos))
      pos[i.dup] <- pos[i.dup] + 1
    }
    stopifnot(min(pos) >= 1,
              max(pos) <= chrom.len,
              anyDuplicated(pos) == 0)
    snp.coords$pos[idx] <- pos
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
##' @return matrix with genotypes in rows and SNPs in columns
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

  if(is.null(snp.ids))
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

##' Genotypes
##'
##' Simulate SNP genotypes as allele dose additively encoded, i.e. 0,1,2, using correlated allele frequencies to mimick genetic structure.
##' @param nb.genos number of genotypes (i.e. individuals) per population
##' @param nb.snps number of SNPs
##' @param div.pops matrix of divergence among populations, with a diagonal of 1's; the closer off-diagonal values are from 1; the weaker the divergence; the further, the stronger
##' @param geno.ids vector of genotype identifiers (if NULL, will be "geno001", etc)
##' @param snp.ids vector of SNP identifiers (if NULL, will be "snp001", etc)
##' @return matrix with genotypes in rows and SNPs in columns
##' @author Andres Legarra [aut], Timothee Flutre [ctb]
##' @export
##' @examples
##' \dontrun{
##' ## weak divergences among populations:
##' weak.div.pops <- diag(3)
##' weak.div.pops[upper.tri(weak.div.pops)] <- 0.9
##' weak.div.pops[lower.tri(weak.div.pops)] <- weak.div.pops[upper.tri(weak.div.pops)]
##' weak.div.pops
##'
##' ## strong divergences among populations:
##' strong.div.pops <- diag(3)
##' strong.div.pops[upper.tri(strong.div.pops)] <- 0.5
##' strong.div.pops[lower.tri(strong.div.pops)] <- strong.div.pops[upper.tri(strong.div.pops)]
##' strong.div.pops
##'
##' X <- simulGenosDoseStruct(div.pops=weak.div.pops)
##' A <- estimGenRel(X)
##' imageWithScale(A, "Weak divergence")
##'
##' X <- simulGenosDoseStruct(div.pops=strong.div.pops)
##' A <- estimGenRel(X)
##' imageWithScale(A, "Strong divergence")
##' }
simulGenosDoseStruct <- function(nb.genos=c(100, 120, 80),
                                 nb.snps=1000,
                                 div.pops=diag(3) * 0.5 + 0.5,
                                 geno.ids=NULL, snp.ids=NULL){
  stopifnot(requireNamespace("MASS", quietly=TRUE),
            length(nb.genos) > 1,
            all(nb.genos %% 1 == 0),
            length(nb.genos) == nrow(div.pops),
            ncol(div.pops) == nrow(div.pops),
            all(sapply(diag(div.pops), all.equal, 1)),
            all(c(div.pops) >= 0),
            all(c(div.pops) <= 1))
  if(! is.null(geno.ids))
    stopifnot(length(geno.ids) == sum(nb.genos))
  if(! is.null(snp.ids))
    stopifnot(length(snp.ids) == nb.snps)

  nb.pops <- length(nb.genos)
  idx.pops <- matrix(nrow=nb.pops, ncol=2)
  idx.pops[1,] <- c(1, nb.genos[1])
  for(j in 2:nb.pops)
    idx.pops[j,] <- c(sum(nb.genos[1:(j-1)]) + 1,
                      sum(nb.genos[1:(j-1)]) + nb.genos[j])

  ## draw allele frequencies from "not-too-different populations" with correlated allele frequencies
  all.freqs <- MASS::mvrnorm(n=nb.snps, mu=rep(0,nb.pops), Sigma=div.pops)
  all.freqs <- qbeta(pnorm(all.freqs), 2, 2)

  ## sample SNP genotypes and assign {0,1,2} coding
  X <- lapply(1:nb.snps, function(i){
    tmp <- rep(NA, sum(nb.genos))
    for(j in 1:nb.pops){
      idx <- idx.pops[j,1]:idx.pops[j,2]
      tmp[idx] <- sample(0:1, nb.genos[j], replace=TRUE,
                         prob=c(1 - all.freqs[i,j], all.freqs[i,j])) +
        sample(0:1, nb.genos[j], replace=TRUE,
               prob=c(1 - all.freqs[i,j], all.freqs[i,j]))
    }
    tmp
  })
  X <- do.call(cbind, X)

  if(is.null(geno.ids))
    geno.ids <- sprintf(fmt=paste0("geno%0", floor(log10(nb.genos))+1, "i"),
                        1:sum(nb.genos))
  rownames(X) <- geno.ids

  if(is.null(snp.ids))
    snp.ids <- sprintf(fmt=paste0("snp%0", floor(log10(nb.snps))+1, "i"),
                       1:nb.snps)
  colnames(X) <- snp.ids

  return(X)
}

##' SNP alleles
##'
##' Simulate alleles for bi-allelic SNPs.
##' @param nb.snps number of SNPs
##' @param snp.ids SNP identifiers (if \code{NULL}, will be generated by default)
##' @param colnames names of the columns for the output data frame
##' @param verbose verbosity level (0/1)
##' @return data frame with SNPs in rows
##' @author Timothee Flutre
##' @seealso \code{\link{simulCoalescent}}, \code{\link{simulGenosDose}}
##' @export
simulRefAltSnpAlleles <- function(nb.snps=NULL, snp.ids=NULL,
                                  colnames=c("ref", "alt"), verbose=1){
  stopifnot(any(! is.null(nb.snps), ! is.null(snp.ids)),
            length(colnames) == 2)
  if(! is.null(nb.snps))
    stopifnot(nb.snps > 0)
  if(! is.null(snp.ids))
    stopifnot(is.character(snp.ids))
  if(all(! is.null(nb.snps), ! is.null(snp.ids)))
    stopifnot(nb.snps == length(snp.ids))

  if(is.null(nb.snps))
    nb.snps <- length(snp.ids)
  if(is.null(snp.ids))
    snp.ids <- sprintf(fmt=paste0("snp%0", floor(log10(nb.snps))+1, "i"),
                       1:nb.snps)

  if(verbose > 0){
    msg <- paste0("simulate a data frame of ", nb.snps, " SNP alleles ...")
    write(msg, stdout())
  }

  tmp <- replicate(nb.snps, sample(x=c("A","T","G","C"), size=2,
                                   replace=FALSE))
  stopifnot(all(tmp[1,] != tmp[2,]))

  alleles <- as.data.frame(t(tmp), row.names=snp.ids,
                           stringsAsFactors=FALSE)
  colnames(alleles) <- colnames

  return(alleles)
}

##' Coalescent with recombination
##'
##' Simulate haplotypes according to an approximation to the coalescent with recombination named the Sequential Coalescent with Recombination Model. Requires the scrm package (Staab et al, 2014).
##' @param nb.inds diploids (thus nb of haplotypes is 2 * nb.inds)
##' @param ind.ids vector of identifiers (one per genotype)
##' @param nb.reps number of independent loci that will be produced (could be seen as distinct chromosomes)
##' @param pop.mut.rate theta = 4 N0 mu
##' @param pop.recomb.rate rho = 4 N0 r
##' @param chrom.len in bp
##' @param other character vector of length 1 with other parameters to the simulator (e.g. time-specific parameters such as "-G 6.93 -eG 0.2 0.0 -eN 0.3 0.5")
##' @param nb.pops number of populations (\code{\link{kmeans}} will then be used to pair haplotypes into diploid genotypes)
##' @param mig.rate migration rate = 4 N0 m (will be symmetric)
##' @param get.trees get gene genealogies in the Newick format
##' @param get.tmrca get time to most recent common ancestor and local tree lengths
##' @param get.alleles get fake alleles sampled in {A,T,G,C}
##' @param permute.alleles if TRUE, the reference alleles are randomly chosen between ancestral and derived alleles
##' @param verbose verbosity level (0/1/2)
##' @return list with haplotypes (list), genotypes as allele doses (matrix) and SNP coordinates (data.frame)
##' @author Timothee Flutre
##' @seealso \code{\link{segSites2snpCoords}}, \code{\link{permuteAllelesInHaplosNum}}, \code{\link{segSites2allDoses}}, \code{\link{simulRefAltSnpAlleles}}, \code{\link{makeCrosses}}
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
simulCoalescent <- function(nb.inds=500,
                            ind.ids=NULL,
                            nb.reps=10,
                            pop.mut.rate=40,
                            pop.recomb.rate=40,
                            chrom.len=5*10^5,
                            other=NULL,
                            nb.pops=1,
                            mig.rate=5,
                            get.trees=FALSE,
                            get.tmrca=FALSE,
                            get.alleles=FALSE,
                            permute.alleles=TRUE,
                            verbose=1){
  requireNamespace("scrm", quietly=TRUE)
  stopifnot(nb.inds > nb.pops,
            is.logical(permute.alleles))
  if(! is.null(other))
    stopifnot(is.character(other),
              length(other) == 1)

  out <- list()

  if(is.null(ind.ids))
    ind.ids <- sprintf(fmt=paste0("ind%0", floor(log10(nb.inds))+1, "i"),
                       1:nb.inds)

  if(verbose > 0){
    msg <- "simulate according to the SCRM ..."
    write(msg, stdout())
  }
  nb.samples <- nb.inds * 2 # e.g. 2 chr1 in ind1, 2 chr1 in ind2, etc
  cmd <- paste0(nb.samples, " ", nb.reps)
  cmd <- paste0(cmd, " -t ", pop.mut.rate)
  cmd <- paste0(cmd, " -r ", pop.recomb.rate, " ", chrom.len)
  if(nb.pops > 1){
    cmd <- paste0(cmd, " -I ", nb.pops)
    nb.inds.per.pop <- rep(0, nb.pops)
    for(p in 1:(nb.pops-1)){
      nb.inds.per.pop[p] <- floor(nb.inds / nb.pops)
      cmd <- paste0(cmd, " ", 2 * nb.inds.per.pop[p])
    }
    nb.inds.per.pop[nb.pops] <- nb.inds - sum(nb.inds.per.pop)
    cmd <- paste0(cmd, " ", 2 * nb.inds.per.pop[nb.pops])
    ## cmd <- paste0(cmd, " ", mig.rate)
    if(is.null(other)){
      cmd <- paste0(cmd, " ", mig.rate)
    } else if(! grepl("-m|-ma", other))
      cmd <- paste0(cmd, " ", mig.rate)
  }
  if(! is.null(other))
    cmd <- paste0(cmd, " ", other)
  if(get.trees)
    cmd <- paste0(cmd, " -T")
  if(get.tmrca)
    cmd <- paste0(cmd, " -L")
  cmd <- paste0(cmd, " -SC abs") # absolute seq positions in bp
  cmd <- paste0(cmd, " -oSFS") # print site freq spectrum, requires -t
  if(verbose > 0){
    msg <- paste0("scrm ", cmd)
    write(msg, stdout())
  }
  out$cmd <- cmd
  sum.stats <- scrm::scrm(cmd)
  if(verbose > 1)
    print(utils::str(sum.stats))

  prefix <- "chr"
  names(sum.stats$seg_sites) <- paste0(prefix, 1:nb.reps)
  nb.snps.per.chr <- sapply(sum.stats$seg_sites, ncol)
  nb.snps <- sum(nb.snps.per.chr)
  if(verbose > 0){
    msg <- paste0("nb of SNPs: ", nb.snps)
    write(msg, stdout())
    print(sapply(sum.stats$seg_sites, ncol))
  }

  snp.ids <- sprintf(fmt=paste0("snp%0", floor(log10(nb.snps))+1, "i"),
                     1:nb.snps)
  snp.coords <- segSites2snpCoords(seg.sites=sum.stats$seg_sites,
                                   snp.ids=snp.ids, chrom.len=chrom.len,
                                   prefix=prefix, verbose=verbose)
  out[["snp.coords"]] <- snp.coords

  if(verbose > 0){
    msg <- "randomize haplotypes to make diploid genotypes ..."
    write(msg, stdout())
  }
  out[["haplos"]] <- list()
  idx <- sample.int(nb.samples)
  if(nb.pops > 1){
    H <- haplosList2Matrix(sum.stats$seg_sites)
    kmH <- stats::kmeans(x=H, centers=nb.pops)
    ## table(kmH$cluster) # to debug
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

  X <- segSites2allDoses(seg.sites=out$haplos, ind.ids=ind.ids,
                         snp.ids=snp.ids, verbose=verbose)
  out[["genos"]] <- X

  if(permute.alleles){
    mafs <- estimSnpMaf(X=out$genos)
    snps.toperm <- names(mafs)[sample.int(n=nb.snps, size=floor(0.5*nb.snps),
                                          replace=FALSE, prob=1 - mafs)]
    out$haplos <- permuteAllelesInHaplosNum(haplos=out$haplos,
                                            snps.toperm=snps.toperm,
                                            verbose=verbose)
    out$genos <- permuteAllelesInGenosDose(X=out$genos,
                                           snps.toperm=snps.toperm,
                                           verbose=verbose)
  }

  if(get.alleles){
    alleles <- simulRefAltSnpAlleles(snp.ids=snp.ids,
                                     colnames=c("first", "second"),
                                     verbose=verbose)
    out[["alleles"]] <- alleles
  }

  if(get.trees)
    out[["trees"]] <- sum.stats$trees
  if(get.tmrca)
    out[["tmrca"]] <- sum.stats$tmrca

  return(out)
}

##' Genomes
##'
##' From a set of genomes simulated together, split them into a training (for estimation) and testing (for prediction) sets.
##' @param genomes list, e.g. returned by \code{\link{simulCoalescent}}
##' @param nb.inds.pred number of genotypes for whom prediction will be performed
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
##' @param haplos matrix of haplotypes, with haplotypes (2 per genotype) in rows and sites in columns
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
  if(nrow(haplos) <= 10){
    nb.intervals <- 2 * nrow(haplos) - 2
    at <- seq(from=0, to=1, length.out=nb.intervals + 1)
    graphics::axis(side=2, at=rev(at[seq(from=1, to=length(at), by=2)]),
                   labels=rownames(haplos),
                   las=1, padj=0)
  } else
    graphics::axis(side=2, at=rev(seq(1, nrow(haplos), length.out=10) / nrow(haplos)),
                   labels=rownames(haplos)[seq(1, nrow(haplos), length.out=10)],
                   las=1, padj=0)

  on.exit(graphics::par(opar))
}

##' Genotype names
##'
##' Return the identifiers of all genotypes
##' @param haplos list of matrices (one per chromosome, with genotypes in rows)
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

##' Haplotypes of an genotype
##'
##' Retrieve both haplotypes for all chromosomes of a given genotype.
##' @param haplos list with one matrix per chromosome, with haplotype in rows and SNPs in columns, as returned by \code{\link{simulCoalescent}}
##' @param ind.name identifier of the genotype to retrieve
##' @return list similar to "haplos"
##' @author Timothee Flutre
##' @export
getHaplosInd <- function(haplos, ind.name){
  stopifnot(is.list(haplos),
            is.character(ind.name))

  ## for each chromosome, retrieve both haplotypes of the given genotype
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

##' Haplotypes of several genotypes
##'
##' Retrieve both haplotypes for all chromosomes of certain genotypes.
##' @param haplos list with one matrix per chromosome, with haplotype in rows and SNPs in columns, as returned by \code{\link{simulCoalescent}}
##' @param ind.names identifier of the genotypes to retrieve
##' @return list similar to "haplos"
##' @author Timothee Flutre
##' @export
getHaplosInds <- function(haplos, ind.names){
  stopifnot(is.list(haplos),
            all(sapply(haplos, methods::is, "matrix")),
            length(unique(sapply(haplos, nrow))) == 1, # same nb of genotypes
            all(sapply(haplos, function(x){! is.null(rownames(x))})), # has row names
            is.character(ind.names))

  ## for each chromosome, retrieve both haplotypes of the given genotypes
  haplos.inds <- lapply(haplos, function(haplos.chr){
    idx <- do.call(c, lapply(ind.names, grep, rownames(haplos.chr)))
    if(length(idx) == 0){
      msg <- paste0("can't find haplotypes for the given genotypes")
      stop(msg)
    } else{
      haplos.chr[idx, , drop=FALSE]
    }
  })

  return(haplos.inds)
}

##' Gamete
##'
##' Make a gamete for a given chromosome by recombining two parental haplotypes.
##' @param haplos.par.chr matrix containing both haplotypes of a parent for a given chromosome (must have dimnames, such as "ind37_h1" in rows and "snp00265" in columns)
##' @param loc.crossovers positions of the crossing overs in terms of SNP indices and not physical coordinates (the coordinate of the first nucleotide is assumed to be 1, and crossing-overs are assumed to occur right after the given positions)
##' @param start.haplo identifier of the haplotype with which to start the gamete (1/2)
##' @return vector
##' @author Timothee Flutre
##' @seealso \code{\link{makeGameteSingleInd}}, \code{\link{drawLocCrossovers}}
##' @export
makeGameteSingleIndSingleChrom <- function(haplos.par.chr, loc.crossovers,
                                           start.haplo=1){
  stopifnot(is.matrix(haplos.par.chr),
            ! is.null(dimnames(haplos.par.chr)),
            is.vector(loc.crossovers),
            min(loc.crossovers) >= 0,
            max(loc.crossovers) < ncol(haplos.par.chr),
            start.haplo %in% c(1, 2))

  nb.snps <- ncol(haplos.par.chr)
  loc.crossovers <- sort(loc.crossovers)

  ind.name <- strsplit(rownames(haplos.par.chr)[1],
                       "_")[[1]]
  ind.name <- paste(ind.name[1:(length(ind.name)-1)], collapse="_")
  gam.chr <- matrix(rep(NA, nb.snps), nrow=1,
                    dimnames=list(ind.name, colnames(haplos.par.chr)))

  ## fill the gamete from the left to the right
  if(loc.crossovers[1] != 0)
    loc.crossovers <- c(0, loc.crossovers)
  for(l in 2:length(loc.crossovers))
    gam.chr[(loc.crossovers[l-1]+1):loc.crossovers[l]] <-
      haplos.par.chr[ifelse((l + start.haplo) %% 2 == 1, 1, 2),
      (loc.crossovers[l-1]+1):loc.crossovers[l]]

  ## fill the gamete until the end
  l <- length(loc.crossovers)
  gam.chr[(loc.crossovers[l]+1):nb.snps] <-
    haplos.par.chr[ifelse((l+1 + start.haplo) %% 2 == 1, 1, 2),
    (loc.crossovers[l]+1):nb.snps]

  stopifnot(sum(is.na(gam.chr)) == 0)

  return(gam.chr)
}

##' Gamete
##'
##' Make a gamete for all chromosomes by recombining two parental haplotypes.
##' @param haplos.par list of matrices containing both haplotypes of a parent for each chromosome
##' @param loc.crossovers list of vectors with positions of the crossing overs (the coordinate of the first nucleotide is assumed to be 1, and crossing-overs are assumed to occur right after the given positions)
##' @param start.haplos list of identifiers (1/2) of the haplotype with which to start the gamete for each chromosome (if NULL, will be sampled)
##' @return list of vectors (one per chromosome)
##' @author Timothee Flutre
##' @seealso \code{\link{makeCross}}, \code{\link{fecundation}}, \code{\link{makeGameteSingleIndSingleChrom}}
##' @export
makeGameteSingleInd <- function(haplos.par, loc.crossovers,
                                start.haplos=as.list(rep(1, length(haplos.par)))){
  stopifnot(is.list(haplos.par),
            is.list(loc.crossovers),
            length(loc.crossovers) == length(haplos.par))
  if(! is.null(start.haplos))
    stopifnot(is.list(start.haplos),
              length(start.haplos) == length(haplos.par),
              all(do.call(c, start.haplos) %in% c(1,2)))

  nb.chrs <- length(loc.crossovers)
  if(is.null(start.haplos))
    start.haplos <- as.list(sample.int(n=2, size=nb.chrs, replace=TRUE))

  gam <- Map(makeGameteSingleIndSingleChrom, haplos.par, loc.crossovers,
             start.haplos)

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
##' Make a cross. If two different genotypes are given, then a fecundation is made with one gamete from each genotype. If the same genotype is given twice, an autofecondation is made with two different gametes from this genotype. If a single genotype is given, a haplodiploidization is made with a single gamete from this genotype.
##' @param haplos.par1 list of matrices (one per chromosome) for the first parent
##' @param loc.crossovers.par1 list of vectors (one per chromosome) for the first parent
##' @param haplos.par2 list of matrices (one per chromosome) for the second parent
##' @param loc.crossovers.par2 list of vectors (one per chromosome) for the second parent
##' @param child.name identifier of the child (if NULL, <parent1>"-x-"<parent2> or <parent1>"-hd")
##' @param howto.start.haplo if 0, haplotypes with which to start the gametes will be chosen at random; but one can also specify 1 (always the first haplotype) or 2 (always the second haplotype)
##' @param verbose verbosity level (0/1)
##' @return list of matrices (one per chromosome) for the child
##' @author Timothee Flutre
##' @seealso \code{\link{makeCrosses}}, \code{\link{drawLocCrossovers}}, \code{\link{makeGameteSingleInd}}, \code{\link{fecundation}}
##' @export
makeCross <- function(haplos.par1,
                      loc.crossovers.par1,
                      haplos.par2=NULL,
                      loc.crossovers.par2=NULL,
                      child.name=NULL,
                      howto.start.haplo=0,
                      verbose=1){
  stopifnot(is.list(haplos.par1),
            is.list(loc.crossovers.par1),
            all(names(haplos.par1) == names(loc.crossovers.par1)),
            howto.start.haplo %in% c(0, 1, 2))
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
  start.haplos <- as.list(rep(howto.start.haplo, length(haplos.par1)))
  if(howto.start.haplo == 0)
    start.haplos <- NULL
  gam.par1 <- makeGameteSingleInd(haplos.par1, loc.crossovers.par1,
                                  start.haplos)
  if(cross.type == "haplodiploidization"){
    if(is.null(child.name))
      child.name <- paste0(name.par1, "-hd")
    haplos.child <- fecundation(gam.par1, gam.par1, child.name)
  } else{
    gam.par2 <- makeGameteSingleInd(haplos.par2, loc.crossovers.par2,
                                    start.haplos)
    if(is.null(child.name))
      child.name <- paste0(name.par1, "-x-", name.par2)
    haplos.child <- fecundation(gam.par1, gam.par2, child.name)
  }

  return(haplos.child)
}

##' Crossing-overs
##'
##' Draw the number and location of crossing-overs per gamete, in terms of SNP indices and not physical coordinates.
##' @param crosses data.frame with three columns, parent1, parent2, child; if parent 1 and 2 are the same, it will be an autofecondation; if parent2 is NA, it will be a haplodiploidization
##' @param nb.snps vector with the nb of SNPs per chromosome, which names are chromosome names
##' @param lambda mean number of crossing-overs (parameter of a Poisson)
##' @param simplistic if TRUE, the value of \code{lambda} is ignored, and a single crossing over per gamete chromosome is assumed, which location is drawn uniformly
##' @param verbose verbosity level (0/1)
##' @return list of lists (one per cross, then one per parent, then one per chromosome) whose names are crosses$child, in the same order
##' @author Timothee Flutre
##' @seealso \code{\link{makeGameteSingleIndSingleChrom}}
##' @examples
##' \dontrun{## make inputs
##' nb.offs <- 100
##' names.offs <- paste0("par1-par2-",
##'                      sprintf(fmt=paste0("%0", floor(log10(nb.offs))+1, "i"),
##'                              1:nb.offs))
##' crosses <- data.frame(parent1=rep("par1", nb.offs),
##'                       parent2=rep("par2", nb.offs),
##'                       child=names.offs,
##'                       stringsAsFactors=FALSE)
##'
##' ## draw locations of crossing-overs
##' loc.crossovers <- drawLocCrossovers(crosses=crosses,
##'                                     nb.snps=setNames(c(2739, 2811),
##'                                                      c("chr1", "chr2")),
##'                                     verbose=1)
##' }
##' @export
drawLocCrossovers <- function(crosses, nb.snps, lambda=2, simplistic=FALSE,
                              verbose=0){
  stopifnot(is.data.frame(crosses),
            ncol(crosses) >= 3,
            all(c("parent1", "parent2", "child") %in% colnames(crosses)),
            sum(is.na(crosses$parent1)) == 0,
            sum(is.na(crosses$child)) == 0,
            is.vector(nb.snps),
            ! is.null(names(nb.snps)),
            is.logical(simplistic))

  loc.crossovers <- list()

  ## draw the number of crossing-overs for each gamete
  parent.names <- c(crosses$parent1, crosses$parent2)
  parent.names <- parent.names[! is.na(parent.names)]
  nb.gametes <- length(parent.names)
  nb.chroms <- length(nb.snps)
  if(simplistic){
    nb.crossovers <- rep(1, nb.gametes * nb.chroms)
  } else
    nb.crossovers <- stats::rpois(nb.gametes * nb.chroms, lambda)
  nb.crossovers[nb.crossovers == 0] <- 1 # at least 1 per chromosome

  ## draw the location of each crossing-over
  drawLocCrossoversOneChr <- function(nb.snps, nb.crossovers){
    sort(sample.int(n=nb.snps, size=nb.crossovers, replace=FALSE, prob=NULL,
                    useHash=FALSE))
  }
  nb.crosses <- nrow(crosses)
  chrom.names <- names(nb.snps)
  cross.idx <- 0
  gam.idx <- 0
  while(cross.idx < nb.crosses){
    cross.idx <- cross.idx + 1
    child <- list()
    gam.idx <- gam.idx + 1
    gam.par1 <- lapply(1:nb.chroms, function(ci){
      drawLocCrossoversOneChr(nb.snps[ci] - 1,
                              nb.crossovers[(gam.idx-1)*nb.chroms + ci])
    })
    names(gam.par1) <- chrom.names
    child[[crosses$parent1[cross.idx]]] <- gam.par1
    if(! is.na(crosses$parent2[cross.idx])){ # (auto)fecondation
      gam.idx <- gam.idx + 1
      gam.par2 <- lapply(1:nb.chroms, function(c){
        drawLocCrossoversOneChr(nb.snps[c] - 1,
                                nb.crossovers[(gam.idx-1)*nb.chroms + c])
      })
      names(gam.par2) <- chrom.names
      child[[crosses$parent2[cross.idx]]] <- gam.par2
    }
    loc.crossovers[[crosses$child[cross.idx]]] <- child
  }

  if(verbose > 0){
    msg <- "distribution of the nb of crossing-overs per chromosome:"
    write(msg, stdout())
    tmp <- table(do.call(c, lapply(loc.crossovers, function(lc.off){
      do.call(c, lapply(lc.off, function(lc.par){
        sapply(lc.par, length)
      }))
    })))
    print(tmp)
  }

  return(loc.crossovers)
}

##' Crosses
##'
##' Make crosses (fecundation, autofecondation, haplodiploidization).
##' @param haplos list of matrices (one per chromosome)
##' @param crosses data.frame with three columns, parent1, parent2, child (no duplicate); if parent 1 and 2 are the same, it will be an autofecondation; if parent2 is NA, it will be a haplodiploidization
##' @param loc.crossovers list of lists (one per cross, then one per parent, then one per chromosome) whose names are crosses$child, in the same order; if NULL, draw many crossing-overs positions at once (as Poisson with parameter 2, assuming all chromosomes have the same length)
##' @param howto.start.haplo if 0, haplotypes with which to start the gametes will be chosen at random; but one can also specify 1 (always the first haplotype) or 2 (always the second haplotype)
##' @param nb.cores the number of cores to use, i.e. at most how many child processes will be run simultaneously (not on Windows)
##' @param verbose verbosity level (0/1)
##' @return list of matrices (one per chromosome) with child haplotypes in rows and SNPs in columns
##' @author Timothee Flutre
##' @seealso \code{\link{makeCross}}, \code{\link{simulCoalescent}}, \code{\link{getHaplosInds}}, \code{\link{drawLocCrossovers}}
##' @examples
##' \dontrun{set.seed(1859)
##' if(require(scrm)){
##'   ## simulate haplotypes
##'   nb.genos <- 2*10^2
##'   nb.chroms <- 10
##'   Ne <- 10^5
##'   chrom.len <- 10^5
##'   mu <- 10^(-8)
##'   c.rec <- 10^(-8)
##'   genomes <- simulCoalescent(nb.inds=nb.genos,
##'                              nb.reps=nb.chroms,
##'                              pop.mut.rate=4 * Ne * mu * chrom.len,
##'                              pop.recomb.rate=4 * Ne * c.rec * chrom.len,
##'                              chrom.len=chrom.len,
##'                              get.alleles=TRUE,
##'                              permute.alleles=TRUE)
##'
##'   ## pick 2 individuals at random as parents
##'   (idx.parents <- sample.int(n=nb.genos, size=2))
##'   genos.parents <- genomes$genos[idx.parents,]
##'   names.parents <- rownames(genos.parents)
##'   haplos.parents <- getHaplosInds(haplos=genomes$haplos,
##'                                  ind.names=names.parents)
##'
##'   ## cross them several times to make offsprings
##'   nb.offs <- 100
##'   names.offs <- paste0(names.parents[1], "-", names.parents[2], "-",
##'                        sprintf(fmt=paste0("%0", floor(log10(nb.offs))+1, "i"),
##'                                1:nb.offs))
##'   crosses <- data.frame(parent1=rep(names.parents[1], nb.offs),
##'                         parent2=rep(names.parents[2], nb.offs),
##'                         child=names.offs,
##'                         stringsAsFactors=FALSE)
##'   loc.crossovers <- drawLocCrossovers(crosses=crosses,
##'                                       nb.snps=sapply(haplos.parents, ncol))
##'   haplos.offs <- makeCrosses(haplos=haplos.parents, crosses=crosses,
##'                              loc.crossovers=loc.crossovers)
##'   genos.offs <- segSites2allDoses(seg.sites=haplos.offs,
##'                                   ind.ids=getIndNamesFromHaplos(haplos.offs),
##'                                   snp.ids=rownames(genomes$snp.coords))
##'
##'   ## look at the first crossing-over in the first offspring
##'   (snps.co <- colnames(haplos.parents$chr1)[loc.crossovers[[1]][[1]]$chr1])
##'   tmp <- subsetDiffHaplosWithinParent(haplos.chr=haplos.parents$chr1[1:2,],
##'                                       snps.tokeep=snps.co)
##'   (idx1 <- which(colnames(tmp) == snps.co[1]))
##'   (snps.tokeep <- colnames(tmp)[(idx1-2):(idx1+2)])
##'   tmp[, snps.tokeep]
##'   haplos.offs$chr1[1:2, snps.tokeep]
##' }
##' }
##' @export
makeCrosses <- function(haplos, crosses, loc.crossovers=NULL,
                        howto.start.haplo=0,
                        nb.cores=1, verbose=1){
  stopifnot(is.list(haplos),
            all(sapply(haplos, methods::is, "matrix")),
            length(unique(sapply(haplos, nrow))) == 1, # same nb of genotypes
            is.data.frame(crosses),
            ncol(crosses) >= 3,
            all(c("parent1", "parent2", "child") %in% colnames(crosses)),
            sum(is.na(crosses$parent1)) == 0,
            sum(is.na(crosses$child)) == 0,
            anyDuplicated(crosses$child) == 0)
  haplos.ind.names <- getIndNamesFromHaplos(haplos)
  idx <- sapply(crosses, is.factor)
  crosses[idx] <- lapply(crosses[idx], as.character)
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
                howto.start.haplo=howto.start.haplo,
                verbose=verbose-1)
    } else #(auto)fecondation
      makeCross(haplos.par1=getHaplosInd(haplos, crosses$parent1[i]),
                loc.crossovers.par1=loc.crossovers[[i]][[crosses$parent1[i]]],
                haplos.par2=getHaplosInd(haplos, crosses$parent2[i]),
                loc.crossovers.par2=loc.crossovers[[i]][[crosses$parent2[i]]],
                child.name=crosses$child[i],
                howto.start.haplo=howto.start.haplo,
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

##' Haplotypes
##'
##' Return a subset of SNPs for both haplotypes within a given parent in order to easily visualize a given crossing-over.
##' @param haplos.chr matrix containing both haplotypes of a given parent, with haplotypes in rows and SNPs in columns
##' @param snps.tokeep vector of SNPs to keep in the subset (typically the ones at which the crossing-overs occurred)
##' @return matrix with both haplotypes in rows and a subset of SNPs in columns
##' @author Timothee Flutre
##' @export
subsetDiffHaplosWithinParent <- function(haplos.chr, snps.tokeep){
  stopifnot(is.matrix(haplos.chr),
            nrow(haplos.chr) == 2,
            all(haplos.chr %in% c(0,1)),
            ! is.null(colnames(haplos.chr)))
  if(! is.null(snps.tokeep))
    stopifnot(is.vector(snps.tokeep),
              all(snps.tokeep %in% colnames(haplos.chr)))

  nb.snps <- ncol(haplos.chr)
  out <- haplos.chr
  idx <- apply(haplos.chr, 2, function(x){
    any(x[1] != x[2])
  })
  if(! is.null(snps.tokeep))
    idx[snps.tokeep] <- TRUE
  out <- out[, idx]

  return(out)
}

.isValidSnpCoords <- function(snp.coords){
  all(is.data.frame(snp.coords),
      ! is.null(rownames(snp.coords)),
      ncol(snp.coords) >= 2,
      "chr" %in% colnames(snp.coords),
      "coord" %in% colnames(snp.coords) | "pos" %in% colnames(snp.coords))
}

##' SNP coordinates from data frame to GRanges
##'
##' Convert a data frame of SNP coordinates to GRanges.
##' @param x data frame of SNP coordinates with at least two columns, "chr" and "coord" (or "pos")
##' @param si output from the "seqinfo" function from the GenomeInfoDb package
##' @return GRanges
##' @author Timothee Flutre
##' @export
snpCoordsDf2Gr <- function(x, si=NULL){
  stopifnot(.isValidSnpCoords(x))

  if("pos" %in% colnames(x))
    colnames(x)[colnames(x) == "pos"] <- "coord"

  out <- df2gr(x=x, seq="chr", start="coord", end="coord", strand=NULL, si=si)

  return(out)
}

##' Distance between SNP pairs
##'
##' For each SNP pair, return the number of "blocks" (i.e. nucleotides) between both SNPs via the \code{\link[GenomicRanges]{distance}} function.
##' @param snp.pairs data.frame with two columns "loc1" and "loc2"
##' @param snp.coords data.frame with SNP identifiers as row names, and two columns, "chr" and "coord" (or "pos")
##' @param nb.cores the number of cores to use
##' @param verbose verbosity level (0/1)
##' @return vector
##' @author Timothee Flutre
##' @seealso \code{\link{estimLd}}, \code{\link{distConsecutiveSnps}}
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
  snp.granges <- snpCoordsDf2Gr(snp.coords)

  if(verbose > 0)
    message("calculate pairwise distances ...")
  dist.loc <- GenomicRanges::distance(x=snp.granges[snp.pairs$loc1],
                                      y=snp.granges[snp.pairs$loc2])

  return(dist.loc)
}

##' SNP genotypes
##'
##' Recode SNP genotypes from additive to dominant.
##' @param X matrix of bi-allelic SNP genotypes encoded in allele doses in {0,1,2}, with genotypes in rows and SNPs in columns; missing values should be encoded as NA
##' @param simplify.imputed if TRUE, imputed genotypes will be considered as homozygotes if less than 0.66 or more than 1.33, and heterozygotes otherwise
##' @return matrix with genotypes in {0,1}
##' @author Timothee Flutre
##' @export
recodeIntoDominant <- function(X, simplify.imputed=FALSE){
  stopIfNotValidGenosDose(X=X, check.hasColNames=FALSE,
                          check.hasRowNames=FALSE, check.noNA=FALSE,
                          check.notImputed=! simplify.imputed)
  X.D <- X
  if(simplify.imputed){
    idx <- indexGenoDoses(X)
    X.D[idx[,"is.1"]] <- 1
    X.D[! idx[,"is.1"]] <- 0
  } else
    X.D[X.D != 1] <- 0

  is.NA <- is.na(X)
  if(any(is.NA))
    X[is.NA] <- NA

  return(X.D)
}

##' Genomic relatedness
##'
##' Estimate genetic relationships between genotypes from their SNP genotypes.
##' Note that this function estimates "relationships" and not "coancestries".
##' See \href{http://dx.doi.org/10.1186/1297-9686-43-27}{Toro et al (2011)}: "for diploid individuals, twice the coancestry coefficient is the additive relationship coefficient, which describes the ratio between the genetic covariance between individuals and the genetic variance of the base population".
##' @param X matrix of bi-allelic SNP genotypes encoded in allele doses in {0,1,2}, with genotypes in rows and SNPs in columns; missing values should be encoded as NA
##' @param afs vector of allele frequencies, corresponding to the alleles whose copies are counted in \code{X} (if NULL, will be calculated with \code{\link{estimSnpAf}})
##' @param thresh threshold on minor allele frequencies below which SNPs are ignored (e.g. 0.01; NULL to skip this step)
##' @param relationships relationship to estimate (additive/dominant/gaussian) where "gaussian" corresponds to the Gaussian kernel from \href{http://dx.doi.org/10.3835/plantgenome2011.08.0024}{Endelman (2011)}
##' @param method \itemize{
##' \item if additive relationships, can be "noia" (see \href{https://doi.org/10.1534/genetics.116.199406}{Vitezica et al, 2017}), "vanraden1" (first method in \href{http://dx.doi.org/10.3168/jds.2007-0980}{VanRaden, 2008}), "toro2011_eq10" (equation 10 using molecular covariance from \href{http://dx.doi.org/10.1186/1297-9686-43-27}{Toro et al, 2011}), "habier" (similar to "vanraden1" but without centering; from \href{http://dx.doi.org/10.1534/genetics.107.081190}{Habier et al, 2007}), "astle-balding" (two times equation 2.2 in \href{http://dx.doi.org/10.1214/09-sts307}{Astle & Balding, 2009}), "yang" (similar to 'astle-balding' but without ignoring sampling error per SNP; from \href{http://dx.doi.org/10.1038/ng.608}{Yang et al, 2010}), "zhou" (centering the genotypes with \code{\link{scale}} and not assuming that rare variants have larger effects; from \href{http://dx.doi.org/10.1371/journal.pgen.1003264}{Zhou et al, 2013}) or "center-std";
##' \item if dominant relationships, can be "noia" (see \href{https://doi.org/10.1534/genetics.116.199406}{Vitezica et al, 2017}), "vitezica" (classical/statistical parametrization from \href{http://dx.doi.org/10.1534/genetics.113.155176}{Vitezica et al, 2013}) or "su" (from \href{http://dx.doi.org/10.1371/journal.pone.0045293}{Su et al, 2012})
##' }
##' @param theta smoothing parameter for "gauss"
##' @param alpha if not NULL, weight to make the matrix invertible as it often is singular because two genotypes are clones or because the centered-coding of SNP genotypes makes the last row predictable from the other ones (see \href{https://doi.org/10.1186/1297-9686-43-25}{Stranden and Christensen, 2011}); otherwise use \code{\link[Matrix]{nearPD}}
##' @param verbose verbosity level (0/1)
##' @return symmetric matrix with the number of SNPs used as an attribute
##' @author Timothee Flutre
##' @examples
##' \dontrun{set.seed(1859)
##' nb.genos <- 200
##' Ne <- 10^4
##' chrom.len <- 10^5
##' mu <- 10^(-8)
##' c <- 10^(-8)
##' genomes <- simulCoalescent(nb.inds=nb.genos,
##'                            pop.mut.rate=4 * Ne * mu * chrom.len,
##'                            pop.recomb.rate=4 * Ne * c * chrom.len,
##'                            chrom.len=chrom.len)
##' X <- genomes$genos
##'
##' A.vr <- estimGenRel(X=X, relationships="additive", method="vanraden1")
##' mean(diag(A.vr)) # should be 1 in a base population at HWE
##' mean(A.vr[upper.tri(A.vr)]) # should be 0 in a base population at HWE
##'
##' D.v <- estimGenRel(X=X, relationships="dominant", method="vitezica")
##' mean(diag(D.v)) # should be 1 in a base population at HWE
##' mean(D.v[upper.tri(D.v)]) # should be 0 in a base population at HWE
##' }
##' @export
estimGenRel <- function(X, afs=NULL, thresh=NULL, relationships="additive",
                        method="vanraden1", theta=0.5, alpha=NULL, verbose=1){
  stopIfNotValidGenosDose(X, check.noNA=FALSE)
  stopifnot(relationships %in% c("additive", "dominant", "gaussian"))
  if(relationships == "additive")
    stopifnot(method %in% c("noia", "vanraden1", "toro2011_eq10", "habier",
                            "astle-balding", "yang",
                            "zhou", "center-std"))
  if(relationships == "dominant")
    stopifnot(method %in% c("noia", "vitezica", "su"))
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
    stopifnot(is.vector(afs),
              is.numeric(afs),
              ! is.null(names(afs)),
              all(colnames(X) %in% names(afs)),
              all(names(afs) %in% colnames(X)))
    afs <- afs[colnames(X)] # put in same order
  }
  if(! is.null(alpha)){
    if(alpha > 0.01){
      msg <- paste0("you want alpha=", alpha, ", but it is recommended to keep it below 0.01")
      warning(msg, immediate.=TRUE)
    }
  }

  gen.rel <- NULL # to be filled and returned

  N <- nrow(X) # nb of genotypes
  P <- ncol(X) # nb of SNPs
  P_init <- P

  if(any(is.na(X))){
    X <- discardMarkersMissGenos(X=X, verbose=verbose)
    P <- ncol(X)
    if(! is.null(afs))
      afs <- afs[colnames(X)]
  }

  if(is.null(afs)){
    if(verbose > 0){
      msg <- "estimate allele frequencies ..."
      write(msg, stdout())
    }
    afs <- estimSnpAf(X=X)
  }

  if(all(relationships == "additive", method == "center-std",
         is.null(thresh))){
    if(any(afs == 0, afs == 1)){
      msg <- paste0("some SNPs have fixed alleles, ",
                    "a threshold of 1% on MAFs will hence be used")
      write(msg, stdout())
      thresh <- 0.01
    }
  }
  if(! is.null(thresh)){
    mafs <- estimSnpMaf(afs=afs)
    if(any(mafs < thresh)){
      X <- discardSnpsLowMaf(X=X, mafs=mafs, thresh=thresh, verbose=verbose)
      P <- ncol(X)
      afs <- afs[colnames(X)]
    }
  }
  if(method == "center-std" & any(afs == 0.5)){
    idx <- which(afs != 0.5)
    X <- X[, idx]
    P <- ncol(X)
    afs <- afs[colnames(X)]
  }

  if(verbose > 0){
    msg <- paste0("estimate relationships with ", P, " SNPs",
                  ifelse(P < P_init, paste0(" out of ", P_init), ""),
                  " ...")
    write(msg, stdout())
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

    } else if(method == "noia"){
      ## compute genotype frequencies at each SNP
      freqG <- t(apply(X, 2, function(x){
        c("A1A1"=sum(x == 0), "A1A2"=sum(x == 1), "A2A2"=sum(x == 2)) / length(x)
      }))

      boundaries <- seq(from=0, to=2, length.out=4)
      is.0 <- (X <= boundaries[2]) # homozygotes for the first allele
      is.1 <- (X > boundaries[2] & X <= boundaries[3]) # heterozygotes
      is.2 <- (X > boundaries[3]) # homozygotes for the second allele

      ## eq 1 of Vitezica et al (2017)
      Z <- matrix(NA, nrow=N, ncol=P, dimnames=dimnames(X))
      for(i in 1:N){
        Z[i, is.0[i,]] <- -(-freqG[is.0[i,],"A1A2"] - 2 * freqG[is.0[i,],"A2A2"])
        Z[i, is.1[i,]] <- -(1 - freqG[is.1[i,],"A1A2"] - 2 * freqG[is.1[i,],"A2A2"])
        Z[i, is.2[i,]] <- -(2 - freqG[is.2[i,],"A1A2"] - 2 * freqG[is.2[i,],"A2A2"])
      }

      ZZt <- tcrossprod(Z, Z)
      gen.rel <- ZZt / (sum(diag(ZZt)) / N)

    } else if(method == "toro2011_eq10"){
      var.afs <- stats::var(afs)
      gen.rel <- 2 * (stats::cov(t(X) / 2) - var.afs) /
        (mean(afs) * mean(1 - afs) - var.afs)

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
      boundaries <- seq(from=0, to=2, length.out=4)
      is.0 <- (X <= boundaries[2]) # homozygotes for the first allele
      is.1 <- (X > boundaries[2] & X <= boundaries[3]) # heterozygotes
      is.2 <- (X > boundaries[3]) # homozygotes for the second allele
      W <- matrix(NA, nrow=N, ncol=P, dimnames=dimnames(X))
      for(i in 1:N){
        W[i, is.0[i,]] <- - 2 * afs[is.0[i,]]^2
        W[i, is.1[i,]] <- 2 * afs[is.1[i,]] * (1 - afs[is.1[i,]])
        W[i, is.2[i,]] <- - 2 * (1 - afs[is.2[i,]])^2
      }
      gen.rel <- tcrossprod(W, W) / sum((2 * afs * (1 - afs))^2)

    } else if(method == "noia"){
      stop("dominance noia not yet implemented")

      ## compute genotype frequencies at each SNP
      freqG <- t(apply(X, 2, function(x){
        c("A1A1"=sum(x == 0), "A1A2"=sum(x == 1), "A2A2"=sum(x == 2)) / length(x)
      }))

      boundaries <- seq(from=0, to=2, length.out=4)
      is.0 <- (X <= boundaries[2]) # homozygotes for the first allele
      is.1 <- (X > boundaries[2] & X <= boundaries[3]) # heterozygotes
      is.2 <- (X > boundaries[3]) # homozygotes for the second allele

      ## eq 2 of Vitezica et al (2017)
      ## TODO
      W <- matrix(NA, nrow=N, ncol=P, dimnames=dimnames(X))
      for(i in 1:N){
        W[i, is.0[i,]] <- NA
        W[i, is.1[i,]] <- NA
        W[i, is.2[i,]] <- NA
      }

      WWt <- tcrossprod(W, W)
      gen.rel <- WWt / (sum(diag(WWt)) / N)

    } else if(method == "su"){
      H <- recodeIntoDominant(X=X)
      H <- sweep(x=H, MARGIN=2, STATS=2*afs*(1-afs), FUN="-")
      gen.rel <- tcrossprod(H, H) /
        (2 * sum(afs * (1 - afs) * (1 - 2 * afs * (1 - afs))))
    }

  } else if(relationships == "gaussian"){
    M <- X - 1 # recode genotypes as {-1,0,1}
    gen.dist <- as.matrix(stats::dist(x=M, method="euclidean")) / (2 * sqrt(P))
    gen.rel <- exp(-(gen.dist / theta)^2)
  }

  ## force to be perfectly symmetric (beyond machine precision)
  gen.rel[lower.tri(gen.rel)] <- t(gen.rel)[lower.tri(t(gen.rel))]

  ## make invertible
  if(! is.null(alpha)){
    gen.rel <- (1 - alpha) * gen.rel + alpha * gen.rel
  }

  ## keep the nb of SNPs used in an attribute
  attr(gen.rel, "nbSnps") <- P

  return(gen.rel)
}

##' Pairwise linkage disequilibrium
##'
##' Estimates linkage disequilibrium between pairs of SNPs when the observations are the genotypes of genotypes, not their gametes (i.e. the gametic phases are unknown).
##' When ignoring kinship and population structure, the estimator of Rogers and Huff (Genetics, 2009) can be used.
##' When kinship and/or population structure are controlled for, the estimator of Mangin et al (Heredity, 2012) is used via their LDcorSV package.
##' @param X matrix of bi-allelic SNP genotypes encoded in allele doses in {0,1,2}, with genotypes in rows and SNPs in columns; missing values should be encoded as NA
##' @param snp.coords data.frame with SNP identifiers as row names, and two columns, "chr" and "pos"
##' @param K matrix of "kinship" (additive genetic relationships)
##' @param pops vector of characters indicating the population of each genotype
##' @param only.chr identifier of a given chromosome
##' @param only.pop identifier of a given population
##' @param only.snp identifier of a given SNP; compatible with neither use.ldcorsv nor use.snpStats (yet?)
##' @param use.ldcorsv required if K and/or pops are not NULL; otherwise use the square of \code{\link{cor}}
##' @param use.snpStats if TRUE, the \code{ld} function of the snpStats package is used (see \href{https://dx.doi.org/10.1159/000101422}{Clayton and Leung, 2007})
##' @param as.cor if TRUE, the square root of the LD estimates is returned
##' @param as.symmat if TRUE, LD values are returned as a symmetric matrix (not with \code{use.ldcorsv} set to TRUE)
##' @param verbose verbosity level (0/1)
##' @return data frame with at least three columns, "loc1", "loc2" and the LD values
##' @author Timothee Flutre
##' @seealso \code{\link{plotLd}}
##' @examples \dontrun{## make fake data
##' library(scrm)
##' set.seed(1859)
##' nb.genos <- 100
##' Ne <- 10^4
##' nb.chrs <- 1
##' chrom.len <- 10^5
##' mu <- 10^(-8)
##' c.rec <- 10^(-8)
##' genomes <- simulCoalescent(nb.inds=nb.genos, nb.reps=nb.chrs,
##'                            pop.mut.rate=4 * Ne * mu * chrom.len,
##'                            pop.recomb.rate=4 * Ne * c.rec * chrom.len,
##'                            chrom.len=chrom.len)
##'
##' ## checks
##' afs <- estimSnpAf(X=genomes$genos)
##' summary(afs)
##' plotHistAllelFreq(afs=afs)
##' mafs <- estimSnpMaf(afs=afs)
##' plotHistMinAllelFreq(maf=mafs)
##' plotHaplosMatrix(haplos=genomes$haplos$chr1)
##'
##' ## subset SNPs
##' min.maf <- 0.15
##' length(snps.tokeep <- rownames(genomes$snp.coords[mafs >= min.maf,]))
##'
##' ## LD estimator of Rogers and Huff
##' system.time(ld <- estimLd(X=genomes$genos[,snps.tokeep],
##'                           snp.coords=genomes$snp.coords[snps.tokeep,]))
##' dim(ld)
##' head(ld)
##' summary(ld$cor2)
##'
##' ## LD estimator of Mangin et al
##' system.time(ld2 <- estimLd(X=genomes$genos[,snps.tokeep],
##'                            snp.coords=genomes$snp.coords[snps.tokeep,],
##'                            use.ldcorsv=TRUE))
##' dim(ld2)
##' head(ld2)
##'
##' ## physical distance between SNP pairs for which LD was computed
##' dis <- distSnpPairs(snp.pairs=ld[, c("loc1","loc2")],
##'                     snp.coords=genomes$snp.coords[snps.tokeep,])
##'
##' ## plot LD
##' plotLd(x=dis, y=sqrt(ld$cor2), estim="r",
##'        main=paste0(length(snps.tokeep), " SNPs with MAF >= ", min.maf),
##'        sample.size=2*nb.genos, add.ohta.kimura=TRUE, Ne=Ne, c=c.rec)
##'
##' ## physical distance between consecutive SNPs
##' tmp <- distConsecutiveSnps(snp.coords=genomes$snp.coords)
##' hist(tmp[["chr1"]], breaks="FD", xlab="in bp",
##'      las=1, col="grey", border="white",
##'      main="Distances between consecutive SNPs")
##' }
##' @export
estimLd <- function(X, snp.coords, K=NULL, pops=NULL,
                    only.chr=NULL, only.pop=NULL, only.snp=NULL,
                    use.ldcorsv=FALSE, use.snpStats=FALSE,
                    as.cor=FALSE, as.symmat=FALSE, verbose=1){
  if(use.ldcorsv)
    stopifnot(requireNamespace("LDcorSV", quietly=TRUE))
  if(use.snpStats)
    stopifnot(requireNamespace("snpStats", quietly=TRUE))
  stopIfNotValidGenosDose(X)
  stopifnot(.isValidSnpCoords(snp.coords))
  if(! is.null(K))
    stopifnot(use.ldcorsv,
              is.matrix(K),
              nrow(K) == ncol(K),
              nrow(K) == nrow(X),
              ! is.null(dimnames(K)),
              all(rownames(K) == colnames(K)),
              all(rownames(K) == rownames(X)),
              ! as.symmat)
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
  if(! is.null(only.snp)){
    stopifnot(only.snp %in% colnames(X),
              only.snp %in% rownames(snp.coords))
    if(use.ldcorsv)
      stop("use.ldcorsv not available with only.snp")
    if(use.snpStats)
      stop("use.snpStats not available with only.snp")
  }

  ld <- NULL # output

  ## subset the data
  subset.snps <- 1:ncol(X)
  subset.inds <- 1:nrow(X)
  if(! is.null(only.chr) | ! is.null(only.pop)){
    if(! is.null(only.chr)){
      if(verbose > 0)
        write("extract subset of SNPs...", stdout())
      subset.snps <- which(snp.coords$chr == only.chr)
      if(! is.null(only.snp))
        if(! only.snp %in% rownames(snp.coords)[subset.snps])
          subset.snps <- c(which(rownames(snp.coords) == only.snp),
                           subset.snps)
    }
    if(! is.null(only.pop)){
      if(verbose > 0)
        write("extract subset of genotypes...", stdout())
      subset.inds <- which(pops == only.pop)
    }
  }
  X <- X[subset.inds, subset.snps]
  if(! is.null(K))
    K <- K[subset.inds, subset.inds]

  ## internal function to convert a symmetric matrix into a data frame
  symmat2longdf <- function(m){
    data.frame(t(utils::combn(colnames(m), 2)),
               m[lower.tri(m)],
               stringsAsFactors=TRUE)
  }

  ## estimate LD
  if(verbose > 0)
    write("estimate pairwise LD...", stdout())
  if(is.null(K)){
    if(is.null(pops)){
      if(use.ldcorsv){
        ld <- LDcorSV::LD.Measures(donnees=X,
                                   V=NA,
                                   S=NA,
                                   data="G", supinfo=FALSE, na.presence=FALSE)
      } else if(use.snpStats){
        X <- methods::as(X, "SnpMatrix")
        LDmat <- snpStats::ld(X, depth=ncol(X) - 1,
                              stats=ifelse(as.cor, "R", "R2"))
        LDmat <- as.matrix(LDmat)
        LDmat[lower.tri(LDmat)] <- t(LDmat)[lower.tri(LDmat)]
        diag(LDmat) <- 1
        isSup1 <- (LDmat > 1)
        if(any(isSup1)){
          LDmat[which(isSup1)] <- 1
        }
        if(as.symmat){
          ld <- LDmat
        } else{
          ld <- symmat2longdf(LDmat)
          colnames(ld) <- c("loc1", "loc2", ifelse(as.cor, "cor", "cor2"))
        }
      } else{
        if(is.null(only.snp)){
          LDmat <- stats::cor(X)
        } else{
          LDmat <- stats::cor(x=X[, only.snp],
                              y=X[, -which(colnames(X) == only.snp)])
          rownames(LDmat) <- only.snp
        }
        if(! as.cor){
          LDmat <- LDmat^2
        }
        if(as.symmat){
          ld <- LDmat
        } else{
          if(is.null(only.snp)){
            ld <- symmat2longdf(LDmat)
            colnames(ld) <- c("loc1", "loc2", ifelse(as.cor, "cor", "cor2"))
          } else{
            ld <- data.frame(loc1=only.snp,
                             loc2=colnames(LDmat))
            ld[[ifelse(as.cor, "cor", "cor2")]] <- LDmat[1,]
          }
        }
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
  if(! as.symmat){
    if(! is.factor(ld$loc1))
      ld$loc1 <- as.factor(ld$loc1)
    if(! is.factor(ld$loc2))
      ld$loc2 <- as.factor(ld$loc2)
    if(use.ldcorsv & as.cor){
      idx <- grep("^r2", colnames(ld))
      for(i in idx){
        ld[,i] <- sqrt(ld[,i])
        colnames(ld)[i] <- sub("2", "", colnames(ld)[i])
      }
    }
  }

  return(ld)
}

##' Pairwise linkage disequilibrium
##'
##' Estimates linkage disequilibrium between pairs of SNPs belonging to the same chromosome when the observations are the genotypes of genotypes, not their gametes (i.e. the gametic phases are unknown).
##' When ignoring kinship and population structure, the estimator of Rogers and Huff (Genetics, 2009) can be used.
##' When kinship and/or population structure are controlled for, the estimator of Mangin et al (Heredity, 2012) is used via their LDcorSV package.
##' @param X matrix of bi-allelic SNP genotypes encoded in allele doses in {0,1,2}, with genotypes in rows and SNPs in columns; missing values should be encoded as NA
##' @param snp.coords data.frame with SNP identifiers as row names, and two columns, "chr" and "pos"
##' @param K matrix of "kinship" (additive genetic relationships)
##' @param pops vector of characters indicating the population of each genotype
##' @param only.pop identifier of a given population
##' @param use.ldcorsv required if K and/or pops are not NULL; otherwise use the square of \code{\link{cor}}
##' @param use.snpStats if TRUE, the \code{ld} function of the snpStats package is used (see \href{https://dx.doi.org/10.1159/000101422}{Clayton and Leung, 2007})
##' @param as.cor if TRUE, the square root of the LD estimates is returned
##' @param as.symmat if TRUE, LD values are returned as a symmetric matrix (not with \code{use.ldcorsv} set to TRUE)
##' @param nb.cores number of cores to estimate LD for each chromosome in parallel
##' @param verbose verbosity level (0/1)
##' @return data frame with at least three columns, "loc1", "loc2" and the LD values
##' @author Timothee Flutre
##' @seealso \code{\link{estimLd}}, \code{\link{plotLd}}
##' @export
estimLdPerChr <- function(X, snp.coords, K=NULL, pops=NULL, only.pop=NULL,
                          use.ldcorsv=FALSE, use.snpStats=FALSE,
                          as.cor=FALSE, as.symmat=FALSE,
                          nb.cores=1, verbose=1){
  stopifnot(.isValidSnpCoords(snp.coords))

  chrs <- sort(unique(as.character(snp.coords$chr)))

  out <- parallel::mclapply(seq_along(chrs), function(i){
    chr <- chrs[i]
    if(all(verbose > 0, nb.cores <= 1))
      write(chr, stdout())
    estimLd(X=X, snp.coords=snp.coords, K=K, pops=pops,
            only.chr=chr, only.pop=only.pop,
            use.ldcorsv=use.ldcorsv, use.snpStats=use.snpStats,
            as.cor=as.cor, as.symmat=as.symmat,
            verbose=ifelse(all(verbose > 0, nb.cores <= 1), verbose, 0))
  }, mc.cores=nb.cores)

  if(as.symmat){
    names(out) <- chrs
  } else{
    out <- do.call(rbind, out)
    rownames(out) <- NULL
  }

  return(out)
}

##' Pairwise linkage disequilibrium
##'
##' Summarize estimates of linkage disequilibrium between pairs of SNPs belonging to the same chromosome over non-overlapping bins.
##' @param ld data frame returned by \code{\link{estimLd}}
##' @param coln.var name of the column in \code{ld} which contains the LD estimates
##' @param coln.dist name of the column in \code{ld} which contains the physical distance between SNPs belonging to the same pair
##' @param bin.width width of each successive bin (they won't overlap)
##' @param max.phy.len maximum physical length to consider
##' @param span the parameter alpha which controls the degree of loess smoothing
##' @param nb.cores the number of cores to use, i.e. at most how many child processes will be run simultaneously (not on Windows)
##' @return list with the LD summarized over bins as well as loess fits of LD vs mid bin
##' @author Timothee Flutre
##' @seealso \code{\link{estimLd}}, \code{\link{estimLdPerChr}}, \code{\link{fitPhyDistVsLd}}
##' @export
summarizeLd <- function(ld, coln.var="cor2", coln.dist="dist.locs",
                        bin.width=500, max.phy.len=10^6, span=0.75,
                        nb.cores=1){
  stopifnot(is.data.frame(ld),
            all(c(coln.dist, coln.var) %in% colnames(ld)),
            max.phy.len > bin.width)

  out <- list()

  ## summarize LD over non-overlapping bins
  bin.starts <- seq(0, max.phy.len, by=bin.width)
  ld.sry <- parallel::mclapply(seq_along(bin.starts)[-1], function(i){
    bin.start <- bin.starts[i-1]
    bin.end <- bin.starts[i]
    idx <- which(ld[[coln.dist]] >= bin.start & ld[[coln.dist]] < bin.end)
    suppressWarnings(betterSummary(ld[[coln.var]][idx]))
  }, mc.cores=nb.cores)
  ld.sry <- do.call(rbind, ld.sry)
  ld.sry <- cbind(ld.sry,
                  bin.start=bin.starts[-length(bin.starts)],
                  bin.mid=bin.starts[-length(bin.starts)] +
                    round(bin.width/2))
  ld.sry <- as.data.frame(ld.sry)

  ## get fitted values from loess of LD summary w.r.t mid bin
  fit.mean <- stats::loess(mean ~ bin.mid, data=ld.sry, span=span,
                           na.action=stats::na.exclude)
  ld.sry$loess.mean <- stats::fitted(fit.mean)
  fit.med <- stats::loess(med ~ bin.mid, data=ld.sry, span=span,
                          na.action=stats::na.exclude)
  ld.sry$loess.med <- stats::fitted(fit.med)
  out$ld.sry <- ld.sry

  return(out)
}

##' Pairwise linkage disequilibrium
##'
##' Summarize estimates of linkage disequilibrium between pairs of SNPs belonging to the same chromosome over non-overlapping bins per chromosome.
##' @param ld data frame returned by \code{\link{estimLd}}
##' @param coln.var name of the column in \code{ld} which contains the LD estimates
##' @param coln.dist name of the column in \code{ld} which contains the physical distance between SNPs belonging to the same pair
##' @param coln.chr name of the column in \code{ld} which contains the chromosome of each SNP pair
##' @param bin.width width of each successive bin (they won't overlap)
##' @param max.phy.len maximum physical length to consider
##' @param span the parameter alpha which controls the degree of loess smoothing
##' @param nb.cores the number of cores to use, i.e. at most how many child processes will be run simultaneously (not on Windows)
##' @param verbose verbosity level (0/1)
##' @return list with a list per chromosome with the LD summarized over bins as well as loess fits of LD vs mid bin
##' @author Timothee Flutre
##' @seealso \code{\link{summarizeLd}}
##' @export
summarizeLdPerChr <- function(ld, coln.var="cor2", coln.dist="dist.locs",
                              coln.chr="chr",
                              bin.width=500, max.phy.len=10^6, span=0.75,
                              nb.cores=1, verbose=0){
  stopifnot(is.data.frame(ld),
            all(c(coln.dist, coln.var) %in% colnames(ld)))

  if(! is.factor(ld[[coln.chr]]))
    ld[[coln.chr]] <- as.factor(ld[[coln.chr]])
  chrs <- levels(ld[[coln.chr]])

  out <- lapply(seq_along(chrs), function(i){
    chr <- chrs[i]
    if(verbose > 0)
      write(chr, stdout())
    summarizeLd(ld=droplevels(ld[ld[[coln.chr]] == chr,]),
                coln.var=coln.var, coln.dist=coln.dist,
                bin.width=bin.width, max.phy.len=max.phy.len, span=span,
                nb.cores=nb.cores)
  })
  names(out) <- chrs

  return(out)
}

##' Pairwise linkage disequilibrium
##'
##' Fit a model of physical distance versus LD (or reciprocally).
##' @param x numeric vector of predictors (LD values or physical distances)
##' @param y numeric vector of responses (LD values or physical distances)
##' @param method loess/lm/loglm/asyreg/biexp
##' @param span the parameter alpha which controls the degree of loess smoothing
##' @param newdata an optional data frame in which to look for variables with which to predict
##' @return list with data, model fit and prediction(s)
##' @author Timothee Flutre
##' @seealso \code{\link{summarizeLd}}
##' @export
fitPhyDistVsLd <- function(x, y, method="loess", span=0.2,
                           newdata){
  stopifnot(length(x) == length(y),
            method %in% c("loess", "lm", "loglm", "asyreg", "biexp"))

  out <- list(fit=NULL, pred=NULL)

  dat <- data.frame(x=x, y=y)
  out$dat <- dat

  if(method == "loess"){
    fit <- try(
        stats::loess(y ~ x,
                     data=dat,
                     na.action=stats::na.exclude,
                     span=span))
  } else if(method == "lm"){
    fit <- try(
        stats::lm(y ~ x,
                  data=dat,
                  na.action=stats::na.exclude))
  } else if(method == "loglm"){
    dat$log.y <- log(dat$y)
    idx <- which(is.nan(dat$log.y) | is.infinite(dat$log.y))
    if(length(idx) > 0)
      dat$log.y[idx] <- NA
    fit <- try(
        stats::lm(log.y ~ x,
                  data=dat,
                  na.action=stats::na.exclude))
  } else if(method == "asyreg"){
    fit <- try(
        stats::nls(y ~ SSasymp(x, Asym, R0, lrc),
                   data=dat,
                   na.action=stats::na.exclude))
  } else if(method == "biexp"){
    fit <- try(
        stats::nls(y ~ SSbiexp(x, A1, lrc1, A2, lrc2),
                   data=dat,
                   na.action=stats::na.exclude))
  }
  out$fit <- fit

  if(! missing(newdata)){
    if(! methods::is(fit, "try-error")){
      out$pred <- stats::predict(fit, newdata=newdata)
      if(method == "loglm")
        out$pred <- exp(out$pred)
    }
  }

  return(out)
}

##' Pairwise linkage disequilibrium
##'
##' Plots the linkage disequilibrium between pairs of SNPs, as a blue density or black points, with a red loess.
##' Possibility to add two analytical approximations of E[r^2] at equilibrium (see McVean, Handbook of Stat Gen, 2007): 1 / (1 + 4 Ne c x) by Sved (1971) and (10 + 4 Ne c x) / (22 + 13 * 4 Ne c x + (4 Ne c x)^2) by Ohta and Kimura (1971).
##' @param x vector of distances between SNPs (see \code{\link{distSnpPairs}})
##' @param y vector of LD estimates (see \code{\link{estimLd}})
##' @param main main title
##' @param estim estimator of pairwise LD corresponding to the values in y (r2/r)
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
##' @param xlim numeric vector of length 2 specifying the x-axis limit (optional)
##' @return invisible list
##' @author Timothee Flutre
##' @seealso \code{\link{estimLd}}, \code{\link{distSnpPairs}}
##' @export
plotLd <- function(x, y, main="", estim="r2",
                   use.density=TRUE,
                   xlab="Physical distance (bp)",
                   ylab=paste0("Linkage disequilibrium (", estim, ")"),
                   span=1/10, degree=1, evaluation=50,
                   sample.size=NULL,
                   add.ohta.kimura=FALSE, add.sved=FALSE, Ne=NULL, c=NULL,
                   xlim){
  stopifnot(is.vector(x),
            is.vector(y),
            estim %in% c("r2","r"))
  if(! is.null(sample.size))
    stopifnot(is.numeric(sample.size))
  if(any(c(add.ohta.kimura, add.sved)))
    stopifnot(! is.null(Ne),
              ! is.null(c))

  out <- list()

  ## plot the pairwise estimates of LD
  lpars <- list(col="red", cex=2)
  if(use.density){
    graphics::smoothScatter(x, y,
                  main=main,
                  xlab=xlab,
                  ylab=ylab,
                  las=1,
                  xlim=xlim)
    pred <- stats::loess.smooth(x, y, span=span, degree=degree,
                                evaluation=evaluation)
    do.call(graphics::lines, c(list(pred), lpars))
    out$loess <- pred
  } else{
    stats::scatter.smooth(x, y, lpars=lpars,
                          main=main, xlab=xlab, ylab=ylab, las=1,
                          span=span, degree=degree, evaluation=evaluation,
                          xlim=xlim)
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
    out$X2 <- X2
  }

  ## add analytical approximations
  if(any(c(add.ohta.kimura, add.sved)))
    scaled.dist <- 4 * Ne * c * x
  if(add.ohta.kimura){
    ok <- (10 + scaled.dist) / (22 + 13 * scaled.dist + scaled.dist^2)
    if(estim == "r")
      ok <- sqrt(ok)
    graphics::points(x, ok, pch=".", col="purple", cex=1.2)
    out$ohta.kimura <- ok
  }
  if(add.sved){
    sved <- 1 / (1 + scaled.dist)
    if(estim == "r")
      sved <- sqrt(sved)
    graphics::points(x, sved, pch=".", col="green", cex=1.2)
    out$sved <- sved
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

  invisible(out)
}

##' Pairwise linkage disequilibrium
##'
##' Plots a summary of linkage disequilibrium between pairs of SNPs per non-overlapping bin (mean, median and quartiles 1 and 3).
##' @param x vector of distances between SNPs (see \code{\link{distSnpPairs}})
##' @param y vector of LD estimates (see \code{\link{estimLd}})
##' @param bin.len length of each bin
##' @param max.dist maximum distance between SNPs to consider
##' @param main main title
##' @param xlab label for the x axis
##' @param ylab label for the y axis
##' @return invisible data frame with the summary per bin
##' @author Timothee Flutre
##' @seealso \code{\link{estimLd}}, \code{\link{distSnpPairs}}, \code{\link{plotLd}}
##' @export
plotLdSry <- function(x, y, bin.len=10^3, max.dist=10^5,
                      main="", xlab="Physical distance (bp)",
                      ylab="Summarized LD"){
  stopifnot(length(x) == length(y))

  idx <- which(x <= max.dist)
  if(length(idx) == 0){
    msg <- paste0("no SNP pair with distance below threshold (", max.dist, ")")
    write(msg, stdout())
    return()
  }
  x <- x[idx]
  y <- y[idx]

  ## summarize LD
  out <- data.frame(bin.start=seq(0, max.dist, by=bin.len),
                    mid.bin=NA,
                    ld.mean=NA,
                    ld.q1=NA,
                    ld.med=NA,
                    ld.q3=NA)
  out <- out[-nrow(out),]
  for(i in 1:nrow(out)){
    idx.bin <- which(x >= out$bin.start[i] &
                     x < out$bin.start[i] + bin.len)
    out$mid.bin[i] <- out$bin.start[i] + bin.len / 2
    out$ld.mean[i] <- mean(y[idx.bin])
    out$ld.q1[i] <- stats::quantile(y[idx.bin], probs=0.25)
    out$ld.med[i] <- stats::median(y[idx.bin])
    out$ld.q3[i] <- stats::quantile(y[idx.bin], probs=0.75)
  }

  ## plot summarized LD
  graphics::plot(x=out$bin.start, y=out$ld.mean,
                 ylim=c(0, 0.5),
                 xlab=xlab, ylab=ylab, main=main,
                 las=1, col="red")
  graphics::points(out$bin.start, out$ld.med, col="black")
  graphics::segments(x0=out$bin.start, y0=out$ld.q1,
                     x1=out$bin.start, y1=out$ld.q3)
  graphics::legend("topright", legend=c("mean", "median", "quartiles 1-3"), bty="n",
                   col=c("red", "black", "black"), pch=c(1, 1, NA),
                   lty=c(NA, NA, 1))

  invisible(out)
}

##' Pairwise linkage disequilibrium
##'
##' Plots a summary of linkage disequilibrium between pairs of SNPs versus physical distance.
##' @param x vector of pysical coordinates of mid bins (see output of \code{\link{summarizeLd}})
##' @param y vector of LD estimates (see output of \code{\link{summarizeLd}})
##' @param xlab label for the x axis
##' @param ylab label for the y axis
##' @param main main title
##' @param fitted.loess see the output of \code{\link{fitPhyDistVsLd}}
##' @param fitted.lm see the output of \code{\link{fitPhyDistVsLd}}
##' @param fitted.loglm see the output of \code{\link{fitPhyDistVsLd}}
##' @param fitted.asyreg see the output of \code{\link{fitPhyDistVsLd}}
##' @param fitted.biexp see the output of \code{\link{fitPhyDistVsLd}}
##' @param ... arguments passed to \code{\link{plot}}
##' @return nothing
##' @seealso \code{\link{summarizeLd}}, \code{\link{fitPhyDistVsLd}}
##' @author Timothee Flutre
##' @export
plotPhyDistVsLdSry <- function(x, y, xlab, ylab, main="",
                               fitted.loess=NULL, fitted.lm=NULL, fitted.loglm=NULL,
                               fitted.asyreg=NULL, fitted.biexp=NULL, ...){
  graphics::plot(x=x, y=y, xlab=xlab, ylab=ylab, main=main, col="darkgrey", pch=1, ...)
  graphics::abline(h=0)
  lgds <- "LD summary"
  cols <- "darkgrey"
  pchs <- 1
  ltys <- NA
  lwds <- NA
  if(! is.null(fitted.loess)){
    graphics::points(x=x, y=fitted.loess, col="blue", pch=1)
    lgds <- c(lgds, "loess")
    cols <- c(cols, "blue")
    pchs <- c(pchs, 1)
    ltys <- c(ltys, NA)
    lwds <- c(lwds, NA)
  }
  if(! is.null(fitted.lm)){
    graphics::lines(x=x, y=fitted.lm, col="red", lty=1, lwd=2)
    lgds <- c(lgds, "LM")
    cols <- c(cols, "red")
    pchs <- c(pchs, NA)
    ltys <- c(ltys, 1)
    lwds <- c(lwds, 2)
  }
  if(! is.null(fitted.loglm)){
    graphics::points(x=x, y=fitted.loglm, col="red", pch=2)
    lgds <- c(lgds, "log LM")
    cols <- c(cols, "red")
    pchs <- c(pchs, 2)
    ltys <- c(ltys, NA)
    lwds <- c(lwds, NA)
  }
  if(! is.null(fitted.asyreg)){
    graphics::points(x=x, y=fitted.asyreg, col="green", pch=3)
    lgds <- c(lgds, "asymptotic regression")
    cols <- c(cols, "green")
    pchs <- c(pchs, 3)
    ltys <- c(ltys, NA)
    lwds <- c(lwds, NA)
  }
  if(! is.null(fitted.biexp)){
    graphics::points(x=x, y=fitted.biexp, col="orange", pch=4)
    lgds <- c(lgds, "biexponential")
    cols <- c(cols, "orange")
    pchs <- c(pchs, 4)
    ltys <- c(ltys, NA)
    lwds <- c(lwds, NA)
  }
  graphics::legend("topright", bty="n", legend=lgds, col=cols, pch=pchs, lty=ltys, lwd=lwds)
}

##' Distance between consecutive SNPs
##'
##' For each pair of consecutive SNPs, return the number of "blocks" (i.e. nucleotides) between both SNPs (to be coherent with \code{\link{distSnpPairs}}).
##' @param snp.coords data.frame with SNP identifiers as row names, and with two columns "chr" and "coord" or "pos"
##' @param only.chr identifier of a given chromosome
##' @param nb.cores the number of cores to use
##' @return list with one component per chromosome
##' @author Timothee Flutre
##' @seealso \code{\link{estimLd}}
##' @examples \dontrun{## make fake data
##' snp.coords <- data.frame(chr=c("chr1","chr1","chr1","chr2"),
##'                          pos=c(150, 131, 171, 17))
##' rownames(snp.coords) <- paste0("snp", 1:nrow(snp.coords))
##'
##' distConsecutiveSnps(snp.coords)
##' }
##' @export
distConsecutiveSnps <- function(snp.coords, only.chr=NULL, nb.cores=1){
  requireNamespace("parallel", quietly=TRUE)
  stopifnot(.isValidSnpCoords(snp.coords))
  if(! "coord" %in% colnames(snp.coords))
    colnames(snp.coords)[colnames(snp.coords) == "pos"] <- "coord"

  snp.coords$chr <- as.character(snp.coords$chr)
  chr.ids <- unique(snp.coords$chr)
  if(! is.null(only.chr)){
    stopifnot(only.chr %in% chr.ids)
    chr.ids <- only.chr
  } else{
    if(requireNamespace("gtools", quietly=TRUE)){
      chr.ids <- gtools::mixedsort(chr.ids)
    } else
      chr.ids <- sort(chr.ids)
  }

  snp.dists <- parallel::mclapply(chr.ids, function(chr.id){
    coords <- snp.coords$coord[snp.coords$chr == chr.id]
    if(length(coords) > 1){
      names(coords) <- rownames(snp.coords)[snp.coords$chr == chr.id]
      coords <- sort(coords)
      dis <- coords[2:length(coords)] - coords[1:(length(coords)-1)] - 1
      names(dis) <- paste(names(dis), names(coords)[-length(coords)], sep="-")
      dis
    } else
      NA
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
##' @seealso \code{\link{pruneSnpsLd}}
##' @export
thinSnps <- function(method, threshold, snp.coords, only.chr=NULL){
  requireNamespace("GenomicRanges", quietly=TRUE)
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
      tmp.gr <- snpCoordsDf2Gr(tmp)
      ovl <- GenomicRanges::findOverlaps(tiles, tmp.gr)
      idx <- sapply(as.list(ovl), `[`, 1)
      names(tmp.gr[idx[! is.na(idx)]])
    }
  }))

  return(out.snp.ids)
}

##' Prune SNPs based on LD
##'
##' Prune SNPs based on their pairwise linkage disequilibriums via sliding windows using the "snpgdsLDpruning" function from the "SNPRelate" package..
##' @param X matrix of bi-allelic SNP genotypes encoded, for each SNP, in number of copies of its second allele, i.e. as allele doses in [0,2], with genotypes in rows and SNPs in columns; the "second" allele is arbitrary, it can correspond to the minor (least frequent) or the major (most frequent) allele; will be transformed into a "gds" object (see next argument); if some values were imputed, they will be automatically thresholded using \code{\link{convertImputedTo012}} (required by the SNPRelate function)
##' @param snp.coords data frame which row names are SNP identifiers, the first column should contain chromosomes as integers (otherwise \code{\link{chromNames2integers}} will be used), and the second column should contain positions; compulsory if the X argument is specified
##' @param gds object of class "SNPGDSFileClass" from the SNPRelate package
##' @param ld.threshold the LD threshold below which SNPs are discarded; will be passed to "snpgdsLDpruning"
##' @param remove.monosnp will be passed to "snpgdsLDpruning"
##' @param maf will be passed to "snpgdsLDpruning"
##' @param missing.rate will be passed to "snpgdsLDpruning"
##' @param method will be passed to "snpgdsLDpruning"
##' @param slide.max.bp will be passed to "snpgdsLDpruning"
##' @param slide.max.n will be passed to "snpgdsLDpruning"
##' @param seed seed for the pseudo-random number generator
##' @param nb.cores will be passed to "snpgdsLDpruning"
##' @param verbose verbosity level (0/1)
##' @return vector of SNP identifiers
##' @author Timothee Flutre
##' @seealso \code{\link{thinSnps}}, \code{\link[SNPRelate]{snpgdsLDpruning}}
##' @export
pruneSnpsLd <- function(X=NULL, snp.coords=NULL, gds=NULL,
                        ld.threshold=0.2,
                        remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
                        method=c("composite", "r", "dprime", "corr"),
                        slide.max.bp=500000, slide.max.n=NA,
                        seed=NULL, nb.cores=1, verbose=1){
  requireNamespace("SNPRelate", quietly=TRUE)
  stopifnot(xor(is.null(X), is.null(gds)))
  if(is.null(gds)){
    if(any(! X %in% c(0,1,2)))
      X <- convertImputedTo012(X)
    stopIfNotValidGenosDose(X=X, check.noNA=FALSE, check.notImputed=TRUE)
    stopifnot(! is.null(snp.coords),
              is.data.frame(snp.coords),
              ! is.null(rownames(snp.coords)),
              all(colnames(X) %in% rownames(snp.coords)),
              ncol(snp.coords) >= 2,
              all(! is.na(snp.coords[,1])),
              all(! is.na(snp.coords[,2])),
              is.numeric(snp.coords[,2]))
    conv.chr2int <- FALSE
    if(! is.numeric(snp.coords[,1])){
      tmp <- suppressWarnings(as.integer(snp.coords[,1]))
      conv.chr2int <- TRUE
    }
  }

  gds.file <- NULL
  if(is.null(gds)){
    gds.file <- tempfile()
    snp.coords <- snp.coords[colnames(X),] # re-order the rows
    if(conv.chr2int){
      cn2i <- chromNames2integers(snp.coords[,1], basic=TRUE)
      snp.coords[,1] <- cn2i$renamed
    }
    SNPRelate::snpgdsCreateGeno(gds.fn=gds.file,
                                genmat=X,
                                sample.id=rownames(X),
                                snp.id=colnames(X),
                                snp.chromosome=snp.coords[,1],
                                snp.position=snp.coords[,2],
                                snpfirstdim=FALSE)
    gds <- SNPRelate::snpgdsOpen(gds.file)
  }

  if(verbose > 1)
    SNPRelate::snpgdsSummary(gds)

  if(! is.null(seed))
    set.seed(seed)
  out.snp.ids <-
    SNPRelate::snpgdsLDpruning(gdsobj=gds,
                               remove.monosnp=remove.monosnp,
                               maf=maf,
                               missing.rate=missing.rate,
                               method=method,
                               slide.max.bp=slide.max.bp,
                               slide.max.n=slide.max.n,
                               ld.threshold=ld.threshold,
                               num.thread=nb.cores,
                               verbose=ifelse(verbose > 0, TRUE, FALSE))
  out.snp.ids <- unlist(out.snp.ids)
  names(out.snp.ids) <- NULL

  if(! is.null(gds.file)){
    SNPRelate::snpgdsClose(gds)
    file.remove(gds.file)
  }

  return(out.snp.ids)
}

##' Genotype imputation
##'
##' Impute missing SNP genotypes with the mean of the non-missing.
##' The code was inspired from the \code{A.mat} function from the \href{https://cran.r-project.org/web/packages/rrBLUP/}{rrBLUP} package.
##' @param X matrix of bi-allelic SNP genotypes encoded in allele doses in {0,1,2}, with genotypes in rows and SNPs in columns; missing values should be encoded as NA
##' @param min.maf minimum minor allele frequency (before imputation) below which SNPs are discarded
##' @param max.miss maximum amount of missing genotypes (before imputation) above which SNPs are discarded
##' @param rm.still.miss if TRUE, remove SNP(s) still with missing genotype(s) after imputation (depending on \code{min.maf} and \code{max.miss})
##' @return matrix with imputed genotypes
##' @author Timothee Flutre
##' @seealso \code{\link{imputeGenosWithMeanPerPop}}
##' @export
imputeGenosWithMean <- function(X, min.maf=0.1, max.miss=0.3, rm.still.miss=TRUE){
  ## requireNamespace("rrBLUP")
  stopIfNotValidGenosDose(X, check.noNA=FALSE, check.isDose=FALSE)

  X.out <- X
  if(any(is.na(X))){
    ## perform imputation
    ## tmp <- rrBLUP::A.mat(X=X - 1, min.MAF=min.maf, max.missing=max.miss,
    ##                      impute.method="mean", return.imputed=TRUE)
    ## X.out <- tmp$imputed + 1
    freq.miss <- calcFreqMissSnpGenosPerSnp(X)
    afs <- estimSnpAf(X=X)
    mafs <- estimSnpMaf(afs=afs)
    tokeep <- which((mafs >= min.maf) & (freq.miss <= max.miss))
    ones <- matrix(data=1, nrow=nrow(X), ncol=1)
    afs.mat <- tcrossprod(ones, matrix(afs[tokeep]))
    idx.na <- which(is.na(X.out[,tokeep]))
    X.out[,tokeep][idx.na] <- 0
    X.out[,tokeep][idx.na] <- (X.out[,tokeep] + 2 * afs.mat)[idx.na]

    ## remove SNPs with remaining missing genotypes
    if(all(rm.still.miss, sum(is.na(X.out)) > 0))
      X.out <- discardMarkersMissGenos(X=X.out, verbose=0)
  }

  return(X.out)
}

##' Genotype imputation
##'
##' Impute missing SNP genotypes with the mean of the non-missing, population by population.
##' @param X matrix of bi-allelic SNP genotypes encoded in allele doses in {0,1,2}, with genotypes in rows and SNPs in columns; missing values should be encoded as NA
##' @param pops data.frame with 2 columns, "ind" and "pop"
##' @param min.maf.pop minimum minor allele frequency per population (before imputation) below which SNPs are discarded
##' @param max.miss.pop maximum amount of missing genotypes per population (before imputation) above which SNPs are discarded
##' @param rm.still.miss if TRUE, remove SNP(s) still with missing genotype(s) in at least one population (depending on \code{min.maf.pop} and \code{max.miss.pop})
##' @param verbose verbosity level (0/1)
##' @return matrix with imputed genotypes
##' @author Timothee Flutre
##' @seealso \code{\link{imputeGenosWithMean}}
##' @export
imputeGenosWithMeanPerPop <- function(X, pops, min.maf.pop=0.1,
                                      max.miss.pop=0.3, rm.still.miss=TRUE,
                                      verbose=1){
  stopIfNotValidGenosDose(X, check.noNA=FALSE, check.isDose=FALSE)
  stopifnot(is.data.frame(pops),
            "ind" %in% colnames(pops),
            "pop" %in% colnames(pops))

  X.out <- X

  pops$ind <- as.character(pops$ind)
  pops$pop <- as.character(pops$pop)
  if(any(is.na(pops$pop)))
    pops <- pops[! is.na(pops$pop),]

  if(verbose > 0){
    msg <- paste0("nb of genotypes: ", nrow(X),
                  "\nnb of SNPs: ", ncol(X),
                  "\nnb of genotypes: ", nrow(X) * ncol(X),
                  "\nnb of missing genotypes: ", sum(is.na(X)),
                  " (", format(100 * sum(is.na(X)) / (nrow(X) * ncol(X)),
                               digits=2), "%)",
                  "\nnb of populations: ", length(unique(pops$pop)))
    write(msg, stdout())
  }

  for(pop in unique(pops$pop)){
    ## identify subset of genotypes
    idxI <- which(rownames(X) %in% pops$ind[pops$pop == pop])
    if(verbose > 0){
      msg <- paste0("impute missing genotypes in ", length(idxI),
                    " genotypes from ", pop, " ...")
      write(msg, stdout())
    }
    ## perform imputation
    X.out[idxI,] <- imputeGenosWithMean(X=X[idxI,], min.maf=min.maf.pop,
                                        max.miss=max.miss.pop,
                                        rm.still.miss=FALSE)
  }
  X.out[X.out < 0 & abs(X.out) < 10^(-6)] <- 0

  if(verbose > 0){
    msg <- paste0("nb of missing genotypes: ", sum(is.na(X.out)),
                  " (", format(100 * sum(is.na(X.out)) / (nrow(X) * ncol(X)),
                               digits=2), "%)")
    write(msg, stdout())
  }

  ## remove SNPs with remaining missing genotypes
  if(all(rm.still.miss, sum(is.na(X.out)) > 0))
    X.out <- discardMarkersMissGenos(X=X.out, verbose=verbose - 1)

  return(X.out)
}

##' Haplotypes
##'
##' Reformat haplotypes from alleles as string into numeric.
##' @param haplos matrix with haplotypes in rows (2 consecutive per genotype) and SNPs in columns
##' @param alleles data.frame with SNPs in rows (names as row names) and alleles in 2 columns
##' @param nb.cores the number of cores to use, i.e. at most how many child processes will be run simultaneously (not on Windows)
##' @return matrix with 0 for the first-column alleles and 1 for the second-column alleles
##' @author Timothee Flutre
##' @seealso \code{\link{segSites2allDoses}}
##' @export
haplosAlleles2num <- function(haplos, alleles, nb.cores=1){
  stopifnot(is.matrix(haplos),
            ! is.null(colnames(haplos)),
            is.data.frame(alleles),
            ! is.null(rownames(alleles)),
            all(colnames(haplos) %in% rownames(alleles)),
            ncol(alleles) == 2)

  out <- matrix(data=NA, nrow=nrow(haplos), ncol=ncol(haplos),
                dimnames=dimnames(haplos))

  alleles <- alleles[colnames(haplos),] # put SNPs in same order

  nb.snps <- ncol(haplos)
  tmp <- parallel::mclapply(1:nb.snps, function(j){
    idx <- which(haplos[,j] == alleles[j, 1])
    if(length(idx) > 0)
      out[idx, j] <<- 0
    idx <- which(haplos[,j] == alleles[j, 2])
    if(length(idx) > 0)
      out[idx, j] <<- 1
  }, mc.cores=nb.cores)

  return(out)
}

##' Compute the inverse of the additive relationship matrix using Henderson algorithm
##'
##' Given the vectors of sire and dam, directly return the additive relationship matrix 'A' (tabular method).
##' Reference:
##' \itemize{
##' \item Henderson, C. R. 1976. Simple Method for Computing the Inverse of a Numerator Relationship Matrix Used in Prediction of Breeding Values. Biometrics 32:69-83.
##' }
##' @param s vector of sires
##' @param d vector of dams
##' @return matrix
##' @author Gota Morota
##' @note Unknown parents should be coded as zero. Last modified by Morota on April 1, 2010.
##' @export
createA <- function(s, d){
  stopifnot(is.vector(s),
            is.vector(d),
            length(s) == length(d))

	n <- length(s)
	N <- n + 1
	A <- matrix(0, ncol=N, nrow=N)

	s <- (s == 0)*(N) + s
	d <- (d == 0)*N + d

	for(i in 1:n){
		A[i,i] <- 1 + A[s[i], d[i]]/2
		for(j in (i+1):n){
			if (j > n) break
			A[i,j] <- (A[i, s[j]] + A[i,d[j]]) / 2
      A[j,i] <- A[i,j]
		}
	}

  return(A[1:n, 1:n])
}

##' Compute the inverse of the additive relationship matrix using Quass algorithm
##'
##' Given the vectors of sire and dam, directly return inverse of additive relationship matrix 'A' without creating the 'A' itself.
##' This is a modification of Henderson's method and unlike createAinv.r, this can be used in inbred populations.
##' Reference:
##' \itemize{
##' \item Quass, R. L. 1976. Computing the Diagonal Elements and Inverse of a Large Numerator Relationship Matrix. Biometrics 32:949-953.
##' }
##' @param s vector of sires
##' @param d vector of dams
##' @return matrix
##' @author Gota Morota
##' @note Unknown parents should be coded as zero. Last modified by Morota on April 2, 2010.
##' @export
quass <- function(s, d){
  stopifnot(is.vector(s),
            is.vector(d),
            length(s) == length(d))

  n <- length(s)
  N <- n
  L <- matrix(0.0, ncol=N, nrow=N)

  s <- (s == 0)*(N) + s
  d <- (d == 0)*N + d

  ## construct L
	for(t in 1:n){

		if(s[t] == N && d[t] == N){
			L[t,t] <- 1.0
			if(t != 1){
				for(j in 1:(t-1)){
					L[t,j] <- 0
				}
			}
		}

		if(s[t] != N && d[t] == N){
			for(j in 1:s[t]){
				L[t,j] <- 0.5 * L[s[t], j]
			}
			tmp <- 0.0
			for(j in 1:s[t]){
				tmp <- tmp + L[t,j]^2
			}
			L[t,t] <- sqrt(1- tmp)
		}

		if(s[t] == N && d[t] != N){
			for(j in 1:d[t]){
				L[t,j] <- 0.5 * L[d[t], j]
			}
			tmp <- 0.0
			for(j in 1:d[t]){
				tmp <- tmp + L[t,j]^2
			}
			L[t,t] <- sqrt(1- tmp)
		}

		if(s[t] != N && d[t] != N ){
			if(s[t] < d[t]){
				p <- s[t]
				q <- d[t]
			}
			else{
				p <- d[t]
				q <- s[t]
			}

			for(j in 1:t-1){
				L[t,j] <- 0.5 * (L[p,j] + L[q,j])
			}

			tmp <- 0.0
			for(j in 1:p){
				tmp <- tmp + L[p,j] * L[q,j]
			}

			tmp2 <- 0.0
			for(j in 1:q){
				tmp2 <- tmp2 + L[t,j]^2
			}

			L[t,t] <- 	sqrt(1 + 0.5 * tmp - tmp2)
		}

	}


  ## calculate inv(A) based on L
	A <- matrix(0.0, ncol=N, nrow=N)

	for(i in 1:n){

		tmp <- 1/L[i,i]^2

		if(s[i] != N && d[i] != N ){
			A[i,i] <- A[i,i] + tmp
			A[i, s[i]] <- A[i, s[i]] - tmp/2.0
			A[s[i], i] <- A[s[i], i] - tmp/2.0
			A[i,d[i]] <- A[i, d[i]] - tmp/2.0
			A[d[i], i] <- A[d[i], i] - tmp/2.0

			A[s[i], s[i]] <- A[s[i], s[i]] + tmp/4.0
			A[s[i], d[i]] <- A[s[i], d[i]] + tmp/4.0
			A[d[i], s[i]] <- A[d[i], s[i]] + tmp/4.0
			A[d[i], d[i]] <- A[d[i], d[i]] + tmp/4.0
		}

		if(s[i] != N && d[i] == N){
			A[i,i] <- A[i,i] + tmp
			A[s[i], i] <- A[s[i], i] - tmp/2.0
			A[i, s[i]] <- A[i, s[i]] - tmp/2.0
			A[s[i], s[i]] <- A[s[i], s[i]] + tmp/4.0
		}

		if(s[i] == N && d[i] != N){
			A[i,i] <- A[i,i] + tmp
			A[d[i], i] <- A[d[i], i] - tmp/2.0
			A[i, d[i]] <- A[i, d[i]] - tmp/2.0
			A[d[i], d[i]] <- A[d[i], d[i]] + tmp/4.0
		}

		if(s[i] == N && d[i] == N){
			A[i,i] <- A[i,i] + tmp
		}

	}

  return(A)
}

##' Animal model
##'
##' Given T traits, I genotypes, Q covariates and N=I*Q phenotypes per trait, simulate phenotypes via the following "animal model": \eqn{Y = W C + Z G_A + Z G_D + E}, where Y is N x T; W is N x Q; Z is N x I; G_A ~ Normal_IxT(0, A, V_{G_A}) with A being IxI and V_{G_A} being TxT; G_D ~ Normal_IxT(0, D, V_{G_D}) with D being IxI and V_{G_D} being TxT; E ~ Normal_NxT(0, Id_N, V_E) with Id_N being NxN and V_E being TxT. Note that using the matrix Normal (MN) is equivalent to using the multivariate Normal (MVN) via the vec operator and Kronecker product: \eqn{Y ~ MN(M, U, V) <=> vec(Y) ~ MVN(vec(M), V \otimes U)}. Spatial heterogeneity can be add (see \href{http://dx.doi.org/10.1139/x02-111}{Dutkowski et al (2002)}).
##' @param T number of traits
##' @param Q number of covariates (as "fixed effects", e.g. replicates)
##' @param mu T-vector of overall means (one per trait), i.e. C[1,1:T]
##' @param mean.C mean of the univariate Normal prior on C[2:Q,1:T] (ignored if Q=1)
##' @param sd.C std dev of the univariate Normal prior on C[2:Q,1:T] (ignored if Q=1)
##' @param A IxI matrix of additive genetic relationships between genotypes (see \code{\link{estimGenRel}} with VanRaden's estimator); if A is singular (i.e. has a large condition number), the Choleski decomposition will be performed with pivoting
##' @param V.G.A scalar (if T=1) or TxT matrix of additive genetic variance-covariance between traits (e.g. 15 when T=1)
##' @param scale.hC.G.A scale of the half-Cauchy prior for sqrt{V_{G_A}} (e.g. 5; used if V.G.A=NULL and T=1)
##' @param nu.G.A degrees of freedom of the Wishart prior for V_{G_A} (used if V.G.A=NULL and T>1)
##' @param D IxI matrix of dominant genetic relationships between genotypes (see \code{\link{estimGenRel}} with Vitezica's estimator); if D is singular (i.e. has a large condition number), the Choleski decomposition will be performed with pivoting
##' @param V.G.D scalar (if T=1) or TxT matrix of dominant genetic variance-covariance between traits (e.g. 3 when T=1; used if D!=NULL)
##' @param scale.hC.G.D scale of the half-Cauchy prior for sqrt{V_{G_D}} (e.g. 5; used if D!=NULL, V.G.D=NULL and T=1)
##' @param nu.G.D degrees of freedom of the Wishart prior for V_{G_D} (used if D!=NULL, V.G.D=NULL and T>1)
##' @param V.E scalar (if T=1) or TxT matrix of error variance-covariance between traits (used if T=1 and err.df=Inf)
##' @param scale.hC.E scale of the half-Cauchy prior for sqrt{V_E} (e.g. 5; used if V.E=NULL and T=1 and err.df=Inf)
##' @param nu.E degrees of freedom of the Wishart prior for V_E (used if V.E=NULL and T>1)
##' @param err.df degrees of freedom of the Student's t-distribution of the errors (e.g. 3; Inf means Normal distribution; will be Inf if T>1)
##' @param perc.NA percentage of missing phenotypes, at random
##' @param seed seed for the pseudo-random number generator
##' @param nb.rows number of rows, when assuming data come from a plant species with a regular spatial conformation
##' @param nb.cols number of columns (nb.rows x nb.columns should be equal to the number of rows of A)
##' @param sigma2.ksi variance of the spatially-correlated errors (only if T=1)
##' @param rho.rows correlation between neighbor plants on the same row (only if T=1)
##' @param rho.cols correlation between neighbor plants on the same column (only if T=1)
##' @return list
##' @author Timothee Flutre
##' @examples
##' \dontrun{## simulate genotypes
##' set.seed(1859)
##' X <- simulGenosDose(nb.genos=200, nb.snps=2000)
##'
##' ## estimate the additive genetic relationships
##' A <- estimGenRel(X, relationships="additive", method="vanraden1", verbose=0)
##'
##' # 1: one trait with heritability h2=0.75, no covariate, Normal errors, no NA
##' model <- simulAnimalModel(T=1, Q=1, A=A, V.G.A=15, V.E=5, seed=1859)
##'
##' # 2: one trait with heritability h2=0.75, three covariates, Normal errors, no NA
##' model <- simulAnimalModel(T=1, Q=3, A=A, V.G.A=15, V.E=5, seed=1859)
##'
##' # 3: one trait with heritability h2=0.75, no covariate, Student errors, no NA
##' model <- simulAnimalModel(T=1, Q=1, A=A, V.G.A=15, err.df=3, seed=1859)
##'
##' # 4: one trait with heritability drawn at random, no covariate, Normal errors, no NA
##' model <- simulAnimalModel(T=1, Q=1, A=A, scale.hC.G.A=10, V.E=5, seed=1859)
##'
##' # 5: two traits with heritabilities drawn at random, no covariate, Normal errors, no NA
##' model <- simulAnimalModel(T=2, Q=1, A=A, nu.G.A=5, nu.E=3, seed=1859)
##'
##' # 6: same as scenario 1 but with spatial coordinates
##' model <- simulAnimalModel(T=1, Q=1, A=A, V.G.A=15, V.E=5, seed=1859, nb.rows=20, nb.cols=10)
##' ## if(require(sp))
##' ##   plot(SpatialPoints(model$coords), axes=TRUE, las=1,
##' ##        xlab="columns", ylab="rows", main="Map")
##' }
##' @export
simulAnimalModel <- function(T=1,
                             Q=3, mu=rep(50,T), mean.C=5, sd.C=2,
                             A, V.G.A=NULL, scale.hC.G.A=NULL, nu.G.A=T,
                             D=NULL, V.G.D=NULL, scale.hC.G.D=NULL, nu.G.D=T,
                             V.E=NULL, scale.hC.E=NULL, nu.E=T,
                             err.df=Inf, perc.NA=0, seed=NULL,
                             nb.rows=NULL, nb.cols=NULL,
                             sigma2.ksi=NULL, rho.rows=NULL, rho.cols=NULL){
  requireNamespace("MASS", quietly=TRUE)
  stopifnot(length(mu) == T,
            is.matrix(A),
            nrow(A) == ncol(A),
            ! is.null(rownames(A)),
            ! is.null(colnames(A)),
            rownames(A) == colnames(A))
  if(! is.null(V.G.A)){
    if(T == 1){
      stopifnot(! is.matrix(V.G.A))
    } else
      stopifnot(is.matrix(V.G.A),
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
      if(T == 1){
        stopifnot(! is.matrix(V.G.D))
      } else
        stopifnot(is.matrix(V.G.D),
                  nrow(V.G.D) == ncol(V.G.D))
    } else{
      if(T == 1){
        stopifnot(! is.null(scale.hC.G.D))
      } else{
        stopifnot(! is.null(nu.G.D))
      }
    }
  }
  if(T == 1 & is.infinite(err.df)){
    if(is.null(V.E)){
      stopifnot(! is.null(scale.hC.E))
    } else
      stopifnot(! is.matrix(V.E))
  }
  if(T > 1){
    if(is.null(V.E)){
      stopifnot(! is.null(nu.E))
    } else
      stopifnot(is.matrix(V.E))
  }
  stopifnot(perc.NA >= 0, perc.NA <= 100)
  if(! is.null(seed))
    set.seed(seed)
  stopifnot(! xor(is.null(nb.rows), is.null(nb.cols)))
  if(all(! is.null(nb.rows), ! is.null(nb.cols)))
    stopifnot(nb.rows * nb.cols == nrow(A))
  stopifnot(! xor(is.null(rho.rows), is.null(rho.cols)))
  stopifnot(! xor(is.null(rho.rows), is.null(sigma2.ksi)))
  if(all(! is.null(rho.rows), ! is.null(rho.cols), ! is.null(sigma2.ksi)))
    stopifnot(T == 1,
              all(! is.null(nb.rows), ! is.null(nb.cols)))

  ## determine the number of genotypes and phenotypes
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
  colnames(Z) <- levels.genos
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

  ## spatial coordinates
  coords <- NULL
  if(all(! is.null(nb.rows), ! is.null(nb.cols)))
    coords <- expand.grid(row=paste0("r", 1:nb.rows),
                          col=paste0("c", 1:nb.cols))

  ## spatially-dependent errors (not working yet)
  add.spatial.errors <- FALSE
  Z.s <- NULL
  if(all(! is.null(rho.rows), ! is.null(rho.cols), ! is.null(sigma2.ksi))){
    add.spatial.errors <- TRUE
    Z.s <- stats::model.matrix(~ coords$row + coords$col - 1)
    corr.mat.spatial <- corrMatAR1(nb.cols, rho.cols) %x% corrMatAR1(nb.rows,
                                                                     rho.rows)
    E.spatial <- matrix(MASS::mvrnorm(n=1, mu=rep(0, ncol(Z.s)),
                                      Sigma=sigma2.ksi * corr.mat.spatial))
  }

  ## independent errors
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
  if(add.spatial.errors)
    Y <- Y + Z %*% E.spatial
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
              data=dat,
              coords=coords))
}

##' Make elements of MME
##'
##' Make elements of Henderson's MME.
##' @param y vector of phenotypes (or matrix with one column)
##' @param X incidence matrix of fixed effects
##' @param Z incidence matrix of random variables
##' @return list of matrices
##' @author Timothee Flutre
##' @seealso \code{\link{solveMme}}
##' @export
makeMmeElements <- function(y, X, Z){
  if(is.vector(y))
    y <- matrix(y, nrow=length(y), ncol=1)
  stopifnot(is.matrix(y),
            is.matrix(X),
            is.matrix(Z),
            ncol(y) == 1,
            nrow(y) == nrow(X),
            nrow(X) == nrow(Z))

  ## for the left-hand side
  tX.X <- crossprod(X, X)
  tZ.X <- crossprod(Z, X)
  tX.Z <- crossprod(X, Z)
  tZ.Z <- crossprod(Z, Z)

  ## for the right-hand side
  tX.y <- crossprod(X, y)
  tZ.y <- crossprod(Z, y)

  return(list(tX.X=tX.X, tZ.X=tZ.X, tX.Z=tX.Z, tZ.Z=tZ.Z,
              tX.y=tX.y, tZ.y=tZ.y))
}

##' Make MME's right-hand side
##'
##' Make Henderson's MME's right-hand side.
##' @param tX.y X being the incidence matrix of fixed effects and y the vector of phenotypes
##' @param tZ.y Z being the incidence matrix of random variables
##' @return matrix
##' @author Timothee Flutre
##' @seealso \code{\link{solveMme}}
##' @export
makeMmeRhs <- function(tX.y, tZ.y){
  stopifnot(is.matrix(tX.y),
            ncol(tX.y) == 1,
            is.matrix(tZ.y),
            ncol(tZ.y) == 1)

  rhs <- c(tX.y, tZ.y)

  ## rhs <- matrix(data=NA, nrow=nrow(tX.y) + nrow(tZ.y), ncol=1)
  ## rhs[1:ncol(X), 1] <- tX.y
  ## rhs[(ncol(X)+1):nrow(rhs), 1] <- tZ.y

  return(rhs)
}

##' Make MME's left-hand side
##'
##' Make Henderson's MME's left-hand side.
##' @param tX.X X being the incidence matrix of fixed effects
##' @param tZ.X Z being the incidence matrix of random variables
##' @param tX.Z see above
##' @param tZ.Z see above
##' @param lambda ratio of the variance component of the errors over the variance component of the additive genotypic values (also called breeding values)
##' @param Ainv inverse of A, the matrix of additive genetic relationships
##' @return matrix
##' @author Timothee Flutre
##' @seealso \code{\link{solveMme}}
##' @export
makeMmeLhs <- function(tX.X, tZ.X, tX.Z, tZ.Z, lambda, Ainv){
  stopifnot(is.matrix(tX.X),
            is.matrix(tZ.X),
            is.matrix(tX.Z),
            is.matrix(tZ.Z),
            nrow(tX.X) == nrow(tX.Z),
            nrow(tZ.Z) == nrow(tZ.X),
            is.numeric(lambda),
            length(lambda) == 1,
            is.matrix(Ainv))

  p <- nrow(tX.X)
  q <- nrow(tZ.Z)
  lhs <- matrix(data=NA, nrow=p+q, ncol=p+q)

  lhs[1:p, 1:p] <- tX.X
  lhs[(p+1):(p+q), 1:p] <- tZ.X
  lhs[1:p, (p+1):(p+q)] <- tX.Z
  lhs[(p+1):(p+q), (p+1):(p+q)] <- tZ.Z + lambda * Ainv

  return(lhs)
}

##' Solve MME
##'
##' Given T=1 response, P fixed effects and Q random variables (for a total of N observations), compute the BLUEs and BLUPs of the following linear mixed model by solving Henderson's mixed model equation (MME): y = X beta + Z u + epsilon, where y is N x 1; X is N x P; Z is N x Q; u ~ Normal_Qx1(0, sigma_u^2 A); epsilon ~ Normal_Nx1(0, sigma^2 Id_N).
##' @param y vector of responses (or matrix with one column)
##' @param X incidence matrix of fixed effects
##' @param Z incidence matrix of random variables
##' @param sigma.u2 variance component of the random variables
##' @param Ainv inverse of A, the variance-covariance matrix of the random variables
##' @param sigma2 variance component of the errors
##' @return vector of length PQ containing the BLUEs of beta and the BLUPs of u
##' @author Timothee Flutre
##' @seealso \code{\link{makeMmeElements}}, \code{\link{makeMmeRhs}}, \code{\link{makeMmeLhs}}
##' @examples
##' \dontrun{## simulate genotypes
##' set.seed(1859)
##' X <- simulGenosDose(nb.genos=200, nb.snps=2000)
##'
##' ## estimate the additive genetic relationships
##' A <- estimGenRel(X, relationships="additive", method="vanraden1", verbose=0)
##'
##' ## simulate phenotypes
##' model <- simulAnimalModel(T=1, Q=3, A=A, V.G.A=15, V.E=5, seed=1859)
##'
##' ## calculate BLUEs and BLUPs
##' if(isSingular(A)){
##'    Ainv <- mpInv(A)
##' } else
##'    Ainv <- solve(A)
##' fit <- solveMme(y=model$Y[,1,drop=FALSE], X=model$W, Z=model$Z,
##'                 sigma.u2=model$V.G.A, Ainv=Ainv, sigma2=model$V.E)
##' cbind(model$C, fit[1:3])
##' cor(model$G.A, fit[4:length(fit)])
##' }
##' @export
solveMme <- function(y, X, Z, sigma.u2, Ainv, sigma2){
  if(is.vector(y))
    y <- matrix(y, nrow=length(y), ncol=1)
  stopifnot(is.matrix(y),
            is.matrix(X),
            is.matrix(Z),
            is.matrix(Ainv),
            ncol(y) == 1,
            nrow(y) == nrow(X),
            nrow(X) == nrow(Z),
            nrow(Ainv) == ncol(Z))

  elems <- makeMmeElements(y=y, X=X, Z=Z)

  rhs <- makeMmeRhs(tX.y=elems$tX.y, tZ.y=elems$tZ.y)

  lhs <- makeMmeLhs(tX.X=elems$tX.X, tZ.X=elems$tZ.X,
                    tX.Z=elems$tX.Z, tZ.Z=elems$tZ.Z,
                    lambda=sigma2 / sigma.u2, Ainv=Ainv)

  theta.hat <- solve(lhs, rhs) # faster than solve(lhs) %*% rhs

  return(as.vector(theta.hat))
}

##' Estimate variance components via REML using EM
##'
##' Given the MME components, estimate variance components by the EM algorithm.
##' References are:
##' \itemize{
##' \item Suzuki, M. 2007. Applied Animal Breeding & Genetics. Course Notes. Obihiro University of Agriculture and Veterinary Medicine.
##' \item Mrode, R.A. 2005. Linear Models for the Prediction of Animal Breeding Values. CAB International, Oxon, UK.
##' }
##' @param Ainv inverse of additive relationship matrix
##' @param y vector of responses
##' @param X incidence matrix for fixed effects
##' @param Z incidence matrix for random effects
##' @param initE initial value for the residual variance
##' @param initU initial value for the additive genetic variance
##' @param verbose verbosity level (0/1/2)
##' @return vector
##' @author Goto Morota [aut], Timothee Flutre [ctb]
##' @note Unknown parents should be coded as zero. Last modified by Morota on April 8, 2010.
##' @seealso \code{\link{makeMmeElements}}, \code{\link{makeMmeRhs}}, \code{\link{makeMmeLhs}}
##' @export
emreml <- function(Ainv, y, X, Z, initE, initU, verbose=0){
	n <- length(y)
	rankX <- qr(X, LAPACK=TRUE)$rank
	rankA <- qr(Ainv, LAPACK=TRUE)$rank

  elems <- makeMmeElements(y=y, X=X, Z=Z)
  Xpy <- elems$tX.y
  Zpy <- elems$tZ.y
  rhs <- makeMmeRhs(tX.y=elems$tX.y, tZ.y=elems$tZ.y)
  XpX <- elems$tX.X
  XpZ <- elems$tX.Z
  ZpX <- elems$tZ.X
  ZpZ <- elems$tZ.Z

	nb.rows.lhs <- length(rhs) # nrow(X) + nrow(Z)
	z <- length(rhs) - dim(ZpZ)[1] + 1

	oldE <- initE
	oldU <- initU
	diff1 <- 1
	diff2 <- 1
	i <- 0
	while(diff1 > 10^(-6) & diff2 > 10^(-6)){
		i <- i + 1

		lambda <- as.vector((oldE / oldU))
    lhs <- makeMmeLhs(tX.X=XpX, tZ.X=ZpX, tX.Z=XpZ, tZ.Z=ZpZ,
                      lambda=lambda, Ainv=Ainv)
		B <- solve(lhs, rhs) #faster than solve(lhs) %*% rhs
		e <-  y - cbind(X,Z) %*% B
		varE <- (t(e) %*% y) / (n - rankX)
		c22 <- solve(lhs)[z:nb.rows.lhs, z:nb.rows.lhs]
		u <- B[z:length(B)]
		## sum(Ainv*c22) is same as sum(diag(Ainv%*%c22))
		varU <- (t(u) %*% Ainv %*% u + sum(Ainv*c22)*oldE) / (rankA)

		diff1 <- abs(varE - oldE)
		diff2 <- abs(varU - oldU)
		if(verbose > 0){
      txt <- paste0("iteration ", i, ":",
                    " varE=", format(varE, digits=5),
                    " varU=", format(varU, digits=5))
      write(txt, stdout())
		}
		oldE <- varE
		oldU <- varU
	}

	return(c(varE=varE, varU=varU))
}

##' Estimate variance components via REML using AI
##'
##' Given the MME components, estimate variance components by the AI algorithm.
##' References are:
##' \itemize{
##' \item Tsuruta, S. 2006. Estimation of Variance Components in Animal Breeding. The University of Georgia.
##' \item Mrode, R.A. 2005. Linear Models for the Prediction of Animal Breeding Values. CAB International, Oxon, UK.
##' }
##' @param A additive relationship matrix
##' @param y vector of responses
##' @param X incidence matrix for fixed effects
##' @param Z incidence matrix for random effects
##' @param initE initial value for the residual variance
##' @param initU initial value for the additive genetic variance
##' @param verbose verbosity level (0/1/2)
##' @return vector
##' @author Goto Morota
##' @note Last modified by Morota on April 8, 2010.
##' @export
aireml <- function(A, y, X, Z, initE, initU, verbose=0){
	N <- length(y)
	Ze <- diag(N)
	oldE <- initE
	oldU <- initU
	AI <- matrix(0, ncol=2, nrow=2) # average information matrix
	s <- matrix(0, ncol=1, nrow=2)  # score function (first derivative of the log-likelihood)

	diff1 <- 1
	diff2 <- 1
	i <- 0
	while(diff1 > 10^(-6) & diff2 > 10^(-7)){
    i <- i + 1

    G <- oldU * A
    R <- oldE * Ze %*% t(Ze)
    V <- Z %*% G %*% t(Z) + R
    Vinv <- solve(V)
    P <- Vinv - Vinv %*% X %*% solve(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv
    ## TODO -> compute log-likelihood:
    ## eq 3 in Johnson & Thompson: l propto -(1/2)[log|V| + log|X' V^-1 X| + y' P y]
    ## see also Mrode (2005) page 240

    AI[1,1] <- sum(diag((t(y) %*% P %*% Z %*% t(Z) %*%
                         P %*% Z %*% t(Z) %*% P %*% y)))
    AI[1,2] <-  sum(diag((t(y) %*% P %*% Z %*% t(Z) %*%
                          P %*% Ze %*% t(Ze) %*% P %*% y)))
    AI[2,1] <- AI[1,2]
    AI[2,2] <- sum(diag((t(y) %*% P %*% Ze %*% t(Ze) %*%
                         P %*% Ze %*% t(Ze) %*% P %*% y)))

    ## ~ eq 92 p.252 in Searle's book
    ## ~ eq 4 in Johnson & Thompson
    s[1,1] <- - sum(diag((P %*% Z %*% t(Z)))) +
      (t(y) %*% P %*% Z %*% t(Z) %*% P %*% y)
    s[2,1] <- - sum(diag((P %*% Ze %*% t(Ze)))) +
      (t(y) %*% P %*% Ze %*% t(Ze) %*% P %*% y)

    ## update of Fisher scoring
    oldvarcomps <- c(oldU, oldE)
    newvarcomps <- oldvarcomps + solve(AI) %*% s
    varU <- newvarcomps[1]
    varE <- newvarcomps[2]

    diff1 <- abs(varE - oldE)
    diff2 <- abs(varU - oldU)
    if(verbose > 0){
      txt <- paste0("iteration ", i, ":",
                    " varE=", format(varE, digits=5),
                    " varU=", format(varU, digits=5))
      write(txt, stdout())
		}
    oldE <- varE
    oldU <- varU
	}

  return(c(varE=varE, varU=varU))
}

## not exported: used only in plantTrialLmmFitCompSel
plantTrialLmmFitFixed <- function(glob.form, dat.noNA, ctl=NULL,
                                  saved.file=NULL,
                                  cl=NULL, nb.cores=1,
                                  verbose=0){
  requireNamespace("MuMIn", quietly=TRUE)
  requireNamespace("parallel", quietly=TRUE)
  stopifnot(is.character(glob.form))
  glob.preds <- trimws(strsplit(glob.form, "\\~|\\+")[[1]])[-1]
  any.rand.var <- any(grepl("\\(*\\)", glob.preds))
  if(any.rand.var)
    requireNamespace("lme4", quietly=TRUE)

  ## determine whether to fit or to load already-computed results
  need.to.fit <- is.null(saved.file)
  if(! is.null(saved.file)){
    if(! file.exists(saved.file))
      need.to.fit <- TRUE
  }

  if(need.to.fit){
    if(verbose > 0){
      msg <- paste0("fit all models with ML",
                    ifelse(any.rand.var, " (lme4 and MuMIn)", " (lm)"),
                    " on ", nb.cores, " core",
                    ifelse(nb.cores <= 1, "", "s"), "...")
      write(msg, stdout())
    }

    if(any.rand.var){
      if(is.null(ctl))
        ctl <- lme4::lmerControl()
      globmod.ml <- lme4::lmer(formula=stats::as.formula(glob.form),
                               data=dat.noNA,
                               na.action="na.fail",
                               REML=FALSE, control=ctl)
    } else
      globmod.ml <- stats::lm(formula=stats::as.formula(glob.form),
                              data=dat.noNA,
                              na.action="na.fail")
    ## print(summary(globmod.ml))

    if(nb.cores <= 1){
      st <- system.time(
        allmod.sel <- suppressMessages(
          MuMIn::dredge(global.model=globmod.ml, rank="AICc")))
    } else{
      should.cl.be.created <- FALSE
      if(is.null(cl)){
        should.cl.be.created <- TRUE
        cl <- parallel::makeCluster(nb.cores, "SOCK")
      }
      parallel::clusterExport(cl=cl, varlist="dat.noNA", envir=environment())
      if(any.rand.var)
        tmp <- parallel::clusterEvalQ(cl=cl, expr=library(lme4))
      st <- system.time(
        allmod.sel <- suppressMessages(
          MuMIn::pdredge(global.model=globmod.ml, cluster=cl, rank="AICc")))
      if(should.cl.be.created)
        parallel::stopCluster(cl)
    } # end of nb.cores > 1
    if(verbose > 0)
      print(st)
    if(! is.null(saved.file)){
      if(verbose > 0){
        msg <- paste0("save all model fits with ML in '", saved.file, "':")
        write(msg, stdout())
      }
      save(allmod.sel, file=saved.file)
      md5.allmod.sel <- tools::md5sum(path.expand(saved.file))
    }
  } else{ # if need.to.fit is FALSE
    if(verbose > 0){
      msg <- paste0("load all model fits with ML from '", saved.file, "':")
      write(msg, stdout())
    }
    md5.allmod.sel <- tools::md5sum(path.expand(saved.file))
    load(saved.file)
  }
  if(verbose > 0){
    if(! is.null(saved.file)){
      msg <- paste0("MD5 sum: ", md5.allmod.sel)
      write(msg, stdout())
    }
  }

  return(allmod.sel)
}

##' Model fit, comparison and selection
##'
##' For a plant field trial, fit by maximum likelihood (ML) a global linear (mixed) model with all the specified terms, then fit by ML various sub-models and compare them (AICc) or test each term (F-test for fixed effects and likelihood ratio test for random variables), finally select the best model and, if there are random variables, re-fit it using restricted maximum likelihood (ReML).
##'
##' @param glob.form character containing the formula for the global model
##' @param dat data frame with all the columns required by \code{glob.form}; if the response contains an inline function (e.g., log or sqrt), the untransformed response may still be required
##' @param part.comp.sel part(s) of the model (fixed and/or random); if only "fixed", the lm/lme4 and MuMIn packages will be used, otherwise the lmerTest package will be used
##' @param alpha.fixed for lmerTest, threshold on p values of fixed effects below which they are kept
##' @param alpha.random for lmerTest, threshold on p values of random effects below which they are kept
##' @param ctl if NULL, will be the output of \code{lmerControl} with default arguments
##' @param saved.file name of the file in which to save all model fits or from which to load all model fits (for MuMIn); for lmerTest, only the fit of the global model will be saved
##' @param nb.cores number of cores; ignored with lmerTest
##' @param cl an object of class "cluster"; if NULL, will be created automaticall based on \code{nb.cores}; ignored with lmerTest
##' @param verbose verbosity level (0/1)
##' @return list
##' @export
plantTrialLmmFitCompSel <- function(glob.form, dat, part.comp.sel="fixed",
                                    alpha.fixed=0.05, alpha.random=0.1,
                                    ctl=NULL,
                                    saved.file=NULL,
                                    nb.cores=1, cl=NULL, verbose=1){
  stopifnot(is.character(glob.form),
            is.data.frame(dat),
            all(part.comp.sel %in% c("fixed", "random")))
  glob.preds <- trimws(strsplit(glob.form, "\\~|\\+")[[1]])[-1]
  any.rand.var <- any(grepl("\\(*\\)", glob.preds))
  if(any.rand.var){
    requireNamespace("lme4", quietly=TRUE)
  } else if("random" %in% part.comp.sel){
    part.comp.sel <- part.comp.sel[-grep("random", part.comp.sel)]
    stopifnot(length(part.comp.sel) > 0)
  }
  if("random" %in% part.comp.sel){
    requireNamespace("lmerTest", quietly=TRUE)
  } else{
    requireNamespace("MuMIn", quietly=TRUE)
    requireNamespace("tools", quietly=TRUE)
    if(nb.cores > 1)
      requireNamespace("parallel", quietly=TRUE)
  }
  if(! is.null(cl))
    stopifnot(methods::is(cl, "cluster"))
  inFctResp <- inlineFctForm(glob.form, only.resp=TRUE)
  if(! all(is.na(inFctResp[[1]]))){
    orig.resp <- inFctResp[[1]][1]
    if(! orig.resp %in% colnames(dat)){
      msg <- paste0("the response '", names(inFctResp)[1], "' is transformed,",
                    " i.e. contains an inline function,",
                    " but the untransformed response '", orig.resp,
                    "' isn't present as column names in 'dat'")
      warning(msg, immediate.=TRUE)
    }
  }

  out <- list()

  ## format input data
  dat <- droplevels(dat)
  dat.noNA <- dat
  if(any(is.na(dat.noNA))){
    if(verbose > 0)
      write("exclude NAs...", stdout())
    dat.noNA <- droplevels(stats::na.exclude(dat.noNA))
  }
  out[["dat"]] <- dat
  out[["dat.noNA"]] <- dat.noNA

  ## fit, compare and select models
  if("random" %in% part.comp.sel){ # LMM and select fix+rand terms
    if(verbose > 0){
      msg <- "fit the global model with ML (lmerTest)"
      write(msg, stdout())
    }
    if(is.null(ctl))
      ctl <- lme4::lmerControl()
    globmod.ml <- do.call(lmerTest::lmer,
                          list(formula=stats::as.formula(glob.form),
                               data=dat.noNA,
                               na.action="na.fail",
                               REML=FALSE,
                               control=ctl))
    if(! is.null(saved.file)){
      if(verbose > 0){
        msg <- paste0("save the global model fit with ML in '", saved.file, "':")
        write(msg, stdout())
      }
      save(globmod.ml, file=saved.file)
    }
    step_res <- lmerTest::step(object=globmod.ml,
                               reduce.fixed="fixed" %in% part.comp.sel,
                               alpha.fixed=alpha.fixed,
                               reduce.random=TRUE,
                               alpha.random=alpha.random)
    bestmod.ml <- lmerTest::get_model(step_res)
  } else{ # LM or LMM and select fix terms
    allmod.sel <- plantTrialLmmFitFixed(glob.form, dat.noNA, ctl,
                                        saved.file,
                                        cl, nb.cores, verbose)
    bestmod.ml <- MuMIn::get.models(allmod.sel, subset=1)[[1]]
  }
  if(verbose > 0){
    msg <- "best model with ML:"
    write(msg, stdout())
    print(stats::formula(bestmod.ml, showEnv=FALSE))
  }
  out[["bestmod.ml"]] <- bestmod.ml
  ## best.form <- Reduce(paste, deparse(stats::formula(bestmod.ml)))
  ## all.preds <- trimws(strsplit(best.form, "\\~|\\+")[[1]])[-1]
  best.form <- stats::formula(bestmod.ml)
  best.preds <- attributes(stats::terms(best.form))$term.labels
  if(any(grep("\\|", best.preds))){
    best.preds.fix <- best.preds[- grep("\\|", best.preds)]
  } else
    best.preds.fix <- best.preds
  best.preds.rnd <- best.preds[grep("\\|", best.preds)]
  best.preds.rnd <- trimws(gsub("1|\\|", "", best.preds.rnd))
  out[["best.form"]] <- best.form
  out[["best.preds"]] <- best.preds
  out[["best.preds.fix"]] <- best.preds.fix
  out[["best.preds.rnd"]] <- best.preds.rnd

  out[["bestmod.reml"]] <- NA
  if(length(best.preds.rnd) == 0){
    warning("can't refit with ReML because no random effect was kept")
  } else{
    if(verbose > 0){
      msg <- "re-fit the best model with ReML..."
      write(msg, stdout())
    }
    if(is.null(ctl))
      ctl <- lme4::lmerControl()
    bestmod.reml <- lme4::lmer(formula=stats::formula(bestmod.ml),
                               data=dat.noNA,
                               na.action="na.fail",
                               REML=TRUE,
                               control=ctl)
    out[["bestmod.reml"]] <- bestmod.reml
  }

  return(out)
}

##' Plot residuals between years
##'
##' Plot residuals between years to check their temporal independence.
##' Especially useful as a diagnostic of model fit for statistical analysis of perennial plants.
##' @param df data frame with a column given by \code{colname.res}, a column named "year", a column named "block" if \code{blocks} is specified, as well as any other column specified in \code{cols.uniq.id}
##' @param colname.res name of the column containing the residuals
##' @param years vector with two years as character
##' @param cols.uniq.id vector of column name(s) allowing to identify each plant and pair it between years or years and blocks, e.g. \code{c("geno","block","rank","location")} or \code{"geno"}
##' @param blocks if not NULL, vector of two blocks (e.g. \code{c("A","A")} or \code{c("A","B")})
##' @param las style of axis labels, see \code{\link[graphics]{par}}
##' @param lgd.pos position of the legend
##' @param ... arguments passed on to \code{\link[graphics]{plot}}, such as \code{main}
##' @return invisible data frame of the data used to make the plot
##' @author Timothee Flutre
##' @export
plotResidualsBtwYears <- function(df, colname.res, years, cols.uniq.id,
                                  blocks=NULL, las=1, lgd.pos="topright",
                                  ...){
  stopifnot(is.data.frame(df),
            all(c(colname.res, "year", cols.uniq.id) %in% colnames(df)))
  if(! is.null(blocks))
    stopifnot("block" %in% colnames(df))
  if("year" %in% cols.uniq.id)
    cols.uniq.id <- cols.uniq.id[-grep("year", cols.uniq.id)]
  stopifnot(length(cols.uniq.id) > 0)

  ## extract the points of interest, insuring coherence between the two years
  if(is.null(blocks)){
    idx.x <- which(df$year == years[1])
    idx.y <- which(df$year == years[2])
  } else{
    idx.x <- which(df$year == years[1] & df$block == blocks[1])
    idx.y <- which(df$year == years[2] & df$block == blocks[2])
  }
  names(idx.x) <- do.call(paste, c(df[idx.x, cols.uniq.id, drop=FALSE],
                                   sep="_"))
  names(idx.y) <- do.call(paste, c(df[idx.y, cols.uniq.id, drop=FALSE],
                                   sep="_"))
  common.names <- names(idx.x)[names(idx.x) %in% names(idx.y)]
  stopifnot(length(common.names) > 0)
  idx.x <- idx.x[common.names]
  idx.y <- idx.y[common.names]
  tmp <- data.frame(x=df[idx.x, colname.res],
                    y=df[idx.y, colname.res])
  tmp <- tmp[stats::complete.cases(tmp),]
  if(nrow(tmp) == 0){
    msg <- paste0("no non-missing data on common plants for years ",
                  years[1], " and ", years[2])
    stop(msg)
  }

  ## make the plot
  if(is.null(blocks)){
    xlab <- years[1]
    ylab <- years[2]
  } else{
    xlab <- paste(blocks[1], years[1])
    ylab <- paste(blocks[2], years[2])
  }
  l <- max(abs(c(tmp$x, tmp$y)))
  graphics::plot(formula=y ~ x, data=tmp, las=las,
                 xlim=c(-l,l), ylim=c(-l,l),
                 xlab=xlab, ylab=ylab, ...)
  graphics::abline(v=0, lty=2)
  graphics::abline(h=0, lty=2)

  ## add information about correlation and R squared
  cor.p <- stats::cor(x=tmp$x, y=tmp$y, method="pearson")

  fit.lm <- stats::lm(y ~ x, data=tmp)
  R2 <- summary(fit.lm)$r.squared
  tmp$fitted.lm <- stats::fitted(fit.lm)
  graphics::lines(x=tmp$x[order(tmp$x)],
                  y=tmp$fitted.lm[order(tmp$x)],
                  col="green", lty=1, lwd=2)

  fit.loess <- stats::loess(y ~ x, data=tmp)
  pseudoR2 <- pseudoR2(dat=tmp$y, res=stats::residuals(fit.loess),
                       pred=stats::predict(fit.loess), method="Efron")
  tmp$fitted.loess <- stats::fitted(fit.loess)
  graphics::lines(x=tmp$x[order(tmp$x)],
                  y=tmp$fitted.loess[order(tmp$x)],
                  col="red", lty=1, lwd=2)

  lgd <- c(bquote("Pearson corr." == .(round(cor.p, 2))),
           bquote("lm, " ~ R^2 ==
                    .(round(R2, 2))),
           bquote("loess, pseudo" ~ R^2 ==
                    .(round(pseudoR2, 2))))
  graphics::legend(lgd.pos,
                   legend=sapply(lgd, as.expression),
                   col=c("0", "green", "red"),
                   lty=c(0, 1, 1),
                   lwd=c(0, 2, 2),
                   bty="n")

  invisible(tmp)
}

##' Broad-sense heritability
##'
##' Estimate broad-sense heritability (squared correlation between predicted and true genotypic effects) on an entry-mean basis:
##' \enumerate{
##' \item via the classical formula for balanced data sets (see Falconer and Mackay, or the introduction of \href{http://www.genetics.org/cgi/doi/10.1534/genetics.107.074229}{Piepho and Mohring (2007)}): H2 = var.g / var.p, where var.p = var.g + var.ge / m + var.e / (r x m) with "m" the number of trials and "r" the number of replicates per trial (for unbalanced data sets, the mean number of trials and replicates per trial are used);
##' \item via the formula of \href{http://dx.doi.org/10.1007/s00122-006-0333-z}{Oakey et al (2006)} for unbalanced data sets: H2 = 1 - trace(G^-1 C_zz) / m, with "m" the number of genotypes.
##' }
##' @param dat data frame of input data after missing data have been excluded, with columns named \code{colname.resp}, "geno" and \code{colname.trial}
##' @param colname.resp name of the column containing the response
##' @param colname.trial name of the column identifying the trials (e.g. \code{"year"}, \code{"year_irrigation"}, etc)
##' @param vc data frame of variance components with columns "grp" and "vcov" (i.e. formatted as \code{as.data.frame(VarCorr())} from the "lme4" package); \code{grp="Residual"} for "var.e"
##' @param geno.var.blups vector of variances of empirical BLUPs of the genotypic effects, g, assuming g ~ MVN(0, G) where G = sigma_g^2 I_m; if not provided, the estimator of Oakey et al won't be computed
##' @return list with the mean number of trials, the mean number of replicates per trial, the broad-sense heritability (classical estimator from Falconer and Mackay, as well as optionally the one from Oakey et al), and a function to compute summary statistics whch can be used for estimating confidence intervals by bootstrap
##' @author Timothee Flutre
##' @export
estimH2means <- function(dat, colname.resp, colname.trial="year", vc,
                         geno.var.blups=NULL){
  requireNamespace("lme4", quietly=TRUE)
  stopifnot(is.data.frame(dat),
            all(! is.na(dat)),
            all(c(colname.resp,"geno",colname.trial) %in% colnames(dat)),
            is.data.frame(vc),
            all(c("grp","vcov") %in% colnames(vc)),
            "geno" %in% vc$grp,
            "Residual" %in% vc$grp)
  if(! is.null(geno.var.blups))
    stopifnot(is.vector(geno.var.blups))

  out <- list()

  ## classical estimator from Falconer and Mackay
  reps.geno.trial <- tapply(dat[[colname.resp]],
                           list(dat[["geno"]], dat[[colname.trial]]),
                           length)
  ## out$reps.geno.trial <- reps.geno.trial # debug
  mean.nb.trials <- mean(apply(reps.geno.trial, 1, function(x){
    sum(! is.na(x))
  }))
  out$mean.nb.trials <- mean.nb.trials
  mean.nb.reps.per.trial <- mean(apply(reps.geno.trial, 2, mean, na.rm=TRUE))
  out$mean.nb.reps.per.trial <- mean.nb.reps.per.trial

  var.geno <- vc[vc$grp == "geno", "vcov"]
  var.pheno <- var.geno
  if(paste0("geno:", colname.trial) %in% vc$grp)
    var.pheno <- var.pheno +
      vc[vc$grp == paste0("geno:", colname.trial), "vcov"] /
      mean.nb.trials
  var.pheno <- var.pheno +
    vc[vc$grp == "Residual", "vcov"] /
    (mean.nb.trials * mean.nb.reps.per.trial)
  H2.classic <- var.geno / var.pheno
  out$H2.classic <- H2.classic

  ## estimator from Oakey et al (2006)
  G <- vc[vc$grp == "geno", "vcov"] *
    diag(length(geno.var.blups))
  Ginv <- solve(G)
  C.zz <- diag(geno.var.blups)
  H2.oakey <- 1 - matrixTrace(Ginv %*% C.zz) / length(geno.var.blups)
  out$H2.oakey <- H2.oakey

  sryStat <- function(.){
    geno.blups <- lme4::ranef(., condVar=TRUE, drop=TRUE)$geno
    geno.var.blups <- stats::setNames(attr(geno.blups, "postVar"),
                                      names(geno.blups))
    tmp <- c(ef=lme4::fixef(.),
             sd.err=stats::sigma(.),
             sd=sqrt(unlist(lme4::VarCorr(.))))
    var.geno <- tmp["sd.geno"]^2
    var.pheno <- var.geno
    if(paste0("sd.geno:", colname.trial) %in% names(tmp))
      var.pheno <- var.pheno +
        tmp[paste0("sd.geno:", colname.trial)]^2 /
        mean.nb.trials
    var.pheno <- var.pheno +
      tmp["sd.err"]^2 / (mean.nb.trials * mean.nb.reps.per.trial)
    tmp <- c(tmp,
             var.geno / var.pheno)
    names(tmp)[length(tmp)] <- "H2.classic"
    G <- tmp["sd.geno"]^2 * diag(length(geno.var.blups))
    Ginv <- solve(G)
    C.zz <- diag(geno.var.blups)
    tmp <- c(tmp,
             1 - matrixTrace(Ginv %*% C.zz) / length(geno.var.blups))
    names(tmp)[length(tmp)] <- "H2.oakey"
    tmp <- c(tmp,
             tmp["sd.geno"] / abs(tmp["ef.(Intercept)"]))
    names(tmp)[length(tmp)] <- "CV.geno"
    return(tmp)
  }
  out$sryStat <- sryStat

  return(out)
}

##' Animal model
##'
##' Given I genotypes, Q covariates and N=I*Q phenotypes for the trait, fit an "animal model" with the lme4 package via the following likelihood: y = W c + Z g_A + Z g_D + epsilon, where y is Nx1; W is NxQ; Z is NxI; g_A ~ Normal_I(0, sigma_A^2 A) with A the known matrix of additive genetic relationships; g_D ~ Normal_I(0, sigma_D^2 D) with D the known matrix of dominant genetic relationships; epsilon ~ Normal_N(0, sigma^2 Id_N); Cov(g_A,g_D)=0; Cov(g_A,e)=0; Cov(g_D,e)=0.
##' Works also for an incomplete design, i.e., when the number of replicates per genotype is unbalanced across genotypes and nrow(data) < N (thanks to N.O. Rode).
##' @param formula formula (see \code{\link[lme4]{lmer}})
##' @param data data.frame containing the data corresponding to formula and relmat (see \code{\link[lme4]{lmer}}); the additive genotypic effect should be named "geno.add" and the dominance genotypic effect, if any, should be named "geno.dom"
##' @param relmat list containing the matrices of genetic relationships (A is compulsory but D is optional); the list should use the same names as the colnames in data (i.e., \code{"geno.add"} and \code{"geno.dom"}) to compute heritability properly; the matrices can be in the "matrix" class (base) or the "dsCMatrix" class (Matrix package); see \code{\link{estimGenRel}}
##' @param REML default is TRUE (use FALSE to compare models with different fixed effects)
##' @param na.action a function that indicates what should happen when the data contain \code{NA}s (e.g. to compute the BLUp of unobserved genotype, use na.action=NULL, see \code{\link[lme4]{lmer}})
##' @param nrep a number or a vector that indicates the number of replicates per genotype if line phenotypic means are provided in data (used to compute individual-level heritability, default = NULL)
##' @param ci.meth method to compute confidence intervals (profile/boot)
##' @param ci.lev level to compute confidence intervals
##' @param nb.boots number of bootstrap replicates; used only if \code{ci.meth="boot"}
##' @param parallel the type of parallel operation to be used, if any (no/multicore/snow)
##' @param ncpus number of processes to be used in parallel operation
##' @param cl an optional object of class \code{"cluster"} returned by \code{\link[parallel]{makeCluster}}; if not supplied, a cluster on the local machine is created for the duration of the call
##' @param verbose verbosity level (0/1)
##' @return list with a \code{\link[lme4]{merMod}} object, a vector of point estimates of variance components, a point estimate of h2, a \code{thpr} object if \code{ci.meth="profile"}, a \code{boot} object if \code{ci.meth="profile"}, a data frame with confidence intervals (if \code{ci.meth} is not NULL)
##' @author Timothee Flutre (inspired by Ben Bolker at http://stackoverflow.com/q/19327088/597069)
##' @note If A is not positive definite, an error will be raised (via \code{\link[base]{chol}}); in such cases, using the \code{nearPD} function from the Matrix package can be useful.
##' @seealso \code{\link{inlaAM}}, \code{\link{jagsAM}}, \code{\link{stanAM}}
##' @examples
##' \dontrun{## simulate genotypes
##' set.seed(1859)
##' X <- simulGenosDose(nb.genos=200, nb.snps=2000)
##'
##' ## simulate phenotypes with only additive part of genotypic values
##' A <- estimGenRel(X, relationships="additive", method="vanraden1", verbose=0)
##' kappa(A)
##' A <- as.matrix(nearPD(A)$mat)
##' kappa(A)
##' modelA <- simulAnimalModel(T=1, Q=3, A=A, V.G.A=15, V.E=5, seed=1859)
##'
##' ## infer with lme4
##' modelA$data$geno.add <- modelA$data$geno
##' fitA <- lmerAM(formula=response1 ~ year + (1|geno.add),
##'                data=modelA$data,
##'                relmat=list(geno.add=A), verbose=0)
##' summary(fitA$merMod)
##' REMLcrit(fitA$merMod)
##' extractAIC(fitA$merMod)
##' summary(residuals(fitA$merMod)) # "deviance residuals"
##' summary(residuals(fitA$merMod) / sigma(fitA$merMod)) # "scaled/Pearson residuals"
##' c(modelA$C); modelA$V.G.A; modelA$V.E
##' fixef(fitA$merMod)
##' coefficients(summary(fitA$merMod))[, "Std. Error"]
##' vc <- as.data.frame(VarCorr(fitA$merMod))
##' c(vc[vc$grp == "geno.add", "vcov"], vc[vc$grp == "Residual", "vcov"])
##' blups.geno <- ranef(fitA$merMod, condVar=TRUE, drop=TRUE)$geno.add
##' var.blups.geno <- setNames(attr(blups.geno, "postVar"), names(blups.geno))
##'
##' ## compute confidence intervals in parallel
##' library(parallel)
##' (nb.cores <- max(1, detectCores() - 1))
##' cl <- makeCluster(spec=nb.cores, type="PSOCK")
##' fitA <- lmerAM(formula=response1 ~ year + (1|geno.add),
##'                data=modelA$data,
##'                relmat=list(geno.add=A), verbose=1,
##'                ci.meth="profile", parallel="snow", ncpus=nb.cores, cl=cl)
##' stopCluster(cl)
##'
##' ## simulate phenotypes with additive and dominant parts of genotypic values
##' D <- estimGenRel(X, relationships="dominant", method="vitezica", verbose=0)
##' kappa(D)
##' modelAD <- simulAnimalModel(T=1, Q=3, A=A, V.G.A=15, V.E=5,
##'                             D=D, V.G.D=3, seed=1859)
##'
##' ## infer with lme4
##' modelAD$data$geno.add <- modelAD$data$geno
##' modelAD$data$geno.dom <- modelAD$data$geno; modelAD$data$geno <- NULL
##' fitAD <- lmerAM(formula=response1 ~ year + (1|geno.add) + (1|geno.dom),
##'                 data=modelAD$data, relmat=list(geno.add=A, geno.dom=D),
##'                 verbose=0)
##' summary(fitAD$merMod)
##' c(modelAD$C); modelAD$V.G.A; modelAD$V.E; modelAD$V.G.D
##' fixef(fitAD$merMod)
##' vc <- as.data.frame(VarCorr(fitAD$merMod))
##' c(vc[vc$grp == "geno.add", "vcov"], vc[vc$grp == "Residual", "vcov"],
##'   vc[vc$grp == "geno.dom", "vcov"])
##' }
##' @export
lmerAM <- function(formula, data, relmat, REML=TRUE, na.action=stats::na.exclude,
		   nrep=NULL,
                   ci.meth=NULL, ci.lev=0.95, nb.boots=10^3,
                   parallel="no", ncpus=1, cl=NULL, verbose=1){
  requireNamespaces(c("lme4", "Matrix"))
  stopifnot(is.data.frame(data),
            all(! duplicated(colnames(data))),
            is.list(relmat),
            ! is.null(names(relmat)),
            "geno.add" %in% names(relmat),
            all(names(relmat) %in% colnames(data)),
            is.logical(REML))
  for(i in seq_along(relmat))
    stopifnot(is.matrix(relmat[[i]]) ||
              Matrix::isSymmetric(relmat[[i]]), # if nearPD()
              ! is.null(rownames(relmat[[i]])),
              ! is.null(colnames(relmat[[i]])),
              rownames(relmat[[i]]) == colnames(relmat[[i]]),
	      !is.data.frame(data[,names(relmat)[i]]),
	      all(data[,names(relmat)[i]] %in% rownames(relmat[[i]])))
  if(! is.null(ci.meth)){
    stopifnot(ci.meth %in% c("profile", "boot"))
    if(ci.meth == "boot"){
      requireNamespace("boot", quietly=TRUE)
      requireNamespace("MASS", quietly=TRUE)
    }
  }
  if(! is.null(parallel))
    stopifnot(parallel %in% c("no",  "multicore", "snow"))

  if(verbose > 0)
    write("parse the formula ...", stdout())
  parsedFormula <- lme4::lFormula(formula=formula,
                                  data=data,
                                  control=lme4::lmerControl(
                                      check.nobs.vs.nlev="ignore",
                                      check.nobs.vs.nRE="ignore"),
                                  na.action=na.action,
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
    relfac[[i]] <- chol(relmat[[i]], pivot=isSingular(relmat[[i]]))
    Ztlist[[i]] <- relfac[[i]] %*% Ztlist[[i]]
  }
  parsedFormula$reTrms[["Ztlist"]] <- Ztlist
  parsedFormula$reTrms[["Zt"]] <- do.call(rbind, Ztlist) # Matrix::rBind if R < 3.2.0

  if(verbose > 0)
    write("make the deviance function ...", stdout())
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

  ## point estimates of variance components
  vc <- as.data.frame(lme4::VarCorr(fit))
  vc <- stats::setNames(vc$vcov, vc$grp)

  ## point estimate of h2
  num <- as.numeric(vc["geno.add"])
  denom <- vc["geno.add"] + vc["Residual"]
  if("geno.dom" %in% names(vc))
    denom <- denom + vc["geno.dom"]
  h2 <- num / denom

  ## point estimate of h2 at the individual level
  if(!is.null(nrep)){
	  denom_ind <- vc["geno.add"] + nrep * vc["Residual"]
	  if("geno.dom" %in% names(vc))
	     denom_ind <- denom_ind + vc["geno.dom"]
	  h2_ind <- num / denom_ind
	  }

  out.prof <- NULL
  out.boot <- NULL
  ci <- NULL
  if(! is.null(ci.meth)){
    if(verbose > 0){
      msg <- paste0("compute confidence intervals (", ci.meth, ")...")
      write(msg, stdout())
    }
    if(ci.meth == "profile"){
      out.prof <- stats::profile(fitted=fit, signames=FALSE,
                                 parallel=parallel, ncpus=ncpus, cl=cl)
      ci <- stats::confint(object=out.prof, level=ci.lev)
      ## suppressMessages(ci <- lme4::confint.merMod(fit, level=ci.lev,
      ##                                             method="profile",
      ##                                             oldNames=FALSE,
      ##                                             parallel=parallel,
      ##                                             ncpus=ncpus, cl=cl))
    } else if(ci.meth == "boot"){
      .bootStat <- function(inputs){
        fit.boot <- lmerAM(formula=inputs$formula, data=inputs$data,
                           relmat=inputs$relmat, REML=inputs$REML,
                           na.action=inputs$na.action, ci.meth=NULL,
                           verbose=0, nrep=inputs$nrep)
        stats.boot <- stats::setNames(
                                 c(lme4::fixef(fit.boot$merMod),
                                   sqrt(fit.boot$vc["geno.add"]),
                                   sqrt(fit.boot$vc["Residual"]),
				   variance=unlist(lme4::VarCorr(fit.boot$merMod)),
				   residual=sigma(fit.boot$merMod)^2),
                                 c(names(lme4::fixef(fit.boot$merMod)),
                                   "sd.geno.add",
                                   "sd.err",
				   names(lme4::VarCorr(fit.boot$merMod)),
				  "ResidualVar"))
        if("geno.dom" %in% names(fit.boot$vc)){
          stats.boot <- c(stats.boot,
                          sqrt(fit.boot$vc["geno.dom"]))
          names(stats.boot)[length(stats.boot)] <- "sd.geno.dom"
        }
        stats.boot <- c(stats.boot,
                        fit.boot$h2)
	names(stats.boot)[length(stats.boot)] <- "h2"
	if(!is.null(nrep)){
		stats.boot <- c(stats.boot,
                        fit.boot$h2_ind)
		names(stats.boot)[length(stats.boot)] <- paste0("h2_ind(nrep=", nrep, ")") 
		} 
        return(stats.boot)
      }
      .bootRanGen <- function(inputs, params){
        outputs <- inputs

        ## parse formula
        tmp <- stats::terms(outputs$formula)
        vars <- as.character(attr(tmp, "variables"))[-1]
        idx.resp <- attr(tmp, "response")
        response <- vars[idx.resp]
        vars <- vars[-idx.resp]
        idx.geno.add <- grep("geno.add", vars)
        vars <- vars[-idx.geno.add]
        idx.geno.dom <- grep("geno.dom", vars)
        if(length(idx.geno.dom) == 1){
          vars <- vars[-idx.geno.dom]
        } else
          idx.geno.dom <- NULL
        if(length(vars) == 0){
          fix <- NULL
        } else
          fix <- vars

        ## get and check model dimensions
        outputs$data$geno.add <- as.factor(outputs$data$geno.add)
        I <- nlevels(outputs$data$geno.add)
        Q <- length(params$fix.eff)
        N <- I * Q
        if(nrow(outputs$data) < N){
          print("note that the design is incomplete")
        } else
          print("note that the design is complete")

        ## make design matrices
        fm <- paste(c("~ 1", fix), collapse="+")
        X <- stats::model.matrix(stats::as.formula(fm),
                                 data=outputs$data)
        stopifnot(ncol(X) == length(params$fix.eff))
        Z <- stats::model.matrix(~ -1 + geno.add,
                                 data=outputs$data)
        stopifnot(ncol(Z) == I)

        ## draw random variables and generate responses
        e <- stats::rnorm(n=nrow(outputs$data), mean=0, sd=params$sd.err)
        g.a <- MASS::mvrnorm(n=1, mu=rep(0, I),
                             Sigma=params$sd.geno.add^2 * relmat$geno.add)
        y <- X %*% params$fix.eff + Z %*% g.a + e
        if(! is.null(idx.geno.dom)){
          g.d <- MASS::mvrnorm(n=1, mu=rep(0, I),
                               Sigma=params$sd.geno.dom^2 * relmat$geno.dom)
          y <- y + g.d
        }
        outputs$data[,response] <- y

        return(outputs)
      }
      inputs <- list(formula=formula, data=data, relmat=relmat,
                     REML=REML, na.action=na.action, nrep=nrep)
      ## .bootStat(inputs) # to debug
      params <- list(sd.err=sqrt(as.numeric(vc["Residual"])),
                     sd.geno.add=sqrt(as.numeric(vc["geno.add"])),
                     fix.eff=lme4::fixef(fit))
      if("geno.dom" %in% names(vc))
        params$sd.geno.dom <- sqrt(as.numeric(vc["geno.dom"]))
      ## .bootRanGen(inputs, params) # to debug
      out.boot <- boot::boot(data=inputs,
                             statistic=.bootStat,
                             R=nb.boots, sim="parametric",
                             ran.gen=.bootRanGen,
                             mle=params,
                             parallel=parallel, ncpus=ncpus, cl=cl)
      ci <- matrix(nrow=length(out.boot$t0),
                   ncol=2,
                   dimnames=list(names(out.boot$t0),
                                 paste0(100 * c((1-ci.lev)/2,
                                                1-(1-ci.lev)/2), "%")))
      for(i in 1:nrow(ci)){
        tmp <- boot::boot.ci(out.boot, conf=ci.lev, type="perc",
                             index=i)$percent
        ci[i,] <- tmp[1, c((ncol(tmp)-1):ncol(tmp))]
      }
    }
  }
  if(is.null(nrep)){
	  return(list(merMod=fit, vc=vc, h2=h2, prof=out.prof, boot=out.boot, ci=ci))
	  }else{
	  return(list(merMod=fit, vc=vc, h2=h2, h2_ind=h2_ind, prof=out.prof, boot=out.boot, ci=ci))
	  }
 
}

##' Animal model
##'
##' Given I genotypes, Q covariates and N=I*Q phenotypes for the trait, fit an "animal model" with the INLA package via the following likelihood: y = W c + Z g_A + Z g_D + epsilon, where y is Nx1; W is NxQ; Z is NxI; g_A ~ Normal_I(0, sigma_A^2 A) with A the known matrix of additive genetic relationships; g_D ~ Normal_I(0, sigma_D^2 D) with D the known matrix of dominant genetic relationships; epsilon ~ Normal_N(0, sigma^2 Id_N); Cov(g_A,g_D)=0; Cov(g_A,e)=0; Cov(g_D,e)=0.
##' @param data data.frame containing the data corresponding to relmat; should have a column grep-able for "response" as well as a column "geno.add" used with matrix A; if a column "geno.dom" exists, it will be used with matrix D; any other column will be interpreted as corresponding to "fixed effects"
##' @param relmat list containing the matrices of genetic relationships (see \code{\link{estimGenRel}}); additive relationships (matrix A) are compulsory, with name "geno.add"; dominant relationships (matrix D) are optional, with name "geno.dom"; can be in the "matrix" class (base) or the "dsCMatrix" class (Matrix package)
##' @param family a string indicating the likelihood family (see \code{\link[INLA]{inla}})
##' @param nb.threads maximum number of threads the inla-program will use (see \code{\link[INLA]{inla}})
##' @param verbose verbosity level (0/1/2)
##' @param silent if equal to TRUE, then the inla-program would be silent (see \code{\link[INLA]{inla}})
##' @return \code{\link[INLA]{inla}} object
##' @author Timothee Flutre
##' @seealso \code{\link{lmerAM}}, \code{\link{jagsAM}}, \code{\link{stanAM}}
##' @examples
##' \dontrun{## simulate genotypes
##' set.seed(1859)
##' X <- simulGenosDose(nb.genos=200, nb.snps=2000)
##'
##' ## simulate phenotypes with only additive part of genotypic values
##' A <- estimGenRel(X, relationships="additive", method="vanraden1", verbose=0)
##' modelA <- simulAnimalModel(T=1, Q=3, A=A, V.G.A=15, V.E=5, seed=1859)
##'
##' ## infer with INLA
##' library(INLA)
##' modelA$data$geno.add <- modelA$data$geno
##' modelA$data$geno <- NULL
##' fitA <- inlaAM(data=modelA$data, relmat=list(geno.add=A), verbose=1)
##' summary(fitA)
##' c(modelA$C); 1/modelA$V.E; 1/modelA$V.G.A
##'
##' ## simulate phenotypes with additive and dominant parts of genotypic values
##' D <- estimGenRel(X, relationships="dominant", method="vitezica", verbose=0)
##' modelAD <- simulAnimalModel(T=1, Q=3, A=A, V.G.A=15, V.E=5,
##'                             D=D, V.G.D=3, seed=1859)
##'
##' ## infer with INLA
##' modelAD$data$geno.add <- modelAD$data$geno
##' modelAD$data$geno.dom <- modelAD$data$geno
##' modelAD$data$geno <- NULL
##' fitAD <- inlaAM(data=modelAD$data, relmat=list(geno.add=A, geno.dom=D),
##'                 verbose=1)
##' summary(fitAD)
##' c(modelAD$C); 1/modelAD$V.E; 1/modelAD$V.G.A; 1/modelAD$V.G.D
##' }
##' @export
inlaAM <- function(data, relmat, family="gaussian",
                   nb.threads=1, verbose=0, silent=TRUE){
  requireNamespace("INLA", quietly=TRUE)
  stopifnot(is.data.frame(data),
            sum(grepl("response", colnames(data))) == 1,
            "geno.add" %in% colnames(data),
            is.list(relmat),
            ! is.null(names(relmat)),
            "geno.add" %in% names(relmat))

  ## make formula with response, intercept and additive genotypic value
  formula <- paste0(colnames(data)[grepl("response", colnames(data))],
                    " ~ 1",
                    " + f(geno.add, model=\"generic0\", constr=TRUE",
                    ", hyper=list(theta=list(param=c(0.5, 0.5), fixed=FALSE))")
  if(! isSingular(relmat[["geno.add"]])){
    formula <- paste0(formula,
                      ", Cmatrix=solve(relmat[[\"geno.add\"]])")
  } else
    formula <- paste0(formula,
                      ", Cmatrix=mpInv(relmat[[\"geno.add\"]])")
  formula <- paste0(formula,
                    ", values=levels(data$geno.add))")

  ## add dominant genotypic value, if present
  if("geno.dom" %in% names(relmat)){
    formula <- paste0(formula,
                      " + f(geno.dom, model=\"generic0\", constr=TRUE",
                      ", hyper=list(theta=list(param=c(0.5, 0.5), fixed=FALSE))")
    if(! isSingular(relmat[["geno.dom"]])){
      formula <- paste0(formula,
                        ", Cmatrix=solve(relmat[[\"geno.dom\"]])")
    } else
      formula <- paste0(formula,
                        ", Cmatrix=mpInv(relmat[[\"geno.dom\"]])")
    formula <- paste0(formula,
                      ", values=levels(data$geno.dom))")
  }

  ## add "fixed effects", if any
  for(cn in colnames(data))
    if(! grepl("response", cn) & cn != "geno.add" & cn != "geno.dom")
      formula <- paste0(formula, " + ", cn)

  ## finalize formula
  formula <- stats::as.formula(formula)
  if(verbose > 0)
    print(formula)

  ## fit the model with INLA
  fit <- INLA::inla(formula=formula,
                    family=family,
                    data=data,
                    control.compute=list(dic=TRUE),
                    num.threads=nb.threads,
                    verbose=ifelse(verbose <= 1, FALSE, TRUE),
                    silent=silent)

  return(fit)
}

##' Animal model
##'
##' Given I genotypes, Q covariates and N=I*Q phenotypes for the trait, fit an "animal model" with the rjags package via the following likelihood: y = W c + Z g_A + Z g_D + epsilon, where y is Nx1; W is NxQ; Z is NxI; g_A ~ Normal_I(0, sigma_A^2 A) with A the known matrix of additive genetic relationships; g_D ~ Normal_I(0, sigma_D^2 D) with D the known matrix of dominant genetic relationships; epsilon ~ Normal_N(0, sigma^2 Id_N); Cov(g_A,g_D)=0; Cov(g_A,e)=0; Cov(g_D,e)=0.
##' @param data data.frame containing the data corresponding to relmat; should have a column grep-able for "response" as well as a column "geno.add" used with matrix A; if a column "geno.dom" exists, it will be used with matrix D; any other column will be interpreted as corresponding to "fixed effects"
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
##' modelA <- simulAnimalModel(T=1, Q=3, A=A, V.G.A=15, V.E=5, seed=1859)
##'
##' ## infer with rjags
##' modelA$data$geno.add <- modelA$data$geno; modelA$data$geno <- NULL
##' fitA <- jagsAM(data=modelA$data, relmat=list(geno.add=A))
##' plotMcmcChain(fitA[[1]], "V.G.A", 1:4, modelA$V.G.A)
##' cbind(truth=c(c(modelA$C), modelA$V.G.A, modelA$V.E),
##'       summaryMcmcChain(fitA[[1]], c("c[1]", "c[2]", "c[3]", "sigma.A2", "V.E")))
##' }
##' @export
jagsAM <- function(data, relmat, inits=NULL,
                   priors=list(fix=list(dist="c", par=5),
                               vc=list(dist="hc", par=5)),
                   nb.chains=1, nb.adapt=10^3, burnin=10^2,
                   nb.iters=10^3, thin=10,
                   progress.bar=NULL, rm.jags.file=TRUE, verbose=0){
  requireNamespace("rjags", quietly=TRUE)
  stopifnot(is.data.frame(data),
            sum(grepl("response", colnames(data))) == 1,
            "geno.add" %in% colnames(data),
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
  y <- data[, grepl("response", colnames(data))]
  W <- matrix(1, nrow=nrow(data), ncol=1)
  colnames(W) <- "mu"
  for(j in 1:ncol(data))
    if(! grepl("response", colnames(data)[j]) & colnames(data)[j] != "geno.add" &
       colnames(data)[j] != "geno.dom"){
      W <- cbind(W, stats::model.matrix(~ data[,j])[,-1])
    }
  colnames(W) <- gsub("data\\[, j\\]", "", colnames(W))
  data.list <- list(N=nrow(data), Q=ncol(W),
                    y=y, W=W, Z=stats::model.matrix(~ data[,"geno.add"] - 1),
                    mean.mu=mean(y), par.c=priors$fix$par,
                    par.vc=priors$vc$par, A=relmat[["geno.add"]],
                    mean.g.A=rep(0, nlevels(data[,"geno.add"])))
  if("geno.dom" %in% names(relmat)){
    data.list$D <- relmat[["geno.dom"]]
    data.list$mean.g.D <- rep(0, nlevels(data[,"geno.add"]))
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

##' Animal model (multivariate)
##'
##' Given T traits, I genotypes, Q covariates and N=I*Q phenotypes per trait, fit an "animal model" with the rstan package via the following likelihood: \eqn{Y = W C + Z G_A + Z G_D + E}, where Y is NxT; W is NxQ; Z is NxI; G_A ~ Normal_IxT(0, A, V_G_A) with A the known matrix of additive genetic relationships; G_D ~ Normal_IxT(0, D, V_G_D) with D the known matrix of dominant genetic relationships; E ~ Normal_NxT(0, Id_N, V_E); Missing phenotypes are jointly imputed with the other unknown variables, and the errors can follow a Student's t distribution to handle outliers.
##' @param data data.frame containing the data corresponding to relmat; should have columns grep-able for "response" as well as a column "geno.add" used with matrix A; if a column "geno.dom" exists, it will be used with matrix D; any other column will be interpreted as corresponding to "fixed effects"
##' @param relmat list containing the matrices of genetic relationships (see \code{\link{estimGenRel}}); additive relationships (matrix A) are compulsory, with name "geno.add"; dominant relationships (matrix D) are optional, with name "geno.dom"; can be in the "matrix" class (base) or the "dsCMatrix" class (Matrix package)
##' @param inits list of initial values (possible to use 1 sub-list per chain, see \code{\link[rjags]{jags.model}}); if NULL, JAGS will choose typical values (usually mean, median, or mode of the prior)
##' @param nb.chains number of independent chains to run (see \code{\link[rjags]{jags.model}})
##' @param nb.adapt number of iterations for adaptation (see \code{\link[rjags]{jags.model}})
##' @param burnin number of initial iterations to discard (see the update function of the rjags package)
##' @param nb.iters number of iterations to monitor (see \code{\link[rjags]{coda.samples}})
##' @param thin thinning interval for monitored iterations (see \code{\link[rjags]{coda.samples}})
##' @param progress.bar type of progress bar (text/gui/none or NULL)
##' @param rm.jags.file remove the file specifying the model written in the JAGS-dialect of the BUGS language
##' @param verbose verbosity level (0/1)
##' @author Timothee Flutre, inspired by code from Najla Saad Elhezzani (arXiv:1507.08638)
##' @seealso \code{\link{jagsAM}}
##' @export
jagsAMmv <- function(data, relmat, inits=NULL,
                     nb.chains=1, nb.adapt=10^3, burnin=10^2,
                     nb.iters=10^3, thin=10,
                     progress.bar=NULL, rm.jags.file=TRUE, verbose=0){
  requireNamespace("rjags", quietly=TRUE)
  if(! is.null(inits))
    stopifnot(is.list(inits))

  ## define model in JAGS dialect of BUGS language in separate file
  st <- Sys.time()
  jags.file <- tempfile(paste0("jagsAMmv_", format(st, "%Y-%m-%d-%H-%M-%S"), "_"),
                        getwd(), ".jags")
  model.txt <- paste0(
      "# Goal: fit an \"animal model\" with rjags",
      "\n# Author: Timothee Flutre (INRA)",
      "\n# Source: rutilstimflutre ", utils::packageVersion("rutilstimflutre"),
      "\n# Date: ", format(st, "%Y-%m-%d %H:%M:%S"))
  cat(model.txt, file=jags.file)
  ## TODO: completing the rest!

  ## read model file and create object
  data.list <- list()
  ## TODO: completing the rest!
  jags <- rjags::jags.model(file=jags.file, data=data.list, inits=inits,
                            n.chains=nb.chains,
                            n.adapt=nb.adapt,
                            quiet=ifelse(verbose > 0, FALSE, TRUE))
  if(rm.jags.file)
    file.remove(jags.file)

  ## update model for burn-in period
  stats::update(jags, n.iter=burnin, progress.bar=progress.bar)

  ## extract samples from model
  vn <- c() # TODO
  fit <- rjags::coda.samples(model=jags,
                             variable.names=vn,
                             n.iter=nb.iters,
                             thin=thin,
                             progress.bar=progress.bar)

  return(fit)
}

##' Animal model (univariate)
##'
##'
##' @param stan.file path to a file; if NULL, a temporary one will be created
##' @param include.dom include a variance component for dominant genetic relationships
##' @param errors.Student use a Student's t distribution for the errors (useful to handle outliers)
##' @param missing.phenos indicating if there is any missing phenotypes to impute jointly with the other unknown variables
##' @return invisible path to stan.file
##' @author Timothee Flutre
##' @seealso \code{\link{stanAM}}
##' @export
stanAMwriteModel <- function(stan.file=NULL,
                             include.dom=FALSE,
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
# persons: Timothee Flutre [cre,aut]"
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
  if(include.dom)
    model.code <- paste0(model.code, "
  matrix[I,I] D;          // dominant genetic relationships matrix")
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
  matrix[I,I] CA;")
  if(include.dom)
    model.code <- paste0(model.code, "
  matrix[I,I] CD;")
  model.code <- paste0(model.code, "
  CA = cholesky_decompose(A);")
  if(include.dom)
    model.code <- paste0(model.code, "
  CD = cholesky_decompose(D);")
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
  if(include.dom)
    model.code <- paste0(model.code, "
  vector[I] g_D_z;          // primitive of g_D
  real<lower=0> sigma_g_D;")
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
  vector[I] g_A;")
  if(include.dom)
    model.code <- paste0(model.code, "
  vector[I] g_D;")
  model.code <- paste0(model.code, "
  c[1] ~ cauchy(loc_mu, scale_mu); // intercept
  for(q in 2:Q)
    c[q] ~ cauchy(0, 5);
  sigma_E ~ cauchy(0, 5);  // implicit half-Cauchy
  sigma_g_A ~ cauchy(0, 5);
  g_A_z ~ normal(0, 1);
  g_A = sigma_g_A * (CA * g_A_z); // implies g_A ~ multi_normal(0, sigma_g_A^2 * A)")
  if(include.dom)
    model.code <- paste0(model.code, "
  sigma_g_D ~ cauchy(0, 5);
  g_D_z ~ normal(0, 1);
  g_D = sigma_g_D * (CD * g_D_z);")
  if(errors.Student){
    if(missing.phenos){
      model.code <- paste0(model.code, "
  y_obs ~ student_t(nu, W_obs * c + Z_obs * g_A")
      if(include.dom)
        model.code <- paste0(model.code, " + Z_obs * g_D")
      model.code <- paste0(model.code, ", sigma_E);
  y_mis ~ student_t(nu, W_mis * c + Z_mis * g_A")
      if(include.dom)
        model.code <- paste0(model.code, " + Z_mis * g_D")
      model.code <- paste0(model.code, ", sigma_E);
")
    } else{
      model.code <- paste0(model.code, "
  y ~ student_t(nu, W * c + Z * g_A")
      if(include.dom)
        model.code <- paste0(model.code, " + Z * g_D")
      model.code <- paste0(model.code, ", sigma_E);
")
    }
  } else{
    if(missing.phenos){
      model.code <- paste0(model.code, "
  y_obs ~ normal(W_obs * c + Z_obs * g_A")
      if(include.dom)
        model.code <- paste0(model.code, " + Z_obs * g_D")
      model.code <- paste0(model.code, ", sigma_E);
  y_mis ~ normal(W_mis * c + Z_mis * g_A")
      if(include.dom)
        model.code <- paste0(model.code, " + Z_mis * g_D")
      model.code <- paste0(model.code, ", sigma_E);
")
    } else{
      model.code <- paste0(model.code, "
  y ~ normal(W * c + Z * g_A")
      if(include.dom)
        model.code <- paste0(model.code, " + Z * g_D")
      model.code <- paste0(model.code, ", sigma_E);
")
    }
  }
  model.code <- paste0(model.code, "}")

  model.code <- paste0(model.code, "

generated quantities {
  vector[I] g_A;")
  if(include.dom)
    model.code <- paste0(model.code, "
  vector[I] g_D;")
  model.code <- paste0(model.code, "
  g_A = sigma_g_A * (CA * g_A_z);")
  if(include.dom)
    model.code <- paste0(model.code, "
  g_D = sigma_g_D * (CD * g_D_z);
")
  model.code <- paste0(model.code, "
}")

  model.code <- paste0(model.code, "
")

  cat(model.code, file=stan.file)

  invisible(stan.file)
}

##' Animal model (univariate)
##'
##' Given I genotypes, Q covariates and N=I*Q phenotypes for the trait, fit an "animal model" with the rstan package via the following likelihood: y = W c + Z g_A + Z g_D + epsilon, where y is Nx1; W is NxQ; Z is NxI; g_A ~ Normal_I(0, sigma_A^2 A) with A the known matrix of additive genetic relationships; g_D ~ Normal_I(0, sigma_D^2 D) with D the known matrix of dominant genetic relationships; epsilon ~ Normal_N(0, sigma^2 Id_N); Cov(g_A,g_D)=0; Cov(g_A,e)=0; Cov(g_D,e)=0. Missing phenotypes are jointly imputed with the other unknown variables, and the errors can follow a Student's t distribution to handle outliers.
##' @param data data.frame containing the data corresponding to relmat; should have a column grep-able for "response" as well as a column "geno.add" used with matrix A; if a column "geno.dom" exists, it will be used with matrix D; any other column will be interpreted as corresponding to "fixed effects"
##' @param relmat list containing the matrices of genetic relationships (see \code{\link{estimGenRel}}); additive relationships (matrix A) are compulsory, with name "geno.add"; dominant relationships (matrix D) are optional, with name "geno.dom"; can be in the "matrix" class (base) or the "dsCMatrix" class (Matrix package)
##' @param errors.Student use a Student's t distribution for the errors (useful to handle outliers)
##' @param nb.chains number of independent chains to run
##' @param nb.iters number of iterations to monitor
##' @param burnin number of initial iterations to discard
##' @param thin thinning interval for monitored iterations
##' @param out.dir working directory in which the stan and sm files will be written
##' @param task.id identifier of the task (used in temporary file names)
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
##' modelA <- simulAnimalModel(T=1, Q=3, A=A, V.G.A=15, V.E=5, seed=1859)
##'
##' ## infer with rstan
##' modelA$data$geno.add <- modelA$data$geno; modelA$data$geno <- NULL
##' fitA <- stanAM(data=modelA$data, relmat=list(geno.add=A))
##' fitA <- rstan::As.mcmc.list(fitA)
##' plotMcmcChain(fitA[[1]], "sigma_g_A", 1:4, sqrt(modelA$V.G.A))
##' cbind(truth=c(c(modelA$C), sqrt(modelA$V.G.A), sqrt(modelA$V.E)),
##'       summaryMcmcChain(fitA[[1]], c("c[1]", "c[2]", "c[3]", "sigma_g_A", "sigma_E")))
##'
##' ## simulate phenotypes with additive and dominant parts of genotypic values
##' D <- estimGenRel(X, relationships="dominant", method="vitezica", verbose=0)
##' modelAD <- simulAnimalModel(T=1, Q=3, A=A, V.G.A=15, V.E=5,
##'                             D=D, V.G.D=3, seed=1859)
##'
##' ## infer with rstan
##' modelAD$data$geno.add <- modelAD$data$geno; modelAD$data$geno <- NULL
##' fitAD <- stanAM(data=modelAD$data, relmat=list(geno.add=A, geno.dom=D))
##' fitAD <- rstan::As.mcmc.list(fitAD)
##' plotMcmcChain(fitAD[[1]], "sigma_g_A", 1:4, sqrt(modelAD$V.G.A))
##' plotMcmcChain(fitAD[[1]], "sigma_g_D", 1:4, sqrt(modelAD$V.G.D))
##' cbind(truth=c(c(modelAD$C), sqrt(modelAD$V.G.A), sqrt(modelAD$V.G.D), sqrt(modelAD$V.E)),
##'       summaryMcmcChain(fitAD[[1]], c("c[1]", "c[2]", "c[3]", "sigma_g_A", "sigma_g_D", "sigma_E")))
##' plot(as.vector(fitAD[[1]][,"sigma_g_A"]), as.vector(fitAD[[1]][,"sigma_g_D"]))
##' }
##' @export
stanAM <- function(data, relmat, errors.Student=FALSE,
                   nb.chains=1, nb.iters=10^3, burnin=10^2, thin=10,
                   out.dir=getwd(),
                   task.id=format(Sys.time(), "%Y-%m-%d-%H-%M-%S"),
                   compile.only=FALSE, rm.stan.file=FALSE, rm.sm.file=FALSE,
                   verbose=0){
  requireNamespace("rstan", quietly=TRUE)
  stopifnot(is.data.frame(data),
            sum(grepl("response", colnames(data))) == 1,
            "geno.add" %in% colnames(data),
            is.list(relmat),
            ! is.null(names(relmat)),
            "geno.add" %in% names(relmat),
            nrow(relmat[["geno.add"]]) == ncol(relmat[["geno.add"]]),
            dir.exists(out.dir))
  if("geno.dom" %in% names(relmat))
    stopifnot(nrow(relmat[["geno.add"]]) == nrow(relmat[["geno.dom"]]))

  ## make the input matrices from the input data.frame
  y <- data[, grepl("response", colnames(data))]
  if(all(is.matrix(y), ncol(y) == 1))
    y <- y[,1]
  is_y_obs <- (! is.na(y))
  W <- matrix(1, nrow=nrow(data), ncol=1)
  colnames(W) <- "mu"
  for(j in 1:ncol(data))
    if(! grepl("response", colnames(data)[j]) & colnames(data)[j] != "geno.add" &
       colnames(data)[j] != "geno.dom"){
      W <- cbind(W, stats::model.matrix(~ data[,j])[,-1])
    }
  colnames(W) <- gsub("data\\[, j\\]", "", colnames(W))
  Z <- stats::model.matrix(~ data[,"geno.add"] - 1)
  stopifnot(ncol(Z) == nrow(relmat[["geno.add"]]))
  N <- nrow(W)
  Q <- ncol(W)
  I <- ncol(Z)

  ## define model in STAN language in separate file
  stan.file <- paste0(out.dir, "/stanAM_", task.id, ".stan")
  stanAMwriteModel(stan.file, "geno.dom" %in% names(relmat), errors.Student,
                   sum(is_y_obs) != N)

  ## compile, or make the input list and run
  sm.file <- paste0(out.dir, "/stanAM_", task.id, "_sm.RData")
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
    if("geno.dom" %in% names(relmat))
      ldat$D <- relmat[["geno.dom"]]
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

##' Animal model (multivariate)
##'
##'
##' @param stan.file path to a file; if NULL, a temporary one will be created
##' @param verbose verbosity level (0/1)
##' @return invisible path to stan.file
##' @author Timothee Flutre
##' @seealso \code{\link{stanAMmv}}
##' @export
stanAMmvwriteModel <- function(stan.file=NULL, verbose=0){
  st <- Sys.time()
  if(is.null(stan.file))
    stan.file <- tempfile(pattern=paste0("stanAMmv_",
                                         format(st, "%Y-%m-%d-%H-%M-%S"),
                                         "_"),
                          tmpdir=getwd(), fileext=".stan")
  if(verbose > 0){
    msg <- paste0("write model in file '", stan.file, "'...")
    write(msg, stdout())
  }

  ## meta-data
  model.code <- "// model: stanAMmv"
  model.code <- paste0(model.code, "
// ref: https://discourse.mc-stan.org/t/multivariate-animal-model/631/13")
  model.code <- paste0(model.code, "
// copyright: ?
// license: ?
// persons: Diogo Melo [cre,aut], Timothee Flutre [ctb]")
  model.code <- paste0(model.code, "
// date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "
")

  ## block: functions
  model.code <- paste0(model.code, "
functions {
  matrix as_matrix(vector X, int N, int K) {
    matrix[N, K] Y;
    for (i in 1:N) {
      Y[i] = to_row_vector(X[((i - 1) * K + 1):(i * K)]);
    }
    return Y;
  }
  vector chol_kronecker_product(matrix LA, matrix LG, vector a) {
    vector[num_elements(a)] new_a;
    new_a = rep_vector(0, num_elements(a));
    for(iA in 1:cols(LA)){
      for(jA in 1:iA){
        if(LA[iA, jA] > 1e-10){ // avoid calculating products between unrelated individuals
          for(iG in 1:cols(LG)){
            for(jG in 1:iG){
              new_a[(cols(LG)*(iA-1))+iG] = new_a[(cols(LG)*(iA-1))+iG] +
                                            LA[iA, jA] * LG[iG, jG] * a[(cols(LG)*(jA-1))+jG];
            }
          }
        }
      }
    }
    return new_a;
  }
}
")

  ## block: data
  model.code <- paste0(model.code, "
data {
  int<lower=1>    K; // number of traits
  int<lower=1>    J; // number of fixed effects
  int<lower=0>    I; // number of genotypes
  int<lower=0>    N; // number of phenotypes
  vector[J]    X[N]; // design matrix for fixed effects
  vector[I]    Z[N]; // design matrix for random effects
  vector[K]    Y[N]; // response variable
  matrix[I, I]    A; // relationship matrix
}
")

  ## block: transformed data
  model.code <- paste0(model.code, "
transformed data {
  matrix[I, I] LA;
  LA = cholesky_decompose(A);
}
")

  ## block: parameters
  model.code <- paste0(model.code, "
parameters {
  matrix[K,J]    beta; // fixed effects
  vector[I*K] a_tilde; // breeding values

  cholesky_factor_corr[K] L_Omega_G;
  vector<lower=0>[K] L_sigma_G;

  cholesky_factor_corr[K] L_Omega_R;
  vector<lower=0>[K] L_sigma_R;
}
")

  ## block: transformed parameters
  model.code <- paste0(model.code, "
transformed parameters {
  matrix[I, K] a;
  a = as_matrix(chol_kronecker_product(LA, diag_pre_multiply(L_sigma_G, L_Omega_G), a_tilde), I, K);
}
")

  ## block: model
  model.code <- paste0(model.code, "
model {
    vector[K] mu[N];
    matrix[K, K] L_Sigma_R;

    L_Sigma_R = diag_pre_multiply(L_sigma_R, L_Omega_R);

    for(n in 1:N)
      mu[n] = beta * X[n] + a * Z[n];

    Y ~ multi_normal_cholesky(mu, L_Sigma_R);

    to_vector(beta) ~ normal(0, 1);
    a_tilde   ~ normal(0, 1);
    L_Omega_G ~ lkj_corr_cholesky(4);
    L_sigma_G ~ cauchy(0, 5);
    L_Omega_R ~ lkj_corr_cholesky(4);
    L_sigma_R ~ cauchy(0, 5);
}
")


  ## block: generated quantities
  model.code <- paste0(model.code, "
generated quantities {
    cov_matrix[K] P;
    cov_matrix[K] G;
    cov_matrix[K] E;
    corr_matrix[K] corrG;
    corr_matrix[K] corrE;

    G = multiply_lower_tri_self_transpose(diag_pre_multiply(L_sigma_G, L_Omega_G));
    E = multiply_lower_tri_self_transpose(diag_pre_multiply(L_sigma_R, L_Omega_R));
    P = G + E;

    corrG = multiply_lower_tri_self_transpose(L_Omega_G);
    corrE = multiply_lower_tri_self_transpose(L_Omega_R);
}
")

  model.code <- paste0(model.code, "
")

  cat(model.code, file=stan.file)

  invisible(stan.file)
}

##' Animal model (multivariate)
##'
##' Given T traits, I genotypes, Q covariates and N=I*Q phenotypes per trait, fit an "animal model" with the rstan package via the following likelihood: \eqn{Y = W C + Z G_A + Z G_D + E}, where Y is NxT; W is NxQ; Z is NxI; G_A ~ Normal_IxT(0, A, V_G_A) with A the known matrix of additive genetic relationships; G_D ~ Normal_IxT(0, D, V_G_D) with D the known matrix of dominant genetic relationships; E ~ Normal_NxT(0, Id_N, V_E); Missing phenotypes are jointly imputed with the other unknown variables, and the errors can follow a Student's t distribution to handle outliers.
##' @param data data.frame containing the data corresponding to relmat; should have columns grep-able for "response" as well as a column "geno.add" used with matrix A; if a column "geno.dom" exists, it will be used with matrix D; any other column will be interpreted as corresponding to "fixed effects"
##' @param relmat list containing the matrices of genetic relationships (see \code{\link{estimGenRel}}); additive relationships (matrix A) are compulsory, with name "geno.add"; dominant relationships (matrix D) are optional, with name "geno.dom"; can be in the "matrix" class (base) or the "dsCMatrix" class (Matrix package)
##' @param errors.Student use a Student's t distribution for the errors (useful to handle outliers)
##' @param nb.chains number of independent chains to run
##' @param nb.iters number of iterations to monitor
##' @param burnin number of initial iterations to discard
##' @param thin thinning interval for monitored iterations
##' @param out.dir working directory in which the stan and sm files will be written
##' @param task.id identifier of the task (used in temporary file names)
##' @param compile.only only compile the model, don't run it
##' @param rm.stan.file remove the file specifying the model written in the STAN language
##' @param rm.sm.file remove the file corresponding to the compiled model
##' @param verbose verbosity level (0/1)
##' @return path to compiled file if \code{compile.only=TRUE}, object of class \code{\link[rstan]{stanfit-class}} otherwise
##' @author Timothee Flutre
##' @seealso \code{\link{stanAM}}
##' @examples
##' \dontrun{## simulate genotypes
##' set.seed(1859)
##' X <- simulGenosDose(nb.genos=200, nb.snps=2000)
##'
##' ## simulate phenotypes with only additive part of genotypic values
##' A <- estimGenRel(X, relationships="additive", method="vanraden1", verbose=0)
##' modelA <- simulAnimalModel(T=2, Q=3, A=A, nu.G.A=5, nu.E=3, seed=1859)
##'
##' ## infer with rstan
##' modelA$data$geno.add <- modelA$data$geno; modelA$data$geno <- NULL
##' fitA <- stanAMmv(data=modelA$data, relmat=list(geno.add=A))
##' fitA <- rstan::As.mcmc.list(fitA)
##' plotMcmcChain(fitA[[1]], "sigma_g_A", 1:4, sqrt(modelA$V.G.A))
##' cbind(truth=c(c(modelA$C), sqrt(modelA$V.G.A), sqrt(modelA$V.E)),
##'       summaryMcmcChain(fitA[[1]], c("c[1]", "c[2]", "c[3]", "sigma_g_A", "sigma_E")))
##'
##' ## simulate phenotypes with additive and dominant parts of genotypic values
##' D <- estimGenRel(X, relationships="dominant", method="vitezica", verbose=0)
##' modelAD <- simulAnimalModel(T=2, Q=3, A=A, nu.G.A=15, nu.E=5,
##'                             D=D, nu.G.D=3, seed=1859)
##'
##' ## infer with rstan
##' modelAD$data$geno.add <- modelAD$data$geno; modelAD$data$geno <- NULL
##' fitAD <- stanAMmv(data=modelAD$data, relmat=list(geno.add=A, geno.dom=D))
##' fitAD <- rstan::As.mcmc.list(fitAD)
##' plotMcmcChain(fitAD[[1]], "sigma_g_A", 1:4, sqrt(modelAD$V.G.A))
##' plotMcmcChain(fitAD[[1]], "sigma_g_D", 1:4, sqrt(modelAD$V.G.D))
##' cbind(truth=c(c(modelAD$C), sqrt(modelAD$V.G.A), sqrt(modelAD$V.G.D), sqrt(modelAD$V.E)),
##'       summaryMcmcChain(fitAD[[1]], c("c[1]", "c[2]", "c[3]", "sigma_g_A", "sigma_g_D", "sigma_E")))
##' plot(as.vector(fitAD[[1]][,"sigma_g_A"]), as.vector(fitAD[[1]][,"sigma_g_D"]))
##' }
##' @export
stanAMmv <- function(data, relmat, errors.Student=FALSE,
                     nb.chains=1, nb.iters=10^3, burnin=10^2, thin=10,
                     out.dir=getwd(),
                     task.id=format(Sys.time(), "%Y-%m-%d-%H-%M-%S"),
                     compile.only=FALSE, rm.stan.file=FALSE, rm.sm.file=FALSE,
                     verbose=0){
  requireNamespace("rstan", quietly=TRUE)
  stopifnot(is.data.frame(data),
            sum(grepl("response", colnames(data))) >= 1,
            "geno.add" %in% colnames(data),
            is.list(relmat),
            ! is.null(names(relmat)),
            "geno.add" %in% names(relmat),
            nrow(relmat[["geno.add"]]) == ncol(relmat[["geno.add"]]),
            dir.exists(out.dir))
  if("geno.dom" %in% names(relmat))
    stopifnot(nrow(relmat[["geno.add"]]) == nrow(relmat[["geno.dom"]]))

  ## make the input matrices from the input data.frame
  Y <- as.matrix(data[, grepl("response", colnames(data))])
  W <- matrix(1, nrow=nrow(data), ncol=1)
  colnames(W) <- "mu"
  for(j in 1:ncol(data))
    if(! grepl("response", colnames(data)[j]) &
       colnames(data)[j] != "geno.add" &
       colnames(data)[j] != "geno.dom"){
      W <- cbind(W, stats::model.matrix(~ data[,j])[,-1])
    }
  colnames(W) <- gsub("data\\[, j\\]", "", colnames(W))
  Z <- stats::model.matrix(~ data[,"geno.add"] - 1)
  stopifnot(ncol(Z) == nrow(relmat[["geno.add"]]))
  T <- ncol(Y)
  N <- nrow(W)
  Q <- ncol(W)
  I <- ncol(Z)
  if(verbose > 0){
    msg <- paste0("dimensions: ", T, " traits",
                  ", ", Q, " fixed effect", ifelse(Q > 1, "s", ""),
                  ", ", I, " individuals")
    write(msg, stdout())
  }

  ## define model in STAN language in separate file
  stan.file <- paste0(out.dir, "/stanAMmv_", task.id, ".stan")
  stanAMmvwriteModel(stan.file, verbose)

  ## compile, or make the input list and run
  sm.file <- paste0(out.dir, "/stanAMmv_", task.id, "_sm.RData")
  if(compile.only){
    if(verbose > 0)
      write(paste0("compile..."), stdout())
    rt <- rstan::stanc(file=stan.file, model_name="stanAMmv")
    sm <- rstan::stan_model(stanc_ret=rt,
                            verbose=ifelse(verbose > 0, TRUE, FALSE))
    save(sm, file=sm.file)
    return(sm.file)
  } else{ # compile (if necessary) and run
    ldat <- list(K=T, # nb of traits
                 J=Q, # nb of fixed effects
                 N=I, # nb of individuals
                 X=W, # design matrix of fixed effects
                 Z=Z, # design matrix of random effects
                 Y=Y, # matrix of responses
                 A=relmat[["geno.add"]])
    if(file.exists(sm.file)){ # run only
      if(verbose > 0)
        write(paste0("run..."), stdout())
      load(sm.file)
      fit <- rstan::sampling(sm,
                             data=ldat,
                             iter=nb.iters + burnin,
                             warmup=burnin,
                             thin=thin,
                             chains=nb.chains)
    } else # compile and run
      if(verbose > 0)
        write(paste0("compile and run..."), stdout())
      fit <- rstan::stan(file=stan.file,
                         data=ldat,
                         iter=nb.iters + burnin,
                         warmup=burnin,
                         thin=thin,
                         chains=nb.chains)
  }
  if(rm.stan.file){
    if(file.exists(stan.file))
      file.remove(stan.file)
  }
  if(rm.sm.file){
    if(file.exists(sm.file))
      file.remove(sm.file)
  }

  return(fit)
}

##' BVSR
##'
##' Simulate phenotypes according to the following model: Y = W c + Z X_A a + Z X_D d + epsilon where Y is N x 1; W is N x Q; c is Q x 1; Z is N x I; X_A/D is I x P and epsilon is N x 1 with epsilon ~ Normal_N(0, sigma^2 Id) and c ~ Normal(mean_a, sd_a) so that sd_a is large ("fixed effect").
##' For SNP p, gamma_p indicates if it is causal, i.e. non-zero additive and/or dominant effect, where Prob(gamma_p=1) is named pi.
##' For the case where pi is small, see Guan & Stephens (2011), Carbonetto & Stephens (2012), Peltola et al (2012), Verzelen (2012).
##' Causal SNP p can have an additive effect, a_p | gamma_p=1 ~ Normal_1(0, sigma_a^2), a dominant effect, d_p | gamma_p=1 ~ Normal_1(0, sigma_d^2), or both.
##' @param Q number of fixed effects, i.e. the intercept plus the number of years during which genotypes are phenotyped (starting in 2010)
##' @param mu overall mean
##' @param mean.c mean of the prior on c[2:Q]
##' @param sd.c std dev of the prior on c[2:Q]
##' @param X matrix of bi-allelic SNP genotypes encoded in allele doses in [0,2], with genotypes in rows and SNPs in columns (SNPs with missing values or low MAF should be discarded beforehand)
##' @param m1 if TRUE, 1 will be subtracted from each entry of X to make X_A (before optional centering and scaling)
##' @param ctr if TRUE, the columns of the X matrix will be centered to make X_A and X_D
##' @param std if TRUE, the columns of the X matrix will also be standardized to make X_A and X_D (note that Habier et al (2007) don't standardize the genotypes, but Janss et al (2012) do)
##' @param pi proportion of marker effects (a) that are non-zero; setting pi at 1 means simulating from the additive infinitesimal model (equivalent to ridge regression)
##' @param pve proportion of phenotypic variance explained by SNPs with non-zero effects ("heritability"); PVE = V[g] / V[y] where y = g + e and g = g_a + g_d (no epistasis); the magnitude of g_a (resp. g_d) depends on whether or not \code{sigma.a2} (resp. \code{sigma.d2}) is set to zero; a value for sigma^2 is then chosen
##' @param sigma.a2 prior variance of the non-zero additive effects
##' @param sigma.d2 prior variance of the non-zero dominant effects (if non-null, a reasonable choice is half of \code{sigma.a2}, as in Servin & Stephens (2007) with their prior D2)
##' @param min.maf minimum minor allele frequency below which SNPs can't have any effect (neither additive nor dominant)
##' @param perc.NA percentage of missing phenotypes, at random
##' @param err.df degrees of freedom of errors' Student's t-distribution
##' @param seed seed for the pseudo-random number generator
##' @param verbose verbosity level (0/1)
##' @return list
##' @author Timothee Flutre
##' @examples
##' \dontrun{## simulate genotypes
##' set.seed(1859)
##' I <- 200
##' P <- 2000
##' X <- simulGenosDose(nb.genos=I, nb.snps=P)
##'
##' ## estimate genetic relationships
##' A <- estimGenRel(X=X, relationships="additive", method="vanraden1")
##' D <- estimGenRel(X=X, relationships="dominant", method="vitezica")
##'
##' ## additive sparse genetic architecture
##' ## choose pi so that sum(gamma * (1 + log(P / sum(gamma)))) < I
##' Q <- 3
##' modelA <- simulBvsr(Q=Q, X=X, pi=0.01, pve=0.7, sigma.a2=1)
##' modelA$sigma2
##'
##' library(lme4)
##' fit1 <- lmer(formula=response1 ~ year + (1|geno), data=modelA$dat)
##' cbind(modelA$c, blues <- fixef(fit1))
##' blups <- ranef(fit1, drop=TRUE)$geno[rownames(X)]
##' cor(modelA$g.A, blups, method="pearson")
##' cor(modelA$g.A, blups, method="spearman")
##' (vc1 <- as.data.frame(VarCorr(fit1)))
##' (pve1.hat <- vc1$vcov[1] / (vc1$vcov[1] + vc1$vcov[2]))
##'
##' fit1A <- lmerAM(formula=response1 ~ year + (1|geno), data=modelA$dat,
##'                 relmat=list(geno=A))
##' ## maybe necessary to use A2 <- nearPD(A)$mat
##' cbind(modelA$c, blues <- fixef(fit1A$merMod))
##' blupsA <- ranef(fit1A$merMod, drop=TRUE)$geno[rownames(X)]
##' cor(modelA$g.A, blupsA, method="pearson")
##' cor(modelA$g.A, blupsA, method="spearman")
##' (vc1A <- as.data.frame(VarCorr(fit1A$merMod)))
##' (pve1A.hat <- vc1A$vcov[1] / (vc1A$vcov[1] + vc1A$vcov[2]))
##'
##' plot(x=modelA$g.A, y=blups, col="blue", asp=1, las=1)
##' points(x=modelA$g.A, y=blupsA, col="red")
##' segments(x0=modelA$g.A, y0=blups, x1=modelA$g.A, y1=blupsA)
##' abline(h=0, v=0, a=0, b=1, lty=2)
##' legend("bottomright", legend=c("without A", "with A"), col=c("blue", "red"),
##'        pch=1, bty="n")
##'
##' library(varbvs)
##' fit2 <- varbvs(X=modelA$X.A, Z=NULL, y=blups, verbose=FALSE)
##' print(fit2s <- summary(fit2))
##' names(sort(modelA$a[modelA$gamma == 1], decreasing=TRUE))
##'
##' (pi.hat <- 10^(fit2s$logodds$x0) / (1 + 10^(fit2s$logodds$x0)))
##' (pi.hat.low <- 10^(fit2s$logodds$a) / (1 + 10^(fit2s$logodds$a)))
##' (pi.hat.high <- 10^(fit2s$logodds$b) / (1 + 10^(fit2s$logodds$b)))
##'
##' w <- c(normalizelogweights(fit2$logw))
##' pips <- c(fit2$alpha %*% w)
##' cols <- rep("black", P)
##' cols[modelA$gamma != 0] <- "red"
##' plot(x=1:P, y=pips, col=cols, las=1, xlab="SNPs", ylab="PIP",
##'      main="Posterior inclusion probabilities")
##'
##' y.pred <- predict(fit2, X=modelA$X.A, Z=NULL)
##' cor(blups, y.pred)
##'
##' ## dominant sparse genetic architecture
##' Q <- 3
##' modelAD <- simulBvsr(Q=Q, X=X, pi=0.01, pve=0.7, sigma.a2=0, sigma.d2=1)
##' modelAD$sigma2
##'
##' library(lme4)
##' dat <- modelAD$dat
##' colnames(dat)[colnames(dat) == "geno"] <- "geno.add"
##' dat$geno.dom <- dat$geno.add
##' fit3A <- lmerAM(formula=response1 ~ year + (1|geno.add), data=dat,
##'                 relmat=list(geno.add=A))
##' (vc3A <- as.data.frame(VarCorr(fit3A$merMod)))
##' (pve3A.hat <- vc3A$vcov[1] / (vc3A$vcov[1] + vc3A$vcov[2]))
##'
##' fit3D <- lmerAM(formula=response1 ~ year + (1|geno.dom), data=dat,
##'                 relmat=list(geno.dom=D))
##' (vc3D <- as.data.frame(VarCorr(fit3D$merMod)))
##' (pve3D.hat <- vc3D$vcov[1] / (vc3D$vcov[1] + vc3D$vcov[2]))
##'
##' fit3AD <- lmerAM(formula=response1 ~ year + (1|geno.add) + (1|geno.dom),
##'                  data=dat, relmat=list(geno.add=A, geno.dom=D))
##' (vc3AD <- as.data.frame(VarCorr(fit3AD$merMod)))
##' (pve3AD.hat <- (vc3AD$vcov[1] + vc3AD$vcov[2]) /
##'                  (vc3AD$vcov[1] + vc3AD$vcov[2] + vc3AD$vcov[3]))
##'
##' ## additive infinitesimal genetic architecture
##' Q <- 3
##' modelA <- simulBvsr(Q=Q, X=X, pi=1, pve=0.7, sigma.a2=1)
##' modelA$sigma2
##'
##' library(rrBLUP)
##' fit.RR <- mixed.solve(y=modelA$Y[,1], Z=modelA$Z %*% modelA$X.A,
##'                       K=NULL, X=modelA$W)
##' fit.RR$Ve
##' fit.RR$Vu
##' fit.AM <- mixed.solve(y=modelA$Y[,1], Z=modelA$Z, K=A, X=modelA$W)
##' fit.AM$Ve
##' fit.AM$Vu
##' afs <- estimSnpAf(X=X)
##' fit.AM$Vu / (2 * sum(afs * (1 - afs))) # see Habier et al (2007)
##' fit.AM$Vu / (fit.AM$Vu + fit.AM$Ve) # same as the specified pve
##' }
##' @export
simulBvsr <- function(Q=3, mu=50, mean.c=5, sd.c=2,
                      X, m1=TRUE, ctr=TRUE, std=FALSE,
                      pi=1, pve=0.7, sigma.a2=1, sigma.d2=0,
                      min.maf=0, perc.NA=0, err.df=Inf, seed=NULL,
                      verbose=1){
  stopIfNotValidGenosDose(X)
  stopifnot(sd.c >= 0,
            is.logical(m1),
            is.logical(ctr),
            is.logical(std),
            pi >= 0,
            pi <= 1,
            pve >= 0,
            pve <= 1,
            sigma.a2 >= 0,
            sigma.d2 >= 0,
            min.maf < 1,
            perc.NA >= 0,
            perc.NA <= 100)
  if(! is.null(seed))
    set.seed(seed)

  if(min.maf > 0)
    X <- discardSnpsLowMaf(X=X, thresh=min.maf, verbose=verbose)

  ## determine the dimensions
  I <- nrow(X)
  P <- ncol(X)
  if(Q > 1){
    N <- Q * I
  } else
    N <- I

  ## incidence matrix of the non-genetic predictors having "fixed effects"
  levels.years <- as.character(seq(from=2010, to=2010+Q-1))
  if(N %% Q == 0){
    years <- rep(levels.years, each=N / Q)
  } else
    years <- sort(sample(x=levels.years, size=N, replace=TRUE))
  years <- as.factor(years)
  if(Q == 1){
    W <- matrix(data=1, nrow=N, ncol=Q)
    rownames(W) <- rownames(X)
  } else
    W <- stats::model.matrix(~ years)
  dat <- data.frame(year=years)

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
  if(Q == 1)
    rownames(Z) <- rownames(X)
  if(m1){
    X.A <- scale(x=X - 1, center=ctr, scale=std)
  } else
    X.A <- scale(x=X, center=ctr, scale=std)
  X.D <- matrix(data=0, nrow=nrow(X.A), ncol=ncol(X.A))
  if(sigma.d2 > 0){
    X.D <- scale(x=recodeIntoDominant(X), center=ctr, scale=std)
    if(std & any(is.na(X.D)))
      stop("use min.maf to avoid NA with std=TRUE for X.D")
  }
  dat$geno <- inds

  ## incidence vector of the causal genetic predictors
  gamma <- stats::setNames(object=stats::rbinom(n=P, size=1, prob=pi),
                           nm=colnames(X))

  ## additive "random effects"
  a <- stats::setNames(object=stats::rnorm(n=P, mean=0, sd=sqrt(sigma.a2)),
                       nm=colnames(X))
  a[gamma == 0] <- 0
  g.A <- X.A %*% a

  ## dominant "random effects"
  d <- stats::setNames(object=stats::rnorm(n=P, mean=0, sd=sqrt(sigma.d2)),
                       nm=colnames(X))
  d[gamma == 0] <- 0
  g.D <- X.D %*% d

  ## errors
  sigma2 <- ((1 - pve) / pve) * stats::var(g.A + g.D)
  if(is.infinite(err.df)){
    epsilon <- matrix(stats::rnorm(n=N, mean=0, sd=sqrt(sigma2)))
  } else
    epsilon <- matrix(stats::rt(n=N, df=err.df, ncp=0))

  ## phenotypes
  Y <- W %*% c + Z %*% g.A + Z %*% g.D + epsilon
  if(perc.NA > 0){
    idx <- sample.int(n=N, size=floor(perc.NA/100 * N))
    Y[idx, 1] <- NA
  }
  if(Q == 1)
    rownames(Y) <- rownames(X)
  t <- 1
  dat[[paste0("response", t)]] <- Y[,t]

  return(list(Y=Y,
              W=W, c=c,
              Z=Z, X.A=X.A, X.D=X.D,
              pi=pi, gamma=gamma,
              sigma.a2=sigma.a2, sigma.d2=sigma.d2,
              pve=pve, sigma2=sigma2,
              a=a, g.A=g.A, d=d, g.D=g.D,
              dat=dat))
}

##' BSLMM
##'
##' Simulate phenotypes according to the Bayesian sparse linear mixed model (Zhou, Carbonetto & Stephens, 2013)): y = W alpha + Z X_c beta-tilde + Z u + epsilon where y is N x 1; W is N x Q; alpha is Q x 1; Z is N x I; X_c is I x P; u is I x 1 and epsilon is N x 1. For all p, beta-tilde_p ~ pi Norm_1(0, sigma_betat^2/sigma^2) + (1 - pi) delta_0; u ~ Norm_I(0, (sigma_u^2/sigma^2) K) and epsilon ~ Norm_N(0, sigma^2 I).
##' @param Q number of fixed effects, i.e. the intercept plus the number of years during which genotypes are phenotyped (starting in 2010)
##' @param mu overall mean
##' @param mean.a mean of the prior on alpha[2:Q]
##' @param sd.a std dev of the prior on alpha[2:Q]
##' @param X matrix of bi-allelic SNP genotypes encoded in allele doses in [0,2], with genotypes in rows and SNPs in columns (SNPs with missing values or low MAF should be discarded beforehand).
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
##' \dontrun{## simulate genotypes
##' set.seed(1859)
##' X <- simulGenosDose(nb.genos=200, nb.snps=2000)
##' afs <- estimSnpAf(X)
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
##' library(rrBLUP)
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
  requireNamespace("MASS", quietly=TRUE)
  stopIfNotValidGenosDose(X)
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

##' Plant association genetics
##'
##' Subset and sort inputs necessary to perform an analysis of plant association genetics, that is, the subset of cultivars with genotypes and phenotypes, and the subset of markers having genotypes and coordinates.
##' @param ids data.frame of identifiers, with (at least) column names \code{cultivar.code} and \code{accession.code}; the outputs will be sorted according to this option
##' @param y vector, matrix or data.frame which row names should be present in \code{ids$cultivar.code} or \code{ids$accession.code}
##' @param X matrix of bi-allelic SNP genotypes encoded in allele doses in [0,2], with genotypes in rows and SNPs in columns; the "second" allele is arbitrary, it corresponds to the second column of \code{alleles}, which can be the minor or the major allele; row names should be present in \code{ids$accession.code} or \code{ids$cultivar.code}, and column names in \code{rownames(snp.coords)}
##' @param snp.coords data.frame with 2 columns \code{coord} and \code{chr}, and SNP identifiers as row names
##' @param alleles data.frame with SNPs in rows (names as row names) and alleles in columns (exactly 2 columns are required); the second column should correspond to the allele which number of copies is counted at each SNP in \code{X}
##' @param rename.chr.prefix if not NULL, the value of this parameter will be passed to \code{\link{chromNames2integers}} so that chromosome names from snp.coords will be renamed as integers
##' @param verbose verbosity level (0/1)
##' @return list with inputs after subsetting and sorting
##' @author Timothee Flutre
##' @seealso \code{\link{chromNames2integers}}
##' @export
rearrangeInputsForAssoGenet <- function(ids, y, X, snp.coords, alleles,
                                        rename.chr.prefix=NULL,
                                        verbose=1){
  ## check inputs separately from each other
  stopifnot(is.data.frame(ids),
            all(c("cultivar.code", "accession.code") %in% colnames(ids)),
            sum(is.na(ids)) == 0)
  ids <- convertFactorColumnsToCharacter(ids)
  ids$cultivar.code <- as.character(ids$cultivar.code)
  ids$accession.code <- as.character(ids$accession.code)
  if(is.vector(y)){
    stopifnot(! is.null(names(y)))
    y <- as.data.frame(y,
                       row.names=names(y))
  } else if(is.matrix(y)){
    stopifnot(! is.null(rownames(y)))
    y <- as.data.frame(y)
  }
  stopifnot(is.data.frame(y),
            ! is.null(rownames(y)))
  y <- convertFactorColumnsToCharacter(y)
  stopIfNotValidGenosDose(X, check.noNA=FALSE)
  stopifnot(.isValidSnpCoords(snp.coords),
            is.data.frame(alleles),
            ! is.null(row.names(alleles)),
            ncol(alleles) == 2)
  alleles <- convertFactorColumnsToCharacter(alleles)
  if(! is.null(rename.chr.prefix))
    stopifnot(is.character(rename.chr.prefix))

  ## determine if the genotypes of y are accession or cultivar codes
  in.geno.y <- NA
  if(all(rownames(y) %in% ids$accession.code)){
    in.geno.y <- "accession"
  } else if(all(rownames(y) %in% ids$cultivar.code)){
    in.geno.y <- "cultivar"
  }
  if(is.na(in.geno.y)){
    msg <- "rownames(y) should be all accession codes or all cultivar codes"
    stop(msg)
  }
  if(in.geno.y == "accession")
    rownames(y) <- ids$cultivar.code[match(rownames(y), ids$accession.code)]

  ## determine if the genotypes of X are accession or cultivar codes
  in.geno.X <- NA
  if(all(rownames(X) %in% ids$accession.code)){
    in.geno.X <- "accession"
  } else if(all(rownames(X) %in% ids$cultivar.code)){
    in.geno.X <- "cultivar"
  }
  if(is.na(in.geno.X)){
    msg <- "rownames(X) should be all accession codes or all cultivar codes"
    stop(msg)
  }
  if(in.geno.X == "accession")
    rownames(X) <- ids$cultivar.code[match(rownames(X), ids$accession.code)]

  ## determine the subset of genotypes in y and X (as cultivar codes)
  cultivars.tokeep <- intersect(rownames(y), rownames(X))
  if(verbose > 0){
    msg <- paste0("nb of genotypes: ", length(cultivars.tokeep))
    write(msg, stdout())
  }
  ## sort them according to "ids"
  cultivars.tokeep <- cultivars.tokeep[order(match(cultivars.tokeep,
                                                   ids$cultivar.code))]
  ids <- ids[ids$cultivar.code %in% cultivars.tokeep,]
  y <- droplevels(y[cultivars.tokeep, , drop=FALSE])
  X <- X[cultivars.tokeep, , drop=FALSE]

  ## determine the SNPs for which genotypes and coordinates are available
  stopifnot(all(colnames(X) %in% rownames(snp.coords)),
            all(rownames(alleles) %in% rownames(snp.coords)))
  snps.tokeep <- intersect(colnames(X), rownames(alleles))
  if(verbose > 0){
    msg <- paste0("nb of SNPs: ", length(snps.tokeep))
    write(msg, stdout())
  }
  snp.coords <- droplevels(snp.coords[snps.tokeep,])
  X <- X[, snps.tokeep, drop=FALSE]
  alleles <- droplevels(alleles[snps.tokeep,])

  ## rename chromosome names as integers
  cn2i <- NULL
  if(! is.null(rename.chr.prefix)){
    cn2i <- chromNames2integers(
        x=stats::setNames(snp.coords$chr, rownames(snp.coords)),
        prefix=rename.chr.prefix)
    snp.coords$chr <- cn2i$renamed
  }

  return(list(ids=ids, y=y, X=X, snp.coords=snp.coords, alleles=alleles,
              cn2i=cn2i))
}

##' QTLRel per chromosome
##'
##' Given I genotypes, P SNPs, Q covariates (e.g. replicates) and N=I*Q phenotypes, the whole likelihood is: y = W alpha + Z X beta + epsilon, where y is Nx1, W is NxQ, Z is NxI, X is IxP and u = X beta.
##' The QTLRel package (Cheng et al, BMC Genetics, 2011) decomposes the inference into an \emph{ad hoc} procedure of two steps.
##' First, the variance components and fixed effects are estimated: y = W alpha + Z u + epsilon.
##' Second the allele effects are tested SNP per SNP: for all p in {1,..,P}, y = W alpha + Z x_p beta_p + Z u + epsilon, using the fixed effects and variance components estimated in the first step.
##' As the SNPs can be in linkage disequilibrium, it is advised (Yang et al, Nature Genetics, 2014) to perform the procedure chromosome per chromosome, which is the goal of this function; but I advise to also run it with chr.ids=NULL.
##' @param y vector or one-column matrix of phenotypes (the order is important and should be in agreement with the other arguments X, W and Z)
##' @param X matrix of bi-allelic SNP genotypes encoded in allele doses in {0,1,2}, with genotypes in rows and SNPs in columns
##' @param snp.coords data.frame with SNP identifiers as row names and two columns named "chr" and "coord" (or "pos")
##' @param thresh threshold on minor allele frequencies below which SNPs are ignored via \code{\link{discardSnpsLowMaf}} (default=0.01; NULL to skip this step; SNPs are ignored for testing, but all are still used to calculate the matrix of additive genetic relationships as rare SNPs are informative in this purpose)
##' @param chr.ids vector of chromosome identifiers to analyze (if NULL, the regular QTLRel procedure is launched, i.e. all chromosomes are used to estimate the variance components)
##' @param W incidence matrix of covariates (should not contain the column of 1's for the intercept; use \code{\link[stats]{model.matrix}} if necessary)
##' @param Z incidence matrix relating phenotypes to genotypes (if nrow(y) and nrow(X) are different, diagonal otherwise; use \code{\link[stats]{model.matrix}} if necessary)
##' @param method.A method to estimate the additive relationships (see \code{\link{estimGenRel}})
##' @param verbose verbosity level (0/1)
##' @return a list with three data.frames as components, variance.components, fixed.effects and scan
##' @author Timothee Flutre
##' @examples
##' \dontrun{## simulate data
##' set.seed(1859)
##' I <- 200
##' P <- 2000
##' Q <- 3
##' N <- Q * I
##' dat <- data.frame(ind=as.factor(rep(paste0("ind", 1:I), times=Q)),
##'                   year=as.factor(rep(paste0(2003:(2003+Q-1)), each=I)))
##' W <- stats::model.matrix(~ year, dat)
##' alpha <- rnorm(n=Q, mean=50, sd=30)
##' X <- simulGenosDose(nb.genos=I, nb.snps=P)
##' beta <- rnorm(n=P, mean=0, sd=2)
##' Z <- stats::model.matrix(~ ind - 1, dat)
##' dat$response <- as.vector(W %*% alpha + Z %*% X %*% beta + rnorm(N))
##' snp.coords <- data.frame(chr=sample(paste0("chr",1:10), P, TRUE),
##'                          coord=sample.int(10^6, P))
##' rownames(snp.coords) <- colnames(X)
##'
##' ## perform the inference
##' library(QTLRel)
##' res <- qtlrelPerChr(dat$response, X, snp.coords, 0.01, "chr1", W=W[,-1], Z=Z)
##' ## res <- qtlrelPerChr(dat$response, X, snp.coords, 0.01, NULL, W=W[,-1], Z=Z)
##' }
##' @export
qtlrelPerChr <- function(y, X, snp.coords, thresh=0.01, chr.ids=NULL, W=NULL, Z=NULL,
                         method.A="vanraden1", verbose=1){
  requireNamespace("QTLRel", quietly=TRUE)
  if(is.matrix(y))
    stopifnot(ncol(y) == 1)
  y <- as.vector(y)
  stopIfNotValidGenosDose(X, check.hasRowNames=FALSE)
  stopifnot(.isValidSnpCoords(snp.coords),
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

##' Genomic control
##'
##' Estimates lambda (\href{http://dx.doi.org/10.1006/tpbi.2001.1542}{Devlin et al, 2001}).
##' @param pvalues vector of raw p values (missing values will be omitted)
##' @return numeric
##' @author Timothee Flutre
##' @examples
##' \dontrun{## simulate SNP genotypes
##' set.seed(1859)
##' library(scrm)
##' nb.genos <- 200
##' Ne <- 10^4
##' chrom.len <- 10^5
##' mu <- 10^(-8)
##' c <- 10^(-8)
##' genomes <- simulCoalescent(nb.inds=nb.genos,
##'                            pop.mut.rate=4 * Ne * mu * chrom.len,
##'                            pop.recomb.rate=4 * Ne * c * chrom.len,
##'                            chrom.len=chrom.len)
##' X <- genomes$genos
##'
##' ## simulate SNP effects
##' Q <- 1
##' modelA <- simulBvsr(Q=Q, X=X, pi=0.1, pve=0.7, sigma.a2=1)
##'
##' ## filter out SNPs with low MAF
##' mafs <- estimSnpMaf(X)
##' maf_threshold <- 0.05
##' subset_snps <- names(mafs[which(mafs >= maf_threshold)])
##' length(subset_snps)
##'
##' ## perform SNP-by-SNP GWAS *without* kinship matrix
##' library(MM4LMM)
##' out_mmest <- MMEst(Y=modelA$Y[,1], X=X[,subset_snps],
##'                    VarList=list(Error=diag(nrow(modelA$Y))),
##'                                 Verbose=FALSE)
##' out_anovatest <- AnovaTest(out_mmest, Type="TypeI")
##' out_gwas <- t(sapply(out_anovatest, function(x){x["Xeffect",]}))
##' head(out_gwas)
##'
##' ## look at p values
##' qqplotPval(out_gwas[,"pval"])
##' genomicControl(out_gwas[,"pval"])
##'
##' ## estimate genetic relatedness
##' A_vr <- estimGenRel(X, relationships="additive", method="vanraden1")
##'
##' ## perform SNP-by-SNP GWAS *with* kinship matrix
##' library(MM4LMM)
##' out_mmest <- MMEst(Y=modelA$Y[,1], X=X[,subset_snps],
##'                    VarList=list(Additive=A_vr,
##'                                 Error=diag(nrow(modelA$Y))),
##'                                 Verbose=FALSE)
##' out_anovatest <- AnovaTest(out_mmest, Type="TypeI")
##' out_gwas <- t(sapply(out_anovatest, function(x){x["Xeffect",]}))
##'
##' ## look at p values
##' qqplotPval(out_gwas[,"pval"])
##' genomicControl(out_gwas[,"pval"])
##' }
##' @export
genomicControl <- function(pvalues){
  if(any(is.na(pvalues)))
    pvalues <- pvalues[! is.na(pvalues)]
  ## https://www.biostars.org/p/43328/#102996
  zscores <- qnorm(pvalues / 2)
  expValH0 <- qchisq(0.5, 1) # ~ 0.455
  lambda <- median(zscores^2) / expValH0
  return(lambda)
}

##' Gibbs sampler from Janss et al (2012)
##'
##' Run the Gibbs sampler from Janss et al (2012).
##' @param y vector of responses of length n
##' @param a used when y is censored
##' @param b used when y is censored
##' @param family "gaussian" or "bernoulli"
##' @param K list of lists (one per prior)
##' @param XF n x pF incidence matrix for fixed effects
##' @param df0 degrees of freedom for the prior variances
##' @param S0 scale for the prior variances
##' @param nIter number of iterations
##' @param burnIn burn-in
##' @param thin thinning
##' @param saveAt if not NULL, results will be saved in this file
##' @param weights optional vector
##' @param verbose verbosity level (0/1/2)
##' @return list
##' @author Luc Janss [aut], Gustavo de los Campos [aut], Timothee Flutre [ctb]
##' @note To reduce dependency, the R version of the "rtrun" function from the "bayesm" package (version 2.2-5 implemented by Peter Rossi and available under the GPL) was copied. Moreover, because the Janss supplement lacked a function, the "dtruncnorm" function from the "truncnorm" package was adapted.
##' @examples
##' \dontrun{## load the wheat data set
##' c("X","Y") %in% ls()
##' library(BLR)
##' data(wheat)
##' c("X","Y") %in% ls()
##' dim(X) # lines x markers
##' dim(Y) # lines x traits
##'
##' ## factorize the genotypes
##' W <- scale(X, center=TRUE, scale=TRUE)
##' G <- tcrossprod(W) / ncol(W)
##' EVD <- eigen(G)
##'
##' ## look at the G matrix
##' imageWithScale(G)
##' library(seriation)
##' order <- seriate(G)
##' imageWithScale(G[order[[1]], order[[2]]])
##'
##' ## look at the factorization
##' plot(EVD$vectors[,2], EVD$vectors[,1], main="Same as figure 3 left",
##'      xlab="Eigenvector 2", ylab="Eigenvector 1", pch=18)
##' plot(x=0:nrow(X), y=c(0, cumsum(EVD$values) / nrow(X)), type="l", las=1,
##'      main="Same as figure 4 left", ylab="", xlab="Number of eigenvalues")
##' abline(h=1, a=0, b=1/nrow(X), lty=2)
##'
##' ## fit the full model
##' fit.full <- gibbsJanss2012(y=Y[,1], K=list(list(V=EVD$vectors, d=EVD$values,
##'                                                 df0=5, S0=3/2)),
##'                            nIter=20000, burnIn=2000, saveAt=tempfile())
##' fit.full$mu # ~= 0
##' fit.full$varE # ~= 0.54
##' fit.full$K[[1]]$varU # ~= 0.52
##' fit.full$K[[1]]$varU / (fit.full$K[[1]]$varU +
##'                         fit.full$varE) # ~= 0.49
##'
##' ## load posterior samples
##' library(coda)
##' post.full <- mcmc(data=read.table(fit.full$saveAt, header=TRUE))
##' post.full <- mcmc(cbind(post.full,
##'                         h2=post.full[,"varU1"] /
##'                           (post.full[,"varU1"] +
##'                            post.full[,"varE"])))
##' plot(post.full)
##' summary(post.full)
##' HPDinterval(post.full)
##'
##' ## fit models with PCs as fixed effects
##' fit.PC <- list()
##' post.PC <- list()
##' PCs.toaccountfor <- c(1, 5, 10, 15, 20)
##' for(i in seq_along(PCs.toaccountfor)){
##'   nb.PCs <- PCs.toaccountfor[i]
##'   fit.PC[[i]] <- gibbsJanss2012(y=Y[,1], XF=as.matrix(EVD$vectors[,1:nb.PCs]),
##'                                 K=list(list(V=EVD$vectors, d=EVD$values,
##'                                             df0=5, S0=3/2)),
##'                                 nIter=20000, burnIn=2000, saveAt=tempfile())
##'   post.PC[[i]] <- mcmc(data=read.table(fit.PC[[i]]$saveAt, header=TRUE))
##'   post.PC[[i]] <- mcmc(cbind(post.full,
##'                              h2=post.full[,"varU1"] /
##'                                (post.full[,"varU1"] +
##'                                 post.full[,"varE"])))
##' }
##'
##' ## reformat results
##' H2W <- matrix(data=NA, nrow=1 + length(PCs.toaccountfor), ncol=2)
##' rownames(H2W) <- c()
##' colnames(H2W) <- c("with_PCs", "without_PCs")
##' H2W[, "without_PCs"] <- VARU_W[,"without_PCs"] / (VARU_W[,"without_PCs"] +
##'                                                   fit.full$varE)
##' H2W[, "with_PCs"] <- VARU_W[,"with_PCs"] / (VARU_W[,"with_PCs"] +
##'                                             fit.full$varE)
##'
##' ## reproduce figure 6
##' plot(H2W[,2] ~ I(0:20), pch=15, col="red", main="Same as figure 6",
##'      xlab="Number of eigenvectors", ylab="Within-group heritability",
##'      type="o", ylim=c(.95,1.05) * range(as.vector(H2W)))
##' abline(h=range(H2W[,2]), lty=2)
##' abline(v=0)
##' points(x=0:20, y=H2W[,1], type="b", pch=15, col="blue")
##'
##' ## clean
##' file.remove(fit.full$saveAt)
##' }
##' @export
gibbsJanss2012 <- function(y, a=NULL, b=NULL, family="gaussian", K, XF=NULL,
                           df0=0, S0=0, nIter=110, burnIn=10, thin=10,
                           saveAt=NULL, weights=NULL, verbose=0){
  stopifnot(is.numeric(y),
            family %in% c("gaussian", "bernoulli"),
            is.list(K))
  if(! is.null(saveAt))
    if(file.exists(saveAt))
      file.remove(saveAt)

  rtrun <- function(mu, sigma, a, b){
    ## the original implementation (supplement of Janss et al, 2012) uses
    ##  the rtrun() function from the R package "bayesm" v2.2-5,
    ##  so I copied it here to avoid depending on the whole package
    ## see https://en.wikipedia.org/wiki/Truncated_normal_distribution#Simulating
    FA <- stats::pnorm(((a - mu) / sigma))
    FB <- stats::pnorm(((b - mu) / sigma))
    return(mu + sigma * stats::qnorm(stats::runif(length(mu)) * (FB-FA) + FA))
  }
  logLikTruncNorm <- function(param=c(0,1), y, a=-Inf, b=+Inf){
    ## the original implementation mentions this function, but it was absent
    ##  from the supplement, so I looked at the "truncnorm" package
    ## see https://en.wikipedia.org/wiki/Truncated_normal_distribution#Definition
    mu <- param[1]
    sigma <- param[2]
    ksi <- (y - mu) / sigma
    alpha <- (a - mu) / sigma
    beta <- (b - mu) / sigma
    Z <- stats::pnorm(beta) - stats::pnorm(alpha)
    n <- length(y)
    loglik <- rep(NA, n)
    inrange <- which(y >= a & y <= b)
    loglik[inrange] <- stats::dnorm(ksi, log=TRUE) - log(sigma) - log(Z)
    loglik[! inrange] <- 0
    return(sum(loglik))
  }
  ## check it, wihout the last sum
  ## all.equal(logLikTruncNorm(y=0.5), dnorm(0.5, log=TRUE))
  ## all.equal(logLikTruncNorm(c(1,2), y=0.5), dnorm(0.5, 1, 2, log=TRUE))
  ## library(truncnorm)
  ## all.equal(exp(logLikTruncNorm(y=0.5)), dtruncnorm(0.5))
  ## all.equal(exp(logLikTruncNorm(c(1,2), y=0.5)), dtruncnorm(0.5, mean=1, sd=2))
  ## all.equal(exp(logLikTruncNorm(y=0.5, a=-1, b=1)), dtruncnorm(0.5, -1, 1))
  ## all.equal(exp(logLikTruncNorm(c(1,2), y=0.5, a=-1, b=1)),
  ##           dtruncnorm(0.5, -1, 1, 1, 2))
  ## all.equal(exp(logLikTruncNorm(c(1,2), y=c(0.5,0.3), a=-1, b=1)),
  ##           dtruncnorm(c(0.5,0.3), -1, 1, 1, 2))

  if(verbose > 0)
    write("prepare inputs...", stdout())
  iter <- 0
  n <- length(y)
  if ((family != "gaussian") & (family != "bernoulli")) {
    stop("Only Gaussian and Bernoulli Outcomes are allowed")
  }
  isYCensored <- FALSE
  if (family == "gaussian") {
    if ((!is.null(a)) | (!is.null(b))) {
      isYCensored <- TRUE
    }
  }
  whichNa <- which(is.na(y))
  nNa <- length(whichNa)
  hasWeights <- !is.null(weights)
  if (!hasWeights) {
    weights <- rep(1, n)
  }
  else {
    if (isYCensored | (family == "bernoulli")) {
      stop("Weights are implemented only for Non-Censored Gaussian")
    }
  }
  sumW2 <- sum(weights^2)
  if (family == "gaussian") {
    if (!isYCensored) {
      mu <- stats::weighted.mean(x = y, w = weights, na.rm = TRUE)
      yStar <- y * weights
      yStar[whichNa] <- mu * weights[whichNa]
    }
    else {
      tmp <- c(mean(y, na.rm = TRUE), stats::sd(y, na.rm = TRUE))
      tmp <- stats::optim(fn = logLikTruncNorm, par = tmp, y = y,
                          a = a, b = b, control = list(fnscale = -1))
      mu <- tmp$par[1]
      varE <- (tmp$par[2]^2)/2
      yStar <- y
      tmp <- rtrun(sigma = sqrt(varE), mu = mu, a = a[whichNa],
                   b = b[whichNa])
      yStar[whichNa] <- tmp
    }
  }
  if (family == "bernoulli") {
    threshold <- 500
    pSuccess <- stats::weighted.mean(x = y, w = weights, na.rm = TRUE)
    mu <- stats::qnorm(sd = 1, mean = 0, p = pSuccess)
    tmpY <- y
    tmpY[whichNa] <- ifelse(stats::rnorm(n = nNa, sd = 1, mean = mu) >
                            0, 1, 0)
    a <- ifelse(tmpY == 0, -threshold, 0)
    b <- ifelse(tmpY == 0, 0, threshold)
    yStar <- rtrun(sigma = 1, mu = mu, a = a, b = b)
  }
  e <- (yStar - weights * mu)
  if (family == "bernoulli") {
    varE <- 1
  }
  else {
    varE <- as.numeric(stats::var(e, na.rm = TRUE)/2)
  }
  sdE <- sqrt(varE)
  nK <- length(K)
  postMu <- 0
  postVarE <- 0
  postYHat <- rep(0, n)
  postYHat2 <- rep(0, n)
  postLogLik <- 0
  hasXF <- !is.null(XF)
  if (hasXF) {
    for (i in 1:n) {
      XF[i, ] <- weights[i] * XF[i, ]
    }
    SVD.XF <- svd(XF)
    SVD.XF$Vt <- t(SVD.XF$v)
    SVD.XF <- SVD.XF[-3]
    pF0 <- length(SVD.XF$d)
    pF <- ncol(XF)
    bF0 <- rep(0, pF0)
    bF <- rep(0, pF)
    namesBF <- colnames(XF)
    post_bF <- bF
    post_bF2 <- bF
    rm(XF)
  }
  for (k in 1:nK) {
    if (hasWeights & is.null(K[[k]]$K)) {
      stop("If weights are provided K is needed")
    }
    if (hasWeights) {
      T <- diag(weights)
      K[[k]]$K <- T %*% K[[k]]$K %*% T
    }
    if (is.null(K[[k]]$V)) {
      K[[k]]$K <- as.matrix(K[[k]]$K)
      tmp <- eigen(K[[k]]$K)
      K[[k]]$V <- tmp$vectors
      K[[k]]$d <- tmp$values
    }
    if (is.null(K[[k]]$tolD)) {
      K[[k]]$tolD <- 1e-12
    }

    tmp <- K[[k]]$d > K[[k]]$tolD
    K[[k]]$levelsU <- sum(tmp)
    K[[k]]$d <- K[[k]]$d[tmp]
    K[[k]]$V <- K[[k]]$V[, tmp]
    if (is.null(K[[k]]$df0)) {
      K[[k]]$df0 = 1
    }
    if (is.null(K[[k]]$S0)) {
      K[[k]]$S0 = 1e-15
    }
    K[[k]]$varU <- varE/nK/2
    K[[k]]$u <- rep(0, n)
    K[[k]]$uStar <- rep(0, K[[k]]$levelsU)
    K[[k]]$postVarU <- 0
    K[[k]]$postUStar <- rep(0, K[[k]]$levelsU)
    K[[k]]$postU <- rep(0, n)
    K[[k]]$postCumMSa <- rep(0, K[[k]]$levelsU)
    K[[k]]$postH1 <- rep(0, K[[k]]$levelsU)
  }

  if(verbose > 0)
    write("perform inference...", stdout())
  if(! is.null(saveAt)){
    tmp <- c("logLik", "mu", "varE", paste0("varU", 1:nK))
    if(hasXF)
      tmp <- c(tmp, paste0("bF", 1:pF))
    write(tmp, ncolumns=length(tmp), file=saveAt, append=FALSE,  sep=" ")
  }
  if(verbose > 0){
    tmpOut <- paste0("iter: ", 0,
                     "; time: ", round(0, 4),
                     "; varE: ", round(varE, 3),
                     "; mu: ", round(mu, 3), "\n")
    cat(tmpOut)
  }
  time <- proc.time()[3]
  for (i in 1:nIter) {
    if (hasXF) {
      sol <- (crossprod(SVD.XF$u, e) + bF0)
      tmp <- sol + stats::rnorm(n = pF0, sd = sqrt(varE))
      bF <- crossprod(SVD.XF$Vt, tmp/SVD.XF$d)
      e <- e + SVD.XF$u %*% (bF0 - tmp)
      bF0 <- tmp
    }
    for (k in 1:nK) {
      e <- e + K[[k]]$u
      rhs <- crossprod(K[[k]]$V, e)/varE
      varU <- K[[k]]$varU * K[[k]]$d
      C <- as.numeric(1/varU + 1/varE)
      SD <- 1/sqrt(C)
      sol <- rhs/C
      tmp <- stats::rnorm(n = K[[k]]$levelsU, sd = SD, mean = sol)
      K[[k]]$uStar <- tmp
      K[[k]]$u <- as.vector(K[[k]]$V %*% tmp)
      e <- e - K[[k]]$u
      tmp <- K[[k]]$uStar/sqrt(K[[k]]$d)
      S <- as.numeric(crossprod(tmp)) + K[[k]]$S0
      df <- K[[k]]$levelsU + K[[k]]$df0
      K[[k]]$varU <- S/stats::rchisq(n = 1, df = df)
    }
    e <- e + weights * mu
    rhs <- sum(weights * e)/varE
    C <- sumW2/varE
    sol <- rhs/C
    mu <- stats::rnorm(n = 1, sd = sqrt(1/C), mean = sol)
    e <- e - weights * mu
    if (family != "bernoulli") {
      df <- n + df0
      SS <- as.numeric(crossprod(e)) + S0
      varE <- SS/stats::rchisq(n = 1, df = df)
    }
    yHat <- yStar - e
    if (family == "bernoulli") {
      yStar <- rtrun(mu = yHat, a = a, b = b, sigma = 1)
      e <- yStar - yHat
    }
    if (nNa > 0) {
      if (family == "gaussian") {
        if (isYCensored) {
          yStar[whichNa] <- rtrun(mu = yHat[whichNa],  a = a[whichNa], b = b[whichNa], sigma = sdE)
        }
        else {
          yStar[whichNa] <- yHat[whichNa] + stats::rnorm(n = nNa,  sd = sdE)
        }
        e[whichNa] <- yStar[whichNa] - yHat[whichNa]
      }
      if (family == "bernoulli") {
        pSuccess <- stats::pnorm(mean = 0, q = yHat[whichNa], sd = 1)
        tmp <- stats::rbinom(size = 1, n = nNa, prob = pSuccess)
        a[whichNa] <- ifelse(tmp == 0, -threshold, 0)
        b[whichNa] <- ifelse(tmp == 0, 0, threshold)
      }
    }
    if (family == "gaussian") {
      tmpE <- e/weights
      tmpSD <- sqrt(varE)/weights
      if (nNa > 0) {
        tmpE <- tmpE[-whichNa]
        tmpSD <- tmpSD[-whichNa]
      }
      logLik <- sum(stats::dnorm(tmpE, sd = tmpSD, log = TRUE))
      if (isYCensored) {
        cdfA <- stats::pnorm(q = a[whichNa], sd = sdE, mean = yHat[whichNa])
        cdfB <- stats::pnorm(q = b[whichNa], sd = sdE, mean = yHat[whichNa])
        logLik <- logLik + sum(log(cdfB - cdfA))
      }
    }
    if (family == "bernoulli") {
      pSuccess <- stats::pnorm(mean = 0, sd = 1, q = yHat[-whichNa])
      tmp <- y[-whichNa]
      logLik <- sum(log(ifelse(y[-whichNa] == 0, (1 - pSuccess), pSuccess)))
    }
    if ((i > burnIn) & (i%%thin == 0)) {
      iter <- iter + 1
      constant <- (iter - 1)/(iter)
      postMu <- postMu * constant + mu/iter
      postVarE <- postVarE * constant + varE/iter
      postYHat <- postYHat * constant + yHat/iter
      for (k in 1:nK) {
        K[[k]]$postVarU <- K[[k]]$postVarU * constant +  K[[k]]$varU/iter
        K[[k]]$postUStar <- K[[k]]$postUStar * constant + K[[k]]$uStar/iter
        K[[k]]$postU <- K[[k]]$postU * constant + K[[k]]$u/iter
        tmp <- cumsum(K[[k]]$uStar^2) / K[[k]]$levelsU # see eq 20
        K[[k]]$postCumMSa <- K[[k]]$postCumMSa * constant + tmp/iter
        tmp <- K[[k]]$uStar^2 > mean(K[[k]]$uStar^2) # see eq 23
        K[[k]]$postH1 <- K[[k]]$postH1 * constant + tmp/iter
      }
      postLogLik <- postLogLik * constant + logLik/iter
      if (hasXF) {
        post_bF <- post_bF * constant + bF/iter
        post_bF2 <- post_bF2 * constant + (bF^2)/iter
      }
    }
    tmp <- i%%thin == 0
    if (! is.null(saveAt) & i > burnIn & (tmp)) {
      tmp <- c(logLik, mu, varE)
      for(k in 1:nK)
        tmp <- c(tmp, K[[k]]$varU)
      if(hasXF)
        tmp <- bF
      write(tmp, ncolumns=length(tmp), file=saveAt, append=TRUE,  sep=" ")
    }
    tmp <- proc.time()[3]
    if(verbose > 0){
      if(verbose == 1 & i %% thin != 0)
        next
      tmpOut <- paste0("iter: ", i,
                       "; time: ", round(tmp - time, 4),
                       "; varE: ", round(varE, 3), "\n")
      cat(tmpOut)
    }
    time <- tmp
  }

  if(verbose > 0)
    write("prepare the output...", stdout())
  out <- list(mu = postMu, fit = list(), varE = postVarE,
              yHat = as.numeric(postYHat/weights),
              weights = weights, K = list(), whichNa = whichNa, df0 = df0,
              S0 = S0, nIter = nIter, burnIn = burnIn, saveAt = saveAt)
  out$fit$postMeanLogLik <- postLogLik
  if (family == "gaussian") {
    tmpE <- (yStar - postYHat)/weights
    tmpSD <- sqrt(postVarE)/weights
    if (nNa > 0) {
      tmpE <- tmpE[-whichNa]
      tmpSD <- tmpSD[-whichNa]
    }
    out$fit$logLikAtPostMean <- sum(stats::dnorm(tmpE, sd = tmpSD,
                                                 log = TRUE))
    if (isYCensored) {
      cdfA <- stats::pnorm(q = a[whichNa], sd = sqrt(postVarE),
                    mean = postYHat[whichNa])
      cdfB <- stats::pnorm(q = b[whichNa], sd = sqrt(postVarE),
                    mean = postYHat[whichNa])
      out$fit$logLikAtPostMean <- out$fit$logLikAtPostMean +
        sum(log(cdfB - cdfA))
    }
  }
  if (family == "bernoulli") {
    pSuccess <- stats::pnorm(mean = 0, sd = 1, q = postYHat[-whichNa])
    tmp <- y[-whichNa]
    out$fit$logLikAtPostMean <- sum(log(ifelse(y[-whichNa] ==
                                               0, (1 - pSuccess), pSuccess)))
  }
  out$fit$pD <- -2 * (out$fit$postMeanLogLik - out$fit$logLikAtPostMean)
  out$fit$DIC <- out$fit$pD - 2 * out$fit$postMeanLogLik
  if (hasXF) {
    out$bF <- as.vector(post_bF)
    out$SD.bF <- as.vector(sqrt(post_bF2 - post_bF^2))
    names(out$bF) <- namesBF
    names(out$SD.bF) <- namesBF
  }
  for (k in 1:nK) {
    out$K[[k]] <- list()
    out$K[[k]]$u <- K[[k]]$postU
    out$K[[k]]$uStar <- K[[k]]$postUStar
    out$K[[k]]$varU <- K[[k]]$postVarU
    out$K[[k]]$cumMSa <- K[[k]]$postCumMSa
    out$K[[k]]$probH1 <- K[[k]]$postH1
    out$K[[k]]$dfo <- K[[k]]$df0
    out$K[[k]]$S0 <- K[[k]]$S0
    out$K[[k]]$tolD <- K[[k]]$tolD
  }
  return(out)
}

##' Falconer's parameterization
##'
##' Plot Falconer's parameterization for a bi-allelic locus with allele A1 and A2.
##' @param a additive effect (gene action) of allele A2
##' @param d dominance effect
##' @param f frequency of allele A2
##' @param version higher version means more details; it is pedagogically relevant to show the first plot and then the others
##' @return invisible list
##' @author Timothee Flutre
##' @examples
##' \dontrun{
##' plotFalconer(a=1, d=0.75, f=0.25, version=1)
##' plotFalconer(a=1, d=0.75, f=0.25, version=2)
##' }
##' @export
plotFalconer <- function(a=1, d=0.75, f=0.25, version=2){
  stopifnot(version %in% c(1, 2))

  ## statistical/populational effects (substitution and deviation)
  (alpha <- a + (1 - 2*f) * d)
  (alpha1 <- -f * alpha)
  stopifnot(all.equal(-f * a - f * (1 - 2*f) * d,
                      alpha1))
  (alpha2 <- (1 - f) * alpha)
  stopifnot(all.equal((1 - f) * a + (1 - f) * (1 - 2*f) * d,
                      alpha2))
  EGij <- (2*f - 1) * a + 2*f*(1 - f) * d
  (mu <- EGij + 2 * alpha1)
  GijA <- stats::setNames(c((0 - 2*f) * alpha,
                          (1 - 2*f) * alpha,
                          (2 - 2*f) * alpha),
                          c("A1A1","A1A2","A2A2"))
  GijD <- stats::setNames(c(-2*f^2*d,
                            2*f*(1-f)*d,
                            -2*(1-f)^2*d),
                          c("A1A1","A1A2","A2A2"))
  (Gij <- data.frame(N2=c(0, 1, 2),
                     y=c(EGij + GijA["A1A1"] + GijD["A1A1"],
                         EGij + GijA["A1A2"] + GijD["A1A2"],
                         EGij + GijA["A2A2"] + GijD["A2A2"]),
                     row.names=c("A1A1","A1A2","A2A2")))
  stopifnot(all.equal(Gij["A1A1","y"], -a),
            all.equal(Gij["A1A2","y"], d),
            all.equal(Gij["A2A2","y"], a))
  stopifnot(all.equal(EGij + GijA["A1A1"] + GijD["A1A1"],
                      mu + alpha * 0 + -2*f^2*d, check.attributes=FALSE),
            all.equal(EGij + GijA["A1A2"] + GijD["A1A2"],
                      mu + alpha * 1 + 2*f*(1-f)*d, check.attributes=FALSE),
            all.equal(EGij + GijA["A2A2"] + GijD["A2A2"],
                      mu + alpha * 2 + -2*(1-f)^2*d, check.attributes=FALSE))
  fit <- stats::lm(y ~ N2, data=Gij)
  (BVs <- stats::fitted(fit))

  ## plot
  graphics::par(mar=c(4, 2, 3, 4) + 0.1)
  graphics::plot(x=Gij$N2, y=Gij$y, las=1, pch=19, cex=1.8,
                 xlim=c(-0.2, 2.2), ylim=c(-1.2*a, 1.4*a),
                 xaxt="n", yaxt="n", xlab="", ylab="",
                 main=paste0("a=", a, ", d=", d, ", f=", f))
  graphics::axis(side=1, at=Gij$N2,
                 labels=c(expression("N"[2]=="0"),
                          expression("N"[2]=="1"),
                          expression("N"[2]=="2")))
  graphics::axis(side=2, at=c(-a, 0, d, a), las=1, labels=c("-a", "0", "d", "+a"))
  graphics::mtext(expression("A"[1]*"A"[1]), side=1, line=2, outer=FALSE, at=0)
  graphics::mtext(expression((1-f)^2), side=1, line=3, outer=FALSE, at=0)
  graphics::mtext(expression("A"[1]*"A"[2]), side=1, line=2, outer=FALSE, at=1)
  graphics::mtext(expression(2*f(1-f)), side=1, line=3, outer=FALSE, at=1)
  graphics::mtext(expression("A"[2]*"A"[2]), side=1, line=2, outer=FALSE, at=2)
  graphics::mtext(expression(f^2), side=1, line=3, outer=FALSE, at=2)
  graphics::text(x=Gij$N2, y=Gij$y, pos=c(1, 3, 1), offset=1,
                 labels=c(expression("G"[11]), expression("G"[12]),
                          expression("G"[22])))
  if(version == 1){
    graphics::abline(h=0, lty=1)
    graphics::segments(x0=-3, y0=Gij["A1A1","y"],
                       x1=Gij["A1A1","N2"], y1=Gij["A1A1","y"],
                       lty=3)
    graphics::segments(x0=-3, y0=Gij["A1A2","y"],
                       x1=Gij["A1A2","N2"], y1=Gij["A1A2","y"],
                       lty=3)
    graphics::segments(x0=-3, y0=Gij["A2A2","y"],
                       x1=Gij["A2A2","N2"], y1=Gij["A2A2","y"],
                       lty=3)
    ## E[G_ij]:
    graphics::abline(h=EGij, lty=2)
    graphics::axis(side=4, at=EGij, las=1,
                   labels=expression(E*"["*G[ij]*"]"))
  }
  if(version >= 2){
    graphics::abline(fit, col="red")
    ## add breeding values:
    graphics::points(x=Gij$N2, y=BVs, pch=21, bg="white", cex=1.8)
    graphics::text(x=Gij$N2, y=BVs, pos=c(2, 4, 2), offset=1,
                   labels=c(expression("G"[11*",A"]),
                            expression("G"[12*",A"]),
                            expression("G"[22*",A"])))
    graphics::axis(side=4, at=c(BVs, 0), las=1,
                   labels=c(expression("(0-2f)"*alpha),
                            expression("(1-2f)"*alpha),
                            expression("(2-2f)"*alpha),
                            "pop.\nmean"))
    ## dominance deviations:
    graphics::segments(x0=c(0, 1, 2), y0=Gij$y,
                       x1=c(0, 1, 2), y1=BVs, lty=3)
    graphics::text(x=Gij$N2, y=BVs + (residuals(fit) / 2),
                   pos=c(4, 2, 4), offset=0.2,
                   labels=c(expression(delta[11]),
                            expression(delta[12]),
                            expression(delta[22])))
    ## population mean:
    x_popmean <- -stats::coef(fit)["(Intercept)"] / stats::coef(fit)["N2"]
    graphics::points(x=x_popmean,
                     y=0, pch=3, cex=2)
    graphics::segments(x0=x_popmean, y0=0,
                       x1=3, y1=0,
                       lty=4)
    ## add slope = alpha:
    graphics::segments(x0=0, y0=BVs["A1A1"], x1=1, y1=BVs["A1A1"], lty=2)
    graphics::arrows(x0=1, y0=BVs["A1A1"], x1=1, y1=BVs["A1A2"], code=3)
    graphics::text(x=1+0.1, y=BVs["A1A1"]+(BVs["A1A2"]-BVs["A1A1"])/2,
                   labels=expression(alpha), cex=1.2)
  }

  invisible(list(a=a, d=d,
                 alpha=alpha, alpha1=alpha1, alpha2=alpha2,
                 EGij=EGij, mu=mu,
                 GijA=GijA, GijD=GijD, Gij=Gij,
                 fit=fit, BVs=BVs))
}

##' Genotype frequencies
##'
##' Calculate the genotype frequencies at a bi-allelic SNP, from its minor allele frequency, assuming Hardy-Weinberg equilibrium.
##' @param maf frequency of the minor allele, a
##' @return vector of genotype frequencies
##' @author Timothee Flutre
##' @examples
##' \dontrun{set.seed(1859)
##' genos <- sample(x=c(0,1,2), size=n, replace=TRUE, prob=maf2genoFreq(maf))
##' table(genos)
##' }
##' @export
maf2genoFreq <- function(maf){
  stopifnot(is.numeric(maf), length(maf) == 1, maf >= 0, maf <= 0.5)
  geno.freq <- c((1 - maf)^2,
                 2 * (1 - maf) * maf,
                 maf^2)
  names(geno.freq) <- c("AA", "Aa", "aa")
  return(geno.freq)
}

##' Proportion of variance explained
##'
##' Computes the additive effect of a bi-allelic SNP (beta) given its PVE, MAF and error standard deviation for the simple linear regression model: for all i in {1,...,n}, y_i = mu + beta * x_i + epsilon_i with epsilon_i ~ N(0, sigma^2).
##' Indeed, for this model: var(y) = beta^2 var(x) + sigma^2.
##' Assuming Hardy-Weinberg equilibrium: x ~ Binomial(2, maf); hence: var(x) = 2 f (1 - f).
##' Moreover: PVE = var(beta x) / var(y).
##' As a consequence, by fixing the PVE, the MAF and sigma, we can deduce a value for beta.
##' @param pve proportion of variance explained
##' @param maf minor allele frequency
##' @param sigma error standard deviation
##' @return numeric
##' @author Timothee Flutre
##' @examples
##' \dontrun{## compare two different PVEs
##' pve2beta(pve=0.7, maf=0.3, sigma=1)
##' pve2beta(pve=0.2, maf=0.3, sigma=1)
##'
##' ## compare two different MAFs, depending on the PVE
##' pve2beta(pve=0.7, maf=0.4, sigma=1)
##' pve2beta(pve=0.7, maf=0.1, sigma=1)
##' pve2beta(pve=0.2, maf=0.4, sigma=1)
##' pve2beta(pve=0.2, maf=0.1, sigma=1)
##' }
##' @export
pve2beta <- function(pve=0.7, maf=0.3, sigma=1){
  var.genos <- maf * (1 - maf)
  beta <- sigma * sqrt(pve / ((1 - pve) * var.genos))
  return(beta)
}

##' Asymptotic Bayes factor
##'
##' Calculate the asymptotic Bayes factor proposed by Wakefield in Genetic Epidemiology 33:79-86 (2009, \url{http://dx.doi.org/10.1002/gepi.20359}).
##' Can also be averaged over a grid of values of W, as done in various papers from Matthew Stephens' lab.
##' @param theta.hat vector of MLE(s) of the additive genetic effect(s)
##' @param V vector of the corresponding variance(s) of \code{theta.hat}
##' @param W vector of variance(s) of the prior on theta (if several values, the ABF will be averaged over them); the vector of default values comes from the single-SNP R implementation of BLMM by Wen for his 2015 article (see https://github.com/xqwen/blmm)
##' @param weights weights used to average over the grid of \code{W} (all equal by default)
##' @param log10 return the log10 of the ABF
##' @return numeric
##' @author Timothee Flutre
##' @export
calcAsymptoticBayesFactorWakefield <- function(theta.hat, V,
                                               W=c(0.1, 0.2, 0.4, 0.8, 1.6),
                                               weights=NULL,
                                               log10=TRUE){
  stopifnot(length(V) == length(theta.hat))
  if(! is.null(names(theta.hat)) & ! is.null(names(V)))
    stopifnot(all(names(theta.hat) == names(V)))

  z2 <- theta.hat^2 / V # Wald statistic

  if(length(theta.hat) == 1){
    log10.ABF <- 0.5 * log10(V) - 0.5 * log10(V + W) +
      (0.5 * z2 * W / (V + W)) / log(10)
  } else{
    log10.ABF <- sapply(seq_along(theta.hat), function(i){
      tmp <- sapply(W, function(Wj){
        calcAsymptoticBayesFactorWakefield(theta.hat[i], V[i], Wj)
      })
      log10WeightedSum(x=tmp, weights=weights)
    })
  }

  if(! is.null(names(theta.hat)))
    names(log10.ABF) <- names(theta.hat)

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
##' @examples
##' \dontrun{## make fake data
##' set.seed(1859)
##' n <- 200
##' mu <- 50
##' maf <- 0.3
##' genos <- sample(x=c(0,1,2), size=n, replace=TRUE, prob=maf2genoFreq(maf))
##' (beta <- pve2beta(pve=0.4, n=n, maf=maf, sigma=1))
##' phenos <- mu + beta * genos + rnorm(n=n, mean=0, sd=1)
##'
##' boxplot(phenos ~ genos, las=1, xlab="genotypes", ylab="phenotypes", at=0:2)
##'
##' ## perform inference via maximum likelihood
##' (fit <- lm(phenos ~ genos))
##' abline(fit)
##'
##' ## compute the Bayes factors
##' (BF <- calcExactBayesFactorServinStephens(G=genos, Y=phenos, sigma.a=0.5,
##'                                           sigma.d=NULL))
##' (aBF <- calcAsymptoticBayesFactorWakefield(theta.hat=coef(fit)[2],
##'                                            V=diag(vcov(fit))["genos"],
##'                                            W=0.5^2))
##' }
##' @export
calcExactBayesFactorServinStephens <- function(G, Y, sigma.a, sigma.d,
                                               log10=TRUE){
  stopifnot(is.vector(G), is.vector(Y))

  subset <- stats::complete.cases(Y) & stats::complete.cases(G)
  Y <- Y[subset]
  G <- G[subset]
  stopifnot(length(Y) == length(G))

  N <- length(G)
  if(is.null(sigma.d)){
    X <- cbind(rep(1,N), G)
    inv.Sigma.B <- diag(c(0, 1/sigma.a^2))
  } else{
    X <- cbind(rep(1,N), G, G == 1)
    inv.Sigma.B <- diag(c(0, 1/sigma.a^2, 1/sigma.d^2))
  }
  inv.Omega <- inv.Sigma.B + t(X) %*% X
  inv.Omega0 <- N
  tY.Y <- t(Y) %*% Y
  log10.BF <- as.numeric(0.5 * log10(inv.Omega0) -
                         0.5 * log10(det(inv.Omega)) -
                         ifelse(is.null(sigma.d),
                                log10(sigma.a),
                                log10(sigma.a) + log10(sigma.d)) -
                         (N/2) * (log10(tY.Y - t(Y) %*% X %*% solve(inv.Omega)
                                        %*% t(X) %*% cbind(Y)) -
                                  log10(tY.Y - N*mean(Y)^2)))

  if(log10)
    return(log10.BF)
  else
    return(10^log10.BF)
}

##' Grid for Bayes Factors
##'
##' Makes the grid of prior variances used to compute Bayes factors as in Wen and Stephens (Annals of Applied Statistics, 2014).
##' @param grid.type "general" indicates the meta-analysis grid (large), otherwise the configuration grid (small) is returned
##' @param no.het if TRUE, the grid is built without heterogeneity
##' @return matrix with two columns named "phi2" and "oma2"
##' @author Timothee Flutre
##' @seealso \code{\link{calcL10ApproximateBayesFactorWenStephens}}
##' @export
makeGridWenStephens <- function(grid.type="general", no.het=FALSE){
  stopifnot(is.character(grid.type),
            is.logical(no.het))

  oma2.plus.phi2 <- c(0.1^2, 0.2^2, 0.4^2, 0.8^2, 1.6^2) # avg eff size
  oma2.over.oma2.plus.phi2 <- c(0, 1/4, 1/2, 3/4, 1) # homogeneity

  if(grid.type != "general"){
    if(no.het){
      oma2.over.oma2.plus.phi2 <- c(1)
    } else
      oma2.over.oma2.plus.phi2 <- c(3/4, 1)
  }

  grid <- matrix(data=NA,
                 nrow=length(oma2.plus.phi2) *
                   length(oma2.over.oma2.plus.phi2),
                 ncol=2)
  colnames(grid) <- c("phi2", "oma2")
  i <- 1
  for(aes in oma2.plus.phi2){
    for(hom in oma2.over.oma2.plus.phi2){
      grid[i,"phi2"] <- aes * (1 - hom)
      grid[i,"oma2"] <- aes * hom
      i <- i + 1
    }
  }

  return(grid)
}

##' Approximate Bayes factor
##'
##' Calculate the log10(ABF) of Wen & Stephens in Annals of Applied Statistics (2014, \url{http://dx.doi.org/10.1214/13-AOAS695}) according to the "exchangeable standardized effects" model.
##' @param sstats matrix of summary statistics with one row per subgroup and three columns, "bhat", "sebhat" and "t"
##' @param phi2 prior variance of the \eqn{b_s} given \eqn{\bar{b}}; controls the prior expected degree of heterogeneity among subgroups
##' @param oma2 prior variance of \eqn{\bar{b}}; controls the prior expected size of the average effect across subgroups
##' @return numeric
##' @author Xiaoquan Wen [aut], Timothee Flutre [ctb,cre]
##' @seealso \code{\link{makeGridWenStephens}}, \code{\link{calcL10ApproximateBayesFactorWen}}
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

##' Approximate Bayes factor
##'
##' Calculate the log10(ABF) under the SSMR model (equation 11 with s=1) of \href{http://dx.doi.org/10.1111/biom.12112}{Wen (Biometrics, 2013)}.
##' @param Y n x r phenotype matrix
##' @param Xg n x p genotype matrix
##' @param Xc n x q design matrix of covariates
##' @param Wg r x r prior var-covar matrix on beta_g (if NULL, requires phi and oma)
##' @param phi used to make Wg
##' @param oma used to make Wg
##' @param alpha mixing proportion
##' @param H r x r matrix parameter of Wishart prior
##' @param nu positive scalar parameter of Wishart prior
##' @param ES if TRUE, use the scale-invariant prior formulation
##' @return numeric
##' @author Xiaoquan Wen [aut], Timothee Flutre [ctb,cre]
##' @seealso \code{\link{calcL10ApproximateBayesFactorWenStephens}}
##' @export
calcL10ApproximateBayesFactorWen <- function(Y, Xg, Xc,
                                             Wg=NULL, phi=NULL, oma=NULL,
                                             alpha=0, H=0, nu=0, ES=TRUE){
  stopifnot(all(is.matrix(Y), is.matrix(Xg), is.matrix(Xc)),
            ! is.null(Wg) || (! is.null(phi) && ! is.null(oma)),
            is.logical(ES))

  l10abf <- NA

  n <- dim(Y)[1]
  r <- dim(Y)[2]
  q <- dim(Xc)[2]
  p <- dim(Xg)[2]
  ## message(paste0("n=",n, " r=",r, " q=",q, " p=",p))

  make_Wg <- function(phi, omega, p, r){
    phi2 <- phi^2
    omg2 <- omega^2
    W <- matrix(ncol=r, rep(omg2,r*r)) + diag(rep(phi2,r))
    Wg <- diag(rep(1,p)) %x% W
    return(Wg)
  }
  if(is.null(Wg))
    Wg <- make_Wg(phi=phi, omega=oma, p=p, r=r)

  X <- cbind(Xc,Xg)

  Sigma_hat <- Sigma_tilde <- diag(rep(0,r))

  if(H == 0){
    H <- diag(rep(1e-8,r))
  }

  if(alpha > 0){
    Sigma_hat <- (t(Y)%*%(diag(rep(1,n)) - X%*%mpInv(t(X)%*%X)%*%t(X))%*%Y)/n
  }

  if(alpha<1){
    Sigma_tilde <- (t(Y)%*%(diag(rep(1,n)) - Xc%*%solve(t(Xc)%*%Xc)%*%t(Xc))%*%Y)/n
  }

  Sigma <- (nu/(nu+n))*H + (n/(n+nu))*(alpha*Sigma_hat + (1-alpha)*Sigma_tilde)
  Sigma_inv <- solve(Sigma)

  if(ES){
    S <- diag(rep(1,p))%x%diag(sqrt(diag(Sigma)))
    Wg <- S%*%Wg%*%S
  }

  Vg_inv <- (t(Xg)%*%Xg - t(Xg)%*%Xc%*%solve(t(Xc)%*%Xc)%*%t(Xc)%*%Xg)%x%Sigma_inv
  vec <- matrix(ncol=1, as.vector(t(Y-Xc%*%solve(t(Xc)%*%Xc)%*%t(Xc)%*%Y)))
  bVi <- t(vec)%*%(Xg%x%Sigma_inv)
  ivw <- diag(rep(1,p*r))+Vg_inv%*%Wg
  l10abf <- (.5*bVi%*%Wg%*%solve(ivw)%*%t(bVi)-0.5*determinant(ivw)$modulus[[1]])/log(10)

  return(as.numeric(l10abf))
}

##' Boxplot of QTL
##'
##' Make a boxplot of a candidate QTL.
##' @param y vector of phenotypes with genotype names
##' @param X matrix of bi-allelic SNP genotypes encoded in allele doses in [0,2], with genotypes in rows and SNPs in columns; missing values should be encoded as NA
##' @param snp character with SNP name corresponding to the candidate QTL to plot
##' @param simplify.imputed if TRUE, imputed genotypes will be transformed back to {0,1,2}
##' @param xlab label of the x-axis
##' @param maf.xlab if TRUE, the minor allele frequency will appear in the label of the x-axis
##' @param ylab label of the y-axis
##' @param main.title main title
##' @param show.points if TRUE, individual points will be shown, with \code{\link{jitter}}, especially useful if some genotypic classes have very low counts
##' @param jit.fact jitter factor used if \code{show.points} is TRUE
##' @param varwidth if TRUE, the boxes are drawn with widths proportional to the square-roots of the number of observations in the groups
##' @param notch if TRUE, a notch is drawn in each side of the boxes (see \code{\link[graphics]{boxplot}})
##' @param suppress.warnings if TRUE, \code{\link{suppressWarnings}} is used for \code{\link[graphics]{boxplot}}
##' @param regline.intercept intercept of the regression line
##' @param regline.slope slope of the regression line
##' @param regline.col color of the regression line
##' @param regline.lty style of the regression line
##' @param regline.legend legend of the regression line
##' @param alleles vector of characters of length 2, the second element being the allele which copies are counted in \code{X}
##' @param title.line line at which the main title should appear
##' @param counts.xticks if TRUE, the sample size of each genotypic class will be added to the x-axis ticks
##' @param mtext.y.line line at which the y-axis label should appear
##' @param mtext.y.cex cex of the y-axis label
##' @param verbose verbosity level (0/1)
##' @param ... other arguments to \code{\link[graphics]{boxplot}}
##' @return invisible list with \code{y} and \code{x} used to make the boxplot (same order, with no missing data)
##' @author Timothee Flutre
##' @examples
##' \dontrun{## simulate genotypes
##' set.seed(1859)
##' I <- 200
##' P <- 2000
##' X <- simulGenosDose(nb.genos=I, nb.snps=P)
##'
##' ## make fake SNP coordinates
##' snp.coords <- data.frame(coord=1:ncol(X),
##'                          chr="chr1",
##'                          row.names=colnames(X),
##'                          stringsAsFactors=FALSE)
##'
##' ## simulate phenotypes (only additive effects)
##' modelA <- simulBvsr(Q=1, X=X, pi=0.01, pve=0.7, sigma.a2=1)
##'
##' ## test SNPs one by one with the univariate LMM
##' fit.u <- gemma(model="ulmm", y=modelA$Y[,1], X=X, snp.coords,
##'                W=modelA$W, out.dir=tempdir(), clean="all")
##'
##' ## diagnostic plots
##' plotHistPval(pvalues=fit.u$tests$p_wald)
##' cols <- rep("black",ncol(X)); cols[modelA$gamma==1] <- "red"
##' pvadj.AA <- qqplotPval(fit.u$tests$p_wald, col=cols, ctl.fdr.bh=TRUE,
##'                        plot.signif=TRUE)
##'
##' ## look at the best candidate QTL
##' snp <- rownames(fit.u$tests[order(fit.u$tests$p_wald),])[1]
##' boxplotCandidateQtl(y=modelA$Y[,1], X=X, snp=snp, main=snp, notch=FALSE,
##'                     show.points=TRUE)
##' fit <- lm(modelA$Y[,1] ~ X[,snp])
##' abline(fit, col="red")
##' abline(a=fit.u$global.mean["beta.hat"] - mean(X[snp]) * fit.u$tests[snp,"beta"],
##'        b=fit.u$tests[snp,"beta"])
##' legend("topright", legend=c("lm","gemma"), col=c("red","black"), lty=1, bty="n")
##' }
##' @export
boxplotCandidateQtl <- function(y, X, snp, simplify.imputed=TRUE,
                                xlab="SNP genotypes", maf.xlab=TRUE,
                                ylab="Phenotypes",
                                main.title=NULL,
                                show.points=FALSE, jit.fact=1,
                                varwidth=TRUE, notch=TRUE, suppress.warnings=TRUE,
                                regline.intercept=NA, regline.slope=NA,
                                regline.col="red", regline.lty=1,
                                regline.legend=NULL,
                                alleles=NULL,
                                title.line=NA,
                                counts.xticks=FALSE,
                                mtext.y.line=3,
                                mtext.y.cex=1,
                                verbose=1, ...){
  if(! is.vector(y)){
    if(is.matrix(y)){
      stopifnot(ncol(y) == 1,
                ! is.null(rownames(y)))
      y <- stats::setNames(y[,1], rownames(y))
    } else if(is.numeric(y) & is.atomic(y)){ # e.g. from ranef()
      stopifnot(! is.null(names(y)))
      y <- stats::setNames(as.vector(y), names(y))
    }
  }
  stopifnot(is.vector(y),
            ! is.null(names(y)),
            is.character(snp),
            snp %in% colnames(X))
  X.snp <- X[names(y), snp, drop=FALSE]
  stopIfNotValidGenosDose(X.snp, check.noNA=FALSE,
                          check.notImputed=ifelse(simplify.imputed, FALSE, TRUE))
  if(! is.null(alleles))
    stopifnot(is.vector(alleles),
              is.character(alleles),
              length(alleles) == 2,
              all(alleles %in% c("A","T","G","C")))

  ## reformat the inputs
  x <- stats::setNames(as.vector(X.snp), rownames(X.snp))
  x <- x[! is.na(x)]
  y <- y[! is.na(y)]
  ind.names <- intersect(names(y), names(x))
  x <- x[ind.names]
  y <- y[ind.names]
  if(simplify.imputed){
    boundaries <- seq(from=0, to=2, length.out=4)
    x[x <= boundaries[2]] <- 0
    x[x > boundaries[2] & x <= boundaries[3]] <- 1
    x[x > boundaries[3]] <- 2
  }

  counts <- table(x)
  af <- mean(x) / 2
  maf <- ifelse(af <= 0.5, af, 1 - af)
  if(verbose > 0){
    msg <- paste0("marker '", snp, "'")
    msg <- paste0(msg, "\ngenotypic classes:")
    for(i in seq_along(counts))
      msg <- paste0(msg, " ", names(counts)[i], "=", counts[i])
    msg <- paste0(msg, " (total=", sum(counts), ")")
    write(msg, stdout())

    msg <- sprintf("minor allele frequency: maf = %.3f", maf)
    write(msg, stdout())

    print(do.call(rbind, tapply(y, factor(x), function(tmp){
      c(mean.y=mean(tmp), med.y=stats::median(tmp), sd.y=stats::sd(tmp))
    })), digits=3)
  }

  ## make boxplot
  if(maf.xlab)
    xlab <- paste0(xlab, " (MAF=", format(maf, digits=3), ")")
  x <- data.frame(dose.num=x,
                  row.names=names(x))
  x$dose.fact <- factor(as.character(x$dose.num), levels=c(0, 1, 2))
  x$final <- x$dose.fact
  x$all <- as.character(x$dose.num)
  if(! is.null(alleles)){
    all.hom.0 <- paste(rep(alleles[1], 2), collapse="")
    all.het <- paste(c(alleles[1], alleles[2]), collapse="")
    all.hom.2 <- paste(rep(alleles[2], 2), collapse="")
    x$all[x$dose.num == 0] <- all.hom.0
    x$all[x$dose.num == 1] <- all.het
    x$all[x$dose.num == 2] <- all.hom.2
    x$all <- factor(x$all, levels=c(all.hom.0, all.het, all.hom.2))
    x$final <- x$all
  }
  if(counts.xticks){
    counts <- table(as.character(x$final))[levels(x$final)]
    for(i in seq_along(counts))
      levels(x$final)[levels(x$final) == names(counts)[i]] <-
        paste0(names(counts)[i], " (", counts[i], ")")
  }
  if(suppress.warnings){
    bp <- suppressWarnings(
        graphics::boxplot(y ~ x$final, las=1, varwidth=varwidth, notch=notch,
                          xlab=xlab, ylab="", at=c(0,1,2), ...))
  } else
    bp <- graphics::boxplot(y ~ x$final, las=1, varwidth=varwidth, notch=notch,
                            xlab=xlab, ylab="", at=c(0,1,2), ...)
  graphics::mtext(text=ylab, side=2, line=mtext.y.line, cex=mtext.y.cex)
  if(! is.null(main.title))
    graphics::title(main=main.title, line=title.line)

  if(show.points){
    for(i in 1:nlevels(x$final)){
      tmp <- y[x$final == levels(x$final)[i]]
      graphics::points(x=jitter(x=rep(i-1, length(tmp)), factor=jit.fact),
                       y=tmp)
    }
  }

  if(all(! is.na(regline.intercept), ! is.na(regline.slope))){
    graphics::abline(a=regline.intercept, b=regline.slope,
                     col=regline.col, lty=regline.lty)
    if(! is.null(regline.legend)){
      legend.location <- "bottomright"
      if(regline.slope < 0)
        legend.location <- "topright"
      graphics::legend(legend.location, legend=regline.legend,
                       col=regline.col, lty=regline.lty, bty="n")
    }
  }

  invisible(list(y=y, x=x))
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
    txt <- paste0(txt, "\n\tnb of genotypes: ", gmap$popSize)
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

##' Subset a pedigree
##'
##' From a pedigree, extract a subset of individuals.
##' @param ped data frame containing the input pedigree with three columns named "parent1", "parent2" and "child"; if "parent1" is different than "parent2", it is an allofecundation; if "parent1" is the same as "parent2" it is an autofecundation; if "parent2" is missing it is a haplodiploidization
##' @param inds vector of individual names (present in the "child" column of "ped") for which their ancestors will be retrieved, up to the founders
##' @return data frame
##' @author Julien Diot [cre,aut], Timothee Flutre [ctb]
##' @seealso \code{\link{plotPedigree}}
##' @export
subsetPedigree <- function(ped, inds){
  if(is.factor(inds))
    inds <- as.character(inds)
  stopifnot(is.data.frame(ped),
            all(c("child","parent1","parent2") %in% colnames(ped)),
            all(! is.na(ped$parent1)),
            all(! is.na(ped$child)),
            all(! duplicated(ped$child)),
            is.character(inds),
            all(inds %in% ped$child))
  if(any(ped$parent2 == "NA"))
    ped$parent2[ped$parent2 == "NA"] <- NA
  ped <- ped[, c("parent1", "parent2", "child")]
  ped <- convertFactorColumnsToCharacter(ped)

  getParents <- function(ped, inds, gen){
    out <- data.frame(parent1=rep(NA, length(inds)),
                      parent2=NA,
                      child=inds,
                      generation=gen,
                      stringsAsFactors=FALSE)
    for(i in seq_along(inds)){
      idx <- which(ped$child == inds[i])
      if(length(idx) > 0)
        out[i, 1:2] <- ped[idx, 1:2]
    }
    return(out)
  }

  ## initialisation
  gen <- 0
  out.ped <- getParents(ped, inds, gen)

  ## ascend genealogy
  inds <- unique(c(out.ped$parent1, out.ped$parent2))
  inds <- inds[! is.na(inds)]
  while(length(inds) > 0){
    gen <- gen + 1
    out.ped <- rbind(out.ped,
                     getParents(ped, inds, gen))
    inds <- unique(c(out.ped$parent1[out.ped$generation == gen],
                     out.ped$parent2[out.ped$generation == gen]))
    inds <- inds[! is.na(inds)]
  }

  ## recalculate generation, and order
  out.ped$generation <- max(out.ped$generation) - out.ped$generation
  out.ped <- out.ped[order(out.ped$generation, out.ped$parent1,
                           out.ped$parent2, out.ped$child),]

  return(out.ped)
}

##' Plot pedigree
##'
##' Plot a pedigree using the "igraph" package.
##' This function was inspired by plot.pedigree() from the "synbreed" package (under GPL-3).
##' It add options for monoecious species and auto-fecondation.
##' @param inds identifiers of the genotypes
##' @param mothers identifiers of their mother; can be NA (case of the founders and haplodiploidization)
##' @param fathers identifiers of their father; can be NA (case of the founders and haplodiploidization)
##' @param generations should start at 0
##' @param sexes "F" for female (circle), "M" for male (square) and "H" for hermaphrodite (triangle); can also be NA (no shape)
##' @param plot.it if TRUE, the pedigree will be plotted
##' @param verbose verbosity level (0/1/2)
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
##' @param xmin see \code{norm_coords} in the "igraph" package
##' @param xmax see \code{norm_coords} in the "igraph" package
##' @param ymin see \code{norm_coords} in the "igraph" package
##' @param ymax see \code{norm_coords} in the "igraph" package
##' @param ... other plotting options; see ?plot.igraph and ?igraph.plotting
##' @return invisible list with objects required to plot the pedigree
##' @author Timothee Flutre
##' @seealso \code{\link{subsetPedigree}}
##' @export
plotPedigree <- function(inds, mothers, fathers, generations, sexes=NULL,
                         plot.it=TRUE, verbose=1,
                         edge.col.mother="black", edge.col.father="darkgrey",
                         vertex.label.color="darkblue", vertex.color="white",
                         vertex.size=20, vertex.shape="none",
                         vertex.label.family="Helvetica", mult.edge.curve=0.25,
                         edge.arrow.width=0.75, edge.arrow.size=0.75,
                         xmin=-1, xmax=1, ymin=-1, ymax=1,
                         ...){
  requireNamespace("igraph", quietly=TRUE)
  if(is.factor(inds))
    inds <- as.vector(inds)
  if(is.factor(mothers))
    mothers <- as.vector(mothers)
  if(is.factor(fathers))
    fathers <- as.vector(fathers)
  stopifnot(is.vector(inds),
            is.vector(mothers),
            is.vector(fathers),
            is.vector(generations))
  if(! is.null(sexes))
    stopifnot(is.vector(sexes))

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

  ## check inds
  inds <- as.character(inds)
  if(verbose > 0){
    msg <- paste0("nb of individuals: ", length(inds))
    write(msg, stdout())
  }
  if(length(unique(inds)) != length(inds)){
    msg <- paste0("only ", length(unique(inds)), " individuals are unique")
    warning(msg)
  }

  ## check mothers
  mothers <- as.character(mothers)
  stopifnot(length(mothers) == length(inds),
            all(mothers[! is.na(mothers)] %in% inds))
  if(verbose > 0){
    uniq.mothers <- unique(mothers[! is.na(mothers)])
    msg <- paste0("nb of unique mothers: ", length(uniq.mothers))
    write(msg, stdout())
  }

  ## check fathers
  fathers <- as.character(fathers)
  stopifnot(length(fathers) == length(inds),
            all(fathers[! is.na(fathers)] %in% inds))
  if(verbose > 0){
    uniq.fathers <- unique(fathers[! is.na(fathers)])
    msg <- paste0("nb of unique fathers: ", length(uniq.fathers))
    write(msg, stdout())
  }

  ## check generations
  generations <- as.numeric(generations)
  stopifnot(length(generations) == length(inds))
  if(verbose > 0){
    msg <- paste0("nb of generations: ", length(unique(generations)))
    write(msg, stdout())
    if(verbose > 1)
      print(table(generations))
  }

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

  ## set up a data frame containing the whole pedigree
  ped.df <- data.frame(ind=inds,
                       par1=mothers,
                       par2=fathers,
                       gen=generations,
                       reprod=NA,
                       founder=FALSE,
                       stringsAsFactors=FALSE)

  ## categorize the modes of reproduction
  idx <- which(! is.na(ped.df$par1) & ! is.na(ped.df$par2) &
               ped.df$par1 != ped.df$par2)
  ped.df$reprod[idx] <- "allofecundation"
  idx <- which(! is.na(ped.df$par1) & ! is.na(ped.df$par2) &
               ped.df$par1 == ped.df$par2)
  ped.df$reprod[idx] <- "autofecundation"
  idx <- which(xor(is.na(ped.df$par1), is.na(ped.df$par2)))
  ped.df$reprod[idx] <- "haplodiploidization"
  if(verbose > 0){
    msg <- "nb of occurrences per reproduction mode:"
    cat(msg)
    print(table(ped.df$reprod, useNA="no"))
  }

  ## identify the founders
  ped.df$founder <- is.na(ped.df$par1) & is.na(ped.df$par2)
  if(verbose > 0){
    msg <- paste0("nb of founders: ", sum(ped.df$founder))
    write(msg, stdout())
    msg <- paste0("nb of non-founding generations: ",
                  length(unique(ped.df$gen[! ped.df$founder])))
    write(msg, stdout())
  }

  ## set up edges (from "ped.df")
  ## need to (1) discard founders and (2) use unique identifiers
  ## "gen_|_ind" because a given individual can be present at multiple
  ## generations, thus also need to (3) update generations
  ped.df.nofnd <- ped.df[! ped.df$founder,
                         c("ind","par1","par2","gen")]
  ped.df.nofnd$gen <- ped.df.nofnd$gen - min(ped.df.nofnd$gen) + 1
  sep <- "_|_"
  ped.df.nofnd$ind.id <- paste0("gen", ped.df.nofnd$gen,
                                sep, ped.df.nofnd$ind)
  stopifnot(all(! duplicated(ped.df.nofnd$ind.id)))
  ped.df.nofnd$par1.id <- paste0("gen", ped.df.nofnd$gen - 1,
                                 sep, ped.df.nofnd$par1)
  ped.df.nofnd$par2.id <- paste0("gen", ped.df.nofnd$gen - 1,
                                 sep, ped.df.nofnd$par2)
  relations <- data.frame(from=c(ped.df.nofnd$par1.id,
                                 ped.df.nofnd$par2.id),
                          to=rep(ped.df.nofnd$ind.id, 2),
                          type=rep(c("maternal","paternal"),
                                   each=nrow(ped.df.nofnd)),
                          stringsAsFactors=FALSE)
  relations <- relations[order(relations$from),]
  rownames(relations) <- NULL
  if(verbose > 0){
    msg <- paste0("nb of relations: ", nrow(relations))
    write(msg, stdout())
  }

  ## set up vertices (from "relations")
  tmp <- unique(c(relations$from, relations$to))
  vertices <- data.frame(name=tmp,
                         generation=as.numeric(sub("gen", "",
                                                   sapply(strsplit(tmp, "_\\|_"),
                                                          `[`, 1))),
                         stringsAsFactors=FALSE)
  if(is.null(sexes)){
    vertices[["label"]] <- sapply(strsplit(tmp, "_\\|_"), `[`, 2)
  } else{
    vertices[["sex"]] <- sexes
    vertices[["label"]] <- paste0(sapply(strsplit(tmp, "_\\|_"), `[`, 2),
                                  "\n(", sexes, ")")
  }
  vertices[["shape"]] <- vertex.shape
  vertices[["size"]] <- vertex.size
  vertices[["color"]] <- vertex.color
  vertices[["label.color"]] <- vertex.label.color
  vertices[["label.family"]] <- vertex.label.family

  ## make "igraph" object
  ped.graph <- igraph::graph_from_data_frame(d=relations,
                                             directed=TRUE,
                                             vertices=vertices)

  ## check multiplicity corresponding to auto-fecundation
  has.autof <- FALSE
  if(igraph::has.multiple(ped.graph)){
    has.autof <- TRUE
    stopifnot(all(igraph::count_multiple(ped.graph) %in% c(1,2)))
  }

  ## set plot coordinates for vertices
  vertices <- vertices[order(vertices$generation),]
  coords <- matrix(data=NA, nrow=nrow(vertices), ncol=2,
                   dimnames=list(vertices$name, c("x", "y")))
  coords[, "y"] <- max(vertices$generation) - vertices$generation
  coords[, "x"] <- order(vertices$generation) -
    cumsum(c(0, table(vertices$generation)))[vertices$generation + 1]
  coords[nrow(coords):1, "x"] <-
    unlist(tapply(coords[, "x"], coords[,"y"], function(x){
      if(length(x) == 1){
        x <- 0
      } else
        x <- rev(scale(x))
      return(x)
    }))
  coords <- igraph::norm_coords(layout=coords,
                                xmin=xmin, xmax=xmax,
                                ymin=ymin, ymax=ymax)

  ## set edge color depending on parental sex
  edge.cols <- relations$type
  edge.cols <- gsub("maternal", edge.col.mother, edge.cols)
  edge.cols <- gsub("paternal", edge.col.father, edge.cols)

  ## tune curvature in case of auto-fecondation
  edge.curvatures <- rep(0, igraph::ecount(ped.graph))
  if(has.autof){
    is.mult <- igraph::count_multiple(ped.graph) > 1
    stopifnot(sum(is.mult) %% 2 == 0)
    edge.curvatures[is.mult] <- mult.edge.curve
    is.mult.dupl <- rep(FALSE, length(is.mult))
    is.mult.dupl[is.mult] <- duplicated(relations[is.mult, "to"])
    edge.curvatures[is.mult.dupl] <- - mult.edge.curve
  }

  ## plot, finally
  if(plot.it)
    igraph::plot.igraph(x=ped.graph,
                        layout=coords, rescale=FALSE,
                        edge.color=edge.cols,
                        edge.curved=edge.curvatures,
                        edge.arrow.width=edge.arrow.width,
                        edge.arrow.size=edge.arrow.size,
                        ...)

  invisible(list(relations=relations,
                 vertices=vertices,
                 graph=ped.graph,
                 layout=coords,
                 edge.color=edge.cols,
                 edge.curved=edge.curvatures))
}
