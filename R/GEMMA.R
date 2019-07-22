## Contains functions useful to wrap GEMMA.

##' Convert genotypes
##'
##' Convert SNP genotypes to the file formats used by BimBam, as specified in its \href{http://www.haplotype.org/download/bimbam-manual.pdf}{manual}.
##' @param X matrix of bi-allelic SNP genotypes encoded in allele doses in {0,1,2}, with genotypes in rows and SNPs in columns; for each SNP, the allele which copy number is counted, should correspond to the second column of \code{alleles}, whether this column corresponds to minor or major alleles (for GEMMA, it should be the minor allele)
##' @param tX matrix with SNPs in rows and genotypes in columns
##' @param alleles data.frame with SNPs in rows (names as row names) and alleles in columns; the second column should correspond to the allele which number of copies is counted at each SNP in \code{X}
##' @param file write the genotype data to this file if not NULL (for instance 'genotypes_bimbam.txt' or 'genotypes_bimbam.txt.gz', but don't use \code{\link[base]{gzfile}})
##' @param format BimBam's format in which the data should be saved (mean/basic)
##' @return invisible data.frame
##' @author Timothee Flutre
##' @export
genoDoses2bimbam <- function(X=NULL, tX=NULL, alleles, file=NULL, format="mean"){
  stopifnot(xor(is.null(X), is.null(tX)),
            is.data.frame(alleles),
            ncol(alleles) == 2,
            ! is.null(row.names(alleles)),
            format %in% c("mean", "basic"))
  if(! is.null(X)){
    stopIfNotValidGenosDose(X, check.hasColNames=FALSE, check.hasRowNames=FALSE)
    tX <- t(X)
  }
  stopifnot(all(rownames(alleles) %in% rownames(tX)),
            all(rownames(tX) %in% rownames(alleles)))
  alleles <- convertFactorColumnsToCharacter(alleles)
  alleles <- alleles[rownames(tX),] # put in same order

  if(format == "mean"){
    out <- cbind(alleles, tX)
    if(! is.null(file))
      utils::write.table(x=out, file=file, quote=FALSE, sep="\t",
                         row.names=TRUE, col.names=FALSE)
  } else if(format == "basic"){
    out <- genoDoses2genoClasses(tX=tX, alleles=alleles, na.string="??")
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

##' GEMMA
##'
##' See Zhou & Stephens (Nature Genetics, 2012), and Zhou et al (PLoS Genetics, 2013).
##' @param model name of the model to fit (ulmm/bslmm)
##' @param y vector of phenotypes with genotype names
##' @param X matrix of bi-allelic SNP genotypes encoded, for each SNP, in number of copies of its second allele, i.e. as allele doses in {0,1,2}, with genotypes in rows and SNPs in columns; the "second" allele is arbitrary, it corresponds to the second column of \code{alleles}, which can be the minor or the major allele
##' @param snp.coords data.frame with SNPs as row names and two columns named "coord" and "chr"
##' @param recode.genos if TRUE, SNP genotypes in X will be recoded so that the minor allele is counted
##' @param alleles data.frame with SNPs in rows (names as row names) and alleles in columns (exactly 2 columns are required); the second column should correspond to the allele which number of copies is counted at each SNP in \code{X}; if NULL, fake alleles will be generated
##' @param maf SNPs with minor allele frequency strictly below this threshold will be discarded
##' @param K.c kinship matrix; if NULL, will be estimated using X via \code{\link{estimGenRel}} with \code{relationships="additive"} and \code{method="zhou"}
##' @param W matrix of covariates with genotypes in rows (names as row names), a first column of 1 for the intercept and, if needed, a second column of covariates values; if NULL, a column of 1 will be used for the intercept
##' @param weights vector of positive weights with genotype names
##' @param exe path to the executable of GEMMA
##' @param out.dir directory in which the output files will be saved
##' @param task.id identifier of the task (used in temporary and output file names)
##' @param verbose verbosity level (0/1)
##' @param clean remove files: none, some (temporary only), all (temporary and results)
##' @param seed seed for the pseudo-random number generator for GEMMA
##' @param burnin number of iterations to discard as burn-in (if model="bslmm")
##' @param nb.iters number of iterations (if model="bslmm")
##' @param thin thining (if model="bslmm")
##' @return list with the following components: recoded, cmd, log, global mean, pve, tests (and hyperparams and params for BSLMM)
##' @author Timothee Flutre [aut,cre], Dalel Ahmed [ctb]
##' @seealso \code{\link{gemmaUlmmPerChr}}
##' @note P. Carbonetto \href{https://github.com/genetics-statistics/GEMMA/issues/12#issuecomment-327008583}{showed} how the PVE of a single SNP can be obtained, assuming no covariate, x is centered and Cov(x,beta) is zero: \code{pve <- var(x) * (beta^2 + se^2)/var(y)}
##' @examples
##' \dontrun{## simulate genotypes
##' set.seed(1859)
##' I <- 200
##' P <- 2000
##' X <- simulGenosDose(nb.genos=I, nb.snps=P)
##'
##' ## make fake SNP coordinates and alleles
##' nb.snps.per.chr <- ncol(X) / 2
##' snp.coords <- data.frame(coord=rep(1:nb.snps.per.chr, 2),
##'                          chr=rep(c("chr1","chr2"), each=nb.snps.per.chr),
##'                          row.names=colnames(X),
##'                          stringsAsFactors=FALSE)
##' alleles <- data.frame(major=rep("A", ncol(X)),
##'                       minor="a",
##'                       row.names=colnames(X),
##'                       stringsAsFactors=FALSE)
##'
##' ## simulate phenotypes (only additive effects)
##' modelA <- simulBvsr(Q=1, X=X, pi=0.01, pve=0.7, sigma.a2=1)
##' summary(abs(modelA$a[modelA$gamma == 1]))
##'
##' ## test SNPs one by one with the univariate LMM
##' fit.u <- gemma(model="ulmm", y=modelA$Y[,1], X=X, snp.coords=snp.coords,
##'                alleles=alleles, W=modelA$W, out.dir=tempdir(), clean="all")
##' fit.u$global.mean
##' fit.u$pve
##' cor(modelA$a[modelA$gamma == 1], fit.u$tests$beta[modelA$gamma == 1])
##' cols <- rep("black",ncol(X)); cols[modelA$gamma==1] <- "red"
##' pvadj.AA <- qqplotPval(fit.u$tests$p_wald, col=cols, ctl.fdr.bh=TRUE,
##'                        plot.signif=TRUE)
##' t(binaryClassif(known.nulls=modelA$gamma == 0,
##'                 called.nulls=pvadj.AA$pv.bh > 0.05))
##' qtl <- names(sort(modelA$a[modelA$gamma == 1], decreasing=TRUE))[1]
##' boxplotCandidateQtl(y=modelA$Y[,1], X=X, snp=qtl, show.points=TRUE, main=qtl)
##' slope <- fit.u$tests[qtl,"beta"] *
##'   ifelse(fit.u$recode[qtl], -1, 1)
##' abline(a=fit.u$global.mean["beta.hat"] - mean(X[,qtl]) * slope,
##'        b=slope, col="red")
##'
##' ## compute the PVE of a single SNP (assuming no covariate)
##' (pve.qtl <- var(X[,qtl]) * (fit.u$tests[qtl,"beta"]^2 + fit.u$tests[qtl,"se"]^2) / var(modelA$Y[,1]))
##'
##' ## same but per chrom
##' ## fit.u2 <- gemmaUlmmPerChr(y=modelA$Y[,1], X=X, snp.coords=snp.coords,
##' ##                           alleles=alleles, W=modelA$W,
##' ##                           out.dir=tempdir(), clean="all")
##'
##' ## fit all SNPs jointly with the BSLMM
##' burnin <- 10^3
##' nb.iters <- 10^4
##' thin <- 10^2
##' fit.bs <- gemma(model="bslmm", y=modelA$Y[,1], X=X, snp.coords=snp.coords,
##'                 alleles=alleles, W=modelA$W, out.dir=tempdir(), clean="all",
##'                 burnin=burnin, nb.iters=nb.iters, thin=thin)
##' posterior.samples <- coda::mcmc(data=fit.bs$hyperparams, start=burnin + 1,
##'                                 end=burnin + nb.iters, thin=thin)
##' summary(posterior.samples)
##'
##' ## simulate phenotypes (only dominant effects)
##' set.seed(1859)
##' modelD <- simulBvsr(Q=1, X=X, pi=0.01, pve=0.7, sigma.a2=0, sigma.d2=1)
##'
##' ## test SNPs as "additive" one by one with the univariate LMM
##' A.z <- estimGenRel(X=X, relationships="additive", method="zhou")
##' fit.u.DA <- gemma(model="ulmm", y=modelD$Y[,1], X=X,
##'                   snp.coords, alleles,
##'                   K.c=A.z, W=modelD$W, out.dir=tempdir(), clean="all")
##' cols <- rep("black",ncol(X)); cols[modelD$gamma==1] <- "red"
##' pvadj.DA <- qqplotPval(fit.u.DA$tests$p_wald, col=cols, ctl.fdr.bh=TRUE,
##'                        plot.signif=TRUE)
##' t(binaryClassif(known.nulls=modelD$gamma == 0,
##'                 called.nulls=pvadj.DA$pv.bh > 0.05))
##'
##' ## test SNPs as "dominant" one by one with the univariate LMM
##' X.D <- recodeIntoDominant(X=X)
##' fit.u.DD <- gemma(model="ulmm", y=modelD$Y[,1], X=X.D,
##'                   snp.coords, alleles,
##'                   K.c=A.z, W=modelD$W, out.dir=tempdir(), clean="all")
##' cols <- rep("black",ncol(X)); cols[modelD$gamma==1] <- "red"
##' pvadj.DD <- qqplotPval(fit.u.DD$tests$p_wald, col=cols, ctl.fdr.bh=TRUE,
##'                        plot.signif=TRUE)
##' t(binaryClassif(known.nulls=modelD$gamma == 0,
##'                 called.nulls=pvadj.DD$pv.bh > 0.05))
##' }
##' @export
gemma <- function(model="ulmm", y, X, snp.coords, recode.genos=TRUE,
                  alleles=NULL, maf=0.01,
                  K.c=NULL, W=NULL, weights=NULL,
                  exe=Sys.which("gemma"),
                  out.dir=getwd(), task.id="gemma", verbose=1, clean="none",
                  seed=1859, burnin=1000, nb.iters=7000, thin=10){
  stopifnot(file.exists(exe))
  stopifnot(system(paste0(exe, " > /dev/null")) == 0)
  stopifnot(model %in% c("ulmm", "bslmm"))
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
  stopIfNotValidGenosDose(X)
  stopifnot(is.vector(y),
            ! is.null(names(y)),
            length(y) == nrow(X),
            all(names(y) == rownames(X)),
            .isValidSnpCoords(snp.coords),
            all(colnames(X) %in% rownames(snp.coords)),
            is.numeric(maf),
            length(maf) == 1,
            maf >= 0,
            maf <= 1)
  if(! is.null(K.c)){
    stopifnot(is.matrix(K.c))
  }
  if(! is.null(W))
    stopifnot(is.matrix(W),
              all(W[,1] == 1),
              nrow(W) == length(y))
  stopifnot(file.exists(out.dir),
            is.character(task.id),
            length(task.id) == 1,
            clean %in% c("none", "some", "all"))
  if(! is.null(alleles))
    stopifnot(is.data.frame(alleles),
              ! is.null(rownames(alleles)),
              ncol(alleles) == 2,
              all(colnames(X) %in% rownames(alleles)))
  if(! is.null(weights))
    stopifnot(is.vector(weights),
              is.numeric(weights),
              ! is.null(names(weights)),
              length(weights) == length(y),
              all(names(weights) == names(y)),
              all(weights[! is.na(weights)] >= 0))

  out <- list()

  ## generate or reformat alleles
  if(is.null(alleles)){
    alleles <- data.frame(first=rep("A", ncol(X)),
                          second=rep("B", ncol(X)),
                          row.names=colnames(X),
                          stringsAsFactors=FALSE)
  } else
    alleles <- alleles[colnames(X),]

  ## for GEMMA, recode SNP genotypes in terms of minor alleles, if necessary
  afs <- estimSnpAf(X=X)
  recoded <- stats::setNames(rep(FALSE, ncol(X)), colnames(X))
  if(recode.genos){
    if(any(afs > 0.5)){
      msg <- paste0(sum(afs > 0.5), " SNPs in X are not encoded",
                    " in terms of their minor allele")
      write(msg, stdout())
      out.recode <- recodeGenosMinorSnpAllele(X=X, alleles=alleles,
                                              verbose=verbose)
      X <- out.recode$X
      alleles <- out.recode$alleles
      recoded <- out.recode$recoded
    }
  }

  ## discard SNPs with low MAF
  mafs <- estimSnpMaf(afs=afs)
  X <- discardSnpsLowMaf(X=X, mafs=mafs, thresh=maf, verbose=verbose)
  alleles <- alleles[colnames(X),]
  recoded <- recoded[colnames(X)]
  out[["recoded"]] <- recoded

  ## reformat snp.coords
  snp.coords <- snp.coords[colnames(X),]
  if(! "coord" %in% colnames(snp.coords))
    colnames(snp.coords)[colnames(snp.coords) == "pos"] <- "coord"
  sc.gemma <- data.frame(snp=rownames(snp.coords),
                         coord=snp.coords$coord,
                         chr=snp.coords$chr,
                         stringsAsFactors=FALSE)

  ## make design matrix, if necessary
  if(is.null(W))
    W <- matrix(data=1, nrow=nrow(X), ncol=1)

  ## prepare input files
  tmp.files <- c()
  tmp.files <- c(tmp.files,
                 genos=paste0(out.dir, "/genos_", task.id, ".txt.gz"))
  genoDoses2bimbam(X=X, alleles=alleles,
                   file=gzfile(tmp.files["genos"]))
  tmp.files <- c(tmp.files,
                 snp.coords=paste0(out.dir, "/snp-coords_", task.id, ".txt"))
  utils::write.table(x=sc.gemma,
                     file=tmp.files["snp.coords"],
                     quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  if(is.null(K.c)){
    K.c <- estimGenRel(X=X, thresh=maf, relationships="additive",
                       method="zhou", verbose=verbose)
  } else
    K.c <- K.c[rownames(X), rownames(X)]
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
  if(! is.null(weights)){
    tmp.files <- c(tmp.files,
                   weights=paste0(out.dir, "/weights_", task.id, ".txt.gz"))
    utils::write.table(x=weights,
                       file=gzfile(tmp.files["weights"]),
                       quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  }

  ## prepare cmd-line and execute it
  cmd <- paste0("cd ", out.dir,
                "; ", exe,
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
  if(! is.null(weights))
    cmd <- paste0(cmd,
                  " -widv ", tmp.files["weights"])
  if(verbose <= 0){
    tmp.files <- c(tmp.files,
                   stdoe=paste0(out.dir, "/stdouterr_gemma_", task.id, ".txt"))
    cmd <- paste0(cmd, " > ", tmp.files["stdoe"], " 2>&1")
  }
  if(verbose > 1)
    write(cmd, stdout())
  system(cmd)
  out[["cmd"]] <- cmd

  ## load output files
  f <- paste0(out.dir, "/results_", task.id, ".log.txt")
  out[["log"]] <- readLines(f)
  if(model == "ulmm"){
    idx <- grep("beta estimate in the null model", out[["log"]])
    if(length(idx) == 1){
      tmp1 <- strsplit(out$log[idx], " ")[[1]]
      tmp2 <- strsplit(out$log[idx+1], " ")[[1]]
      out[["global.mean"]] <- c(beta.hat=as.numeric(tmp1[length(tmp1)]),
                                se.beta.hat=as.numeric(tmp2[length(tmp2)]))
    }
  } else if(model == "bslmm"){
    idx <- grep("estimated mean", out[["log"]])
    if(length(idx) == 1){
      tmp1 <- strsplit(out$log[idx], " ")[[1]]
      out[["global.mean"]] <- c(beta.hat=as.numeric(tmp1[length(tmp1)]),
                                se.beta.hat=NA)
    }
  }
  idx <- grep("pve estimate in the null model", out[["log"]])
  if(length(idx) == 1){
    tmp1 <- strsplit(out$log[idx], " ")[[1]]
    tmp2 <- strsplit(out$log[idx+1], " ")[[1]]
    out[["pve"]] <- c(pve.hat=as.numeric(tmp1[length(tmp1)]),
                      se.pve.hat=as.numeric(tmp2[length(tmp2)]))
  }
  if(clean == "all")
    file.remove(f)
  if(model == "ulmm"){
    f <- paste0(out.dir, "/results_", task.id, ".assoc.txt")
    tmp <- utils::read.table(file=f, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    rownames(tmp) <- tmp$rs
    out[["tests"]] <- tmp
    if(clean == "all")
      file.remove(f)
  } else if(model == "bslmm"){
    f <- paste0(out.dir, "/results_", task.id, ".hyp.txt")
    out[["hyperparams"]] <-
      utils::read.table(file=f, header=FALSE, sep="\t", skip=1,
                        stringsAsFactors=FALSE,
                        col.names=c("h", "pve", "rho", "pge",
                                    "pi", "n_gamma", ""))
    out[["hyperparams"]][7] <- NULL
    if(clean == "all")
      file.remove(f)
    f <- paste0(out.dir, "/results_", task.id, ".param.txt")
    out[["params"]] <-
      utils::read.table(file=f, header=FALSE, sep="\t", skip=1,
                        stringsAsFactors=FALSE,
                        col.names=c("chr", "rs", "ps", "n_miss",
                                    "alpha", "beta", "gamma"))
    rownames(out[["params"]]) <- out[["params"]]$rs
    if(clean == "all")
      file.remove(f)
  }

  ## clean temporary files
  if(clean != "none")
    for(f in tmp.files)
      file.remove(f)

  return(out)
}

##' GEMMA uLMM per chromosome
##'
##' See Zhou & Stephens (Nature Genetics, 2012).
##' @param y vector of phenotypes with genotype names
##' @param X matrix of bi-allelic SNP genotypes encoded, for each SNP, in number of copies of its second allele, i.e. as allele doses in {0,1,2}, with genotypes in rows and SNPs in columns; the "second" allele is arbitrary, it corresponds to the second column of \code{alleles}, which can be the minor or the major allele
##' @param snp.coords data.frame with SNPs as row names and two columns named "coord" and "chr"
##' @param recode.genos if TRUE, SNP genotypes in X will be recoded so that the minor allele is counted
##' @param alleles data.frame with SNPs in rows (names as row names) and alleles in columns (exactly 2 columns are required); the second column should correspond to the allele which number of copies is counted at each SNP in \code{X}; if NULL, fake alleles will be generated
##' @param maf SNPs with minor allele frequency strictly below this threshold will be discarded
##' @param chr.ids set of chromosome identifiers to analyze (optional, all by default)
##' @param W matrix of covariates with genotypes in rows (names as row names), a first column of 1 and a second column of covariates values
##' @param weights vector of positive weights with genotype names
##' @param exe path to the executable of GEMMA
##' @param out.dir directory in which the data files will be found
##' @param task.id identifier of the task (used in output file names)
##' @param clean remove files: none, some (temporary only), all (temporary and results)
##' @param verbose verbosity level (0/1)
##' @return a data.frame with GEMMA's output for all chromosomes
##' @author Timothee Flutre [aut,cre], Dalel Ahmed [ctb]
##' @seealso \code{link{gemma}}, \code{\link{plotHistPval}}, \code{\link{qqplotPval}}
##' @export
gemmaUlmmPerChr <- function(y, X, snp.coords, recode.genos=TRUE,
                            alleles=NULL, maf=0.01,
                            chr.ids=NULL, W, weights=NULL,
                            exe=Sys.which("gemma"),
                            out.dir=getwd(), task.id="gemma",
                            clean="none", verbose=1){
  stopifnot(file.exists(exe))
  stopifnot(system(paste0(exe, " > /dev/null")) == 0)
  if(is.matrix(y)){
    stopifnot(ncol(y) == 1,
              ! is.null(rownames(y)))
    y <- stats::setNames(y[,1], rownames(y))
  }
  stopIfNotValidGenosDose(X)
  stopifnot(is.vector(y),
            ! is.null(names(y)),
            length(y) == nrow(X),
            all(names(y) == rownames(X)),
            .isValidSnpCoords(snp.coords),
            all(colnames(X) %in% rownames(snp.coords)),
            is.numeric(maf),
            length(maf) == 1,
            maf >= 0,
            maf <= 1,
            is.matrix(W),
            all(W[,1] == 1),
            nrow(W) == length(y),
            file.exists(out.dir),
            is.character(task.id),
            length(task.id) == 1,
            clean %in% c("none", "some", "all"))
  if(! is.null(alleles))
    stopifnot(is.data.frame(alleles),
              ! is.null(rownames(alleles)),
              ncol(alleles) == 2,
              all(colnames(X) %in% rownames(alleles)))
  if(! is.null(weights))
    stopifnot(is.vector(weights),
              is.numeric(weights),
              ! is.null(names(weights)),
              length(weights) == length(y),
              all(names(weights) == names(y)),
              all(weights[! is.na(weights)] >= 0))

  out <- list()

  ## generate or reformat alleles
  if(is.null(alleles)){
    alleles <- data.frame(first=rep("A", ncol(X)),
                          second=rep("B", ncol(X)),
                          row.names=colnames(X),
                          stringsAsFactors=FALSE)
  } else
    alleles <- alleles[colnames(X),]

  ## for GEMMA, recode SNP genotypes in terms of minor alleles, if necessary
  afs <- estimSnpAf(X=X)
  recoded <- stats::setNames(rep(FALSE, ncol(X)), colnames(X))
  if(recode.genos){
    if(any(afs > 0.5)){
      msg <- paste0(sum(afs > 0.5), " SNPs in X are not encoded",
                    " in terms of their minor allele")
      write(msg, stdout())
      out.recode <- recodeGenosMinorSnpAllele(X=X, alleles=alleles,
                                              verbose=verbose)
      X <- out.recode$X
      alleles <- out.recode$alleles
      recoded <- out.recode$recoded
    }
  }

  ## discard SNPs with low MAF
  mafs <- estimSnpMaf(afs=afs)
  X <- discardSnpsLowMaf(X=X, mafs=mafs, thresh=maf, verbose=verbose)
  alleles <- alleles[colnames(X),]
  recoded <- recoded[colnames(X)]
  out[["recoded"]] <- recoded

  ## reformat snp.coords
  snp.coords <- snp.coords[colnames(X),]
  if(! "coord" %in% colnames(snp.coords))
    colnames(snp.coords)[colnames(snp.coords) == "pos"] <- "coord"

  if(is.null(chr.ids))
    chr.ids <- sort(unique(as.character(snp.coords$chr)))

  ## launch GEMMA for each chromosome
  out.pc <- list()
  for(chr.id in chr.ids){
    subset.snp.ids <- (snp.coords$chr == chr.id)
    if(verbose > 0)
      write(paste0(chr.id, ": ", sum(subset.snp.ids), " SNPs"), stdout())
    K.c <- estimGenRel(X=X[, ! subset.snp.ids], thresh=maf,
                       relationships="additive", method="zhou",
                       verbose=verbose)
    out.chr.id <- gemma(model="ulmm",
                        y=y,
                        X=X[, subset.snp.ids],
                        snp.coords=snp.coords[subset.snp.ids,],
                        recode.genos=recode.genos,
                        alleles=alleles[subset.snp.ids,],
                        maf=maf, K.c=K.c, W=W, weights=weights,
                        exe=exe,
                        out.dir=out.dir,
                        task.id=paste0(task.id, "-", chr.id),
                        verbose=verbose-1, clean=clean)
    out.pc[[chr.id]] <- out.chr.id
  }

  ## combine the outputs of each chromosome
  out[["cmd"]] <- sapply(out.pc, `[[`, "cmd")
  out[["log"]] <- lapply(out.pc, `[[`, "log")
  out[["global.mean"]] <- t(sapply(out.pc, `[[`, "global.mean"))
  out[["pve"]] <- t(sapply(out.pc, `[[`, "pve"))
  out[["tests"]] <- do.call(rbind, lapply(out.pc, `[[`, "tests"))
  rownames(out$tests) <- out$tests$rs
  out$tests <- out$tests[colnames(X),]

  return(out)
}
