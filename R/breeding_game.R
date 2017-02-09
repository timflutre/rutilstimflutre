## Contains functions useful for the "breeding game"

##' Set up breeding game
##'
##' Set up the directories and names for the breeding game.
##' Already-existing directories are not re-created.
##' @param root.dir path to the root directory
##' @param shared.dir path to the shared directory (e.g. via Dropbox; root.dir if NULL)
##' @param nb.breeders number of breeders (an additional "test" breeder will be created)
##' @param lang language to be used (fr/en)
##' @return list
##' @author Timothee Flutre
##' @export
setUpBreedingGame <- function(root.dir, shared.dir=NULL, nb.breeders=3,
                              lang="fr"){
  stopifnot(is.character(root.dir),
            dir.exists(root.dir),
            lang %in% c("fr", "en"))
  if(! is.null(shared.dir))
    stopifnot(is.character(shared.dir),
              dir.exists(shared.dir))

  out <- list(root.dir=root.dir)

  if(lang == "fr"){
    truth.dir <- paste0(root.dir, "/", "verite")
  } else if(lang == "en"){
    truth.dir <- paste0(root.dir, "/", "truth")
  }
  if(! dir.exists(truth.dir))
    dir.create(truth.dir)
  out$truth.dir <- truth.dir

  if(is.null(shared.dir)){
    if(lang == "fr"){
      shared.dir <- paste0(root.dir, "/", "partage")
    } else if(lang == "en"){
      shared.dir <- paste0(root.dir, "/", "shared")
    }
  }
  if(! dir.exists(shared.dir))
    dir.create(shared.dir)
  out$shared.dir <- shared.dir

  if(lang == "fr"){
    init.dir <- paste0(shared.dir, "/", "donnees_initiales")
  } else if(lang == "en"){
    init.dir <- paste0(shared.dir, "/", "initial_data")
  }
  if(! dir.exists(init.dir))
    dir.create(init.dir)
  out$init.dir <- init.dir

  breeders <- c("test")
  if(lang == "fr"){
    breeders <- c(breeders, paste0("selectionneur", 1:nb.breeders))
  } else if(lang == "en"){
    breeders <- c(breeders, paste0("breeder", 1:nb.breeders))
  }
  breeder.dirs <- c()
  for(breeder in breeders){
    truth.breeder.dir <- paste0(truth.dir, "/", breeder)
    if(! dir.exists(truth.breeder.dir))
      dir.create(truth.breeder.dir)
    breeder.dirs[[breeder]] <- paste0(shared.dir, "/", breeder)
    if(! dir.exists(breeder.dirs[breeder]))
      dir.create(breeder.dirs[breeder])
  }
  out$breeders <- breeders
  out$breeder.dirs <- breeder.dirs

  dbname <- paste0(root.dir, "/breeding-game.sqlite")
  out$dbname <- dbname

  return(out)
}

##' Check breeding game
##'
##' Check if the given breeder is part of the game.
##' @param breeder name of the breeder (e.g. "breeder3" or "selectionneur2", depending on the language)
##' @param root.dir path to the root directory
##' @param lang language to be used (fr/en)
##' @return logical
##' @author Timothee Flutre
##' @export
doesBreederExist <- function(breeder, root.dir, lang="fr"){
  stopifnot(is.character(breeder),
            is.character(root.dir),
            dir.exists(root.dir),
            is.character(lang),
            lang %in% c("fr", "en"))

  path.to.dir <- root.dir
  if(lang == "fr"){
    path.to.dir <- paste0(path.to.dir, "/partage")
  } else if(lang == "en"){
    path.to.dir <- paste0(path.to.dir, "/shared")
  }
  path.to.dir <- paste0(path.to.dir, "/", breeder)

  return(dir.exists(path.to.dir))
}

##' Set up breeding game
##'
##' Retrieve the paths to the directories used for the breeding game.
##' @param root.dir path to the root directory
##' @return list
##' @author Timothee Flutre
##' @export
getBreedingGameSetup <- function(root.dir){
  stopifnot(is.character(root.dir),
            length(root.dir) == 1,
            dir.exists(root.dir))

  out <- list(root.dir=root.dir)

  if(dir.exists(paste0(root.dir, "/verite"))){
    lang <- "fr"
    out$truth.dir <- paste0(root.dir, "/verite")
    out$shared.dir <- paste0(root.dir, "/partage")
    out$init.dir <- paste0(out$shared.dir, "/donnees_initiales")
    tmp <- length(Sys.glob(paste0(out$shared.dir, "/selectionneur*")))
    out$breeders <- c("test", paste0("selectionneur", 1:tmp))

  } else if(dir.exists(paste0(root.dir, "/truth"))){
    lang <- "en"
    out$truth.dir <- paste0(root.dir, "/truth")
    out$shared.dir <- paste0(root.dir, "/shared")
    out$init.dir <- paste0(out$shared.dir, "/initial_data")
    tmp <- length(Sys.glob(paste0(out$shared.dir, "/breeder*")))
    out$breeders <- c("test", paste0("breeder", 1:tmp))

  } else
    stop("can't determine the langage used for the breeding game")

  out$breeder.dirs <- c()
  for(breeder in out$breeders)
    out$breeder.dirs[[breeder]] <- paste0(out$shared.dir, "/", breeder)

  out$dbname <- paste0(root.dir, "/breeding-game.sqlite")

  return(out)
}

##' Simul breeding game
##'
##' Simulate the random part of traits 1 and 2 jointly.
##' @param h2 narrow-sense heritability of both traits
##' @param sigma.beta2 variance of the SNP effects for both traits
##' @param prop.pleio proportion of SNPs having an effect on both traits
##' @param cor.pleio genotypic correlation at these SNPs
##' @param X matrix of SNP genotypes encoded in allele dose in {0,1,2}; will be encoded in {-1,0,1} afterwards
##' @param center.genos if TRUE, the genotypes will be centered per SNP (but for the breeding game, it should be set to FALSE as the SNP effects won't change from generation to generation)
##' @param plot if TRUE, the first-trait genotypic values will be plotted against the second-trait
##' @param verbose verbosity level (0/1)
##' @return matrix of genotypic values with noise
##' @author Timothee Flutre
##' @export
simulTraits12Rnd <- function(h2=c(0.3, 0.4), sigma.beta2=c(10^(-5), 10^(-5)),
                             prop.pleio=0.4, cor.pleio=-0.7,
                             X, center.genos=FALSE,
                             plot=TRUE, verbose=1){
  stopifnot(length(h2) == 2,
            all(h2 >= 0, h2 <= 1),
            all(sigma.beta2 >= 0),
            all(prop.pleio >= 0, prop.pleio <= 1),
            all(cor.pleio >= -1, cor.pleio <= 1),
            is.logical(center.genos),
            is.logical(plot))
  stopIfNotValidGenosDose(X)

  I <- nrow(X)
  P <- ncol(X)

  Sigma.beta.nopleio <- matrix(c(sigma.beta2[1], 0, 0, sigma.beta2[2]),
                               nrow=2, ncol=2)
  Beta <- MASS::mvrnorm(n=P, mu=c(0,0), Sigma=Sigma.beta.nopleio)
  if(verbose > 0)
    message(paste0("cor(Beta[,1], Beta[,2]) = ",
                   format(stats::cor(Beta[,1], Beta[,2]), digits=3)))

  cov.pleio <- cor.pleio * sqrt(sigma.beta2[1] * sigma.beta2[2])
  Sigma.beta.pleio <- matrix(c(sigma.beta2[1], cov.pleio, cov.pleio,
                               sigma.beta2[2]), nrow=2, ncol=2)
  length(idx.pleio <- sample.int(n=P, size=floor(prop.pleio * P)))
  Beta[idx.pleio,] <- MASS::mvrnorm(n=length(idx.pleio), mu=c(0,0),
                                    Sigma=Sigma.beta.pleio)

  if(verbose > 0){
    message(paste0("cor(Beta[idx.pleio,1], Beta[idx.pleio,2]) = ",
                   format(stats::cor(Beta[idx.pleio,1], Beta[idx.pleio,2]),
                          digits=3)))
    message(paste0("cor(Beta[,1], Beta[,2]) = ",
                   format(stats::cor(Beta[,1], Beta[,2]), digits=3)))
  }

  X.tmp <- scale(x=X - 1, center=center.genos, scale=FALSE)
  G.A <- X.tmp %*% Beta
  if(verbose > 0){
    message(paste0("cor(G.A[,1], G.A[,2]) = ",
                   format(stats::cor(G.A[,1], G.A[,2]), digits=3)))
    message(paste0("mean(G.A[,1]) = ", format(mean(G.A[,1]), digits=3)))
    message(paste0("mean(G.A[,2]) = ", format(mean(G.A[,2]), digits=3)))
    message(paste0("var(G.A[,1]) = ", format(stats::var(G.A[,1]), digits=3)))
    message(paste0("var(G.A[,2]) = ", format(stats::var(G.A[,2]), digits=3)))
  }
  if(plot){
    regplot(G.A[,1], G.A[,2], xlab="G.A[,1]", ylab="G.A[,2]")
    graphics::abline(h=0, v=0, lty=2)
  }

  sigma2 <- c(((1 - h2[1]) / h2[1]) * stats::var(G.A[,1]),
              ((1 - h2[2]) / h2[2]) * stats::var(G.A[,2]))
  if(verbose > 0){
    message(paste0("sigma2[1] = ", format(sigma2[1], digits=3)))
    message(paste0("sigma2[2] = ", format(sigma2[2], digits=3)))
  }
  Sigma <- matrix(c(sigma2[1], 0, 0, sigma2[2]), nrow=2, ncol=2)
  Epsilon <- MASS::mvrnorm(n=I, mu=c(0,0), Sigma=Sigma)

  Y <- G.A + Epsilon
  if(verbose > 0){
    message(paste0("var(G.A[,1]) / var(Y[,1]) = ",
                   format(stats::var(G.A[,1]) / stats::var(Y[,1]), digits=3)))
    message(paste0("var(G.A[,2]) / var(Y[,2]) = ",
                   format(stats::var(G.A[,2]) / stats::var(Y[,2]), digits=3)))
  }

  return(list(h2=h2, sigma.beta2=sigma.beta2, prop.pleio=prop.pleio,
              cor.pleio=cor.pleio, X=X, Beta=Beta, G.A=G.A,
              sigma2=sigma2, Y=Y))
}

##' Example for breeding game
##'
##' Make a file with examples of plant material to request (autofecundation, allofecundation, haplodiploidization)
##' @param out.dir path to the directory in which the file will be saved
##' @param lang language to be used (fr/en)
##' @return invisible data.frame
##' @author Timothee Flutre
##' @seealso \code{\link{readCheckBreedPlantFile}}
##' @export
makeExamplePlantFile <- function(out.dir, lang="fr"){
  stopifnot(dir.exists(out.dir),
            lang %in% c("fr", "en"))

  plants <-
    data.frame(parent1=c("Coll0037", "37-8.1", "37-8.1"),
               parent2=c("Coll0008", "37-8.1", NA),
               child=c("37-8.1", "37-8.1.1", "37-8.1.HD1"),
               explanations=c("allofecundation", "autofecundation",
                              "haplodiploidization"),
               stringsAsFactors=FALSE)

  if(lang == "fr"){
    f <- paste0(out.dir, "/exemple_requete_croisements.txt")
  } else if(lang == "en")
    f <- paste0(out.dir, "/example_request_crosses.txt")

  cat("# the table (below) contains the crosses to be made\n",
      file=f, append=FALSE)
  cat("# only columns 'parent1', 'parent2' and 'child' are compulsory\n",
      file=f, append=TRUE)
  cat("# each row corresponds to a child\n",
      file=f, append=TRUE)
  cat("# the 'child' column shouldn't have any duplicate\n",
      file=f, append=TRUE)
  cat("# only the 'parent2' column can be empty\n",
      file=f, append=TRUE)
  cat("# individual names should only use [a-z], [A-Z], [0-9], [_-] (no space, comma, etc)\n",
      file=f, append=TRUE)
  cat("# use write.table(x=crosses, file=\"<breeder>_crosses.txt\", quote=FALSE, sep=\"\\t\", na=\"\", row.names=FALSE, col.names=TRUE)\n",
      file=f, append=TRUE)
  cat("# lines starting with '#' will be ignored\n",
      file=f, append=TRUE)
  suppressWarnings(utils::write.table(x=plants, file=f, append=TRUE, quote=FALSE,
                               sep="\t", na="", row.names=FALSE,
                               col.names=TRUE))

  invisible(plants)
}

##' Example for breeding game
##'
##' Make a file with examples of data to request (phenotypes, genotypes).
##' @param out.dir path to the directory in which the file will be saved
##' @param lang language to be used (fr/en)
##' @return invisible data.frame
##' @author Timothee Flutre
##' @seealso \code{\link{readCheckBreedDataFile}}
##' @export
makeExampleDataFile <- function(out.dir, lang="fr"){
  stopifnot(dir.exists(out.dir),
            lang %in% c("fr", "en"))

  dat <- data.frame(ind=c(paste0("ind", 7:10), "ind7", "ind31"),
                    task=c(rep("pheno", 3), "geno", "geno", "geno"),
                    details=c("2", "1", "5", "hd", "snp15492", "ld"),
                    stringsAsFactors=FALSE)

  if(lang == "fr"){
    f <- paste0(out.dir, "/exemple_requete_donnees.txt")
  } else if(lang == "en")
    f <- paste0(out.dir, "/example_request_data.txt")

  cat("# the list (below) contains individuals to genotype and/or phenotype\n",
      file=f, append=FALSE)
  cat("# columns 'ind', 'task' and 'details' are compulsory\n",
      file=f, append=TRUE)
  cat("# individual names should only use [a-z], [A-Z], [0-9], [_-] (no space, comma, etc)\n",
      file=f, append=TRUE)
  cat("# the 'task' column should contain 'pheno' or 'geno'\n",
      file=f, append=TRUE)
  cat("# if 'task=pheno', the 'details' column should contain the number of plots\n",
      file=f, append=TRUE)
  cat("# if 'task=geno', the 'details' column should contain 'hd', 'ld' or the SNP identifier\n",
      file=f, append=TRUE)
  cat("# individuals should not be duplicated within each task\n",
      file=f, append=TRUE)
  cat("# the total number of requested plots (task=pheno) should not exceed the total available\n",
      file=f, append=TRUE)
  cat("# use write.table(x=inds, file=\"<breeder>_inds.txt\", quote=FALSE, sep=\"\\t\", row.names=FALSE, col.names=TRUE)\n",
      file=f, append=TRUE)
  cat("# lines starting with '#' will be ignored\n",
      file=f, append=TRUE)
  suppressWarnings(utils::write.table(x=dat, file=f, append=TRUE, quote=FALSE,
                               sep="\t", row.names=FALSE, col.names=TRUE))

  invisible(dat)
}

##' Read for breeding game
##'
##' Read and check a file supposed to contain requests about plant material.
##' It should have 3 columns named \code{parent1}, \code{parent2} and \code{child}.
##' @param f path to the input file (columns should be separated by a tabulation)
##' @param df data.frame (if the file was already read)
##' @param max.nb.hd maximum number of haplodiploidization to request in a single file
##' @return invisible data.frame
##' @author Timothee Flutre
##' @seealso \code{\link{makeExamplePlantFile}}
##' @export
readCheckBreedPlantFile <- function(f=NULL, df=NULL, max.nb.hd=1000){
  stopifnot(! is.null(f) || ! is.null(df),
            is.numeric(max.nb.hd),
            length(max.nb.hd) == 1,
            max.nb.hd > 0)

  if(is.null(df)){
    stopifnot(file.exists(f))
    df <- utils::read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  }

  stopifnot(is.data.frame(df),
            ncol(df) >= 3,
            all(c("parent1", "parent2", "child") %in% colnames(df)),
            all(! is.na(df$parent1)),
            all(! is.na(df$child)),
            anyDuplicated(df$child) == 0,
            all(! grepl("[^[:alnum:]._-]", df$child)),
            sum(is.na(df$parent2)) <= max.nb.hd)

  invisible(df)
}

##' Read for breeding game
##'
##' Read and check a file supposed to contain requests about pheno/geno data.
##' It should have 3 columns named \code{ind}, \code{task} and \code{details}.
##' @param f path to the input file (columns should be separated by a tabulation)
##' @param df data.frame (if the file was already read)
##' @param max.nb.plots maximum number of plots
##' @param subset.snps list with two components named "ld" and "hd" containing vector of SNP identifiers
##' @param max.nb.inds maximum number of unique individuals for which at least one request is made
##' @return invisible data.frame
##' @author Timothee Flutre
##' @seealso \code{\link{makeExampleDataFile}}
##' @export
readCheckBreedDataFile <- function(f=NULL, df=NULL, max.nb.plots, subset.snps,
                                   max.nb.inds=1000){
  stopifnot(! is.null(f) || ! is.null(df),
            is.numeric(max.nb.plots),
            length(max.nb.plots) == 1,
            max.nb.plots > 0,
            is.list(subset.snps),
            all(names(subset.snps) %in% c("ld", "hd")))

  if(is.null(df)){
    stopifnot(file.exists(f))
    df <- utils::read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  }

  stopifnot(is.data.frame(df),
            ncol(df) >= 3,
            all(c("ind", "task", "details") %in% colnames(df)),
            all(! is.na(df$ind)),
            length(unique(df$ind)) <= max.nb.inds,
            all(! grepl("[^[:alnum:]._-]", df$ind)),
            all(! is.na(df$task)),
            all(df$task %in% c("pheno", "geno")))

  if("pheno" %in% df$task){
    tmp <- suppressWarnings(as.numeric(df$details[df$task ==
                                                  "pheno"]))
    stopifnot(all(! is.na(tmp)),
              ! anyDuplicated(df$ind[df$task == "pheno"]),
              sum(as.numeric(df$details[df$task == "pheno"])) <=
              max.nb.plots)
  }

  if("geno" %in% df$task){
    stopifnot(all(grepl("hd|ld|snp", df$details[df$task == "geno"])))
    idx.notsnp <- df$task == "geno" & ! grepl("snp", df$details)
    stopifnot(! anyDuplicated(df$ind[idx.notsnp]),
              all(df$details[df$task == "geno" &
                             ! df$details %in% c("ld","hd")]
                  %in% subset.snps[["hd"]]))
  }

  invisible(df)
}

##' Count types
##'
##' Count the types of breeding requests.
##' @param df data.frame from readCheckBreedPlantFile or readCheckBreedDataFile
##' @return named vector
##' @author Timothee Flutre
##' @export
countRequestedBreedTypes <- function(df){
  stopifnot(is.data.frame(df),
            ! is.null(colnames(df)))

  types <- NULL

  if("parent1" %in% colnames(df)){
    types <- stats::setNames(c(sum(df$parent1 != df$parent2 &
                            ! is.na(df$parent2)),
                        sum(df$parent1 == df$parent2 &
                            ! is.na(df$parent2)),
                        sum(is.na(df$parent2))),
                      c("allofecundation", "autofecundation",
                        "haplodiploidization"))
  } else if("ind" %in% colnames(df)){
    types <- stats::setNames(c(ifelse("pheno" %in% df$task,
                               sum(as.numeric(df$details[df$task ==
                                                                "pheno"])),
                               0),
                        sum(df$task == "geno" &
                            df$details == "hd"),
                        sum(df$task == "geno" &
                            df$details == "ld"),
                        sum(df$task == "geno" &
                            ! df$details %in% c("hd", "ld"))),
                      c("pheno", "geno-hd", "geno-ld", "geno-single-snp"))
  }

  return(types)
}
