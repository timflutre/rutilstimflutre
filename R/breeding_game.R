## Contains functions useful for the "breeding game"

##' Set up
##'
##' Set up the directories and names for the breeding game.
##' Already-existing directories are not re-created.
##' @param root.dir path to the root directory
##' @param shared.dir path to the shared directory (e.g. via Dropbox; root.dir if NULL)
##' @param nb.breeders number of breeders (in EN, will be named "breeder<i>"; in FR, will be named "selectionneur<i>"; a "test" will also be created)
##' @param lang language to be used (default=fr/en)
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

  truth.dir <- paste0(root.dir, "/", ifelse(lang == "fr", "verite",
                                            "truth"))
  if(! dir.exists(truth.dir))
    dir.create(truth.dir)
  out$truth.dir <- truth.dir

  if(is.null(shared.dir))
    shared.dir <- paste0(root.dir, "/", ifelse(lang == "fr", "partage",
                                               "shared"))
  if(! dir.exists(shared.dir))
    dir.create(shared.dir)
  out$shared.dir <- shared.dir

  init.dir <- paste0(shared.dir, "/", ifelse(lang == "fr", "donnees_initiales",
                                             "initial_data"))
  if(! dir.exists(init.dir))
    dir.create(init.dir)
  out$init.dir <- init.dir

  breeders <- c("test", paste0(ifelse(lang == "fr", "selectionneur",
                                      "breeder"), 1:nb.breeders))
  breeder.dirs <- list()
  for(breeder in breeders){
    truth.breeder.dir <- paste0(truth.dir, "/", breeder)
    if(! dir.exists(truth.breeder.dir))
      dir.create(truth.breeder.dir)
    breeder.dirs[[breeder]] <- paste0(shared.dir, "/", breeder)
    if(! dir.exists(breeder.dirs[[breeder]]))
      dir.create(breeder.dirs[[breeder]])
  }
  out$breeders <- breeders
  out$breeder.dirs <- breeder.dirs

  dbname <- paste0(shared.dir, "/breeding-game.sqlite")
  out$dbname <- dbname

  return(out)
}

##' Check
##'
##' Check if the given breeder is part of the game.
##' @param breeder name of the breeder (e.g. "breeder3" or "selectionneur2", depending on the language)
##' @param root.dir path to the root directory
##' @param lang language to be used (default=fr/en)
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

##' Set up
##'
##' Retrieve the paths to the directories used for the breeding game.
##' @param root.dir path to the root directory
##' @return list
##' @author Timothee Flutre
##' @export
getBreedingDirs <- function(root.dir){
  stopifnot(is.character(root.dir),
            dir.exists(root.dir))

  out <- list(root.dir=root.dir)

  if(dir.exists(paste0(root.dir, "/verite"))){
    lang <- "fr"
    out$truth.dir <- paste0(root.dir, "/verite")
    out$shared.dir <- paste0(root.dir, "/partage")
    out$init.dir <- paste0(out$shared.dir, "/donnees_initiales")
    nb.breeders <- length(Sys.glob(paste0(root.dir, "/selectionneur")))
    out$breeders <- paste0("selectionneur", 1:nb.breeders)
    out$breeder.dirs <- paste0(out$shared.dir, "/selectionneur", 1:nb.breeders)

  } else if(dir.exists(paste0(root.dir, "/truth"))){
    lang <- "en"
    out$truth.dir <- paste0(root.dir, "/truth")
    out$shared.dir <- paste0(root.dir, "/shared")
    out$init.dir <- paste0(out$shared.dir, "/initial_data")
    nb.breeders <- length(Sys.glob(paste0(root.dir, "/breeder")))
    out$breeders <- paste0("breeder", 1:nb.breeders)
    out$breeder.dirs <- paste0(out$shared.dir, "/breeder", 1:nb.breeders)

  } else
    stop("can't determine the langage used for the breeding game")

  return(out)
}

##' Example for breeding game
##'
##' Make a file with examples of plant material to request (autofecundation, allofecundation, haplodiploidization)
##' @param out.dir path to the directory in which the file will be saved
##' @return invisible data.frame
##' @author Timothee Flutre
##' @export
makeExamplePlantFile <- function(out.dir){
  stopifnot(dir.exists(out.dir))

  plants <-
    data.frame(parent1=c("Coll0037", "37-8.1", "37-8.1"),
               parent2=c("Coll0008", "37-8.1", NA),
               child=c("37-8.1", "37-8.1.1", "37-8.1.HD1"),
               explanations=c("allofecundation", "autofecundation",
                              "haplodiploidization"),
               stringsAsFactors=FALSE)

  f <- paste0(out.dir, "/requested_crosses_example.txt")

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
##' @return invisible data.frame
##' @author Timothee Flutre
##' @export
makeExampleDataFile <- function(out.dir){
  stopifnot(dir.exists(out.dir))

  dat <- data.frame(ind=c(paste0("ind", 7:10), "ind7", "ind31"),
                    task=c(rep("pheno", 3), "geno", "geno", "geno"),
                    details=c("2", "1", "5", "hd", "snp15492", "ld"),
                    stringsAsFactors=FALSE)

  f <- paste0(out.dir, "/requested_data_example.txt")

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
##' @param f path to the input file
##' @param df data.frame (if the file was already read)
##' @param max.nb.hd maximum number of haplodiploidization to request in a single file
##' @return invisible data.frame
##' @author Timothee Flutre
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
##' @param f path to the input file
##' @param df data.frame (if the file was already read)
##' @param max.nb.plots maximum number of plots
##' @param subset.snps list with two components named "ld" and "hd" containing vector of SNP identifiers
##' @param max.nb.inds maximum number of unique individuals for which at least one request is made
##' @return invisible data.frame
##' @author Timothee Flutre
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
