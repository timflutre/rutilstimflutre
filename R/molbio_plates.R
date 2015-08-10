## Contains functions handling plates as often used in molecular biology.

##' Initialize plate(s)
##'
##' Initialize plates as matrices with missing data.
##' @param n number of plates
##' @param nrow vector of number of rows for each plate
##' @param ncol vector of number of columns for each plate
##' @param names vector of names for each plate
##' @return list of matrices, one per plate, in the "wide" format
##' @author Timothee Flutre
init.plates <- function(n, nrow, ncol, names){
  plates <- list()

  for(i in 1:n)
    plates[[names[i]]] <-
      matrix(data=NA, nrow=nrow[i],
             ncol=ncol[i],
             dimnames=list(LETTERS[1:nrow[i]], 1:ncol[i]))

  return(plates)
}

##' Describe plate
##'
##' Describe the given plate.
##' @param plate matrix
##' @param plate.name string
##' @param verbose verbosity level (0/default=1)
##' @return invisible named vector
##' @author Timothee Flutre
desc.plate <- function(plate, plate.name, verbose=1){
  out <- setNames(object=rep(NA, 3),
                  nm=c("nb.wells", "nb.empty.wells", "nb.samples"))

  out["nb.wells"] <- nrow(plate) * ncol(plate)
  out["nb.empty.wells"] <- sum(is.na(plate))
  samples <- c(plate)
  samples <- unique(samples[! is.na(samples)])
  out["nb.samples"] <- length(samples)

  if(verbose > 0)
    write(paste0(plate.name,
                 " ", out["nb.wells"],
                 " ", out["nb.empty.wells"],
                 " ", out["nb.samples"]),
          stdout())

  invisible(out)
}

##' Load plate(s)
##'
##' Read each file into a matrix and gather them into a list.
##' @param files vector of paths to file(s), one per plate, in csv (sep=",") or other format (sep="<tab>")
##' @param verbose verbosity level (0/default=1)
##' @return list of matrices, one per plate, in the "wide" format
##' @author Timothee Flutre
load.plates <- function(files, verbose=1){
  plates <- list()

  if(verbose > 0)
    write("plate nb.wells nb.empty.wells nb.samples", stdout())

  for(i in seq_along(files)){
    plate <- NULL
    plate.name <- strsplit(x=basename(files[i]), split="\\.")[[1]][1]
    file.ext <- rev(strsplit(x=basename(files[i]), split="\\.")[[1]])[1]
    if(file.ext == "csv"){
      plate <- read.csv(file=files[i], stringsAsFactors=FALSE,
                        row.names=1)
    } else
      plate <- read.table(file=files[i], header=TRUE, sep="\t",
                          row.names=1, stringsAsFactors=FALSE)
    if(ncol(plate) == 0)
      stop(paste0(plate.name, ": 0 columns"))
    plate <- as.matrix(plate)
    colnames(plate) <- as.character(1:ncol(plate))
    plate[plate == ""] <- NA
    if(verbose > 0)
      desc.plate(plate=plate, plate.name=plate.name, verbose=1)
    plates[[plate.name]] <- plate
  }

  return(plates)
}

##' Plot plate
##'
##' Plot the arrangement of samples on the given plate.
##' @param plate matrix
##' @param main string containing the text for the main title
##' @return nothing
##' @author Timothee Flutre
plotPlate <- function(plate, main="Plate"){
  stopifnot(is.matrix(plate))

  par(mar=c(3, 3, 5, 1) + 0.1)
  plot(x=0, y=0, type="n",
       xlim=c(0.7, ncol(plate)+0.3),
       ylim=c(0.7, nrow(plate)+0.3),
       main="", xlab="", ylab="", xaxt="n", yaxt="n")
  mtext(text=main, side=3, line=3, cex=2, font=2)
  axis(side=3, at=1:ncol(plate), labels=colnames(plate))
  axis(side=2, at=1:nrow(plate), labels=rev(rownames(plate)), las=1)

  text(x=rep(1:ncol(plate), each=nrow(plate)),
       y=rep(nrow(plate):1, ncol(plate)),
       c(plate))
}

##' Lengthen plate
##'
##' Lengthen a "wide" plate into 3 columns for easier processing.
##' @param plate.w matrix of a plate in the "wide" format
##' @return data.frame of a plate in the "long" format (1 well per row)
##' @author Timothee Flutre
lengthen.plate <- function(plate.w){
  stopifnot(is.matrix(plate.w),
            ! is.null(rownames(plate.w)),
            ! is.null(colnames(plate.w)))

  nb.samples <- nrow(plate.w) * ncol(plate.w)
  plate.l <- data.frame(sample=rep(NA, nb.samples),
                        row=rep(NA, nb.samples),
                        col=rep(NA, nb.samples))

  sample.id <- 1
  for(i in 1:nrow(plate.w)){
    for(j in 1:ncol(plate.w)){
      plate.l$sample[sample.id] <- plate.w[i,j]
      plate.l$row[sample.id] <- rownames(plate.w)[i]
      plate.l$col[sample.id] <- colnames(plate.w)[j]
      sample.id <- sample.id + 1
    }
  }

  return(plate.l)
}

##' Find empty wells
##'
##' Identify empty wells, if any, in the given plate.
##' @param plate.w matrix of a plate in the "wide" format
##' @return 2 column data.frame (row;col) corresponding to empty wells
##' @author Timothee Flutre
empty.wells <- function(plate.w){
  plate.l <- lengthen.plate(plate.w)
  empty.idx <- is.na(plate.l$sample)
  return(plate.l[empty.idx, c("row", "col")])
}

##' Randomize plate
##'
##' Randomize the given plate.
##' @param plate matrix
##' @param rand.scheme string
##' @param seed integer
##' @param file path to file in which to write the randomised plate
##' @return matrix
##' @author Timothee Flutre
rand.plate <- function(plate, rand.scheme="1x96", seed=NULL, file=NULL){
  stopifnot(is.matrix(plate),
            ! is.null(dimnames(plate)),
            nrow(plate) * ncol(plate) == 96,
            rand.scheme %in% c("1x96", "2x48"))

  if(! is.null(seed))
    set.seed(seed)

  out <- matrix(data=NA, nrow=nrow(plate), ncol=ncol(plate), byrow=FALSE,
                dimnames=dimnames(plate))

  if(rand.scheme == "1x96"){
    indices <- seq(from=1, to=96, by=1)
    rand.indices <- sample(x=indices, size=length(indices), replace=FALSE)
    out[indices] <- plate[rand.indices]
  } else if(rand.scheme == "2x48"){
    indices1 <- seq(from=1, to=48, by=1)
    rand.indices1 <- sample(x=indices1, size=length(indices1), replace=FALSE)
    out[indices1] <- plate[rand.indices1]
    indices2 <- seq(from=49, to=96, by=1)
    rand.indices2 <- sample(x=indices2, size=length(indices2), replace=FALSE)
    out[indices2] <- plate[rand.indices2]
  }

  if(! is.null(file)){
    write(x=paste0("#seed=", seed), file=file, append=FALSE)
    write.table(x=t(c("", colnames(out))), file=file, append=TRUE,
                quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    write.table(x=out, file=file, append=TRUE, quote=FALSE, sep="\t",
                row.names=TRUE, col.names=FALSE)
  }

  return(out)
}
