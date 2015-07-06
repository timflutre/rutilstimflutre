## Contains functions handling plates as often used in molecular biology.

##' Initialize plates as matrices with missing data.
##'
##'
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

##' Read each file into a matrix and gather them into a list.
##'
##'
##' @param files vector of paths to csv file, one per plate
##' @param verbose verbosity level (0/default=1)
##' @return list of matrices, one per plate, in the "wide" format
##' @author Timothee Flutre
load.plates <- function(files, verbose=1){
  plates <- list()

  for(i in seq_along(files)){
    plate.name <- strsplit(x=basename(files[i]), split="\\.")[[1]][1]
    plate <- read.csv(file=files[i], stringsAsFactors=FALSE,
                      row.names=1)
    plate <- as.matrix(plate)
    colnames(plate) <- as.character(1:ncol(plate))
    plate[plate == ""] <- NA
    if(verbose > 0){
      write(paste0(plate.name, ": ", sum(is.na(plate)),
                   " missing data"),
            stdout())
    }
    plates[[plate.name]] <- plate
  }

  return(plates)
}

##' Plot the arrangement of samples on a given plate.
##'
##'
##' @param plate matrix
##' @param main string containing the text for the main title
##' @return nothing
##' @author Timothee Flutre
plot.plate <- function(plate, main="Plate"){
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

##' Lengthen a "wide" plate into 3 columns for easier processing.
##'
##'
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

##' Identify empty wells, if any, in a plate.
##'
##'
##' @param plate.w matrix of a plate in the "wide" format
##' @return 2 column data.frame (row;col) corresponding to empty wells
##' @author Timothee Flutre
empty.wells <- function(plate.w){
  plate.l <- lengthen.plate(plate.w)
  empty.idx <- is.na(plate.l$sample)
  return(plate.l[empty.idx, c("row", "col")])
}
