## Contains functions useful for population genetics.

##' STRUCTURE
##'
##' Read the output file from STRUCTURE.
##' @param out.file path to the output file
##' @return list
##' @author Timothee Flutre
##' @export
readOutSTRUCTURE <- function(out.file){
  stopifnot(file.exists(out.file))

  out <- list()

  tmp <- readLines(out.file, warn=FALSE)

  idx.start <- grep("^Inferred ancestry of individuals:$", tmp)
  idx.end <- grep("^Estimated Allele Frequencies in each cluster$", tmp)
  ancestries <- tmp[(idx.start+2):(idx.end-3)]
  ancestries <- as.data.frame(t(sapply(strsplit(ancestries, " "),
                                       function(x){x[! x == ""]})),
                              stringsAsFactors=FALSE)
  stopifnot(all(ancestries[,4] == ":"))
  nb.clusters <- ncol(ancestries) - 4
  colnames(ancestries) <- c("index", "label", "perc.miss", "sep",
                            paste0("cluster", 1:nb.clusters))
  rownames(ancestries) <- ancestries$label
  ancestries$index <- as.integer(ancestries$index)
  for(i in 1:nb.clusters){
    col.idx <- ncol(ancestries) - i + 1
    ancestries[, col.idx] <- as.numeric(ancestries[, col.idx])
  }
  out$ancestries <- ancestries

  return(out)
}

##' STRUCTURE to CLUMPP
##'
##' Gather output files from STRUCTURE (several runs at the same number of clusters, K) into a single file for input to CLUMPP.
##' @param struct.files vector of file names
##' @param clumpp.file file name
##' @param verbose verbosity level (0/1/2)
##' @return invisible list of lists, one per STRUCTURE file
##' @author Timothee Flutre
##' @seealso \code{\link{readOutSTRUCTURE}}
##' @export
structure2clumpp <- function(struct.files, clumpp.file, verbose=1){
  stopifnot(all(file.exists(struct.files)),
            ! clumpp.file %in% struct.files)
  if(file.exists(clumpp.file))
    file.remove(clumpp.file)

  out <- list()

  if(verbose > 0){
    msg <- paste0("gather ", length(struct.files), " output files",
                  " from STRUCTURE",
                  " to use with CLUMPP")
    write(msg, stdout())
  }

  for(i in seq_along(struct.files)){
    tmp <- readOutSTRUCTURE(struct.files[i])
    out[[struct.files[i]]] <- tmp
    tmp <- tmp$ancestries
    tmp <- cbind(tmp$index,           # col 1) ignored by CLUMPP
                 tmp$index,           # col 2) integer identifying the individuals
                 tmp$perc.miss,       # col 3) ignored by CLUMPP
                 rep(1, nrow(tmp)),   # col 4) ignored by CLUMPP
                 tmp$sep,             # col 5) ignored by CLUMPP
                 tmp[,5:ncol(tmp)])   # cols 6...) membership coefficients
    utils::write.table(x=tmp, file=clumpp.file, append=TRUE, quote=FALSE,
                       sep=" ", row.names=FALSE, col.names=FALSE)
  }

  invisible(out)
}
