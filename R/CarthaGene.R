## Contains functions useful to wrap CarthaGene.

## Thanks to Matthieu Falque (INRA Moulon)

##' Open CarthaGene
##'
##' Open FIFOs to interact with CarthaGene.
##' @param file path to the file containing the dataset (will be loaded if not NULL)
##' @return list
##' @author Matthieu Falque [aut], Timothee Flutre [ctb]
##' @export
openCarthagene <- function(file=NULL){
  exe.name <- "carthagene"
  stopifnot(file.exists(Sys.which(exe.name)),
            capabilities()["fifo"])
  if(! is.null(file))
    stopifnot(file.exists(file))

  out <- list()

  name.fifo.in <- "CGfifoIn"
  name.fifo.out <- "CGfifoOut"
  if(file.exists(name.fifo.in))
    file.remove(name.fifo.in)
  if(file.exists(name.fifo.out))
    file.remove(name.fifo.out)
  out$name.fifo.in <- name.fifo.in
  out$name.fifo.out <- name.fifo.out

  ## input FIFO to send commands to CG
  fifo.in <- fifo(name.fifo.in, "w+", blocking=FALSE)
  out$fifo.in <- fifo.in

  ## output FIFO to get results from CG
  cmd <- paste0("mkfifo ", name.fifo.out)
  system(cmd, intern=FALSE, ignore.stderr=FALSE, wait=TRUE, input=NULL)
  fifo.out <- fifo(name.fifo.out, "r+", blocking=FALSE)
  out$fifo.out <- fifo.out

  ## launch CG in the background
  cmd <- paste0(exe.name, " <", name.fifo.in, " >", name.fifo.out)
  system(cmd, intern=FALSE, ignore.stderr=FALSE, wait=FALSE, input=NULL)

  ## load the dataset
  if(! is.null(file)){
    cmd <- paste0("dsload ",file)
    tmp <- runCarthagene(out, cmd)
    tmp <- readLines(fifo.out, n=-1) # empty the buffer
  }

  return(out)
}

##' Run a CarthaGene command
##'
##' Send a command to CarthaGene and wait for the result.
##' @param cg list output from \code{\link{openCarthagene}}
##' @param cmd character specifying the command to run
##' @return invisible vector of lines corresponding to the output from CarthaGene
##' @author Matthieu Falque [aut], Timothee Flutre [ctb]
##' @export
runCarthagene <- function(cg, cmd){
  stopifnot(is.list(cg),
            length(cg) >= 4,
            all(c("fifo.in", "fifo.out") %in% names(cg)),
            is.character(cmd))

  writeLines(cmd, cg$fifo.in)

  ## dummy command, just to know when CG has finished loading the dataset
  dmy.cmd <- "group -u"
  writeLines(dmy.cmd, cg$fifo.in)

  ## wait for answer
  out <- readLines(cg$fifo.out, n=1)
  idx <- grep("group DistThres LODThres", out)
  while(length(idx) == 0){
    out <- c(out, readLines(cg$fifo.out, n=1))
    idx <- grep("group DistThres LODThres", out)
  }

  ## remove the lines corresponding to the dummy command
  out <- out[-((length(out)-1):length(out))]

  invisible(out)
}

##' Parse CarthaGene's output
##'
##' Parse the output of CarthaGene's command "mrkinfo".
##' @param out.mrkinfo vector of character output from \code{\link{runCarthagene}} corresponding to "mrkinfo"
##' @return data frame with one row per marker and columns "id", "name", etc
##' @author Timothee Flutre
##' @export
parseCgMrkinfo <- function(out.mrkinfo){
  stopifnot(is.character(out.mrkinfo))

  idx <- grep("Num                Names : Sets Merges", out.mrkinfo)
  if(length(idx) == 0){
    msg <- paste0("can't parse CarthaGene's output;",
                  " does it correspond to the 'mrkinfo' command?")
    stop(msg)
  }

  tmp <- strsplit(out.mrkinfo[-c(1:idx)], " ")
  tmp <- lapply(tmp, function(x){
    x[! x %in% c("", ":")]
  })
  if(unique(sapply(tmp, length)) != 3){
    msg <- paste0("can't parse CarthaGene's output;",
                  " does it correspond to the 'mrkinfo' command?")
    stop(msg)
  }
  out <- as.data.frame(do.call(rbind, tmp), stringsAsFactors=FALSE)
  colnames(out) <- c("id", "name", "datasets")
  out$id <- as.integer(out$id)

  return(out)
}

##' Parse CarthaGene's output
##'
##' Parse the output of CarthaGene's command "group".
##' @param out.group vector of character output from \code{\link{runCarthagene}} corresponding to "group"
##' @param mrk.info data frame output from \code{\link{parseCgMrkinfo}} (optional; useful to have the marker names instead of only the marker identifiers)
##' @return list with one data frame per linkage group
##' @author Timothee Flutre
##' @export
parseCgGroup <- function(out.group, mrk.info=NULL){
  stopifnot(is.character(out.group))
  if(! is.null(mrk.info))
    stopifnot(is.data.frame(mrk.info),
              ncol(mrk.info) >= 2,
              all(c("id", "name") %in% colnames(mrk.info)))

  out <- list()

  idx <- grep("Group ID : Marker ID List ..", out.group)
  if(length(idx) == 0){
    msg <- paste0("can't parse CarthaGene's output;",
                  " does it correspond to the 'group' command?")
    stop(msg)
  }

  nb.linkgroups <- as.integer(out.group[length(out.group)])
  for(lg in 1:nb.linkgroups){
    tmp <- strsplit(out.group[idx + lg], ":")[[1]]
    if(length(tmp) != 2){
      msg <- paste0("can't parse data for linkage group '", lg, "'")
      stop(msg)
    }
    out[[lg]] <- data.frame(id=as.integer(strsplit(tmp[2], " ")[[1]][-1]),
                            stringsAsFactors=FALSE)
  }

  if(! is.null(mrk.info)){
    out <- lapply(out, function(lg){
      lg$name <- mrk.info$name[match(lg$id, mrk.info$id)]
      lg
    })
  }

  return(out)
}

##' Parse CarthaGene's output
##'
##' Parse the output of CarthaGene's command "heaprint".
##' @param out.heaprint vector of character output from \code{\link{runCarthagene}} corresponding to "heaprint"
##' @return data frame with one row per map in the heap
##' @author Timothee Flutre
##' @export
parseCgHeaprint <- function(out.heaprint){
  stopifnot(is.character(out.heaprint))

  out <- NULL

  ptn <- "^Map[ -]+[0-9] : log10-likelihood ="
  idx <- grep(pattern=ptn, x=out.heaprint, perl=TRUE)
  if(length(idx) == 0){
    msg <- paste0("can't parse CarthaGene's output;",
                  " does it correspond to the 'heaprint' command?")
    stop(msg)
  }

  lines.log10lik <- out.heaprint[idx]
  out <- data.frame(id=rep(NA, length(lines.log10lik)),
                    log10.lik=NA,
                    stringsAsFactors=FALSE)
  ptn <- "^Map[ ]+([-]?[0-9]+).*$"
  out$id <- as.numeric(gsub(pattern=ptn, replacement="\\1",
                            lines.log10lik))
  ptn <- "^.*=[ ]+([-]?[0-9]+[\\.]]?[0-9]+)$"
  out$log10.lik <- as.numeric(gsub(pattern=ptn, replacement="\\1",
                                   lines.log10lik))

  return(out)
}

##' Parse CarthaGene's output
##'
##' Parse the output of CarthaGene's command "maprintd" and "bestmaprintd", and add the cumulative Kosambi distance.
##' @param out.maprintd vector of character output from \code{\link{runCarthagene}} corresponding to "maprintd" or "bestmaprintd"
##' @return data frame corresponding to the given linkage group
##' @author Timothee Flutre
##' @export
parseCgMaprintd <- function(out.maprintd){
  stopifnot(is.character(out.maprintd))

  out <- NULL

  ## remove empty lines
  out.maprintd <- out.maprintd[- which(out.maprintd == "")]

  ## find the index of the first row of the table
  ptn <- "^Pos[ ]+Id[ ]+name"
  idx <- grep(pattern=ptn, x=out.maprintd, perl=TRUE)
  if(length(idx) == 0){
    msg <- paste0("can't parse CarthaGene's output;",
                  " does it correspond to the 'maprintd' command?")
    stop(msg)
  }
  idx.first <- idx + 1

  ## find the index of the last two rows of the table
  ptn <- "^[ ]+[0-9]+[ ]markers, log10-likelihood ="
  idx <- grep(pattern=ptn, x=out.maprintd, perl=TRUE)
  if(length(idx) == 0){
    msg <- paste0("can't parse CarthaGene's output;",
                  " does it correspond to the 'maprintd' command?")
    stop(msg)
  }
  idx.last <- idx - 1
  stopifnot(idx.last > idx.first)

  ## parse each row of the table
  tmp <- strsplit(out.maprintd[idx.last + 1], " ")[[1]]
  tmp <- tmp[- which(tmp == "")]
  nb.markers <- as.numeric(tmp[1])
  out <- data.frame(id=rep(NA, nb.markers),
                    name=NA,
                    dist.haldane=NA,
                    cum.dist.haldane=NA,
                    dist.kosambi=NA,
                    rec.ratio=NA,
                    lod.2pt=NA,
                    stringsAsFactors=FALSE)
  for(i in idx.first:(idx.last - 2)){
    tokens <- strsplit(out.maprintd[i], "\\s")[[1]]
    tokens <- tokens[- which(tokens == "")]
    out[i - idx.first + 1,] <- tokens[c(2, 3, 4, 6, 8, 10, 12)]
  }
  tokens <- strsplit(out.maprintd[idx.last - 1], "\\s")[[1]]
  tokens <- tokens[- which(tokens == "")]
  out[nb.markers, 1:2] <- tokens[2:3]
  tokens <- strsplit(out.maprintd[idx.last], "\\s")[[1]]
  tokens <- tokens[- which(tokens == "")]
  out[nb.markers, 4] <- tokens[1]

  ## convert to int/num when appropriate
  out$id <- as.integer(out$id)
  out$dist.haldane <- as.numeric(out$dist.haldane)
  out$cum.dist.haldane <- as.numeric(out$cum.dist.haldane)
  out$dist.kosambi <- as.numeric(out$dist.kosambi)
  out$rec.ratio <- as.numeric(out$rec.ratio)
  out$lod.2pt <- as.numeric(out$lod.2pt)

  ## add the cumulative Kosambi distance
  out$cum.dist.kosambi <- NA
  out$cum.dist.kosambi[1] <- 0.0
  out$cum.dist.kosambi[2:nrow(out)] <- cumsum(out$dist.kosambi[-nrow(out)])
  out <- out[, c(1,2,3,4,5,8,6,7)] # re-order columns

  return(out)
}

##' Close CarthaGene
##'
##' Close FIFO used to interact with CarthaGene.
##' @param cg list output from \code{\link{openCarthagene}}
##' @param file path to the file containing the dataset (if not NULL, the temporary file ".2pt", written by CarthaGene based on its name, will be removed)
##' @return nothing
##' @author Matthieu Falque [aut], Timothee Flutre [ctb]
##' @export
closeCarthagene <- function(cg, file=NULL){
  stopifnot(is.list(cg),
            length(cg) >= 4,
            all(c("fifo.in", "fifo.out") %in% names(cg)),
            all(c("name.fifo.in", "name.fifo.out") %in% names(cg)))
  if(! is.null(file))
    stopifnot(is.character(file))

  ## close CG
  cmd <- "exit"
  writeLines(cmd, cg$fifo.in)

  ## close FIFOs
  close(cg$fifo.in)
  close(cg$fifo.out)

  ## remove them
  success <- file.remove(cg$name.fifo.in)
  if(! success){
    msg <- paste0("can't remove input FIFO '", cg$name.fifo.in, "'")
    warning(msg, call.=FALSE, immediate.=TRUE)
  }
  success <- file.remove(cg$name.fifo.out)
  if(! success){
    msg <- paste0("can't remove output FIFO '", cg$name.fifo.out, "'")
    warning(msg, call.=FALSE, immediate.=TRUE)
  }

  ## remove the temporary file ".2pt"
  if(! is.null(file)){
    tmp.file <- paste0(file, ".2pt")
    if(file.exists(tmp.file)){
      success <- file.remove(tmp.file)
      if(! success){
        msg <- paste0("can't remove temporary file '", tmp.file, "'")
        warning(msg, call.=FALSE, immediate.=TRUE)
      }
    }
  }
}
