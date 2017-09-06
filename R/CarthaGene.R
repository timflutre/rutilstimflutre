## Contains functions useful to wrap CarthaGene.

## Thanks to Matthieu Falque (INRA Moulon)

##' Open CarthaGene
##'
##' Open FIFOs to interact with CarthaGene.
##' @param file path to the file containing the dataset (will be loaded if not NULL)
##' @param task.id identifier of the task (used in FIFO names)
##' @return list
##' @author Matthieu Falque [aut], Timothee Flutre [ctb]
##' @export
openCarthagene <- function(file=NULL, task.id=NULL){
  exe.name <- "carthagene"
  stopifnot(file.exists(Sys.which(exe.name)),
            capabilities()["fifo"])
  if(! is.null(file)){
    stopifnot(file.exists(file))
    file.2pt <- paste0(file, ".2pt")
    if(file.exists(file.2pt)){
      msg <- paste0("file '", file.2pt, "' created by CarthaGene",
                    " already exists\nyou may want to remove it first",
                    " to let CarthaGene recreate it")
      warning(msg)
    }
  }
  if(! is.null(task.id))
    stopifnot(is.character(task.id))

  out <- list()

  name.fifo.in <- "CGfifoIn"
  if(! is.null(task.id))
    name.fifo.in <- paste0(name.fifo.in, "_", task.id)
  name.fifo.out <- "CGfifoOut"
  if(! is.null(task.id))
    name.fifo.out <- paste0(name.fifo.out, "_", task.id)
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
##' @param verbose verbosity level (0/1)
##' @return data frame with one row per marker and columns "id", "name", etc
##' @author Timothee Flutre
##' @export
parseCgMrkinfo <- function(out.mrkinfo, verbose=1){
  stopifnot(is.character(out.mrkinfo))

  if(verbose > 0){
    msg <- "parse the output of 'mrkinfo'..."
    write(msg, stdout())
    utils::flush.console()
  }

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
  colnames(out) <- c("id", "name", "dataset")
  out$id <- as.integer(out$id)

  if(verbose > 0){
    msg <- paste0("nb of data set(s): ", length(unique(out$dataset)),
                  "\ntotal nb of markers: ", nrow(out))
    write(msg, stdout())
  }

  return(out)
}

##' Parse CarthaGene's output
##'
##' Parse the output of CarthaGene's command "group".
##' @param out.group vector of character output from \code{\link{runCarthagene}} corresponding to "group"
##' @param mrk.info data frame output from \code{\link{parseCgMrkinfo}} (optional; useful to have the marker names instead of only the marker identifiers)
##' @param mrk2chr named vector which values are markers and names are chromosomes; before usage, marker names ending in "_m" (an indication of parental coding) will be edited to remove such a suffix
##' @param verbose verbosity level (0/1)
##' @return data frame with at least two columns, "linkage.group" and "id", another column, "locus", if mrk.info is provided, and another column, "chr", if mrk2chr is provided
##' @author Timothee Flutre
##' @export
parseCgGroup <- function(out.group, mrk.info=NULL, mrk2chr=NULL, verbose=1){
  stopifnot(is.character(out.group))
  if(! is.null(mrk.info))
    stopifnot(is.data.frame(mrk.info),
              ncol(mrk.info) >= 2,
              all(c("id", "name") %in% colnames(mrk.info)))
  if(! is.null(mrk2chr))
    stopifnot(is.character(mrk2chr),
              ! is.null(names(mrk2chr)))

  if(verbose > 0){
    msg <- "parse the output of 'group'..."
    write(msg, stdout())
    utils::flush.console()
  }

  idx <- grep("Group ID : Marker ID List ..", out.group)
  if(length(idx) == 0){
    msg <- paste0("can't parse CarthaGene's output;",
                  " does it correspond to the 'group' command?")
    stop(msg)
  }

  out <- list()
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

  out <- cbind(linkage.group=rep(1:nb.linkgroups,
                                 sapply(out, nrow)),
               do.call(rbind, out))

  if(verbose > 0){
    msg <- paste0("nb of linkage groups: ", length(unique(out$linkage.group)),
                  "\nnb of markers per linkage group:")
    write(msg, stdout())
    print(table(out$linkage.group))
  }

  if(! is.null(mrk.info))
    out$locus <- mrk.info$name[match(out$id, mrk.info$id)]

  if(! is.null(mrk2chr)){
    out$chr <- NA
    out.locus <- gsub("_m$", "", out$locus)
    idx <- which(out.locus %in% names(mrk2chr))
    out$chr[idx] <- mrk2chr[out.locus[idx]]
    if(any(is.na(out$chr))){
      msg <- paste0("nb of markers with no chromosome information: ",
                    sum(is.na(out$chr)))
      warning(msg)
    }
    if(verbose > 0){
      msg <- "nb of markers per chromosome in each linkage group:"
      write(msg, stdout())
      lg2chr <- tapply(1:nrow(out), factor(out$linkage.group), function(idx){
        stats::setNames(out$chr[idx], out$locus[idx])
      })
      table.lg2chr <- lapply(lg2chr, table, useNA="always")
      print(table.lg2chr)

      msg <- "nb of markers per linkage group for each chromosome:"
      write(msg, stdout())
      chr2lg <- list()
      chr.names <- unique(out$chr)
      for(chr in chr.names){
        idx <- which(out$chr == chr)
        chr2lg[[chr]] <- stats::setNames(out$linkage.group[idx],
                                         out$locus[idx])
      }
      table.chr2lg <- lapply(chr2lg, table, useNA="always")
      print(table.chr2lg)
    }
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
                    locus=NA,
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

##' Define linkage groups with CarthaGene
##'
##' Defines linkage groups with CarthaGene ("group" command), so that two markers whose 2-point Haldane distance is below the given distance threshold, and 2-point LOD is above the LOD threshold, will be put in the same linkage group.
##' Choose these thresholds by trial and error until it makes sense given the ploidy of the species of interest.
##' @param file path to the file containing the dataset (used if cg is NULL)
##' @param task.id identifier of the task (used if cg is NULL)
##' @param cg list returned by \code{\link{openCarthagene}} (used if file and task.id are NULL)
##' @param dist.thresh distance threshold
##' @param lod.thresh LOD threshold
##' @param mrk2chr named vector which values are markers and names are chromosomes
##' @param close.cg if TRUE, the connection to CarthaGene will be closed and the temporary file removed
##' @param verbose verbosity level (0/1)
##' @return data frame returned by \code{\link{parseCgGroup}}, with at least three columns, "linkage.group", "id" and "locus", and another column, "chr", if mrk2chr is provided
##' @author Timothee Flutre
##' @seealso \code{\link{openCarthagene}}, \code{\link{runCarthagene}}, \code{\link{parseCgMrkinfo}}, \code{\link{parseCgGroup}}, \code{\link{closeCarthagene}}
##' @export
defLinkgroupsWithCarthagene <- function(file=NULL, task.id=NULL, cg=NULL,
                                        dist.thresh=0.3, lod.thresh=3,
                                        mrk2chr=NULL,
                                        close.cg=FALSE, verbose=1){
  stopifnot(xor(all(is.null(file), is.null(task.id)), is.null(cg)),
            is.numeric(dist.thresh),
            is.numeric(lod.thresh))

  if(is.null(cg))
    cg <- openCarthagene(file=file, task.id=task.id)
  if(FALSE){ # for debugging purposes
    out <- runCarthagene(cg, "cgversion")
    out <- runCarthagene(cg, "dsinfo")
  }

  out.mrkinfo <- runCarthagene(cg, "mrkinfo")
  mrk.info <- parseCgMrkinfo(out.mrkinfo, verbose)

  cmd <- paste("group", dist.thresh, lod.thresh)
  if(verbose > 0){
    msg <- paste0("execute '", cmd, "'...")
    write(msg, stdout())
    utils::flush.console()
  }
  st <- system.time(
      out.group <- runCarthagene(cg, cmd))
  if(verbose > 0)
    print(st)
  linkgroups <- parseCgGroup(out.group, mrk.info, mrk2chr, verbose)

  if(close.cg)
    closeCarthagene(cg, file)

  return(linkgroups)
}

##' Estimate marker order and genetic distances with CarthaGene
##'
##' Estimates marker order and genetic distances for a given linkage group with CarthaGene (commands "buildfw", "flips", polish" and "squeeze").
##' @param cg list returned by \code{\link{openCarthagene}}
##' @param linkgroups data frame returned by \code{\link{defLinkgroupsWithCarthagene}}, towhich a column of logicals named "todrop" may have been added, indicating which marker(s) to drop via "mrkdel" before executing "buildfw"
##' @param lg.id identifier of the linkage group of interest
##' @param keep.thresh the minimum difference in log-likelihood between the best insertion point and the second best insertion point required for the map to be considered in the future
##' @param add.thresh the minimum difference in log-likelihood between the best insertion point and the second best insertion point required for a locus to be insertable; this threshold is also used to filter out the differences in log-likelihood reported by the post-processing option (only differences lower than the threshold will be reported)
##' @param init.order if not NULL, vector containing at least three marker identifiers in a given order to start from
##' @param postprocessing if 0, only framework mapping is performed; if 1, post-processing is also applied (see "add.thresh"); if 2, post-processing is applied only to markers in "mrk.list"
##' @param flips.size size of the sliding window for the "flips" command (cannot exceed 9)
##' @param flips.lod.thresh maximum difference of log-likelihood with the best map so that "flips" reports the permutation
##' @param flips.iter if 1, "flips" will be iterated as long as a new, improved map has been found; specify O, otherwise
##' @param squeeze.thresh distance threshold in cMorgan (Haldane distance) for genetics data or cRay for Radiation Hybrids data used by "squeeze" to expunge non-reliable locus from the map
##' @param mrk2phy named vector which names are the locus and values are the physical distances on the chromosome corresponding to the linkage group of interest
##' @param verbose verbosity level (0/1)
##' @return data frame returned by \code{\link{parseCgMaprintd}} containing the best map for the given linkage group according to CarthaGene, and possibly a column "physical.distance"
##' @author Timothee Flutre
##' @seealso \code{\link{openCarthagene}}, \code{\link{runCarthagene}}, \code{\link{defLinkgroupsWithCarthagene}}, \code{\link{parseCgMaprintd}}, \code{\link{closeCarthagene}}
##' @export
estMrkOrderGenDistsWithCarthagene <- function(cg, linkgroups, lg.id,
                                              keep.thresh=3, add.thresh=3,
                                              init.order=NULL, postprocessing=1,
                                              flips.size=4, flips.lod.thresh=3,
                                              flips.iter=1, squeeze.thresh=20,
                                              mrk2phy=NULL, verbose=1){
  stopifnot(is.list(cg),
            is.data.frame(linkgroups),
            all(c("linkage.group", "id", "locus") %in% colnames(linkgroups)),
            lg.id %in% linkgroups$linkage.group,
            is.numeric(keep.thresh),
            is.numeric(add.thresh),
            is.numeric(flips.size),
            flips.size <= 9,
            is.numeric(flips.lod.thresh),
            flips.iter %in% c(0,1),
            is.numeric(squeeze.thresh))
  if("todel" %in% colnames(linkgroups))
    stopifnot(is.logical(linkgroups$todel))
  if(! is.null(init.order))
    stopifnot(is.vector(init.order),
              length(init.order) >= 3,
              all(init.order %in% linkgroups$id))
  if(! is.null(mrk2phy))
    stopifnot(is.vector(mrk2phy),
              is.numeric(mrk2phy),
              ! is.null(names(mrk2phy)))

  cmd <- paste0("mrkselset [groupget ", lg.id, "]")
  runCarthagene(cg, cmd)
  is.lg.todel <- FALSE
  if("todel" %in% colnames(linkgroups))
    is.lg.todel <- linkgroups$linkage.group == lg.id & linkgroups$todel
  if(any(is.lg.todel)){
    for(mrk.name in linkgroups$locus[is.lg.todel])
      runCarthagene(cg, paste0("mrkdel ", mrk.name))
  }
  if(is.null(init.order)){
    init.order <- "{}"
  } else
    init.order <- paste("{", paste(init.order, collapse=" "), "}")
  cmd <- paste("buildfw", keep.thresh, add.thresh, init.order, postprocessing)
  if(verbose > 0){
    msg <- paste0("execute '", cmd, "'...")
    write(msg, stdout())
    utils::flush.console()
  }
  st <- system.time(
      out.buildfw <- runCarthagene(cg, cmd))
  if(verbose > 0)
    print(st)

  if(verbose > 0){
    out.heaprint <- runCarthagene(cg, "heaprint")
    print(maps.info <- parseCgHeaprint(out.heaprint))
  }

  cmd <- paste("flips", flips.size, flips.lod.thresh, flips.iter)
  if(verbose > 0){
    msg <- paste0("execute '", cmd, "'...")
    write(msg, stdout())
    utils::flush.console()
  }
  st <- system.time(
      out.flips <- runCarthagene(cg, cmd))
  if(verbose > 0)
    print(st)

  if(verbose > 0){
    msg <- "execute 'polish'..."
    write(msg, stdout())
    utils::flush.console()
  }
  st <- system.time(
      out.polish <- runCarthagene(cg, "polish"))
  if(verbose > 0)
    print(st)

  cmd <- paste("squeeze", squeeze.thresh)
  if(verbose > 0){
    msg <- paste0("execute '", cmd, "'...")
    write(msg, stdout())
    utils::flush.console()
  }
  st <- system.time(
      out.squeeze <- runCarthagene(cg, cmd))
  if(verbose > 0)
    print(st)

  if(verbose > 0){
    out.heaprint <- runCarthagene(cg, "heaprint")
    print(maps.info <- parseCgHeaprint(out.heaprint))
  }

  if(verbose > 0){
    msg <- "execute 'bestprintd'..."
    write(msg, stdout())
    utils::flush.console()
  }
  out.bestprintd <- runCarthagene(cg, "bestprintd")
  bestmap <- parseCgMaprintd(out.bestprintd)
  if(! is.null(mrk2phy)){
    bestmap$physical.distance <- NA
    locus.init <- gsub("_m$", "", bestmap$locus)
    has.phy <- locus.init %in% names(mrk2phy)
    if(any(has.phy))
      bestmap$physical.distance[has.phy] <- mrk2phy[locus.init[has.phy]]
  }

  return(bestmap)
}
