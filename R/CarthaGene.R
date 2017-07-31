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
##' @return vector of lines
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
