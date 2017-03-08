##' @import data.table
##' @import lme4
##' @import Matrix
##' @import Rcpp

.onAttach <- function(libname, pkgname) {
  msg <- paste0("package '", pkgname,
                "' (version ", utils::packageVersion(pkgname), ")",
                " is loaded",
                "\ndev at https://github.com/timflutre/rutilstimflutre")
  packageStartupMessage(msg)
}

##' Namespaces
##'
##' Stops if any package's namespace can't be loaded.
##' @param packages vector of package names
##' @return nothing
##' @author Timothee Flutre
##' @export
requireNamespaces <- function(packages){
  stopifnot(is.vector(packages),
            is.character(packages))
  for(pkg in packages)
    if(! requireNamespace(pkg, quietly=TRUE))
      stop(paste0("package '", pkg, "' needed for this function to work."),
           call.=FALSE)
}

##' File
##'
##' Remove the last extension in a file name.
##' @param file character
##' @param fileext character
##' @return character
##' @author Timothee Flutre
##' @examples
##' f <- "data.txt.gz"
##' removeFileExtension(file=f, fileext=".gz")
##' removeFileExtension(file=f, fileext=".txt.gz")
##' @export
removeFileExtension <- function(file, fileext){
  stopifnot(is.character(file),
            strsplit(fileext, "")[[1]][1] == ".")
  out <- strsplit(file, fileext)[[1]][1]
  return(out)
}

##' Data.frame
##'
##' Convert the factor columns into character.
##' @param x data.frame
##' @return data.frame
##' @author Timothee Flutre
##' @export
convertFactorColumnsToCharacter <- function(x){
  stopifnot(is.data.frame(x))
  idx <- sapply(x, is.factor)
  x[idx] <- lapply(x[idx], as.character)
  return(x)
}
