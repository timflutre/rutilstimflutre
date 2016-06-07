##' @import data.table
##' @import lme4
##' @import Matrix

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
##' removeFileExtension(file=f, fileext="gz")
##' @export
removeFileExtension <- function(file, fileext){
  stopifnot(is.character(file))

  out <- file

  tmp <- strsplit(file, "\\.")[[1]]
  if(tmp[length(tmp)] == fileext)
    out <- paste0(tmp[1:(length(tmp)-1)], collapse=".")

  return(out)
}
