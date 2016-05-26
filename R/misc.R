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
