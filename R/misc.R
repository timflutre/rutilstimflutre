##' @import data.table
##' @import lme4
##' @import Matrix
##' @import Rcpp
##' @import stats

## https://stackoverflow.com/questions/66816638
utils::globalVariables(".")

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
##' Remove extension(s) from a file name.
##' @param file character (will be passed to \code{\link{strsplit}})
##' @param fileext character containing a regular expression (will be passed to \code{\link{strsplit}})
##' @return character
##' @author Timothee Flutre
##' @examples
##' f <- "data.txt.gz"
##' removeFileExtension(file=f, fileext="\\.gz")
##' removeFileExtension(file=f, fileext="\\.txt\\.gz")
##' @export
removeFileExtension <- function(file, fileext){
  stopifnot(is.character(file),
            is.character(fileext))
  out <- strsplit(x=file, split=fileext)[[1]][1]
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

##' Inline function in formula
##'
##' Detect the presence of inline function(s) in a formula and parse them.
##' @param form character corresponding to a formula, e.g. \code{"y ~ 1 + x"}
##' @param only.resp if TRUE, the presence of inline function(s) is detected only in the response
##' @return list with one component per term of the formula (first the response and then the predictors, if any), each component containing a vector of characters with the untransformed variable first and then the inline function(s) from the outermost to the innermost, or NA if there is no inline function corresponding to this term
##' @author Timothee Flutre
##' @export
inlineFctForm <- function(form, only.resp=TRUE){
  stopifnot(is.character(form))
  if(! only.resp){
    msg <- "only.resp=FALSE is not implemented (yet)"
    stop(msg)
  }

  out <- list()

  if(only.resp)
    form.terms <- trimws(strsplit(form, "~")[[1]][1])

  for(x in form.terms){
    out[[x]] <- NA
    untransf <- regmatches(x,
                           regexec("\\(([^\\(\\)]*)\\)", x))[[1]]
    if(length(untransf) == 2){
      out[[x]] <- untransf[2]
      inFcts <- regmatches(x,
                           gregexpr("([^\\(]*)\\(", x))[[1]]
      if(length(inFcts) > 0)
        out[[x]] <- c(out[[x]],
                      sub("\\(", "", inFcts))
    }
  }

  return(out)
}

##' Reformat
##'
##' Reformats the upper triangular part (and, optionally, the diagonal) of a matrix from the wide format into the long format.
##' @param mat matrix
##' @param diag if FALSE, the diagonal entries will be ignored
##' @param outColNs names of the columns in the output
##' @param stringsAsFactors if TRUE, character vectors will be converted to factors
##' @return data frame
##' @author Zheyuan Li [aut] (https://stackoverflow.com/a/41745373/597069), Timothee Flutre [ctb]
##' @export
matWide2Long <- function(mat, diag=TRUE, outColNs=c("row","col","val"),
                         stringsAsFactors=FALSE){
  stopifnot(length(outColNs) == 3,
            all(is(outColNs, "character")))

  ind <- which(upper.tri(mat, diag=TRUE), arr.ind=TRUE)
  nn <- dimnames(mat)
  out <- data.frame(row = nn[[1]][ind[, 1]],
                    col = nn[[2]][ind[, 2]],
                    val  = mat[ind],
                    stringsAsFactors=stringsAsFactors)
  colnames(out) <- outColNs

  if(! diag){
    idx <- which(out[,1] != out[,2])
    out <- droplevels(out[idx,])
  }
  rownames(out) <- NULL

  return(out)
}
