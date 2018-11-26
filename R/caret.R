## Contains functions useful for caret.
## https://topepo.github.io/caret/using-your-own-model-in-train.html#model-components

##' Fit with varbvs for caret
##'
##' @param x current predictors used to fit the model
##' @param y current outcome used to fit the model
##' @param wts optional instance weights
##' @param param current tuning parameter values
##' @param lev class levels of the outcome (or NULL in regression)
##' @param last logical for whether the current fit is the final fit
##' @param weights ?
##' @param classProbs logical for whether class probabilities should be computed 
##' @param ... arguments passed on to \code{varbvs::varbvs}
##' @return output of \code{varbvs::varbvs}
##' @author Timothee Flutre
##' @export
caretFitVarbvs <- function(x, y, wts, param, lev, last, weights, classProbs, ...){
  varbvs::varbvs(X=x, y=y, weights=wts, ...)
}

##' Predict with varbvs for caret
##'
##' @param modelFit model produced by \code{\link{caretFitVarbvs}}
##' @param newdata predictor values of the instances being predicted (e.g. out-of-bag samples)
##' @param submodels optional list of tuning parameters only used with the "loop" element
##' @return vector
##' @author Timothee Flutre
##' @export
caretPredictVarbvs <- function(modelFit, newdata, submodels=NULL){
  stats::predict(object=modelFit, X=newdata, Z=NULL)
}

##' Grid with varbvs for caret
##'
##' @param x predictors
##' @param y outcome
##' @param len value of \code{tuneLength} that is potentially passed in through \code{train}
##' @param search either \code{"grid"} or \code{"random"}
##' @return data frame of tuning parameter combinations with a column for each parameter
##' @author Timothee Flutre
##' @export
caretGridVarbvs <- function(x, y, len=NULL, search="grid"){
  data.frame(intercept=TRUE)
}

##' Fit with BGLR for caret
##'
##' The argument named ETA is compulsory for \code{BGLR::BGLR}.
##' Only its first component will be taken into account.
##' The design matrix of its first component will be set to be \code{x}.
##' The model of its first component is used in output files.
##' @param x current predictors used to fit the model
##' @param y current outcome used to fit the model
##' @param wts optional instance weights
##' @param param current tuning parameter values
##' @param lev class levels of the outcome (or NULL in regression)
##' @param last logical for whether the current fit is the final fit
##' @param weights ?
##' @param classProbs logical for whether class probabilities should be computed 
##' @param ETA two-level list used to specify the regression function for \code{BGLR::BGLR}
##' @param saveAt string that may include a path and a pre-fix that will be added to the name of the files that are saved as \code{BGLR::BGLR} runs
##' @param keep.samples logical for whether the samples will be returned (as a \code{coda::mcmc.list})
##' @param nIter number of iterations 
##' @param burnIn number of burn-in
##' @param thin thinning
##' @param ... arguments passed on to \code{BGLR::BGLR}
##' @return output of \code{BGLR::BGLR} and, optionally, the samples 
##' @author Timothee Flutre
##' @export
caretFitBglr <- function(x, y, wts, param, lev, last, weights, classProbs,
                         ETA, saveAt, nIter, burnIn, thin, keep.samples=FALSE,
                         ...){
  ETA[[1]]$X <- x
  if(missing(wts))
    wts <- NULL

  out <- BGLR::BGLR(y=y, ETA=ETA, saveAt=saveAt, weights=wts, nIter=nIter,
                    burnIn=burnIn, thin=thin, ...)

  if(keep.samples){
    out$post.samples <- cbind(mu=utils::read.table(paste0(saveAt, "mu.dat"))[,1],
                              varE=utils::read.table(paste0(saveAt, "varE.dat"))[,1],
                              pi=utils::read.table(paste0(saveAt, "ETA_1_par",
                                                   ETA[[1]]$model, ".dat"))[,1],
                              varB=utils::read.table(paste0(saveAt, "ETA_1_par",
                                                            ETA[[1]]$model, ".dat"))[,2])
    out$post.samples <- coda::mcmc.list(coda::mcmc(out$post.samples, start=thin,
                                       end=nIter, thin=thin))
    out$post.samples <- stats::window(out$post.samples, start=burnIn+1)
    out$ess <- coda::effectiveSize(out$post.samples)
  }

  for(f in c(paste0(saveAt, c("mu.dat", "varE.dat",
                               paste0("ETA_1_par", ETA[[1]]$model, ".dat"),
                               "ETA_1_b.bin"))))
    if(file.exists(f))
      file.remove(f)

  out
}

##' Predict with BGLR for caret
##'
##' @param modelFit model produced by \code{\link{caretFitBglr}}
##' @param newdata predictor values of the instances being predicted (e.g. out-of-bag samples)
##' @param submodels optional list of tuning parameters only used with the "loop" element
##' @return vector
##' @author Timothee Flutre
##' @export
caretPredictBglr <- function(modelFit, newdata, submodels=NULL){
  out <- newdata %*% modelFit$ETA[[1]]$b
  out[,1]
}

##' Grid with BGLR for caret
##'
##' @param x predictors
##' @param y outcome
##' @param len value of \code{tuneLength} that is potentially passed in through \code{train}
##' @param search either \code{"grid"} or \code{"random"}
##' @return data frame of tuning parameter combinations with a column for each parameter
##' @author Timothee Flutre
##' @export
caretGridBglr <- function(x, y, len=NULL, search="grid"){
  data.frame(intercept=TRUE)
}
