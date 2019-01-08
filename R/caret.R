## Contains functions useful for caret.
## https://topepo.github.io/caret/using-your-own-model-in-train.html#model-components

##' Summary for caret
##'
##' Returns various metrics comparing predicted data with the values of the hold-out fold.
##' @param data data frame with a column named \code{"obs"} and another named \code{"pred"}
##' @param lev levels (unused here)
##' @param model model (unused here)
##' @param plot if TRUE, observations regressed on predictions will be plotted
##' @return vector with the root mean square error, variances of observations and predictions (if 0, most other metrics will be NA or NaN), Pearson and Spearman correlations between all data points, the 50\% best and the 25\% best, as well as the intercept, slope and determination coefficient of the simple linear regression \code{lm(obs ~ pred)} (\href{https://doi.org/10.1016/j.ecolmodel.2008.05.006}{Pineiro et al., 2008}), and the modelling efficiency (\href{https://doi.org/10.1016/0304-3800(93)90105-2}{Mayer and Butler, 1993})
##' @author Timothee Flutre
##' @examples
##' \dontrun{set.seed(1859)
##' n <- 100
##'
##' ## without bias and low coef of determination
##' d <- data.frame(pred=rnorm(n=n, mean=30, sd=10))
##' d$obs <- 0 + 1 * d$pred + rnorm(n=n, mean=0, sd=8)
##' (out1 <- caretSummary(data=d, plot=TRUE))
##'
##' ## with bias and high coef of determination
##' d <- data.frame(pred=rnorm(n=n, mean=30, sd=10))
##' d$obs <- 10 + 0.5 * d$pred + rnorm(n=n, mean=0, sd=1)
##' (out2 <- caretSummary(data=d, plot=TRUE))
##' }
##' @export
caretSummary <- function(data, lev=NULL, model=NULL, plot=FALSE){
  stopifnot(is.data.frame(data),
            all(c("obs", "pred") %in% colnames(data)))

  data <- data[order(data$obs, decreasing=TRUE),] # sort best -> worse
  nb.inds <- nrow(data)
  idx.best50p <- floor(0.50 * nb.inds)
  idx.best25p <- floor(0.25 * nb.inds)
  fit <- stats::lm(obs ~ pred, data=data)
  coefOls <- as.numeric(stats::coef(fit))
  R2 <- summary(fit)$r.squared
  R2.adj <- summary(fit)$adj.r.squared
  mod.eff <- 1 - sum((data$obs - data$pred)^2) /
    sum((data$obs - mean(data$obs))^2)

  if(plot){
    graphics::plot(x=data$pred, y=data$obs,
                   main="Assessment of prediction accuracy",
                   xlab="Predictions", ylab="Observations",
                   asp=1, las=1, pch=20)
    graphics::abline(a=0, b=1, lty=2)
    graphics::abline(fit, lty=1)
    graphics::legend("topleft", bty="n",
                     legend=sapply(c(bquote("Adjusted"~R^2==.(round(R2,2))),
                                     bquote(y==.(round(coefOls[1],2))+.(round(coefOls[2],2))~"x")),
                                   as.expression))
  }

  out <- c(rmse=sqrt(mean((data$pred - data$obs)^2)),
           var.obs=stats::var(data$obs),
           var.pred=stats::var(data$pred),
           corP=stats::cor(data$obs, data$pred, method="pearson"),
           corS=stats::cor(data$obs, data$pred, method="spearman"),
           corP.best50p=stats::cor(data$obs[1:idx.best50p],
                                   data$pred[1:idx.best50p],
                                   method="pearson"),
           corS.best50p=stats::cor(data$obs[1:idx.best50p],
                                   data$pred[1:idx.best50p],
                                   method="spearman"),
           corP.best25p=stats::cor(data$obs[1:idx.best25p],
                                   data$pred[1:idx.best25p],
                                   method="pearson"),
           corS.best25p=stats::cor(data$obs[1:idx.best25p],
                                   data$pred[1:idx.best25p],
                                   method="spearman"),
           reg.intercept=coefOls[1],
           reg.slope=coefOls[2],
           reg.R2=R2,
           reg.R2.adj=R2.adj,
           mod.eff=mod.eff)
  return(out)
}

##' Fit with rrBLUP for caret
##'
##' @param x current predictors used to fit the model
##' @param y current outcome used to fit the model
##' @param wts optional instance weights
##' @param param current tuning parameter values
##' @param lev class levels of the outcome (or NULL in regression)
##' @param last logical for whether the current fit is the final fit
##' @param weights ?
##' @param classProbs logical for whether class probabilities should be computed
##' @param ... arguments passed on to \code{rrBLUP::mixed.solve}
##' @return output of \code{rrBLUP::mixed.solve}
##' @author Timothee Flutre
##' @export
caretFitRrblup <- function(x, y, wts, param, lev, last, weights, classProbs, ...){
  rrBLUP::mixed.solve(y=y, Z=x, K=NULL, ...)
}

##' Predict with rrBLUP for caret
##'
##' @param modelFit model produced by \code{\link{caretFitRrblup}}
##' @param newdata predictor values of the instances being predicted (e.g. out-of-bag samples)
##' @param submodels optional list of tuning parameters only used with the "loop" element
##' @return vector
##' @author Timothee Flutre
##' @export
caretPredictRrblup <- function(modelFit, newdata, submodels=NULL){
  newdata %*% modelFit$u
}

##' Grid with rrBLUP for caret
##'
##' @param x predictors
##' @param y outcome
##' @param len value of \code{tuneLength} that is potentially passed in through \code{train}
##' @param search either \code{"grid"} or \code{"random"}
##' @return data frame of tuning parameter combinations with a column for each parameter
##' @author Timothee Flutre
##' @export
caretGridRrblup <- function(x, y, len=NULL, search="grid"){
  data.frame(intercept=TRUE)
}

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

##' Fit with qtl for caret
##'
##' Fits the simple or multiple interval mapping model for usage with caret.
##' First, the genotypes in the JoinMap format are saved into a temporary file.
##' Second, a \code{cross} object is made.
##' Third, SIMQTL or MIMQTL is run (calc.genoprob, scanone or scantwo).
##' Fourth, allelic effects are estimated.
##' @param x current predictors used to fit the model
##' @param y current outcome used to fit the model
##' @param wts optional instance weights
##' @param param current tuning parameter values
##' @param lev class levels of the outcome (or NULL in regression)
##' @param last logical for whether the current fit is the final fit
##' @param weights ?
##' @param classProbs logical for whether class probabilities should be computed
##' @param genmap genetic map
##' @param pop.type population type in the JoinMap format
##' @param nperm number of permutations
##' @param threshold.alpha vector
##' @param plot logical
##' @param phase ?
##' @param QTLmethod character (SIM/MIM)
##' @param QTLposition ?
##' @param p2d path to directory with results of \code{scantwo} (passed to MIMQTL)
##' @param nb.cores number of cores (passed to MIMQTL)
##' @param verbose verbosity level (0/1/2)
##' @return output of \code{\link{SIMQTL}} or \code{\link{MIMQTL}}
##' @author Timothee Flutre
##' @export
caretFitQtl <- function(x, y, wts, param, lev, last, weights, classProbs,
                        genmap, pop.type="CP", nperm=10, threshold.alpha=0.05,
                        plot=FALSE, phase, QTLmethod="SIM",
                        p2d="", nb.cores=1,
                        QTLposition=NULL, verbose=0){
  requireNamespace("qtl")

  stopifnot(pop.type == "CP",
            QTLmethod %in% c("SIM", "MIM"))
  stopifnot(xor(QTLmethod == "MIM", length(plot) == 1),
            xor(QTLmethod == "SIM", length(plot) == 2),
            length(y) == dim(x)[1], length(phase) == dim(x)[2],
            dim(genmap)[1] == dim(x)[2],
            class(y) == "numeric", class(x) == "matrix",
            class(threshold.alpha) == "numeric",
            class(plot) == "logical",
            class(nperm) == "numeric")

  ## here x must have genotypes in rows and locus in columns
  Genotypes <- rownames(x)
  Locus <- colnames(x)
  x <- t((x))

  ## here x must have genotypes in columns and locus in rows
  rownames(x) <- Locus
  x <- as.data.frame(x)
  colnames(x) <- NULL
  tmpf <- tempfile(pattern="genos_joinmap", fileext=".loc")
  writeSegregJoinMap(pop.name="caretFitQtl",
                     pop.type=pop.type,
                     locus=rownames(x),
                     segregs=getJoinMapSegregs(x),
                     genos=as.data.frame(x),
                     phases=phase,
                     file=tmpf,
                     verbose=verbose)

  gendat <- qtl:::read.cross.mq.loc(tmpf)
  gendat <- as.matrix(gendat$genotypes)
  rownames(gendat) <- Genotypes # in gendat, genotypes are in row and locus in columns
  stopifnot(dim(gendat)[1] == dim(x)[2],
            dim(gendat)[2] == dim(x)[1])
  is.rm <- file.remove(tmpf)
  phenos <- data.frame(Genotype=Genotypes, y=y)
  phenos$Genotype <- as.character(phenos$Genotype)
  stopifnot(dim(phenos)[1] == dim(x)[2])
  if(pop.type == "CP")
    cross.type <- "4way"
  cross <- setupQtlCrossObject(gendat=gendat, cross.type="4way",
                               genmap=genmap, phenos=phenos)
  stopifnot(all(class(cross) == c("4way", "cross")))
  ## cross$geno$chr$map is in the wrong format (numeric instead of matrix)
  nb_chr <- length(lapply(cross$geno, attributes))
  for(i in 1:nb_chr){
    cross[["geno"]][[i]][["map"]] <- as.matrix(cross[["geno"]][[i]][["map"]])
    cross[["geno"]][[i]][["map"]] <- rbind(t(cross[["geno"]][[i]][["map"]]), t(cross[["geno"]][[i]][["map"]]))
  }

  ## QTL detection (calc.genoprob, scanone or scantwo with permutations)
  colnames(cross$pheno)[1] <- "indiv"
  if(QTLmethod == "SIM"){
    fit.qtl <- SIMQTL(cross, response.in.cross=TRUE, numeric.chr.format=TRUE,
                      pheno.col="y", geno.joinmap=x, threshold.alpha=threshold.alpha,
                      nperm=nperm, method="em", phase=phase, plot=plot,
                      QTL_position=NULL, verbose=verbose)
  } else if(QTLmethod == "MIM"){
    fit.qtl <- MIMQTL(cross, response.in.cross=TRUE,
                      pheno.col="y", geno.joinmap=x,
                      range_nb_qtl_max=seq(1:4), response=NULL,
                      numeric.chr.format=TRUE, nrun=10, nperm=nperm,
                      additive.only=TRUE, method="hk", phase=phase, # em method is too long
                      p2d=p2d,
                      scan2file="",  type.CI="LOD-1",
                      threshold.alpha=threshold.alpha,
                      QTL_position=NULL,
                      nb.cores=nb.cores,
                      verbose=verbose)
  }

  return(fit.qtl)
}

##' Predict with qtl for caret
##'
##' @param modelFit model produced by \code{\link{caretFitQtl}}; should be a list with a component named \code{"allelic.effects"} which contains a vector in the same order as the columns of newdata and, as values, the estimated allele effects (0 if not significant, non-zero otherwise)
##' @param newdata predictor values of the instances being predicted (e.g. out-of-bag samples); should be in the JoinMap format
##' @param submodels optional list of tuning parameters only used with the "loop" element
##' @return vector
##' @author Charlotte Brault [aut], Timothee Flutre [ctb]
##' @export
caretPredictQtl <- function(modelFit, newdata, submodels=NULL) {
  stopifnot(is.list(modelFit),
            "allelic.effects" %in%names(modelFit),
            ! is.null(rownames(newdata)))

  newdata <- as.data.frame(t(newdata))
  loc.seg.phase.clas <- data.frame(locus=rownames(newdata),
                                   seg=getJoinMapSegregs(newdata),
                                   phase=NA, clas=NA,
                                   stringsAsFactors=FALSE)
  jm <- cbind(loc.seg.phase.clas, newdata)
  X <- joinMap2designMatrix(jm=jm, use.phase=FALSE, rm.col.zeros=FALSE,
                            verbose=0)
  X <- X[, -1] # remove intercept
  allelic.effects <- as.matrix(modelFit$allelic.effects$effect)
  rownames(allelic.effects) <- modelFit$allelic.effects$predictor

  ## remove allelic effects that are not in modelFit$allelic.effects
  X <- subset(X, select=rownames(allelic.effects))

  stopifnot(dim(X)[2] == nrow(allelic.effects))

  predictions <- X %*% allelic.effects
  return(predictions[,1])
}

##' Grid with qtl for caret
##'
##' @param x predictors
##' @param y outcome
##' @param len value of \code{tuneLength} that is potentially passed in through \code{train}
##' @param search either \code{"grid"} or \code{"random"}
##' @return data frame of tuning parameter combinations with a column for each parameter
##' @author Timothee Flutre
##' @export
caretGridQtl <- function(x, y, len=NULL, search="grid"){
  data.frame(intercept=TRUE)
}
