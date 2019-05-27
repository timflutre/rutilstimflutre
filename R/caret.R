## Contains functions useful for caret.
## https://topepo.github.io/caret/using-your-own-model-in-train.html#model-components

##' Summary for caret
##'
##' Returns various metrics comparing predicted data with the values of the hold-out fold.
##' @param data data frame with a column named \code{"obs"} and another named \code{"pred"}
##' @param lev levels (unused here)
##' @param model model (unused here)
##' @param plot if TRUE, observations regressed on predictions will be plotted
##' @return vector with the root mean square error, variances of observations and predictions (if 0, most other metrics will be NA or NaN), Pearson and Spearman correlations between all data points, the 50\% best and the 25\% best, as well as the intercept, slope and determination coefficient of the simple linear regression \code{lm(obs ~ pred)} (\href{https://doi.org/10.1016/j.ecolmodel.2008.05.006}{Pineiro et al., 2008}), the statistic and p value for testing null bias (\href{https://tel.archives-ouvertes.fr/tel-00985747v2}{Baey, 2014, pages 52-53}), and the modelling efficiency (\href{https://doi.org/10.1016/0304-3800(93)90105-2}{Mayer and Butler, 1993})
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
  is.pred.constant <- ifelse(stats::sd(data$pred) == 0, TRUE, FALSE)
  nb.inds <- nrow(data)
  idx.best50p <- floor(0.50 * nb.inds)
  idx.best25p <- floor(0.25 * nb.inds)
  fit <- stats::lm(obs ~ pred, data=data)
  coefOls <- as.numeric(stats::coef(fit))
  R2 <- summary(fit)$r.squared
  R2.adj <- summary(fit)$adj.r.squared
  m2.minus.m1 <- nb.inds * coefOls[1]^2 +
    2 * nb.inds * coefOls[1] * (coefOls[2] - 1) * mean(data$pred) +
    (coefOls[2] - 1)^2 * sum(data$pred^2)
  y.minus.m1 <- sum((data$obs - stats::fitted(fit))^2)
  stat.nobias <- (m2.minus.m1 / 2) /
    (y.minus.m1 / (nb.inds - 2))
  pval.nobias <- stats::pf(q=stat.nobias, df1=2, df2=nb.inds - 2,
                           lower.tail=FALSE)
  mod.eff <- 1 - sum((data$obs - data$pred)^2) /
    sum((data$obs - mean(data$obs))^2)

  if(plot){
    graphics::plot(x=data$pred, y=data$obs,
                   main="Assessment of prediction accuracy",
                   xlab="Predictions", ylab="Observations",
                   asp=1, las=1, pch=20)
    graphics::abline(a=0, b=1, lty=2)
    if(! is.pred.constant)
      graphics::abline(fit, lty=1)
    graphics::legend("topleft", bty="n",
                     legend=sapply(c(bquote("Adjusted"~R^2==.(round(R2,2))),
                                     bquote(y==.(round(coefOls[1],2))+.(round(coefOls[2],2))~"x")),
                                   as.expression))
  }

  out <- c(rmse=sqrt(mean((data$pred - data$obs)^2)),
           var.obs=stats::var(data$obs),
           var.pred=stats::var(data$pred),
           corP=ifelse(is.pred.constant, NA,
                       stats::cor(data$obs, data$pred, method="pearson")),
           corS=ifelse(is.pred.constant, NA,
                       stats::cor(data$obs, data$pred, method="spearman")),
           corP.best50p=ifelse(is.pred.constant, NA,
                               stats::cor(data$obs[1:idx.best50p],
                                          data$pred[1:idx.best50p],
                                          method="pearson")),
           corS.best50p=ifelse(is.pred.constant, NA,
                               stats::cor(data$obs[1:idx.best50p],
                                          data$pred[1:idx.best50p],
                                          method="spearman")),
           corP.best25p=ifelse(is.pred.constant, NA,
                               stats::cor(data$obs[1:idx.best25p],
                                          data$pred[1:idx.best25p],
                                          method="pearson")),
           corS.best25p=ifelse(is.pred.constant, NA,
                               stats::cor(data$obs[1:idx.best25p],
                                          data$pred[1:idx.best25p],
                                          method="spearman")),
           reg.intercept=coefOls[1],
           reg.slope=coefOls[2],
           reg.R2=R2,
           reg.R2.adj=R2.adj,
           stat.nobias=stat.nobias,
           pval.nobias=pval.nobias,
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
##' @param phase vector of character (length nb of markers) in the format '{--}'
##' @param QTLmethod character (SIM/MIM)
##' @param QTL_position data frame with columns linkage.group, genetic.distance and locus to plot abline at simulated QTL position
##' @param p2d path to directory with results of \code{scantwo} (passed to MIMQTL)
##' @param nb.cores number of cores (passed to MIMQTL)
##' @param cross.geno.prob list containing genotypic class probabilities calculated from calc.genoprob.
##' If not passed (NULL), prediction is done with allelic effects at the closest marker.
##' @param verbose verbosity level (0/1/2)
##' @return output of \code{\link{SIMQTL}} or \code{\link{MIMQTL}}
##' @author Timothee Flutre
##' @export
caretFitQtl <- function (x, y, wts, param, lev, last, weights, classProbs,
                         genmap, tuneThreshold=TRUE, alpha = c(0.05, NA),
                         pop.type = "CP", nperm = 10,  plot = FALSE, 
                         phase, QTLmethod = "SIM", p2d = "", 
                         nb.cores =parallel::detectCores()-2,
                         QTL_position = NULL, 
                         cross.geno.prob = NULL, verbose = 0) 
{
  requireNamespace("qtl")
  stopifnot(pop.type == "CP", QTLmethod %in% c("SIM", "MIM"))
  stopifnot(xor(QTLmethod == "MIM", length(plot) == 1), 
            xor(QTLmethod == "SIM", length(plot) == 2),
            length(y) == dim(x)[1], length(phase) == dim(x)[2],
            dim(genmap)[1] == dim(x)[2],
            class(y) == "numeric", 
            class(x) == "matrix", 
            class(plot) == "logical", class(nperm) == "numeric")
  Genotypes <- rownames(x)
  Genotypes <- make.names(Genotypes)
  Locus <- colnames(x)
  x <- t((x))
  rownames(x) <- Locus
  x <- as.data.frame(x)
  colnames(x) <- NULL
  tmpf <- tempfile(pattern = "genos_joinmap", fileext = ".loc")
  writeSegregJoinMap(pop.name = "caretFitQtl", pop.type = pop.type, 
                     locus = rownames(x), segregs = getJoinMapSegregs(x), 
                     genos = as.data.frame(x), phases = phase, file = tmpf, 
                     verbose = verbose)
  gendat <- qtl:::read.cross.mq.loc(tmpf)
  gendat <- as.matrix(gendat$genotypes)
  rownames(gendat) <- Genotypes
  stopifnot(dim(gendat)[1] == dim(x)[2], dim(gendat)[2] == 
              dim(x)[1])
  is.rm <- file.remove(tmpf)
  phenos <- data.frame(Genotype = Genotypes, y = y)
  phenos$Genotype <- as.character(phenos$Genotype)
  stopifnot(dim(phenos)[1] == dim(x)[2])
  if (pop.type == "CP") 
    cross.type <- "4way"
  cross <- setupQtlCrossObject(gendat = gendat, cross.type = "4way", 
                               genmap = genmap, phenos = phenos)
  for (i in 1:length(cross$geno)) {
    cross$geno[[i]]$map <- matrix(data = rep(cross$geno[[i]]$map, 2), nrow = 2, byrow = TRUE)
    colnames(cross$geno[[i]]$map) <- colnames(cross$geno[[i]]$data)
  }
  stopifnot(all(class(cross) == c("4way", "cross")))
  colnames(cross$pheno)[1] <- "indiv"
  if (QTLmethod == "SIM") {
    fit.qtl <- SIMQTL(cross, response.in.cross = TRUE, numeric.chr.format = TRUE, 
                      pheno.col = "y", geno.joinmap = x, 
                      threshold = ifelse(tuneThreshold, param$threshold, NA), 
                      alpha=alpha,
                      nperm = nperm, method = "em", phase = phase, plot = plot, 
                      QTL_position = QTL_position, verbose = verbose)
  }
  else if (QTLmethod == "MIM") {
    fit.qtl <- MIMQTL(cross, response.in.cross = TRUE, pheno.col = "y", 
                      geno.joinmap = x, 
                      threshold = ifelse(tuneThreshold, param$threshold, NA), 
                      alpha = alpha, 
                      range_nb_qtl_max = seq(1:4), response = NULL, 
                      numeric.chr.format = TRUE, nrun = 10, nperm = nperm, 
                      additive.only = TRUE, method = "hk", phase = phase, 
                      p2d = p2d, scan2file = "", type.CI = "LOD-1",
                      QTL_position = QTL_position, nb.cores = nb.cores, 
                      verbose = verbose)
  }
  if (!is.null(cross.geno.prob)) {
    chr_names <- names(lapply(lapply(cross$geno, names), 
                              names))
    list_locus <- 0
    for (chr in chr_names) {
      tmp <- names(attr(cross.geno.prob[[chr]], "map")[1,])
      tmp <- paste0(chr, "@", tmp)
      list_locus <- append(list_locus, tmp)
    }
    list_locus <- list_locus[-1]
    list_locus_clgeno <- 0
    for (i in 1:length(list_locus)) {
      tmp <- c(paste0(list_locus[i], "@AC"), 
               paste0(list_locus[i], "@AD"), 
               paste0(list_locus[i], "@BC"),
               paste0(list_locus[i], "@BD"))
      list_locus_clgeno <- append(list_locus_clgeno, tmp)
    }
    list_locus_clgeno <- list_locus_clgeno[-1]
    stopifnot(length(list_locus) * 4 == length(list_locus_clgeno))
    geno_prob <- matrix(nrow = length(Genotypes), ncol = length(list_locus_clgeno), 
                        dimnames = list(Genotypes, list_locus_clgeno))
    for (i in 1:ncol(geno_prob)) {
      name <- colnames(geno_prob)[i]
      chr <- strsplit(name, "@", fixed = TRUE)[[1]][1]
      loc <- strsplit(name, "@", fixed = TRUE)[[1]][2]
      clgeno <- strsplit(name, "@", fixed = TRUE)[[1]][3]
      rownames(cross.geno.prob[[chr]]) <- make.names(rownames(cross.geno.prob[[chr]]))
      geno_prob[, i] <- cross.geno.prob[[chr]][Genotypes, loc, clgeno]
    }
    qtl.df <- fit.qtl$qtl.df
    select_names <- list()
    if (all(!is.na(qtl.df$linkage.group))) {
      if (verbose > 0) {
        print(paste(nrow(qtl.df), " QTLs found"))
      }
      for (i in 1:nrow(qtl.df)) {
        map <- (attr(cross.geno.prob[[qtl.df$linkage.group[i]]], 
                     "map")[1, ])
        loc <- names(map[map == qtl.df$position[i]])
        select_names[[i]] <- c(paste0(qtl.df$linkage.group[i], "@", loc, "@AC"), 
                               paste0(qtl.df$linkage.group[i], "@", loc, "@BC"),
                               paste0(qtl.df$linkage.group[i], "@", loc, "@AD"), 
                               paste0(qtl.df$linkage.group[i], "@", loc, "@BD"))
      }
      stopifnot(all(unlist(select_names) %in% colnames(geno_prob)))
      if (verbose > 1) {
        print(unlist(select_names))
      }
      X.loc <- geno_prob[, unlist(select_names)]
      fit.lm <- stats::lm(y ~ X.loc)
      coeff <- as.matrix(fit.lm$coefficients[-1])
      rownames(coeff) <- substr(rownames(coeff), start=6, stop=nchar(rownames(coeff)))
      if(!is.na(alpha[2])){
        summary.fit <- summary(fit.lm)$coefficient
        rownames(summary.fit)[-1] <- substr(rownames(summary.fit)[-1], start=6, 
                                            stop=nchar(rownames(summary.fit)[-1]))
        coeff <- as.matrix(summary.fit[,1][summary.fit[,4] < alpha[2]])
      }
      stopifnot(!is.null(rownames(coeff)))
      geno.effects <- data.frame(predictor = colnames(geno_prob), effect = 0)
      for (i in 1:length(coeff)) {
        geno.effects[geno.effects$predictor %in% rownames(coeff)[i], "effect"] <- ifelse(is.na(coeff[i]), 0, coeff[i])
      }
    }
    else {
      geno.effects <- data.frame(predictor = colnames(geno_prob), effect = 0)
      if (verbose > 0) {
        print("No QTL found")
      }
    }
  }
  else {
    geno.effects <- fit.qtl$allelic.effects
  }
  return(list(geno.effects = geno.effects, qtl.df = fit.qtl$qtl.df, 
              cross.geno.prob = cross.geno.prob))
}

##' Predict with qtl for caret
##'
##' @param modelFit model produced by \code{\link{caretFitQtl}}; should be a list with a component named \code{"geno.effects"} which contains a vector in the same order as the columns of newdata and, as values, the estimated allele effects (0 if not significant, non-zero otherwise) or genotypic class effects.
##' @param newdata predictor values of the instances being predicted (e.g. out-of-bag samples); should be in the JoinMap format
##' @param submodels optional list of tuning parameters only used with the "loop" element
##' @return vector
##' @author Charlotte Brault [aut], Timothee Flutre [ctb]
##' @export
caretPredictQtl <- function(modelFit, newdata, submodels = NULL){
  stopifnot(is.list(modelFit),
            all(c("geno.effects","cross.geno.prob") %in% names(modelFit)))

  newdata <- as.data.frame(t(newdata))
  Genotypes <- colnames(newdata) ; Genotypes <- make.names(Genotypes)

  if(!is.null(modelFit$cross.geno.prob)){
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Construct geno_prob matrix for newdata (VS)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cross.geno.prob <- modelFit$cross.geno.prob
    # For each chromosome, extract each locus name and paste the chromosome name
    list_locus <- 0
    for(chr in names(cross.geno.prob)){
      tmp <-  names(attr(cross.geno.prob[[chr]], "map")[1,])
      tmp <- paste0(chr, "@", tmp)
      list_locus <- append(list_locus, tmp)
    }
    list_locus <- list_locus[-1]

    # Add genotypic classes
    list_locus_clgeno <- 0
    for (i in 1:length(list_locus)){
      tmp <- c(paste0(list_locus[i], "@AC"),
               paste0(list_locus[i], "@AD"),
               paste0(list_locus[i], "@BC"),
               paste0(list_locus[i], "@BD"))
      list_locus_clgeno <- append(list_locus_clgeno, tmp)
    }
    list_locus_clgeno <- list_locus_clgeno[-1]
    stopifnot(length(list_locus)*4 ==  length(list_locus_clgeno))

    # Save names of genotypes to create the final matrix geno_prob
    geno_prob <- matrix(nrow=length(Genotypes), ncol=length(list_locus_clgeno),
                        dimnames=list(Genotypes, list_locus_clgeno))
    stopifnot(all(colnames(geno_prob) %in% modelFit$geno.effects$predictor))

    # Extract genotypes probabilities for genotypes in VS only
    for(i in 1:ncol(geno_prob)){
      name <- colnames(geno_prob)[i]
      chr <- strsplit(name, "@", fixed=TRUE)[[1]][1]
      loc <- strsplit(name, "@", fixed=TRUE)[[1]][2]
      clgeno <- strsplit(name, "@", fixed=TRUE)[[1]][3]
      rownames(cross.geno.prob[[chr]]) <- make.names(rownames(cross.geno.prob[[chr]]))
      geno_prob[,i] <- cross.geno.prob[[chr]][Genotypes, loc, clgeno]
    }
    predictions <- geno_prob %*% modelFit$geno.effects$effect
    ### Use allelic effects at marker position
  } else {
    pred <- loc.seg.phase.clas <- data.frame(locus = rownames(newdata),
                                             seg = getJoinMapSegregs(newdata),
                                             phase = NA, clas = NA,
                                             stringsAsFactors = FALSE)
    jm <- cbind(loc.seg.phase.clas, newdata)
    X <- joinMap2designMatrix(jm = jm, use.phase = FALSE, rm.col.zeros = FALSE,
                              verbose = 0)
    X <- X[, -1]
    allelic.effects <- as.matrix(modelFit$geno.effects$effect)
    rownames(allelic.effects) <- modelFit$geno.effects$predictor
    X <- subset(X, select = rownames(allelic.effects))
    stopifnot(dim(X)[2] == nrow(allelic.effects))
    predictions <- X %*% allelic.effects
  }
  return(predictions)[,1]
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
