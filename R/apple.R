##' Biennial bearing index
##'
##' Compute various types of the biennial bearing index (BBI) as defined in \href{https://dx.doi.org/10.1093/jxb/ert297}{Durand et al (2013)}.
##' @param dat data frame
##' @param coln.prod column name for the production variable; missing data will be discarded
##' @param coln.geno column name of the genotype
##' @param coln.rep column name of the tree replicate
##' @param coln.year column name of the year (the column data should be convertible to \code{numeric})
##' @param coln.epsilon column name of the residuals (used only for \code{type="res_norm"})
##' @param type type of BBI (classic/norm/res_norm)
##' @param verbose verbosity level (0/1)
##' @return numeric
##' @author Timothee Flutre
##' @export
BBI <- function(dat, coln.prod="prod", coln.geno="geno",
                coln.rep="tree", coln.year="year", coln.epsilon=NULL,
                type="classic",
                verbose=0){
  stopifnot(is.data.frame(dat),
            coln.prod %in% colnames(dat),
            coln.rep %in% colnames(dat),
            coln.geno %in% colnames(dat),
            coln.year %in% colnames(dat),
            type %in% c("classic", "norm", "res_norm"))
  if(type == "res_norm")
    stopifnot(! is.null(coln.epsilon))

  ## reformat the input data frame
  dat <- dat[! is.na(dat[[coln.prod]]),]
  dat <- droplevels(dat)
  dat[[coln.rep]] <- as.factor(dat[[coln.rep]])
  dat[[coln.year]] <- as.factor(dat[[coln.year]])
  dat[[coln.year]] <- as.numeric(levels(dat[[coln.year]]))[dat[[coln.year]]]
  dat <- dat[order(dat[[coln.year]], dat[[coln.geno]], dat[[coln.rep]]),]

  ## prepare the output vector
  dat[[coln.geno]] <- as.factor(dat[[coln.geno]])
  G <- nlevels(dat[[coln.geno]])
  out <- rep(NA, G)
  names(out) <- levels(dat[[coln.geno]])

  ## compute the index
  for(g in 1:G){
    geno <- levels(dat[[coln.geno]])[g]
    trees.g <- unique(dat[dat[[coln.geno]] == geno, coln.rep])
    R.g <- length(trees.g)
    if(R.g == 0)
      next

    num <- 0
    denom <- 0
    if(type == "classic"){
      for(r in 1:R.g){
        tree <- trees.g[r]
        y.g.r <- dat[dat[[coln.geno]] == geno &
                     dat[[coln.rep]] == tree,
                     coln.prod]
        T.gr <- length(y.g.r)
        denom <- denom + T.gr - 1
        for(t in 2:T.gr)
          if(y.g.r[t-1] != 0 || y.g.r[t] != 0)
            num <- num +
              abs(y.g.r[t] - y.g.r[t-1]) /
              (y.g.r[t-1] + y.g.r[t])
      }
      out[g] <- (2 * num) / denom

    } else if(type %in% c("norm", "res_norm")){
      num.num <- 0
      denom.num <- 0
      num.denom <- 0
      denom.denom <- 0
      for(r in 1:R.g){
        tree <- trees.g[r]
        y.g.r <- dat[dat[[coln.geno]] == geno &
                     dat[[coln.rep]] == tree,
                     coln.prod]
        epsilon.g.r <- NULL
        if(type == "res_norm")
          epsilon.g.r <- dat[dat[[coln.geno]] == geno &
                             dat[[coln.rep]] == tree,
                             coln.epsilon]
        T.gr <- length(y.g.r)
        num.denom <- num.denom + T.gr - 1
        denom.num <- denom.num + y.g.r[1]
        denom.denom <- denom.denom + T.gr
        for(t in 2:T.gr){
          denom.num <- denom.num + y.g.r[t]
          if(type == "norm"){
            num.num <- num.num +
              abs(y.g.r[t] - y.g.r[t-1])
          } else if(type == "res_norm")
            num.num <- num.num +
              abs(epsilon.g.r[t] - epsilon.g.r[t-1])
        }
      }
      out[g] <- (num.num / num.denom) /
        (denom.num / denom.denom)
    }

  }

  return(out)
}

##' Genotypic values for fruit bearing
##'
##' Compute the AR(1) coefficient for each genotype averaged over tree replicates following the model in \href{https://dx.doi.org/10.1093/jxb/ert297}{Durand et al (2013)}.
##' @param dat data frame
##' @param coln.epsilon column name of the residuals
##' @param coln.geno column name of the genotype
##' @param coln.rep column name of the tree replicate
##' @param coln.year column name of the year (the column data should be convertible to \code{numeric})
##' @return vector
##' @author Timothee Flutre
##' @export
genoAr1Coef <- function(dat, coln.epsilon="residual", coln.geno="geno",
                        coln.rep="tree", coln.year="year"){
  stopifnot(is.data.frame(dat),
            coln.epsilon %in% colnames(dat),
            coln.rep %in% colnames(dat),
            coln.geno %in% colnames(dat),
            coln.year %in% colnames(dat))

  dat[[coln.geno]] <- as.factor(dat[[coln.geno]])
  genos <- levels(dat[[coln.geno]])
  G <- nlevels(dat[[coln.geno]])

  fit.ar1.gr <-
    do.call(rbind, lapply(1:G, function(g){
      geno <- genos[g]
      trees.g <- unique(dat[dat[[coln.geno]] == geno, coln.rep])
      R.g <- length(trees.g)
      out <- c()
      for(r in 1:R.g){
        tmp <- droplevels(dat[dat[[coln.geno]] == geno &
                              dat[[coln.rep]] == as.character(r),])
        fit <- stats::arima(x=tmp$epsilonB, order=c(1,0,0),
                            method="ML") ## CSS can return an error
        out <- c(out, stats::coef(fit)["ar1"])
      }
      out
    }))
  rownames(fit.ar1.gr) <- genos

  fit.ar1.g <- apply(fit.ar1.gr, 1, mean)

  return(fit.ar1.g)
}
