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
    trees.g <- as.character(unique(dat[dat[[coln.geno]] == geno, coln.rep]))
    R.g <- length(trees.g)
    if(verbose > 0){
      txt <- paste0(geno, ": ", R.g)
      write(txt, stdout())
    }
    if(R.g == 0)
      next

    if(type == "classic"){
      num <- 0
      denom <- 0
      for(r in 1:R.g){
        tree <- trees.g[r]
        y.g.r <- dat[dat[[coln.geno]] == geno &
                     dat[[coln.rep]] == tree,
                     coln.prod]
        T.gr <- length(y.g.r)
        denom <- denom + T.gr - 1
        if(T.gr > 1){
          for(t in 2:T.gr)
            if(y.g.r[t-1] != 0 || y.g.r[t] != 0)
              num <- num +
                abs(y.g.r[t] - y.g.r[t-1]) /
                (y.g.r[t-1] + y.g.r[t])
        }
      }
      if(denom != 0)
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
        if(T.gr > 1){
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
      }
      if(all(denom.num != 0, denom.denom != 0))
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
##' @param method fitting method passed to \code{arima} (ML/CSS/ML|CSS); if \code{"ML|CSS"}, first "ML" then "CSS" will be tried
##' @param transform.pars if true, the AR parameters are transformed to ensure that they remain in the region of stationarity
##' @param optim.control list of control parameters for \code{optim}
##' @param verbose verbosity level (0/1)
##' @return vector
##' @author Timothee Flutre
##' @export
genoAr1Coef <- function(dat, coln.epsilon="residual", coln.geno="geno",
                        coln.rep="tree", coln.year="year",
                        method="ML|CSS", transform.pars=TRUE,
                        optim.control=list(),
                        verbose=0){
  stopifnot(is.data.frame(dat),
            coln.epsilon %in% colnames(dat),
            coln.rep %in% colnames(dat),
            coln.geno %in% colnames(dat),
            coln.year %in% colnames(dat),
            method %in% c("ML", "CSS", "ML|CSS", "CSS|ML"))

  dat <- droplevels(dat)
  dat[[coln.geno]] <- as.factor(dat[[coln.geno]])
  genos <- levels(dat[[coln.geno]])
  G <- nlevels(dat[[coln.geno]])

  fit.ar1.gr <-
    lapply(1:G, function(g){
      geno <- genos[g]
      trees.g <- as.character(unique(dat[dat[[coln.geno]] == geno, coln.rep]))
      R.g <- length(trees.g)
      out <- rep(NA, R.g)
      for(r in 1:R.g){
        if(verbose > 0){
          txt <- paste0(geno, ": ", r, "/", R.g)
          write(txt, stdout())
        }
        tmp <- droplevels(dat[dat[[coln.geno]] == geno &
                              dat[[coln.rep]] == trees.g[r],])
        if(nrow(tmp) > 1){
          for(m in strsplit(method, "\\|")[[1]]){
            fit <- try(stats::arima(x=tmp[[coln.epsilon]],
                                    order=c(1,0,0),
                                    method=m,
                                    transform.pars=transform.pars,
                                    optim.control=optim.control),
                       silent=ifelse(verbose <= 0, TRUE, FALSE))
            if(! methods::is(fit, "try-error"))
              break
          }
          if(methods::is(fit, "Arima"))
            out[r] <- stats::coef(fit)["ar1"]
        }
      }
      out
    })
  names(fit.ar1.gr) <- genos

  fit.ar1.g <- sapply(fit.ar1.gr, mean)

  return(fit.ar1.g)
}
