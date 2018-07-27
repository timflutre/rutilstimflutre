## Contains functions useful for the "breeding game"
## https://github.com/timflutre/PlantSelBreedGame

##' Set up breeding game
##'
##' Set up the directories and names for the breeding game.
##' Already-existing directories are not re-created.
##' @param root.dir path to the root directory
##' @param shared.dir path to the shared directory (e.g. via Dropbox; root.dir if NULL)
##' @param create.admin if TRUE, a game master named "admin" will be created (by default, only a "tester" named "test" is created)
##' @param nb.breeders number of breeders; by default, the "admin" (which will have the "game master" status in "setup.Rmd") and "test" (which will have the "tester" status) breeders will be created, and other regular players can be created via the interface
##' @return list
##' @author Timothee Flutre
##' @export
setUpBreedingGame <- function(root.dir, shared.dir=NULL,
                              create.admin=FALSE,
                              nb.breeders=0){
  stopifnot(is.character(root.dir),
            dir.exists(root.dir))
  if(! is.null(shared.dir))
    stopifnot(is.character(shared.dir),
              dir.exists(shared.dir))

  out <- list(root.dir=root.dir)

  truth.dir <- paste0(root.dir, "/", "truth")
  if(! dir.exists(truth.dir))
    dir.create(truth.dir)
  out$truth.dir <- truth.dir

  if(is.null(shared.dir))
    shared.dir <- paste0(root.dir, "/", "shared")
  if(! dir.exists(shared.dir))
    dir.create(shared.dir)
  out$shared.dir <- shared.dir

  init.dir <- paste0(shared.dir, "/", "initial_data")
  if(! dir.exists(init.dir))
    dir.create(init.dir)
  out$init.dir <- init.dir

  breeders <- c("test")
  if(create.admin)
    breeders <- c(breeders, "admin")
  if(nb.breeders > 0)
    breeders <- c(breeders, paste0("breeder", 1:nb.breeders))
  breeder.dirs <- c()
  for(breeder in breeders){
    truth.breeder.dir <- paste0(truth.dir, "/", breeder)
    if(! dir.exists(truth.breeder.dir))
      dir.create(truth.breeder.dir)
    breeder.dirs[[breeder]] <- paste0(shared.dir, "/", breeder)
    if(! dir.exists(breeder.dirs[breeder]))
      dir.create(breeder.dirs[breeder])
  }
  out$breeders <- breeders
  out$breeder.dirs <- breeder.dirs

  dbname <- paste0(root.dir, "/breeding-game.sqlite")
  out$dbname <- dbname

  return(out)
}

##' Get the breeding game setup
##'
##' Retrieve the paths to the directories used for the breeding game.
##' @param root.dir path to the root directory
##' @return list
##' @author Timothee Flutre
##' @export
getBreedingGameSetup <- function(root.dir){
  stopifnot(is.character(root.dir),
            length(root.dir) == 1,
            dir.exists(root.dir))

  out <- list(root.dir=root.dir)

  out$truth.dir <- paste0(root.dir, "/truth")
  out$shared.dir <- paste0(root.dir, "/shared")
  out$init.dir <- paste0(out$shared.dir, "/initial_data")
  tmp <- basename(Sys.glob(paste0(out$shared.dir, "/*")))
  for(x in tmp){
    if(x != "initial_data")
      out$breeders <- c(out$breeders, x)
  }

  out$breeder.dirs <- c()
  for(breeder in out$breeders)
    out$breeder.dirs[[breeder]] <- paste0(out$shared.dir, "/", breeder)

  out$dbname <- paste0(root.dir, "/breeding-game.sqlite")

  return(out)
}

##' Get the breeding game constants
##'
##' Retrieve the constants used to parametrized the breeding game from the SQLite database.
##' @param dbname name of the SQLite database (full path)
##' @param table name of the table
##' @return list
##' @author Timothee Flutre
##' @export
getBreedingGameConstants <- function(dbname, table="constants"){
  requireNamespace("DBI")
  requireNamespace("RSQLite")
  stopifnot(file.exists(dbname),
            is.character(table))

  out.list <- list()

  ## retrieve the content of the table
  db <- DBI::dbConnect(RSQLite::SQLite(), dbname=dbname)
  query <- paste0("SELECT *",
                  " FROM ", table)
  out.df <- DBI::dbGetQuery(db, query)
  DBI::dbDisconnect(db)

  ## reformat
  out.list <- as.list(out.df$value)
  names(out.list) <- out.df$item
  for(i in seq_along(out.list))
    out.list[[i]] <- tryCatch({
      as.numeric(out.list[[i]])
    }, warning = function(c){ # ex.: case of 'max.upload.pheno.field'
      out.list[[i]]
    })

  return(out.list)
}

##' Simul breeding game
##'
##' Make the structure of the data.frame that will be given to the players at the beginning of the game.
##' @param nb.lines.per.year number of lines phenotyped each year
##' @param nb.years number of years of phenotyping
##' @param nb.plots.per.line.per.year number of plots per line per year
##' @param first.year numeric of the year of the first phenotyping
##' @param line.ids vector of line identifiers (should be sorted)
##' @return data.frame
##' @author Timothee Flutre
##' @seealso \code{\link{makeDfPhenos}}
##' @export
makeDfInitPhenos <- function(nb.lines.per.year=150, nb.years=10,
                             nb.plots.per.line.per.year=2,
                             first.year=2005, line.ids){
  stopifnot(is.numeric(nb.lines.per.year),
            is.numeric(nb.years),
            is.numeric(nb.plots.per.line.per.year),
            is.numeric(first.year),
            is.character(line.ids),
            ! is.unsorted(order(line.ids)))

  nb.plots <- nb.lines.per.year * nb.plots.per.line.per.year
  latest.year <- first.year + nb.years - 1
  nb.new.lines.per.year <- nb.lines.per.year / nb.plots.per.line.per.year
  stopifnot(length(line.ids) == nb.lines.per.year + nb.new.lines.per.year *
            (nb.years - 1))

  df <- data.frame(ind=NA,
                   year=as.factor(rep(first.year:latest.year, each=nb.plots)),
                   plot=as.factor(rep(1:nb.plots, times=nb.years)),
                   pathogen=NA,
                   trait1.raw=NA,
                   trait1=NA,
                   trait2=NA,
                   trait3=NA)

  ## fill the identifiers column
  df$ind[df$year == levels(df$year)[1]] <-
    rep(line.ids[1:nb.new.lines.per.year], each=nb.plots.per.line.per.year)
  for(j in 1:nb.years){
    year <- levels(df$year)[j]
    idx <- (1 + (j-1) * nb.new.lines.per.year) :
      ((j-1) * nb.new.lines.per.year + nb.lines.per.year)
    df$ind[df$year == year] <- rep(line.ids[idx],
                                   each=nb.plots.per.line.per.year)
  }
  df$ind <- as.factor(df$ind)

  return(df)
}

##' Simul breeding game
##'
##' Make the structure of the data.frame that will be given to the players when they request phenotyping during the game.
##' @param ind.ids vector of genotype identifiers (will be sorted using \code{\link{sort}}))
##' @param nb.plots.per.ind vector with the number of plots at which each genotype should be phenotype (will be sorted in the same order as \code{ind.ids})
##' @param year numeric of the year at which phenotyping occurs
##' @param pathogen if TRUE, the pathogen will be present the given year
##' @return data.frame
##' @author Timothee Flutre
##' @seealso \code{\link{makeDfInitPhenos}}
##' @export
makeDfPhenos <- function(ind.ids, nb.plots.per.ind, year, pathogen){
  stopifnot(is.character(ind.ids),
            is.numeric(nb.plots.per.ind),
            length(nb.plots.per.ind) == length(ind.ids),
            is.numeric(year),
            is.logical(pathogen))

  tmp <- data.frame(ind.id=ind.ids, nb.plots=nb.plots.per.ind,
                    stringsAsFactors=FALSE)
  tmp <- tmp[order(tmp$ind.id),]

  df <- data.frame(ind=rep(tmp$ind.id, tmp$nb.plots),
                   year=as.factor(rep(year, sum(tmp$nb.plots))),
                   plot=as.factor(1:sum(tmp$nb.plots)),
                   pathogen=pathogen,
                   trait1.raw=NA,
                   trait1=NA,
                   trait2=NA,
                   trait3=NA)

  return(df)
}

##' Simul breeding game
##'
##' Simulate the additive SNP effects of traits 1 and 2 jointly, with some degree of pleiotropy.
##' @param snp.ids vector of SNP identifiers for which effects should be simulated
##' @param sigma.beta2 vector of variances of additive SNP effects for both traits
##' @param prop.pleio proportion of SNPs having an effect on both traits
##' @param cor.pleio genotypic correlation at these SNPs
##' @param verbose verbosity level (0/1)
##' @return list with a matrix of SNP effects and a vector identifying which SNPs are pleiotropic
##' @author Timothee Flutre
##' @examples
##' \dontrun{set.seed(1859)
##' X <- simulGenosDose(nb.genos=825, nb.snps=2000)
##' snp.effects <- simulSnpEffectsTraits12(snp.ids=colnames(X))
##' sum(snp.effects$is.pleiotropic)
##' regplot(snp.effects$Beta[,1], snp.effects$Beta[,2])
##' }
##' @export
simulSnpEffectsTraits12 <- function(snp.ids,
                                    sigma.beta2=c(trait1=0.017,
                                                  trait2=10^(-4)),
                                    prop.pleio=0.4,
                                    cor.pleio=-0.7,
                                    verbose=1){
  requireNamespace("MASS")
  stopifnot(is.character(snp.ids),
            length(snp.ids) > 0,
            ! anyDuplicated(snp.ids),
            is.vector(sigma.beta2),
            is.numeric(sigma.beta2),
            length(sigma.beta2) == 2,
            ! is.null(names(sigma.beta2)),
            all(sigma.beta2 >= 0),
            is.numeric(prop.pleio),
            length(prop.pleio) == 1,
            all(prop.pleio >= 0, prop.pleio <= 1),
            is.numeric(cor.pleio),
            length(cor.pleio) == 1,
            all(cor.pleio >= -1, cor.pleio <= 1))

  nb.snps <- length(snp.ids)
  traits <- names(sigma.beta2)

  if(verbose > 0){
    msg <- "first, simulate without pleiotropy..."
    write(msg, stdout())
  }
  Sigma.beta.nopleio <- matrix(c(sigma.beta2[1], 0, 0, sigma.beta2[2]),
                               nrow=2, ncol=2,
                               dimnames=list(traits, traits))
  Beta <- MASS::mvrnorm(n=nb.snps, mu=c(0,0), Sigma=Sigma.beta.nopleio)
  rownames(Beta) <- snp.ids
  colnames(Beta) <- traits
  if(verbose > 0){
    msg <- paste0("var(Beta[,1]) = ",
                  format(stats::var(Beta[,1]), digits=2),
                  "\nvar(Beta[,2]) = ",
                  format(stats::var(Beta[,2]), digits=2),
                  "\ncor(Beta[,1], Beta[,2]) = ",
                  format(stats::cor(Beta[,1], Beta[,2]), digits=3))
    write(msg, stdout())
  }

  if(verbose > 0){
    msg <- "second, add pleiotropy..."
    write(msg, stdout())
  }
  cov.pleio <- cor.pleio * sqrt(sigma.beta2[1] * sigma.beta2[2])
  Sigma.beta.pleio <- matrix(c(sigma.beta2[1], cov.pleio, cov.pleio,
                               sigma.beta2[2]),
                             nrow=2, ncol=2,
                             dimnames=list(traits, traits))
  length(idx.pleio <- sample.int(n=nb.snps,
                                 size=floor(prop.pleio * nb.snps)))
  Beta[idx.pleio,] <- MASS::mvrnorm(n=length(idx.pleio), mu=c(0,0),
                                    Sigma=Sigma.beta.pleio)
  is.pleiotropic <- stats::setNames(rep(FALSE, nb.snps), snp.ids)
  is.pleiotropic[idx.pleio] <- TRUE

  if(verbose > 0){
    msg <- paste0("var(Beta[,1]) = ",
                  format(stats::var(Beta[,1]), digits=2),
                  "\nvar(Beta[,2]) = ",
                  format(stats::var(Beta[,2]), digits=2),
                  "\ncor(Beta[idx.pleio,1], Beta[idx.pleio,2]) = ",
                  format(stats::cor(Beta[idx.pleio,1], Beta[idx.pleio,2]),
                         digits=3))
    msg <- paste0(msg, "\ncor(Beta[,1], Beta[,2]) = ",
                  format(stats::cor(Beta[,1], Beta[,2]), digits=3))
    write(msg, stdout())
  }

  return(list(sigma.beta2=sigma.beta2, prop.pleio=prop.pleio,
              cor.pleio=cor.pleio, is.pleiotropic=is.pleiotropic, Beta=Beta))
}

##' Simul breeding game
##'
##' Simulate phenotypes of traits 1 and 2 jointly: Y = Z.J * Alpha + Z.I * G.A + E.
##' @param dat data.frame specifying the structure of the initial data set of phenotypes, with columns "ind", "year", etc
##' @param mu vector of global means, for each trait
##' @param sigma.alpha2 vector of variance for the "year" effects, for each trait (ignored if \code{Alpha} is not NULL)
##' @param Alpha matrix of "year" effects, for each trait, the years corresponding to those indicated in \code{dat} (if NULL, will be simulated using \code{sigma.alpha2})
##' @param X matrix of bi-allelic SNP genotypes encoded in allele dose in {0,1,2}, with individuals in rows in the same order as the levels of \code{dat$ind}
##' @param afs vector of allele frequencies (to center X so that the genotypic values are also centered) which names are the SNPs
##' @param Beta matrix of additive SNP effects, for each trait
##' @param h2 vector of heritabilities, with the name of each trait, for instance \code{c(trait1=0.3, trait2=0.4)} (if NULL, \code{sigma2} will be used, but do not specify both)
##' @param sigma2 vector of error variances, with the name of each trait, for instance \code{c(trait1=0.467, trait2=0.210)} (if NULL, \code{h2} will be used, but do not specify both)
##' @param cor.E.inter.trait correlation of errors between both traits
##' @param set.neg.to.zero if TRUE, the negative phenotypic values, if any, will be set to zero
##' @param verbose verbosity level (0/1)
##' @return list
##' @author Timothee Flutre
##' @seealso \code{\link{simulSnpEffectsTraits12}}, \code{\link{makeDfInitPhenos}}
##' @examples
##' \dontrun{set.seed(1859)
##' X <- simulGenosDose(nb.genos=825, nb.snps=2000)
##' snp.effects12 <- simulSnpEffectsTraits12(snp.ids=colnames(X))
##' dat <- makeDfInitPhenos(line.ids=rownames(X))
##' phenos <- simulTraits12(dat, X=X, Beta=snp.effects12$Beta)
##' regplot(phenos$G.A[,1], phenos$G.A[,2])
##' }
##' @export
simulTraits12 <- function(dat,
                          mu=c(trait1=100, trait2=15),
                          sigma.alpha2=c(trait1=230, trait2=10),
                          Alpha=NULL,
                          X,
                          afs,
                          Beta,
                          h2=NULL,
                          sigma2=c(trait1=330, trait2=0.6),
                          cor.E.inter.trait=0,
                          set.neg.to.zero=TRUE,
                          verbose=1){
  requireNamespace("MASS")
  if(! is.matrix(X)){
    if(is.vector(X)){
      warning("X should be a matrix; you should use drop=FALSE")
      X <- matrix(data=X, nrow=1)
    }
  }
  stopIfNotValidGenosDose(X)
  stopifnot(is.data.frame(dat),
            ncol(dat) >= 4,
            ! is.null(colnames(dat)),
            is.vector(mu),
            is.numeric(mu),
            all(mu >= 0),
            ! is.null(names(mu)),
            length(mu) == 2,
            all(c("ind","year","plot",names(mu)) %in% colnames(dat)),
            is.factor(dat$year),
            is.factor(dat$ind),
            ! is.unsorted(order(dat$year)),
            any(! is.null(sigma.alpha2), ! is.null(Alpha)),
            rownames(X) == levels(dat$ind),
            is.numeric(afs),
            length(afs) == ncol(X),
            ! is.null(names(afs)),
            all(names(afs) == colnames(X)),
            is.matrix(Beta),
            ncol(Beta) == length(mu),
            ! is.null(colnames(Beta)),
            nrow(Beta) == ncol(X),
            all(rownames(Beta) == colnames(X)),
            xor(is.null(h2), is.null(sigma2)),
            xor(! is.null(h2), ! is.null(sigma2)),
            is.numeric(cor.E.inter.trait),
            all(cor.E.inter.trait >= -1, cor.E.inter.trait <= 1))
  if(is.null(Alpha)){
    stopifnot(is.vector(sigma.alpha2),
              is.numeric(sigma.alpha2),
              all(sigma.alpha2 >= 0),
              length(sigma.alpha2) == length(mu),
              ! is.null(names(sigma.alpha2)),
              all(names(sigma.alpha2) == names(mu)))
  } else
    stopifnot(is.matrix(Alpha),
              ncol(Alpha) == length(mu),
              nrow(Alpha) == nlevels(dat$year),
              ! is.null(colnames(Alpha)),
              colnames(Alpha) == names(mu),
              ! is.null(rownames(Alpha)),
              rownames(Alpha) == levels(dat$year))
  if(is.null(h2)){
    stopifnot(! is.null(sigma2),
              is.vector(sigma2),
              is.numeric(sigma2),
              all(sigma2 >=0),
              length(sigma2) == length(mu),
              ! is.null(names(sigma2)),
              names(sigma2) == names(mu))
  } else
    stopifnot(is.null(sigma2),
              is.vector(h2),
              is.numeric(h2),
              all(h2 >= 0, h2 <= 1),
              length(h2) == length(mu),
              ! is.null(names(h2)),
              names(h2) == names(mu))

  out <- list()

  ## get dimensions
  N <- nrow(dat)          # number of phenotypic observations
  I <- nlevels(dat$ind)   # number of individuals (genotypes)
  J <- nlevels(dat$year)  # number of years
  K <- nlevels(dat$plot)  # number of plots
  T <- length(mu)         # number of traits
  inds <- levels(dat$ind)
  years <- levels(dat$year)
  traits <- names(mu)

  ## create phenotypic matrix
  Y <- matrix(data=NA, nrow=N, ncol=T)
  colnames(Y) <- traits

  ## make the predictors
  ## global intercept
  Z.mu <- matrix(1, nrow=N, ncol=1)

  ## "year" effects
  if(J == 1){
    Z.J <- matrix(data=1, nrow=N, ncol=1)
  } else
    Z.J <- stats::model.matrix(~ year - 1, data=dat)
  colnames(Z.J) <- years
  if(is.null(Alpha)){
    Alpha <- matrix(c(stats::rnorm(n=J, mean=0, sd=sqrt(sigma.alpha2[1])),
                      stats::rnorm(n=J, mean=0, sd=sqrt(sigma.alpha2[2]))),
                    nrow=J, ncol=T,
                    dimnames=list(years, traits))
  }
  if(verbose > 0){
    msg <- paste0("'year' effects: ")
    write(msg, stdout())
    print(Alpha)
  }

  ## genotypic values
  if(nlevels(dat$ind) == 1){
    Z.I <- matrix(data=1, nrow=nrow(dat), ncol=1)
  } else
    Z.I <- stats::model.matrix(~ ind - 1, data=dat)
  colnames(Z.I) <- inds
  ## center X as in Vitezica et al (2013)
  tmp <- matrix(rep(1, I)) %*% (2 * afs)
  X.center <- X - tmp
  G.A <- X.center %*% Beta
  if(verbose > 0){
    msg <- paste0("cor(G.A[,1], G.A[,2]) = ",
                  format(stats::cor(G.A[,1], G.A[,2]), digits=3))
    msg <- paste0(msg, "\nmean(G.A[,1]) = ", format(mean(G.A[,1]), digits=3))
    msg <- paste0(msg, "\nmean(G.A[,2]) = ", format(mean(G.A[,2]), digits=3))
    msg <- paste0(msg, "\nvar(G.A[,1]) = ", format(stats::var(G.A[,1]),
                                                   digits=3))
    msg <- paste0(msg, "\nvar(G.A[,2]) = ", format(stats::var(G.A[,2]),
                                                   digits=3))
    write(msg, stdout())
  }

  M <- Z.mu %*% t(as.matrix(mu)) + Z.J %*% Alpha + Z.I %*% G.A

  ## make the errors
  if(is.null(sigma2)){
    sigma2 <- c(((1 - h2[1]) / h2[1]) * stats::var(G.A[,1]),
                ((1 - h2[2]) / h2[2]) * stats::var(G.A[,2]))
    names(sigma2) <- traits
    if(verbose > 0){
      msg <- paste0("sigma2[1] = ", format(sigma2[1], digits=3))
      msg <- paste0(msg, "\nsigma2[2] = ", format(sigma2[2], digits=3))
      write(msg, stdout())
    }
  }
  cov.E.inter.trait <- cor.E.inter.trait * sqrt(sigma2[1] * sigma2[2])
  Sigma <- matrix(c(sigma2[1], cov.E.inter.trait,
                    cov.E.inter.trait, sigma2[2]),
                  nrow=2, ncol=2,
                  dimnames=list(traits, traits))
  E <- MASS::mvrnorm(n=N, mu=c(0,0), Sigma=Sigma)

  ## make the phenotypes
  Y <- M + E
  if(set.neg.to.zero){
    for(t in 1:T){
      is.neg <- Y[,t] < 0
      if(any(is.neg))
        Y[is.neg, t] <- 0
    }
  }

  ## fill the output
  out$N <- N
  out$I <- I
  out$J <- J
  out$mu <- mu
  out$Z.J <- Z.J
  out$sigma.alpha2 <- sigma.alpha2
  out$Alpha <- Alpha
  out$Z.I <- Z.I
  out$Beta <- Beta
  out$G.A <- G.A
  out$h2 <- h2
  out$sigma2 <- sigma2
  out$cor.E.inter.trait <- cor.E.inter.trait
  out$Y <- Y

  return(out)
}

##' Simul breeding game
##'
##' Simulate phenotypes of trait 3.
##' @param dat data.frame specifying the structure of the initial data set of phenotypes
##' @param X matrix of bi-allelic SNP genotypes encoded in allele dose in {0,1,2}
##' @param afs vector of allele frequencies (ignored if \code{qtn.id} is not NULL)
##' @param subset.snps vector of SNP identifiers to which the QTN should belong (ignored if \code{qtn.id} is not NULL)
##' @param qtn.id identifier of the causal SNP (QTN)
##' @param resist.genos vector of genotypes in allele dose at the causal SNP providing the resistance (ignored if \code{qtn.id} is NULL)
##' @param prob.resist.no.qtl probability for a genotype without the QTL to be resistant a year with the pathogen
##' @param verbose verbosity level (0/1)
##' @return list
##' @author Timothee Flutre
##' @seealso \code{\link{makeDfInitPhenos}}
##' @export
simulTrait3 <- function(dat, X, afs=NULL, subset.snps=NULL,
                        qtn.id=NULL, resist.genos=NULL,
                        prob.resist.no.qtl=0.02, verbose=1){
  stopifnot(is.data.frame(dat),
            ncol(dat) >= 4,
            ! is.null(colnames(dat)),
            c("ind", "pathogen") %in% colnames(dat))
  stopIfNotValidGenosDose(X)
  stopifnot(is.numeric(prob.resist.no.qtl),
            all(prob.resist.no.qtl >= 0, prob.resist.no.qtl <= 1))
  if(is.null(qtn.id)){
    stopifnot(! is.null(afs),
              ! is.null(subset.snps))
  } else
    stopifnot(! is.null(resist.genos))

  out <- list()

  ## choose the causal SNP (QTN)
  if(is.null(qtn.id)){
    mafs <- estimSnpMaf(afs=afs)
    qtn.id <- names(sample(x=which(mafs >= 0.14 & mafs <= 0.16 &
                                   names(mafs) %in% subset.snps),
                           size=1))
    if(verbose > 0){
      msg <- paste0("QTN of trait3: ", qtn.id)
      write(msg, stdout())
      print(table(X[, qtn.id]))
    }
    is.minor.counted <- afs[qtn.id] < 0.5
    if(is.minor.counted){
      resist.genos <- c(1, 2)
    } else{
      resist.genos <- c(0, 1)
    }
  }

  ## identify resistant individuals based on their genotypes at the QTN
  inds.resist <- rownames(X)[X[, qtn.id] %in% resist.genos]

  ## make the phenotypes
  ## 1. no pathogen -> no symptom
  ## 2. pathogen + QTL -> no symptom
  ## 3. pathogen + no QTL -> symptoms (but not always)
  y <- rep(NA, nrow(dat))
  y[! dat$pathogen] <- 0
  idx <- which(dat$pathogen & (dat$ind %in% inds.resist))
  y[idx] <- 0
  idx <- which(dat$pathogen & (! dat$ind %in% inds.resist))
  y[idx] <- stats::rbinom(n=length(idx), size=1, prob=1 - prob.resist.no.qtl)

  ## fill the output
  out$prob.resist.no.qtl <- prob.resist.no.qtl
  out$qtn.id <- qtn.id
  out$resist.genos <- resist.genos
  out$y <- y

  return(out)
}

##' Example for breeding game
##'
##' Make a file with examples of plant material to request (autofecundation, allofecundation, haplodiploidization)
##' @param out.dir path to the directory in which the file will be saved
##' @return invisible data.frame
##' @author Timothee Flutre
##' @export
makeExamplePlantFile <- function(out.dir){
  stopifnot(dir.exists(out.dir))

  plants <-
    data.frame(parent1=c("Coll0037", "37-8.1", "37-8.1"),
               parent2=c("Coll0008", "37-8.1", NA),
               child=c("37-8.1", "37-8.1.1", "37-8.1.HD1"),
               explanations=c("allofecundation", "autofecundation",
                              "haplodiploidization"),
               stringsAsFactors=FALSE)

  f <- paste0(out.dir, "/example_request_plant_material.txt")

  cat("# the table (below) contains the requested plant material\n",
      file=f, append=FALSE)
  cat("# only columns 'parent1', 'parent2' and 'child' are compulsory\n",
      file=f, append=TRUE)
  cat("# each row corresponds to a child\n",
      file=f, append=TRUE)
  cat("# the 'child' column shouldn't have any duplicate\n",
      file=f, append=TRUE)
  cat("# only the 'parent2' column can be empty\n",
      file=f, append=TRUE)
  cat("# individual names should only use [a-z], [A-Z], [0-9], [_-] (no space, comma, etc)\n",
      file=f, append=TRUE)
  cat("# use write.table(x=plants, file=\"<breeder>_<year>_plant_material.txt\", quote=FALSE, sep=\"\\t\", na=\"\", row.names=FALSE, col.names=TRUE)\n",
      file=f, append=TRUE)
  cat("# lines starting with '#' will be ignored\n",
      file=f, append=TRUE)
  suppressWarnings(utils::write.table(x=plants, file=f, append=TRUE,
                                      quote=FALSE, sep="\t", na="",
                                      row.names=FALSE, col.names=TRUE))

  invisible(plants)
}

##' Example for breeding game
##'
##' Make a file with examples of data to request (phenotypes, genotypes).
##' @param out.dir path to the directory in which the file will be saved
##' @return invisible data.frame
##' @author Timothee Flutre
##' @export
makeExampleDataFile <- function(out.dir){
  stopifnot(dir.exists(out.dir))

  dat <- data.frame(ind=c(paste0("ind", 7:10), "ind7", "ind31"),
                    task=c(rep("pheno-field", 2), "pheno-patho", "geno", "geno", "geno"),
                    details=c("2", "5", "1", "hd", "snp15492", "ld"),
                    stringsAsFactors=FALSE)

  f <- paste0(out.dir, "/example_request_data.txt")

  cat("# the list (below) contains individuals to genotype and/or phenotype\n",
      file=f, append=FALSE)
  cat("# columns 'ind', 'task' and 'details' are compulsory\n",
      file=f, append=TRUE)
  cat("# individual names should only use [a-z], [A-Z], [0-9], [_-] (no space, comma, etc)\n",
      file=f, append=TRUE)
  cat("# the 'task' column should contain 'pheno-field', 'pheno-patho' or 'geno'\n",
      file=f, append=TRUE)
  cat("# if 'task=pheno-field', the 'details' column should contain the number of plots\n",
      file=f, append=TRUE)
  cat("# if 'task=pheno-patho', the 'details' column should contain the number of replicates\n",
      file=f, append=TRUE)
  cat("# if 'task=geno', the 'details' column should contain 'hd', 'ld' or the SNP identifier\n",
      file=f, append=TRUE)
  cat("# individuals should not be duplicated within each task\n",
      file=f, append=TRUE)
  cat("# the total number of requested plots (task=pheno-field) should not exceed the total available\n",
      file=f, append=TRUE)
  cat("# use write.table(x=inds, file=\"<breeder>_inds.txt\", quote=FALSE, sep=\"\\t\", row.names=FALSE, col.names=TRUE)\n",
      file=f, append=TRUE)
  cat("# lines starting with '#' will be ignored\n",
      file=f, append=TRUE)
  suppressWarnings(utils::write.table(x=dat, file=f, append=TRUE, quote=FALSE,
                               sep="\t", row.names=FALSE, col.names=TRUE))

  invisible(dat)
}
