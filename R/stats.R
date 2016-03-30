## Contains functions useful for statistics.

##' Summary
##'
##' Print the output of \code{\link{summary}} in a single line.
##' @param x vector of numbers
##' @param spec specifier, see \code{\link{sprintf}}
##' @return nothing
##' @author Timothee Flutre
##' @export
prettyPrintSummary <- function(x, spec="%.2f"){
  stopifnot(is.vector(x))
  tmp <- summary(x)
  fmt <- paste0("min=", spec, " q1=", spec, " med=", spec, " mean=", spec,
                " q3=", spec, " max=", spec)
  txt <- sprintf(fmt, tmp[1], tmp[2], tmp[3], tmp[4], tmp[5], tmp[6])
  print(txt)
}

##' Return the Root Mean Squared Error
##'
##'
##' @title Root Mean Squared Error
##' @param error vector \eqn{\hat{\theta}_i - \theta_i}
##' @return numeric
##' @author Timothee Flutre
##' @export
rmse <- function(error){
  sqrt(mean(error^2))
}

##' Return the Mean Absolute Error
##'
##'
##' @title Mean Absolute Error
##' @param error vector \eqn{\hat{\theta}_i - \theta_i}
##' @return numeric
##' @author Timothee Flutre
##' @export
mae <- function(error){
  mean(abs(error))
}

##' Return the Mean Signed Difference
##'
##'
##' @title Mean Signed Difference
##' @param error vector \eqn{\hat{\theta}_i - \theta_i}
##' @return numeric
##' @author Timothee Flutre
##' @export
msd <- function(error){
  mean(error)
}

##' Return the number of true positives, false positives, true negatives,
##' false negatives, true positive proportion (sensitivity), false positive
##' proportion, accuracy, true negative proportion (specificity), false
##' discovery proportion, false negative proportion and positive predictive
##' value (precision)
##'
##' Both input vectors should be sorted beforehand
##' @title Binary classification
##' @param known.nulls vector of booleans (TRUE if the null is true)
##' @param called.nulls vector of booleans (TRUE if the null is accepted)
##' @return vector with names
##' @author Timothee Flutre
##' @export
binaryClassif <- function(known.nulls, called.nulls){
  ## http://en.wikipedia.org/wiki/Sensitivity_and_specificity
  ##
  ##                                  CALLED
  ##                     Accepted null     Rejected null
  ##
  ##       true null         TN (U)            FP (V)          n0
  ## TRUTH
  ##       false null        FN (T)            TP (S)          n1
  ##
  ##                         a                 r               n
  stopifnot(is.vector(known.nulls), is.vector(called.nulls),
            length(known.nulls) == length(called.nulls),
            sum(! is.logical(known.nulls)) == 0,
            sum(! is.logical(called.nulls)) == 0)

  n <- length(known.nulls) # total number of tests
  n0 <- sum(known.nulls)   # nb of true nulls
  n1 <- n - n0             # nb of "false nulls" (i.e. "true alternatives")
  a <- sum(called.nulls)   # nb of accepted nulls ("called not significant")
  r <- n - a               # nb of rejected nulls ("called significant", "discoveries")

  ## true positive = reject a false null
  tp <- sum(which(! called.nulls) %in% which(! known.nulls))

  ## false positive = reject a true null (type I error, "false alarm")
  fp <- sum(which(! called.nulls) %in% which(known.nulls))

  ## true negatives = accept a true null
  tn <- sum(which(called.nulls) %in% which(known.nulls))

  ## false negatives = accept a false null (type II error, "miss")
  fn <- sum(which(called.nulls) %in% which(! known.nulls))

  tpp <- tp / n1        # true positive prop (sensitivity)
  fpp <- fp / n0        # false positive prop
  acc <- (tp + tn) / n  # accuracy
  tnp <- tn / n0        # true negative prop (specificity), = 1 - fpp
  fdp <- fp / r         # false discovery prop
  fnp <- fn / a         # false negative prop
  ppv <- tp / r         # positive predictive value (precision)

  return(c(n=n, n0=n0, n1=n1, a=a, r=r,
           tp=tp, fp=fp, tn=tn, fn=fn,
           tpp=tpp, fpp=fpp, acc=acc, tnp=tnp, fdp=fdp, fnp=fnp, ppv=ppv))
}

##' Stable computation of \eqn{log_{10}(\sum_i w_i 10^x_i)}
##'
##' Use equal weights if not specified.
##' @title Log of weighted sum
##' @param x vector
##' @param weights weights
##' @return numeric
##' @author Timothee Flutre
##' @export
log10WeightedSum <- function(x, weights=NULL){
  stopifnot(is.numeric(x), is.vector(x))
  if(! is.null(weights)){
    stopifnot(is.vector(weights),
              length(weights) == length(x))
  } else
    weights <- rep(1/length(x), length(x))
  m <- max(x)
  m + log10(sum(weights * 10^(x - m)))
}

##' Return the Moore-Penrose pseudo-inverse of a matrix
##'
##' Golub & Van Loan, Matrix Computations, 3rd edition, ch5, p257
##' @title Pseudo-inverse
##' @param mat matrix
##' @return matrix
##' @author Timothee Flutre
##' @export
mpInv <- function(mat){
  mat.svd <- svd(mat)
  mat.svd$v %*% diag(1/mat.svd$d) %*% t(mat.svd$u)
}

##' Principal component analysis
##'
##' Via the singular value decomposition (SVD): X = U D V^T.
##' This is mostly for teaching purposes, see \code{\link[stats]{prcomp}} otherwise.
##' @param X data matrix with N rows and P columns
##' @param ct use TRUE to center the columns of X (recommended), FALSE otherwise
##' @param sc use TRUE to scale the columns of X (if different units), FALSE otherwise
##' @param plot use TRUE to show a plot of PC1 versus PC2
##' @param main main title of the plot
##' @param cols N-vector of colors
##' @return list with the rotated matrix (= X V) which columns corresponds to "principal components", and with the proportion of variance explained per PC
##' @author Timothee Flutre
##' @export
pca <- function(X, ct=TRUE, sc=FALSE, plot=FALSE, main="PCA", cols=NULL){
  stopifnot(is.matrix(X),
            is.logical(ct),
            is.logical(sc),
            is.logical(plot))
  if(plot & ! is.null(cols))
    stopifnot(is.vector(cols),
              length(cols) == nrow(X))

  X <- scale(x=X, center=ct, scale=sc)

  res.svd <- svd(x=X) # X = U D V^T

  rotation <- X %*% res.svd$v
  colnames(rotation) <- paste0("PC", 1:ncol(rotation))

  prop.vars <- res.svd$d / sqrt(max(1, nrow(X) - 1))
  prop.vars <- prop.vars / sum(prop.vars)
  names(prop.vars) <- colnames(rotation)

  if(plot){
    plot(x=rotation[,1], y=rotation[,2], las=1,
         xlab=paste0("PC1 (", format(100 * prop.vars[1], digits=3), "%)"),
         ylab=paste0("PC2 (", format(100 * prop.vars[2], digits=3), "%)"),
         main=main, type=ifelse(is.null(cols), "p", "n"))
    abline(h=0, lty=2); abline(v=0, lty=2)
    if(! is.null(cols))
      points(x=rotation[,1], y=rotation[,2], col=cols, pch=20)
  }

  return(list(rotation=rotation,
              prop.vars=prop.vars))
}

##' Return the number of PCs that minimizes the average squared partial
##' correlation
##'
##' Shriner (Heredity, 2011)
##' @param X genotype matrix (0,1,2) with SNPs in rows and individuals in columns
##' @return integer
##' @author Daniel Shriner [aut], Timothee Flutre [cre]
##' @export
getNbPCsMinimAvgSqPartCor <- function(X){
  if(nrow(X) < ncol(X))
    warning("input matrix doesn't seem to have genes/snps in rows and samples in columns")
  mu <- apply(X, 1, mean, na.rm=TRUE)
  X <- X - mu
  X2 <- cor(X, use="complete.obs")
  a <- eigen(X2)
  a$values[a$values<0] <- 0
  b <- diag(a$values, nrow=length(a$values))
  loadings <- a$vectors %*% sqrt(b)
  partial <- function(x) {
    c <- loadings[,1:x]
    partcov <- X2 - (c %*% t(c))
    d <- diag(partcov)
    if(any(is.element(NaN,d), is.element(0,d), length(d[d<0])!=0)) {
      map <- 1
    } else {
      d <- 1/(sqrt(d))
      e <- diag(d, nrow=length(d))
      pr <- e %*% partcov %*% e
      map <- (sum(pr^2) - ncol(X2)) / (ncol(X2) * (ncol(X2) - 1))
    }
    return(map)
  }
  fm <- sapply(1:(ncol(X2) - 1), partial)
  fm <- c((sum(X2^2) - ncol(X2))/(ncol(X2) * (ncol(X2) - 1)), fm)
  return(max(1, which.min(fm) - 1))
}

##' Quantile-normalize a vector of numbers to a standard normal distribution.
##'
##' Bolstad et al, Bioinformatics, 2003
##' @param x vector of numeric data
##' @param break.ties.rand break ties randomly (default=TRUE)
##' @param seed see for the pseudo-random number generator (default=1859)
##' @return vector
##' @author Timothee Flutre
##' @export
quantNorm <- function(x, break.ties.rand=TRUE, seed=1859){
  stopifnot(is.vector(x), is.numeric(x), is.logical(break.ties.rand),
            is.numeric(seed))

  out <- setNames(object=rep(NA, length(x)), nm=names(x))

  if(break.ties.rand){
    if(! is.null(seed))
      set.seed(seed)
    idx <- sample.int(n=length(x))
    tmp <- qqnorm(y=x[idx], plot.it=FALSE)$x
    out <- tmp[sort(idx, index.return=TRUE)$ix]
  } else
    out <- qqnorm(y=x, plot.it=FALSE)$x

  return(out)
}

##' Scales a correlation matrix into the corresponding covariance matrix efficiently.
##'
##' Use sweep.
##' See \url{https://en.wikipedia.org/wiki/Covariance_matrix#Correlation_matrix}.
##' @param x correlation matrix
##' @param sd standard deviations
##' @return matrix
##' @author Timothee Flutre
##' @export
cor2cov <- function(x, sd){
  ## D <- diag(sd); return(D %*% x %*% D)
  return(sweep(sweep(x, 1, sd, "*"), 2, sd, "*"))
}

##' The matrix-variate Normal distribution
##'
##' Random generation for the matrix-variate Normal distribution.
##' See \url{https://en.wikipedia.org/wiki/Matrix_normal_distribution}.
##' @param n number of observations
##' @param M mean matrix
##' @param U between-row covariance matrix
##' @param V between-column covariance matrix
##' @return array
##' @author Timothee Flutre
##' @export
rmatnorm <- function(n=1, M, U, V){
  stopifnot(nrow(M) == nrow(U),
            ncol(M) == nrow(V),
            nrow(U) == ncol(U),
            nrow(V) == ncol(V))

  ## chol returns upper triangular factor of Cholesky decomp
  chol.tU <- chol(t(U)) # A
  chol.V <- chol(V) # B

  ## for X ~ MN(M, AA', B'B): draw Z ~ MN(0, I, I), then X = M + A Z B
  tmp <- lapply(1:n, function(i){
    Z <- matrix(data=rnorm(n=nrow(M) * ncol(M), mean=0, sd=1),
                nrow=nrow(M), ncol=ncol(M))
    matrix(data=M + chol.tU %*% Z %*% chol.V,
           nrow=nrow(M), ncol=ncol(M))
  })

  return(array(data=do.call(c, tmp),
               dim=c(nrow(M), ncol(M), n)))
}

## to check rmatnorm above:
## Sigma <- matrix(c(3,2,2,4), nrow=2, ncol=2)
## rho <- Sigma[2,1] / prod(sqrt(diag(Sigma)))
## samples <- rmatnorm(n=100, M=matrix(0, nrow=10^3, ncol=2),
##                     U=diag(10^3), V=Sigma)
## tmp <- t(apply(samples, 3, function(mat){
##   c(var(mat[,1]), var(mat[,2]), cor(mat[,1], mat[,2]))
## }))
## summary(tmp) # corresponds well to Sigma

##' MCMC diagnostics
##'
##' Plot the trace, autocorrelation, running average and density side by side for a single chain.
##' TODO: show dark background for burn-in in traceplot
##' @param res.mcmc object of class "mcmc" (not "mcmc.list"!)
##' @param param.name parameter name in the columns of res.mcmc
##' @param subplots vector indicating which sub-plot(s) to make (1 for trace, 2 for autocorrelation, 3 for running quantiles, 4 for density)
##' @param pe point estimate (optional, but required if "hi" is given)
##' @param hi half-interval (optional)
##' @return nothing
##' @author Timothee Flutre
##' @export
plotMcmcChain <- function(res.mcmc, param.name, subplots=1:4,
                          pe=NULL, hi=NULL){
  if(! requireNamespace("coda", quietly=TRUE))
    stop("Pkg coda needed for this function to work. Please install it.",
         call.=FALSE)
  stopifnot(coda::is.mcmc(res.mcmc),
            is.character(param.name),
            param.name %in% colnames(res.mcmc),
            is.vector(subplots), is.numeric(subplots))
  if(! is.null(hi))
    stopifnot(! is.null(pe))

  changed.par <- FALSE
  if(1 %in% subplots & 2 %in% subplots & 3 %in% subplots & 4 %in% subplots &
     all(par("mfrow") == c(1, 1))){
    par(mfrow=c(2,2))
    changed.par <- TRUE
  }

  if(1 %in% subplots){
    ## ts.plot(res.mcmc[, param.name],
    coda::traceplot(res.mcmc[, param.name],
                    ## xlab="Iterations",
                    ylab="Trace",
                    main=paste0(param.name))
    if(! is.null(pe)){
      abline(h=pe, col="red")
      if(! is.null(hi)){
        abline(h=pe + 2 * hi, col="red", lty=2)
        abline(h=pe - 2 * hi, col="red", lty=2)
      }
    }
  }

  if(2 %in% subplots){
    ## acf(res.mcmc[, param.name],
    coda::autocorr.plot(res.mcmc[, param.name],
                        main=paste0(param.name),
                        ## xlab="Lag",
                        ## ylab="Autocorrelation",
                        auto.layout=FALSE,
                        ask=FALSE)
  }

  if(3 %in% subplots){
    ## ts.plot(cumsum(res.mcmc[, param.name]) / seq(along=res.mcmc[, param.name]),
    coda::cumuplot(res.mcmc[, param.name],
                   probs=c(0.025, 0.5, 0.975),
                   main=paste0(param.name),
                   ## xlab="Iterations",
                   ylab="Running quantiles",
                   auto.layout=FALSE,
                   ask=FALSE)
    if(! is.null(pe)){
      abline(h=pe, col="red")
      if(! is.null(hi)){
        abline(h=pe + 2 * hi, col="red", lty=2)
        abline(h=pe - 2 * hi, col="red", lty=2)
      }
    }
  }

  if(4 %in% subplots){
    ## plot(density(res.mcmc[, param.name]), type="l",
    coda::densplot(res.mcmc[, param.name],
                   ylab="Density",
                   main=paste0(param.name))
    if(! is.null(pe)){
      abline(v=pe, col="red")
      if(! is.null(hi)){
        abline(v=pe + 2 * hi, col="red", lty=2)
        abline(v=pe - 2 * hi, col="red", lty=2)
      }
    }
  }

  if(changed.par)
    par(mfrow=c(1,1))
}
