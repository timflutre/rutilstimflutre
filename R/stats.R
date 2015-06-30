## Contains functions useful for statistics.

##' Return the Root Mean Squared Error
##'
##'
##' @title Root Mean Squared Error
##' @param error vector \eqn{\hat{\theta}_i - \theta_i}
##' @return numeric
##' @author Timothée Flutre
rmse <- function(error){
  sqrt(mean(error^2))
}

##' Return the Mean Absolute Error
##'
##'
##' @title Mean Absolute Error
##' @param error vector \eqn{\hat{\theta}_i - \theta_i}
##' @return numeric
##' @author Timothée Flutre
mae <- function(error){
  mean(abs(error))
}

##' Return the Mean Signed Difference
##'
##'
##' @title Mean Signed Difference
##' @param error vector \eqn{\hat{\theta}_i - \theta_i}
##' @return numeric
##' @author Timothée Flutre
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
##' @author Timothée Flutre
binary.classif <- function(known.nulls, called.nulls){
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
##' Use equal weights if not specified
##' @title Log of weighted sum
##' @param x vector
##' @param weights weights
##' @return numeric
##' @author Timothée Flutre
log10.weighted.sum <- function(x, weights=NULL){
  if(is.null(weights))
    weights <- rep(1/length(x), length(x))
  max <- max(x)
  max + log10(sum(weights * 10^(x - max)))
}

##' Return the Moore-Penrose pseudo-inverse of a matrix
##'
##'
##' @title Pseudo-inverse
##' @param mat matrix
##' @return matrix
##' @author Timothée Flutre
mp.inv <- function(mat){
  mat.svd <- svd(mat)
  mat.svd$v %*% diag(1/mat.svd$d) %*% t(mat.svd$u)
}

##' Return the number of PCs that minimizes the average squared partial
##' correlation
##'
##' Shriner (Heredity, 2011)
##' @param X genotype matrix (0,1,2) with SNPs in rows and individuals in columns
##' @return integer
##' @author Shriner
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
##' TODO: add ref
##' @param x vector of numeric data
##' @param break.ties.rand break ties randomly (default=TRUE)
##' @param seed see for the pseudo-random number generator (default=1859)
##' @return vector
##' @author Timothée Flutre
quant.norm <- function(x, break.ties.rand=TRUE, seed=1859){
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
##' Use sweep. https://en.wikipedia.org/wiki/Covariance_matrix#Correlation_matrix
##' @param x correlation matrix
##' @param sd standard deviations
##' @return matrix
##' @author Timothée Flutre
cor2cov <- function(x, sd){
  ## D <- diag(sd); return(D %*% x %*% D)
  return(sweep(sweep(x, 1, sd, "*"), 2, sd, "*"))
}
