## Contains functions useful for statistics.

##' Summary
##'
##' Print the output of \code{\link{summary}} in a single line, as well as the total number of observations and the number of missing data.
##' @param x vector of numbers
##' @param spec specifier, see \code{\link{sprintf}}
##' @return invisible summary
##' @author Timothee Flutre
##' @export
prettyPrintSummary <- function(x, spec="%.2f"){
  stopifnot(is.vector(x),
            is.numeric(x))

  out <- stats::setNames(rep(NA, 8),
                         c("n", "na", "min", "q1", "med", "mean", "q3", "max"))

  out["n"] <- length(x)
  isNa <- is.na(x)
  if(any(isNa))
    x <- x[! isNa]
  out["na"] <- sum(isNa)
  out["min"] <- min(x)
  out["q1"] <- stats::quantile(x, 0.25)
  out["med"] <- stats::median(x)
  out["mean"] <- mean(x)
  out["q3"] <- stats::quantile(x, 0.75)
  out["max"] <- max(x)

  fmt <- paste(paste0(names(out), "=", spec), collapse=" ")
  txt <- sprintf(fmt=fmt, out["n"], out["na"], out["min"], out["q1"], out["med"],
                 out["mean"], out["q3"], out["max"])
  print(txt)

  invisible(out)
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

##' Moore-Penrose pseudo-inverse
##'
##' Return the Moore-Penrose pseudo-inverse of a matrix (Golub & Van Loan, Matrix Computations, 3rd edition, ch5, p257).
##' @title Pseudo-inverse
##' @param x matrix
##' @return matrix
##' @author Timothee Flutre
##' @export
mpInv <- function(x){
  stopifnot(is.matrix(x))

  mat.svd <- svd(x)
  out <- mat.svd$v %*% diag(1/mat.svd$d) %*% t(mat.svd$u)
  dimnames(out) <- dimnames(x)

  return(out)
}

##' Principal component analysis
##'
##' It is performed via the singular value decomposition (SVD) of the usually-centered data matrix, Xc = U D V^T. This is mostly for teaching purposes, see \code{\link[stats]{prcomp}} otherwise. A good reason to center the data matrix for PCA is given in \href{http://link.springer.com/10.1007/s11063-007-9069-2}{Miranda et al (2008)}.
##' @param X data matrix with N rows and P columns
##' @param ct use TRUE to center the columns of X (recommended), FALSE otherwise
##' @param sc use TRUE to scale the columns of X (if different units), FALSE otherwise
##' @param plot if not NULL, use "points" to show a plot with \code{\link[graphics]{points}} of PC1 versus PC2, and "text" to use \code{\link[graphics]{text}} with row names of \code{X} as labels
##' @param main main title of the plot
##' @param cols N-vector of colors
##' @return list with the rotated matrix (= X V) which rows correspond to the original rows after translation towards the sample mean (if center=TRUE) and rotation onto the "principal components" (eigenvectors of the sample covariance matrix), and with the proportion of variance explained per PC
##' @author Timothee Flutre
##' @seealso \code{\link{plotPca}}
##' @examples
##' \dontrun{set.seed(1859)
##' genomes <- simulCoalescent(nb.inds=200, nb.pops=3, mig.rate=3)
##' X <- genomes$genos
##' A <- estimGenRel(X)
##' imageWithScale(A, main="Additive genetic relationships") # we clearly see 3 clusters
##' out.prcomp <- prcomp(x=X, retx=TRUE, center=TRUE, scale.=FALSE) # uses SVD
##' summary(out.prcomp)$importance[,1:4]
##' out.prcomp$sdev[1:4]
##' (out.prcomp$sdev^2 / sum(out.prcomp$sdev^2))[1:4]
##' head(out.prcomp$rotation[, 1:4]) # first four PCs
##' head(out.prcomp$x[, 1:4]) # rotated data
##' out.princomp <- princomp(x=X) # uses EVD, and requires more units than variables
##' out.pca <- pca(X=X, ct=TRUE, sc=FALSE)
##' out.pca$prop.vars[1:4]
##' head(out.pca$rotation[, 1:4]) # rotated data
##' }
##' @export
pca <- function(X, ct=TRUE, sc=FALSE, plot=NULL, main="PCA",
                cols=rep("black", nrow(X))){
  stopifnot(is.matrix(X),
            is.logical(ct),
            is.logical(sc))
  if(! is.null(plot))
    stopifnot(plot %in% c("points", "text"),
              is.vector(cols),
              length(cols) == nrow(X))

  X <- scale(x=X, center=ct, scale=sc)

  res.svd <- svd(x=X) # X = U D V^T

  rotation <- X %*% res.svd$v
  colnames(rotation) <- paste0("PC", 1:ncol(rotation))

  prop.vars <- res.svd$d / sqrt(max(1, nrow(X) - 1))
  prop.vars <- prop.vars^2 / sum(prop.vars^2)
  names(prop.vars) <- colnames(rotation)

  if(! is.null(plot))
    plotPca(rotation=rotation, prop.vars=prop.vars, plot=plot, main=main,
            cols=cols)

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
  X2 <- stats::cor(X, use="complete.obs")
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

  out <- stats::setNames(object=rep(NA, length(x)), nm=names(x))

  if(break.ties.rand){
    if(! is.null(seed))
      set.seed(seed)
    idx <- sample.int(n=length(x))
    tmp <- stats::qqnorm(y=x[idx], plot.it=FALSE)$x
    out <- tmp[sort(idx, index.return=TRUE)$ix]
  } else
    out <- stats::qqnorm(y=x, plot.it=FALSE)$x

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

##' Singular matrix
##'
##' Assess if a matrix is singular, i.e. not full rank, by comparing its condition number with the relative machine precision.
##' As \href{http://www.stat.wisc.edu/~st849-1/Rnotes/ModelMatrices.html#sec-2_2}{explained by Douglas Bates}, "a matrix is regarded as being numerically singular when its reciprocal condition number is less than the relative machine precision".
##' @param x matrix
##' @return logical
##' @author Timothee Flutre
##' @export
isSingular <- function(x){
  return(1 / kappa(x) <= 100 * .Machine$double.eps)
}

##' Matrix-variate Normal distribution
##'
##' Random generation for the matrix-variate Normal distribution.
##' See \url{https://en.wikipedia.org/wiki/Matrix_normal_distribution}.
##' @param n number of observations
##' @param M mean matrix
##' @param U between-row covariance matrix
##' @param V between-column covariance matrix
##' @param pivot 2-element vector with values TRUE/FALSE/"auto", where TRUE (FALSE) means using pivoting (or not) for Choleski decomposition of U and/or V (see \code{\link[base]{chol}}); useful when U and/or V are singular; with "auto", this will be automatically determined
##' @return array
##' @author Timothee Flutre
##' @examples
##' Sigma <- matrix(c(3,2,2,4), nrow=2, ncol=2)
##' rho <- Sigma[2,1] / prod(sqrt(diag(Sigma)))
##' samples <- rmatnorm(n=100, M=matrix(0, nrow=10^3, ncol=2),
##'                     U=diag(10^3), V=Sigma)
##' tmp <- t(apply(samples, 3, function(mat){
##'   c(var(mat[,1]), var(mat[,2]), cor(mat[,1], mat[,2]))
##' }))
##' summary(tmp) # corresponds well to Sigma
##' @export
rmatnorm <- function(n=1, M, U, V,
                     pivot=c(U=FALSE, V=FALSE)){
  stopifnot(is.matrix(M),
            is.matrix(U),
            is.matrix(V),
            nrow(M) == nrow(U),
            ncol(M) == nrow(V),
            nrow(U) == ncol(U),
            nrow(V) == ncol(V),
            is.vector(pivot),
            all(c("U","V") %in% names(pivot)),
            all(pivot %in% c(TRUE, FALSE, "auto")))

  ## chol() returns upper triangular factor of Cholesky decomp
  if(pivot["U"] == "auto")
    pivot["U"] <- isSingular(U)
  chol.tU <- chol(t(U), pivot=pivot["U"]) # A
  if(pivot["V"] == "auto")
    pivot["V"] <- isSingular(V)
  chol.V <- chol(V, pivot=pivot["V"]) # B

  ## for X ~ MN(M, AA', B'B): draw Z ~ MN(0, I, I), then X = M + A Z B
  tmp <- lapply(1:n, function(i){
    Z <- matrix(data=stats::rnorm(n=nrow(M) * ncol(M), mean=0, sd=1),
                nrow=nrow(M), ncol=ncol(M))
    matrix(data=M + chol.tU %*% Z %*% chol.V,
           nrow=nrow(M), ncol=ncol(M))
  })

  return(array(data=do.call(c, tmp),
               dim=c(nrow(M), ncol(M), n)))
}

##' AR(1)
##'
##' Return a first-order auto-regressive correlation matrix, as noted in equation 4 of \href{http://dx.doi.org/10.1139/x02-111}{Dutkowski et al (2002)}.
##' @param n dimension of the matrix (number of rows and columns)
##' @param rho correlation between successive variables
##' @return matrix
##' @author Timothee Flutre
##' @export
corrMatAR1 <- function(n, rho){
  stopifnot(is.numeric(n),
            n >= 0,
            is.numeric(rho),
            abs(rho) <= 1)

  M <- diag(n)
  M <- row(M) - col(M)
  M <- abs(M)

  return(rho^M)
}

##' AR(1)
##'
##' Return a first-order auto-regressive covariance matrix.
##' @param n dimension of the matrix (number of rows and columns)
##' @param rho correlation between successive variables
##' @param sigma2 variance of the errors (also called "innovations" in the time-series literature)
##' @return matrix
##' @author Timothee Flutre
##' @export
covMatAR1 <- function(n, rho, sigma2){
  stopifnot(is.numeric(n),
            n >= 0,
            is.numeric(rho),
            abs(rho) <= 1,
            is.numeric(sigma2),
            sigma2 > 0)

  corrMat <- corrMatAR1(n=n, rho=rho)

  Sigma <- (sigma2 / (1 - rho^2)) * corrMat

  return(Sigma)
}

##' AR(1)
##'
##' Return a first-order auto-regressive precision matrix, as in section 1.1.1 of Rue and Held (2005).
##' @param n dimension of the matrix (number of rows and columns)
##' @param rho correlation between successive variables
##' @param sigma2 variance of the errors (also called "innovations" in the time-series literature)
##' @return sparse matrix
##' @author Timothee Flutre
##' @export
precMatAR1 <- function(n, rho, sigma2){
  requireNamespace("Matrix")
  stopifnot(is.numeric(n),
            n >= 0,
            is.numeric(rho),
            abs(rho) <= 1,
            is.numeric(sigma2),
            sigma2 > 0)

  tmp <- Matrix::bandSparse(n=n, m=n, k=c(1, -1),
                            diagonals=list(rep(- rho, n - 1),
                                           rep(- rho, n - 1)))

  Q <- (1 / sigma2) * (Matrix::Diagonal(x=c(1, rep(1 + rho^2, n - 2), 1)) +
                       tmp)

  return(Q)
}

##' P values
##'
##' Plot the histogram of p values with freq=FALSE to imitate figure 1 of Storey & Tibshirani (2003).
##' The x-axis goes from 0 to 1 with equidistant bins of width 0.05, and the y-axis is in the density scale.
##' For a given bin, the density of p values in this bin corresponds to the relative frequency in this bin divided by the bin width, where the relative frequency is the absolute frequency (i.e. raw counts in this bin) divided by the sum of all counts (i.e. over all bins).
##' In the density scale, the area of a given bin is equal to its density times the bin width, thus the area of the whole histogram is equal to 1.
##' Assuming all null hypotheses are true, the "density" histogram represents the density of the discrete version of the Uniform distribution, where the height of each bin is 1, so that the whole histogram area remains equal to 1.
##' @param pvalues vector of raw p values (missing values will be omitted)
##' @param breaks vector giving the breakpoints between histogram cells (see \code{\link[graphics]{hist}})
##' @param freq as in \code{\link[graphics]{hist}}; if TRUE, the histogram graphic is a representation of frequencies, the "counts" component of the result; if FALSE, probability densities, component "density", are plotted (so that the histogram has a total area of one)
##' @param main an overall title for the plot (default: "Density of <nb of non-NA p values> p values")
##' @param col a colour to be used to fill the bars (see \code{\link[graphics]{hist}})
##' @param border the color of the border around the bars (see \code{\link[graphics]{hist}})
##' @param pi0 estimate of the proportion of null hypotheses
##' @return invisible output of \code{\link[graphics]{hist}}
##' @seealso \code{\link{qqplotPval}}
##' @author Timothee Flutre
##' @examples
##' set.seed(1859)
##' P <- 1000; P1 <- 100; thresh <- 0.05
##' pvalues.0 <- setNames(runif(n=P-P1, min=0, max=1),
##'                       paste0("null",1:(P-P1)))
##' pvalues.1 <- setNames(rbeta(n=P1, shape1=1, shape2=10^3),
##'                       paste0("alt", 1:P1))
##' pvalues <- c(pvalues.0, pvalues.1)
##' pvalues[c(1,10)] <- NA
##' plotHistPval(pvalues, pi0=(P-P1)/P, freq=TRUE)
##' pval.bonf <- p.adjust(pvalues, method="bonferroni")
##' sum(pval.bonf <= thresh, na.rm=TRUE) # 9
##' pval.bh <- p.adjust(pvalues, method="BH")
##' sum(pval.bh <= thresh, na.rm=TRUE) # 104
##' if(require(qvalue)){
##'   qv <- qvalue(p=pvalues, fdr.level=thresh, pfdr=TRUE)
##'   summary(qv)
##' }
##' @export
plotHistPval <- function(pvalues, breaks=seq(0, 1, 0.05), freq=FALSE,
                         main=NULL, col="grey", border="white", pi0=NULL){
  stopifnot(is.numeric(pvalues),
            is.vector(pvalues),
            is.vector(breaks))

  isna <- is.na(pvalues)
  if(any(isna)){
    warning(paste0(sum(isna), " NA are discarded"))
    pvalues <- pvalues[! is.na(pvalues)]
  }

  if(is.null(main))
    main <- paste0("Histogram of ", length(pvalues), " p values")

  out <- graphics::hist(pvalues, breaks=breaks, xlim=c(0,1), freq=freq, las=1,
              xlab=expression(italic(p)~values), main=main,
              ## ylab=ifelse(freq, "Frequency", "Density"),
              col=col, border=border)

  ## height one would expect if all tested hypotheses were null
  if(freq){
    graphics::abline(h=length(pvalues) / length(breaks), lty=2)
  } else
    graphics::abline(h=1, lty=2)

  legs <- paste0("expected ", ifelse(freq, "frequency", "density"),
                 " if all hypotheses are null")
  cols <- "black"
  ltys <- 2
  lwds <- 1

  ## height of the estimate of the proportion of null hypotheses
  if(! is.null(pi0)){
    if(freq){
      graphics::abline(h=pi0 * length(pvalues) / length(breaks), lty=3)
    } else
      graphics::abline(h=pi0, lty=3)
    legs <- c(legs, paste0("estimated ", ifelse(freq, "frequency", "density"),
                           " of null hypotheses"))
    cols <- c(cols, "black")
    ltys <- c(ltys, 3)
    lwds <- c(lwds, 1)
  }

  graphics::legend("topright", legend=legs, col=cols, lty=ltys, lwd=lwds, bty="n")

  invisible(out)
}

##' Q-Q plot for p values
##'
##' Produce a quantile-quantile plot for p values and display its confidence interval.
##' A quantile is an order statistic, and the j-th order statistic from a Uniform(0,1) sample has a Beta(j,N-j+1) distribution (Casella & Berger, 2001, 2nd edition, p230). Let us assume we have N independent p values, \eqn{\{p_1,\ldots,p_N\}}. Under the null, they are independent and identically uniformly distributed: \eqn{\forall i \; p_i \sim \mathcal{U}_{[0,1]}}. Therefore, the 95\% confidence interval for the j-th quantile of the set of p values can be calculated with: qbeta(0.95, j, N-j+1). See also the qqman package.
##' @param pvalues vector of raw p values (missing values will be omitted)
##' @param plot.conf.int show the confidence interval
##' @param xlab a title for the x axis (see default)
##' @param ylab a title for the x axis (see default)
##' @param main an overall title for the plot (default: "Q-Q plot (<length(pvalues)> p values)")
##' @param col vector of plotting color(s) for the points (default is all points in black)
##' @param ... graphical parameters other than xlim, ylim, xlab, ylab, las and col
##' @return invisible vector of pvalues with NA omitted
##' @seealso \code{\link{plotHistPval}}
##' @author Timothee Flutre (inspired by an anonymous comment at http://gettinggeneticsdone.blogspot.fr/2009/11/qq-plots-of-p-values-in-r-using-ggplot2.html)
##' @examples
##' set.seed(1859)
##' P <- 1000; P1 <- 100; thresh <- 0.05
##' pvalues.0 <- setNames(runif(n=P-P1, min=0, max=1),
##'                       paste0("null",1:(P-P1)))
##' pvalues.1 <- setNames(rbeta(n=P1, shape1=1, shape2=10^3),
##'                       paste0("alt", 1:P1))
##' pvalues <- c(pvalues.0, pvalues.1)
##' pvalues[c(1,10)] <- NA
##' plotHistPval(pvalues, pi0=(P-P1)/P, freq=TRUE)
##' pvalues <- qqplotPval(pvalues)
##' pval.bonf <- p.adjust(pvalues, method="bonferroni")
##' sum(pval.bonf <= thresh) # 9
##' names(pvalues)[pval.bonf <= thresh]
##' pval.bh <- p.adjust(pvalues, method="BH")
##' sum(pval.bh <= thresh) # 104
##' names(pvalues)[pval.bh <= thresh]
##' lim.pval.bonf <- sort(pvalues[pval.bonf <= thresh], decreasing=TRUE)[1]
##' abline(h=-log10(lim.pval.bonf), lty=2)
##' lim.pval.BH <- sort(pvalues[pval.bh <= thresh], decreasing=TRUE)[1]
##' abline(h=-log10(lim.pval.BH), lty=3)
##' legend(x=2.2, y=1.8, legend=c(paste0(sum(pval.bonf <= thresh),
##'  " tests significant at 5%\nwith FWER controlled via Bonferroni"),
##'  paste0(sum(pval.bh <= thresh), " tests significant at 5%",
##'  "\nwith FDR controlled via BH")), lty=c(2,3), bty="n", y.intersp=2)
##' @export
qqplotPval <- function(pvalues, plot.conf.int=TRUE,
                       xlab=expression(Expected~~-log[10](italic(p)~values)),
                       ylab=expression(Observed~~-log[10](italic(p)~values)),
                       main=NULL, col=NULL){
  stopifnot(is.numeric(pvalues),
            is.vector(pvalues))
  if(! is.null(col))
    stopifnot(is.vector(col),
              length(col) == length(pvalues))

  if(is.null(col))
    col <- rep(1, length(pvalues))

  isna <- is.na(pvalues)
  if(any(isna)){
    warning(paste0(sum(isna), " NA are discarded"))
    pvalues <- pvalues[! isna]
    col <- col[! isna]
  }

  N <- length(pvalues)
  expected <- - log10(1:N / N)
  observed <- - log10(pvalues)
  MAX <- max(c(expected, observed))

  if(plot.conf.int){
    c95 <- rep(0, N)
    c05 <- rep(0, N)
    for(j in 1:N){
      c95[j] <- stats::qbeta(0.95, j, N-j+1)
      c05[j] <- stats::qbeta(0.05, j, N-j+1)
    }
    c95 <- - log10(c95)
    c05 <- - log10(c05)
    graphics::plot(expected, c95, ylim=c(0,MAX), xlim=c(0,MAX), type="l",
         axes=FALSE, xlab="", ylab="")
    graphics::par(new=T)
    graphics::plot(expected, c05, ylim=c(0,MAX), xlim=c(0,MAX), type="l",
         axes=FALSE, xlab="", ylab="")
    graphics::par(new=T)
  }

  if(is.null(main))
    main <- paste0("Q-Q plot (", N, " p values)")

  graphics::plot(x=sort(expected), y=sort(observed),
       xlim=c(0,MAX), ylim=c(0,MAX),
       las=1, col=col[order(observed)],
       xlab=xlab, ylab=ylab, main=main)
  graphics::abline(0, 1, col="red")

  invisible(pvalues)
}

##' FDR
##'
##' Estimate pi0 (proba for a null hypothesis to be true) via the EBF procedure (Wen, arXiv:1311.3981).
##' @param log10.bfs vector containing the log10(BF) of each test (NA will be discarded)
##' @param verbose verbosity level (0/1)
##' @return numeric
##' @author Timothee Flutre
##' @seealso \code{\link{estimatePi0WithQbf}}, \code{\link{controlBayesFdr}}
##' @export
estimatePi0WithEbf <- function(log10.bfs, verbose=1){
  stopifnot(is.numeric(log10.bfs),
            is.vector(log10.bfs))

  pi0.ebf <- 1

  isna <- is.na(log10.bfs)
  if(any(isna)){
    warning(paste0(sum(isna), " NA are discarded"))
    log10.bfs <- log10.bfs[! isna]
  }
  if(verbose > 0)
    message(paste0("nb of tests: ", length(log10.bfs)))

  tmp <- log10.bfs[order(log10.bfs)] # sort in increasing order

  d0 <- which(cumsum(10^tmp) / seq_along(10^tmp) >= 1)[1]
  if(length(d0) > 0 & ! is.na(d0)){
    if(verbose > 0)
      message(paste0("cutoff at the ", d0, "-th BF"))

    pi0.ebf <- d0 / length(tmp)
    if(verbose > 0)
      message(paste0("estimate pi0-hat = ",
                     format(x=pi0.ebf, scientific=TRUE, digits=6)))
  }

  return(pi0.ebf)
}

##' FDR
##'
##' Estimate pi0 (proba for a null hypothesis to be true) via the QBF procedure (Wen, arXiv:1311.3981).
##' @param log10.bfs matrix with tests in rows and two columns, the true log10(BF) and the gamma-quantile log10(BF) under the null
##' @param gamma level of the quantile (e.g. 0.5 for the median)
##' @param verbose verbosity level (0/1)
##' @return numeric
##' @seealso \code{\link{estimatePi0WithEbf}}, \code{\link{controlBayesFdr}}
##' @author Timothee Flutre
##' @export
estimatePi0WithQbf <- function(log10.bfs, gamma=0.5, verbose=1){
  stopifnot(is.numeric(log10.bfs),
            is.matrix(log10.bfs),
            ncol(log10.bfs) == 2,
            ! any(is.na(log10.bfs)))
  if(verbose > 0)
    message(paste0("nb of tests: ", nrow(log10.bfs)))

  pi0.qbf <- sum(log10.bfs[,1] <= log10.bfs[,2]) / (nrow(log10.bfs) * gamma)
  if(verbose > 0)
    message(paste0("estimate pi0-hat = ",
                   format(x=pi0.qbf, scientific=TRUE, digits=6)))

  return(pi0.qbf)
}

##' FDR
##'
##' Call significant tests by controlling the Bayesian FDR via the procedure from Newton et al (Biostatistics, 2004) also described in Muller et al (JASA, 2006).
##' @param log10.bfs vector containing the log10(BF) of each test
##' @param pi0 estimate of the proba for a null hypothesis to be true
##' @param fdr.level threshold below which a null is rejected
##' @param verbose verbosity level (0/1)
##' @return logical vector with TRUE if null is rejected (thus called significant)
##' @seealso \code{\link{estimatePi0WithEbf}}, \code{\link{estimatePi0WithQbf}}
##' @author Timothee Flutre
##' @examples
##' set.seed(1859)
##' log10.bfs <- rgamma(n=1000, shape=0.5, rate=1) - 0.5 # fake but looks realistic
##' pi0 <- estimatePi0WithEbf(log10.bfs)
##' signif <- controlBayesFdr(log10.bfs, pi0)
##' @export
controlBayesFdr <- function(log10.bfs, pi0, fdr.level=0.05, verbose=1){
  stopifnot(is.numeric(log10.bfs),
            is.vector(log10.bfs),
            ! is.na(pi0))
  if(verbose > 0)
    message(paste0(length(log10.bfs), " tests and pi0-hat = ",
                   format(x=pi0, scientific=TRUE, digits=6)))

  ## compute the posterior probability of each test being null
  post.null <- pi0 / (pi0 + (1-pi0) * 10^log10.bfs)

  ## find the cutoff for which their cumulative mean is >= fdr.level
  idx <- order(post.null)
  post.null <- post.null[idx]
  for(L in 1:length(post.null))
    if(mean(post.null[1:L]) >= fdr.level)
      break
  log10.bf.L <- log10.bfs[idx[L]]
  if(verbose > 0)
    message(paste0(L, " significant tests, at cutoff log10(BF)=", log10.bf.L))

  significants <- log10.bfs >= log10.bf.L

  return(significants)
}

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
  requireNamespaces("coda")
  stopifnot(coda::is.mcmc(res.mcmc),
            is.character(param.name),
            param.name %in% colnames(res.mcmc),
            is.vector(subplots), is.numeric(subplots))
  if(! is.null(hi))
    stopifnot(! is.null(pe))

  changed.par <- FALSE
  if(1 %in% subplots & 2 %in% subplots & 3 %in% subplots & 4 %in% subplots &
     all(graphics::par("mfrow") == c(1, 1))){
    graphics::par(mfrow=c(2,2))
    changed.par <- TRUE
  }

  if(1 %in% subplots){
    ## ts.plot(res.mcmc[, param.name],
    coda::traceplot(res.mcmc[, param.name],
                    ## xlab="Iterations",
                    ylab="Trace",
                    main=paste0(param.name))
    if(! is.null(pe)){
      graphics::abline(h=pe, col="red")
      if(! is.null(hi)){
        graphics::abline(h=pe + 2 * hi, col="red", lty=2)
        graphics::abline(h=pe - 2 * hi, col="red", lty=2)
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
      graphics::abline(h=pe, col="red")
      if(! is.null(hi)){
        graphics::abline(h=pe + 2 * hi, col="red", lty=2)
        graphics::abline(h=pe - 2 * hi, col="red", lty=2)
      }
    }
  }

  if(4 %in% subplots){
    ## plot(density(res.mcmc[, param.name]), type="l",
    coda::densplot(res.mcmc[, param.name],
                   ylab="Density",
                   main=paste0(param.name))
    if(! is.null(pe)){
      graphics::abline(v=pe, col="red")
      if(! is.null(hi)){
        graphics::abline(v=pe + 2 * hi, col="red", lty=2)
        graphics::abline(v=pe - 2 * hi, col="red", lty=2)
      }
    }
  }

  if(changed.par)
    graphics::par(mfrow=c(1,1))
}

##' MCMC results
##'
##' Return effective sample size, posterior mean, time-series standard error of the posterior mean (i.e. without ignoring autocorrelation; based on an estimate of the spectral density at 0 via \code{\link[coda]{spectrum0.ar}}), posterior standard deviation, and posterior quantiles (at 0.025, 0.5, 0.975).
##' @param res.mcmc object of class "mcmc" (not "mcmc.list"!)
##' @param param.names parameter names in the columns of res.mcmc
##' @return matrix with parameters in rows
##' @author Timothee Flutre
##' @export
summaryMcmcChain <- function(res.mcmc, param.names){
  requireNamespaces("coda")
  stopifnot(coda::is.mcmc(res.mcmc),
            all(is.character(param.names)),
            all(param.names %in% colnames(res.mcmc)))

  ## safespec0 in coda is not exported
  safespectrum0 <- function(x){
    result <- try(coda::spectrum0.ar(x)$spec)
    if(class(result) == "try-error") result <- NA
    result
  }

  out <- t(apply(res.mcmc[, param.names, drop=FALSE], 2, function(samples){
    c(coda::effectiveSize(samples),
      mean(samples),
      sqrt(safespectrum0(samples) / length(samples)),
      stats::sd(samples),
      stats::quantile(samples, c(0.025, 0.5, 0.975)))
  }))
  colnames(out) <- c("ess", "mean", "se.mean", "sd",
                     "0.025q", "0.5q", "0.975q")

  return(out)
}
