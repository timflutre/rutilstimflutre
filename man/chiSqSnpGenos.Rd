% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quantgen.R
\name{chiSqSnpGenos}
\alias{chiSqSnpGenos}
\title{Chi-squared for Hardy-Weinberg}
\usage{
chiSqSnpGenos(
  X,
  mafs = NULL,
  c = 0.5,
  thresh.c = 0.01,
  calc.with.D = TRUE,
  calc.pval = TRUE
)
}
\arguments{
\item{X}{matrix of bi-allelic SNP genotypes encoded as allele doses in \{0,1,2\}, with genotypes in rows and SNPs in columns; missing values should be encoded as NA}

\item{mafs}{vector of minor allele frequencies; if NULL, the frequencies will be estimated by feeding \code{X} to \code{\link{estimSnpMaf}}}

\item{c}{continuity correction (\code{0} means "no correction"; usually, \code{0.5} is used, from \href{http://dx.doi.org/10.2307/2983604}{Yates (1934)})}

\item{thresh.c}{threshold on minor allele frequencies below which the continuity correction isn't applied (used when \code{c > 0}); see \href{https://dx.doi.org/10.1016\%2Fj.ajhg.2009.11.019}{Graffelman (2010)}}

\item{calc.with.D}{calculate the chi2 statistic with D, the deviation from independence for the heterozygote, as in equation 1 from Graffelman & Camarena (2007), which only requires the number of heterozygotes; otherwise, use equation 4}

\item{calc.pval}{calculate the p values associated with the test statistics (Chi-squared distribution with one degree of freedom)}
}
\value{
matrix
}
\description{
Calculate the chi2 statistic to test for Hardy-Weinberg equilibrium.
See \href{http://www.karger.com/?doi=10.1159/000108939}{Graffelman & Camarena (2007)} and \href{https://dx.doi.org/10.1002/gepi.20612}{Shriner (2011)}.
}
\examples{
\dontrun{set.seed(1859)
library(scrm)
nb.genos <- 100
Ne <- 10^4
chrom.len <- 10^5
mu <- 10^(-8)
c <- 10^(-8)
genomes <- simulCoalescent(nb.inds=nb.genos,
                           pop.mut.rate=4 * Ne * mu * chrom.len,
                           pop.recomb.rate=4 * Ne * c * chrom.len,
                           chrom.len=chrom.len)
X <- genomes$genos
out <- chiSqSnpGenos(X)
head(out)
sum(p.adjust(out[,"pvalue"], "BH") <= 0.05)
## library(HardyWeinberg) # available on CRAN
## cts <- countGenotypicClasses(X=X)[, -4]
## colnames(cts) <- c("AA","AB","BB")
## out2 <- HWChisqMat(X=cts, cc=0.5, verbose=FALSE)
## lapply(out2, head)
## HWTernaryPlot(X=cts, region=2)
}
}
\seealso{
\code{\link{estimSnpMaf}}, \code{\link{countGenotypicClasses}}
}
\author{
Timothee Flutre
}
