% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findhap.R
\name{writeGenosFindhap}
\alias{writeGenosFindhap}
\title{Write SNP genotypes for findhap}
\usage{
writeGenosFindhap(X, genos.file)
}
\arguments{
\item{X}{matrix of bi-allelic SNP genotypes encoded in number of copies of the 2nd allele, i.e. as allele doses in \{0,1,2\}, with genotypes in rows and SNPs in columns; missing values should be encoded as NA; the maximum length of genotypes identifiers is 10 characters}

\item{genos.file}{path to the file in which genotypes will be saved}
}
\value{
invisible data frame which is written in \code{genos.file}
}
\description{
Write SNP genotypes into a file formatted for findhap.
}
\examples{
\dontrun{## fake genotypes
set.seed(1)
nb.inds <- 3
nb.snps <- 5
X <- matrix(data=sample(c(0,1,2,NA), size=nb.inds * nb.snps, replace=TRUE),
            nrow=nb.inds, ncol=nb.snps,
            dimnames=list(paste0("ind", 1:nb.inds),
                          paste0("snp", 1:nb.snps)))
f <- tempfile()
out <- writeGenosFindhap(X, f)
}
}
\author{
Timothee Flutre
}
