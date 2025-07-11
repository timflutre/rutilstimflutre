% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quantgen.R
\name{recodeIntoDominant}
\alias{recodeIntoDominant}
\title{SNP genotypes}
\usage{
recodeIntoDominant(X, simplify.imputed = FALSE)
}
\arguments{
\item{X}{matrix of bi-allelic SNP genotypes encoded in allele doses in \{0,1,2\}, with genotypes in rows and SNPs in columns; missing values should be encoded as NA}

\item{simplify.imputed}{if TRUE, imputed genotypes will be considered as homozygotes if less than 0.66 or more than 1.33, and heterozygotes otherwise}
}
\value{
matrix with genotypes in \{0,1\}
}
\description{
Recode SNP genotypes from additive to dominant.
}
\author{
Timothee Flutre
}
