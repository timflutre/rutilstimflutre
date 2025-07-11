% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quantgen.R
\name{genoDoses2ASMap}
\alias{genoDoses2ASMap}
\title{Convert genotypes}
\usage{
genoDoses2ASMap(X = NULL, tX = NULL)
}
\arguments{
\item{X}{matrix of bi-allelic SNP genotypes encoded, for each SNP, in number of copies of its second allele, i.e. as allele doses in \{0,1,2\}, with genotypes in rows and SNPs in columns; the "second" allele is arbitrary, it can correspond to the minor (least frequent) or the major (most frequent) allele}

\item{tX}{matrix with SNPs in rows and genotypes in columns}
}
\value{
list of two data frames, one per parental map
}
\description{
Convert SNP genotypes of a bi-parental cross from "allele doses" for usage in the "ASMap" package.
}
\author{
Timothee Flutre
}
