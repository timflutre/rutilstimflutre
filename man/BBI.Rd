% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apple.R
\name{BBI}
\alias{BBI}
\title{Biennial bearing index}
\usage{
BBI(
  dat,
  coln.prod = "prod",
  coln.geno = "geno",
  coln.rep = "tree",
  coln.year = "year",
  coln.epsilon = NULL,
  type = "classic",
  verbose = 0
)
}
\arguments{
\item{dat}{data frame}

\item{coln.prod}{column name for the production variable; missing data will be discarded}

\item{coln.geno}{column name of the genotype}

\item{coln.rep}{column name of the tree replicate}

\item{coln.year}{column name of the year (the column data should be convertible to \code{numeric})}

\item{coln.epsilon}{column name of the residuals (used only for \code{type="res_norm"})}

\item{type}{type of BBI (classic/norm/res_norm)}

\item{verbose}{verbosity level (0/1)}
}
\value{
numeric
}
\description{
Compute various types of the biennial bearing index (BBI) as defined in \href{https://dx.doi.org/10.1093/jxb/ert297}{Durand et al (2013)}.
}
\author{
Timothee Flutre
}
