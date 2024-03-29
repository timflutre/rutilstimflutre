% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/molbio_plates.R
\name{loadPlates}
\alias{loadPlates}
\title{Load plate(s)}
\usage{
loadPlates(files, sep = NULL, verbose = 1)
}
\arguments{
\item{files}{vector of paths to file(s), one per plate}

\item{sep}{separator of columns; if NULL, read.csv will be used}

\item{verbose}{verbosity level (0/default=1)}
}
\value{
list of matrices, one per plate, in the "wide" format
}
\description{
Read each file into a matrix and gather them into a list.
}
\author{
Timothee Flutre
}
