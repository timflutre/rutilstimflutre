% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CarthaGene.R
\name{closeCarthagene}
\alias{closeCarthagene}
\title{Close CarthaGene}
\usage{
closeCarthagene(cg, file = NULL)
}
\arguments{
\item{cg}{list output from \code{\link{openCarthagene}}}

\item{file}{path to the file containing the dataset (if not NULL, the temporary file ".2pt", written by CarthaGene based on its name, will be removed)}
}
\value{
nothing
}
\description{
Close FIFO used to interact with CarthaGene.
}
\author{
Matthieu Falque [aut], Timothee Flutre [ctb]
}
