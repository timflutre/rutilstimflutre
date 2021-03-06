% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grapevine.R
\name{getSsrGenosFromPlantGrapePdf}
\alias{getSsrGenosFromPlantGrapePdf}
\title{SSR genotypes from PlantGrape}
\usage{
getSsrGenosFromPlantGrapePdf(file, var.name, sep = "/")
}
\arguments{
\item{file}{path to the PDF file}

\item{var.name}{name of the grape variety}

\item{sep}{character to use to separate both alleles of a given SSR}
}
\value{
matrix with one row and as many columns as SSRs, usable by "df2genind" from the "adegenet" package
}
\description{
Extract SSR genotypes from a PDF produced by \href{http://plantgrape.plantnet-project.org/fr/}{PlantGrape} for a given grape variety.
}
\author{
Timothee Flutre
}
