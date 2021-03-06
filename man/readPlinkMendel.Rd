% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bioinfo_sequences.R
\name{readPlinkMendel}
\alias{readPlinkMendel}
\title{Read PLINK}
\usage{
readPlinkMendel(prefix, suffix = "", verbose = 1)
}
\arguments{
\item{prefix}{prefix of the file to read (e.g. \code{"~/work/output_plink"})}

\item{suffix}{suffix of the file to read
\itemize{
\item fmendel: one row per parental pair
\item imendel: one row per individual per parental pair
\item lmendel: one row per variant (SNP)
\item mendel: one row per error (Mendelian violation)
}}

\item{verbose}{verbosity level (0/1)}
}
\value{
data frame
}
\description{
Read output files from the \href{https://www.cog-genomics.org/plink2}{PLINK} software used with the \code{--mendel} option.
For potentially large files ("lmendel" and "mendel"), the \code{fread} function from the \href{https://cran.r-project.org/package=data.table}{data.table} package is used as it is faster than \code{\link[utils]{read.table}}.
}
\examples{
\dontrun{## assuming properly-formatted data, launch PLINK via a system call
cmd <- "plink --mendel --bed input.bed --bim input.bim --fam input.fam --out output"
system(cmd)

pl.mend <- readPlinkMendel(prefix="output", suffix="mendel")
str(pl.mend)

pl.imend <- readPlinkMendel(prefix="output", suffix="imendel")
head(pl.imend[order(pl.imend$N, decreasing=TRUE),])
}
}
\author{
Timothee Flutre
}
