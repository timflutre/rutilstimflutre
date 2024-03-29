% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quantgen.R
\name{genoClasses2JoinMap}
\alias{genoClasses2JoinMap}
\title{Convert genotypes}
\usage{
genoClasses2JoinMap(
  x,
  reformat.input = TRUE,
  na.string = "--",
  thresh.counts = NULL,
  thresh.na = NULL,
  is.F2 = FALSE,
  verbose = 1
)
}
\arguments{
\item{x}{data.frame of bi-allelic SNP genotypes, with SNPs in rows and genotypes in columns; row names should contain SNP identifiers, the first column should contain the SNP genotypes of the first parent (traditionnaly the mother), the second column should contain the SNP genotypes of the second parent (traditionnaly the father), and the remaining columns should contain the SNP genotypes of the offsprings (full siblings); if it is a matrix, it will be silently transformed into a data frame}

\item{reformat.input}{if TRUE, the function \code{\link{reformatGenoClasses}} will be used}

\item{na.string}{a character to be interpreted as NA values by \code{\link{reformatGenoClasses}}}

\item{thresh.counts}{threshold (per SNP) on the number of offsprings having a particular genotypic class, below which counts are converted into \code{NA}}

\item{thresh.na}{threshold (per SNP) on the number of offsprings having \code{NA} (applied after \code{thresh.counts})}

\item{is.F2}{if TRUE, it is assumed that the two outbred parents were crossed once to make a F1 which was then autofecundated multiple times to make the offsprings; note that the SNP genotypes of the F1 should not be given (and they are often unavailable)}

\item{verbose}{verbosity level (0/1)}
}
\value{
data.frame
}
\description{
Convert SNP genotypes of an outbred bi-parental cross from genotypic classes into the \href{https://www.kyazma.nl/index.php/JoinMap/}{JoinMap} format.
For the moment, missing genotypes in parents result in the SNP being ignored, but we could imagine using genotypes in offsprings to impute such cases.
}
\examples{
\dontrun{
nb.snps <- 6
x <- data.frame(par1=c("AA", "GC", "CG", "AT", NA, "AA"),
                par2=c("AT", "GC", "GG", "AT", "AT", "AT"),
                off1=c("AA", "GG", "CG", "AA", "AA", "AT"),
                off2=c("AT", "GG", "CG", "AT", "AT", "AA"),
                off3=c("AT", "GG", "GG", "TT", "TT", NA),
                off4=c(NA, NA, NA, NA, NA, NA),
                row.names=paste0("snp", 1:nb.snps),
                stringsAsFactors=FALSE)
genoClasses2JoinMap(x=x, reformat.input=TRUE, thresh.na=2, verbose=1)
}
}
\seealso{
\code{\link{reformatGenoClasses}}, \code{\link{writeSegregJoinMap}}
}
\author{
Timothee Flutre
}
