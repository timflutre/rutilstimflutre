% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bioinfo_sequences.R
\name{seqIdStartEnd2GRanges}
\alias{seqIdStartEnd2GRanges}
\title{Make GRanges}
\usage{
seqIdStartEnd2GRanges(
  seq.id,
  seq.start = NULL,
  seq.end = NULL,
  dict.file = NULL,
  subseq.name = NULL
)
}
\arguments{
\item{seq.id}{sequence identifier(s) to work on (e.g. \code{"chr2"} or \code{c("chr2","chr7")})}

\item{seq.start}{start(s) of the sequence(s) to work on (if NULL, whole seq; see \code{dict.file})}

\item{seq.end}{end(s) of the sequence(s) to work on (if NULL, whole seq; see \code{dict.file})}

\item{dict.file}{path to the SAM dict file (see \url{https://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary}) if \code{seq.id} is specified with no start/end}

\item{subseq.name}{names of the subsequence(s) (optional; e.g. \code{"gene34"} or \code{"snp005"})}
}
\value{
\code{GRanges} object from the \code{GenomicRanges} package
}
\description{
Return a \code{GRanges} object from sequence identifiers, starts and ends.
}
\seealso{
\code{\link{vcf2dosage}}, \code{\link{vcf2genoClasses}}
}
\author{
Timothee Flutre
}
