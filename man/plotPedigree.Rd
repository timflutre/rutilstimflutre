% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quantgen.R
\name{plotPedigree}
\alias{plotPedigree}
\title{Plot pedigree}
\usage{
plotPedigree(
  inds,
  mothers,
  fathers,
  generations,
  sexes = NULL,
  plot.it = TRUE,
  verbose = 1,
  edge.col.mother = "black",
  edge.col.father = "darkgrey",
  vertex.label.color = "darkblue",
  vertex.color = "white",
  vertex.size = 20,
  vertex.shape = "none",
  vertex.label.family = "Helvetica",
  mult.edge.curve = 0.25,
  edge.arrow.width = 0.75,
  edge.arrow.size = 0.75,
  xmin = -1,
  xmax = 1,
  ymin = -1,
  ymax = 1,
  ...
)
}
\arguments{
\item{inds}{identifiers of the genotypes}

\item{mothers}{identifiers of their mother; can be NA (case of the founders and haplodiploidization)}

\item{fathers}{identifiers of their father; can be NA (case of the founders and haplodiploidization)}

\item{generations}{should start at 0}

\item{sexes}{"F" for female (circle), "M" for male (square) and "H" for hermaphrodite (triangle); can also be NA (no shape)}

\item{plot.it}{if TRUE, the pedigree will be plotted}

\item{verbose}{verbosity level (0/1/2)}

\item{edge.col.mother}{see ?igraph.plotting}

\item{edge.col.father}{see ?igraph.plotting}

\item{vertex.label.color}{see ?igraph.plotting}

\item{vertex.color}{see ?igraph.plotting}

\item{vertex.size}{see ?igraph.plotting}

\item{vertex.shape}{see ?igraph.plotting}

\item{vertex.label.family}{see ?igraph.plotting}

\item{mult.edge.curve}{see ?igraph.plotting}

\item{edge.arrow.width}{see ?igraph.plotting}

\item{edge.arrow.size}{see ?igraph.plotting}

\item{xmin}{see \code{norm_coords} in the "igraph" package}

\item{xmax}{see \code{norm_coords} in the "igraph" package}

\item{ymin}{see \code{norm_coords} in the "igraph" package}

\item{ymax}{see \code{norm_coords} in the "igraph" package}

\item{...}{other plotting options; see ?plot.igraph and ?igraph.plotting}
}
\value{
invisible list with objects required to plot the pedigree
}
\description{
Plot a pedigree using the "igraph" package.
This function was inspired by plot.pedigree() from the "synbreed" package (under GPL-3).
It add options for monoecious species and auto-fecondation.
}
\seealso{
\code{\link{subsetPedigree}}
}
\author{
Timothee Flutre
}
