% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quantgen.R
\name{calcAsymptoticBayesFactorWakefield}
\alias{calcAsymptoticBayesFactorWakefield}
\title{Asymptotic Bayes factor}
\usage{
calcAsymptoticBayesFactorWakefield(
  theta.hat,
  V,
  W = c(0.1, 0.2, 0.4, 0.8, 1.6),
  weights = NULL,
  log10 = TRUE
)
}
\arguments{
\item{theta.hat}{vector of MLE(s) of the additive genetic effect(s)}

\item{V}{vector of the corresponding variance(s) of \code{theta.hat}}

\item{W}{vector of variance(s) of the prior on theta (if several values, the ABF will be averaged over them); the vector of default values comes from the single-SNP R implementation of BLMM by Wen for his 2015 article (see https://github.com/xqwen/blmm)}

\item{weights}{weights used to average over the grid of \code{W} (all equal by default)}

\item{log10}{return the log10 of the ABF}
}
\value{
numeric
}
\description{
Calculate the asymptotic Bayes factor proposed by Wakefield in Genetic Epidemiology 33:79-86 (2009, \url{http://dx.doi.org/10.1002/gepi.20359}).
Can also be averaged over a grid of values of W, as done in various papers from Matthew Stephens' lab.
}
\author{
Timothee Flutre
}
