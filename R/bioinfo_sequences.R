## Contains functions handling sequences often used in bioinformatics.

##' Calculate the GC content of a set of sequences.
##'
##' Requires the Biostrings package.
##' @param x vector of sequences (e.g. "AGGT"), possibly with names
##' @return vector
##' @author Timothée Flutre
gc.content <- function(x){
  stopifnot(is.character(x))
	suppressPackageStartupMessages(library(Biostrings))

  sapply(x, function(xi){
    sum(alphabetFrequency(DNAString(xi),
                          baseOnly=TRUE, as.prob=TRUE)[c("C","G")])
  })
}

##' Align each sequence against each other (and itself).
##'
##' Requires the Biostrings package.
##' @title All pairwise alignments
##' @param x vector of sequences (e.g. "AGGT"), possibly with names
##' @param type type of alignment (default="global", i.e. Needleman-Wunsch)
##' @return list of instances of class PairwiseAlignments
##' @author Timothée Flutre
all.pair.aligns <- function(x, type="global", ...){
	stopifnot(is.character(x))
	suppressPackageStartupMessages(library(Biostrings))

	aligns <- list()

	for(i in 1:(length(x)-1)){
		for(j in i:length(x)){
			aligns[[paste0(i,"-",j)]] <-
				pairwiseAlignment(pattern=x[j],
													subject=x[i],
													type=type,
													substitutionMatrix=
													nucleotideSubstitutionMatrix(match=1,
																											 mismatch=0,
																											 baseOnly=FALSE,
																											 type="DNA"),
													...)
		}
	}

	return(aligns)
}

##' Extract statistics from all pairwise alignments
##'
##' Requires the Biostrings package.
##' @param aligns list of instances of class PairwiseAlignments (see all.pair.aligns())
##' @param nb.sequences number of sequences
##' @return list of matrices
##' @author Timothée Flutre
stats.all.pair.aligns <- function(aligns, nb.sequences){
	stopifnot(is.list(aligns))
	suppressPackageStartupMessages(library(Biostrings))

	scores <- matrix(data=NA, nrow=nb.sequences, ncol=nb.sequences)
	dists <- matrix(data=NA, nrow=nb.sequences, ncol=nb.sequences)
	pids <- matrix(data=NA, nrow=nb.sequences, ncol=nb.sequences)
	nmatchs <- matrix(data=NA, nrow=nb.sequences, ncol=nb.sequences)
	nmismatchs <- matrix(data=NA, nrow=nb.sequences, ncol=nb.sequences)
	ninss <- matrix(data=NA, nrow=nb.sequences, ncol=nb.sequences)
	ndels <- matrix(data=NA, nrow=nb.sequences, ncol=nb.sequences)

	for(i in 1:(nb.sequences-1)){
		for(j in i:nb.sequences){
			## message(paste0(i,"-",j))
			scores[i,j] <- score(aligns[[paste0(i,"-",j)]])
			dists[i,j] <- nedit(aligns[[paste0(i,"-",j)]])
			pids[i,j] <- pid(aligns[[paste0(i,"-",j)]])
			nmatchs[i,j] <- nmatch(aligns[[paste0(i,"-",j)]])
			nmismatchs[i,j] <- nmismatch(aligns[[paste0(i,"-",j)]])
			ninss[i,j] <- sum(nindel(aligns[[paste0(i,"-",j)]])@insertion[,"Length"] != 0) # insertions in pattern wrt subject
			ndels[i,j] <- sum(nindel(aligns[[paste0(i,"-",j)]])@deletion[,"Length"] != 0) # idem
		}
	}

	return(list(scores=scores, dists=dists, pids=pids, nmatchs=nmatchs,
							nmismatchs=nmismatchs, ninss=ninss, ndels=ndels))
}
