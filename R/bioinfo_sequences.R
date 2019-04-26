## Contains functions handling sequences often used in bioinformatics.

##' Rename chromosomes
##'
##' Rename chromosomes into integers, especially useful with FImpute, GEMMA, qqman, SNPRelate, etc.
##' @param x vector of chromosome names (if factor, will be converted to character)
##' @param basic if TRUE, the unique chromosome identifiers will be sorted alphanumerically (which requires gtools), and will be renamed as 1, ..., nb of chromosomes (for SNPRelate)
##' @param prefix characters to be removed at the start of each chromosome name (case will be ignored)
##' @param toreplace2 second set of characters to be removed in chromosome names (start or end), especially useful for "random" chromosomes (case will be ignored)
##' @param prefix2 second prefix to be removed at the start of each chromosome name for those for which "prefix" didn't work (case will be ignored)
##' @param thresh.max.chr.int threshold on the maximum chromosome integer, above which the renaming will simply be the sequence from 1 to the number of chromosomes
##' @return data frame with the same length as \code{x}, and two columns named "original" and "renamed"
##' @author Timothee Flutre
##' @examples
##' \dontrun{## example from grapevine
##' chroms <- c("chr1", "chr1_random", "chr10", "chr10_random", "chrUn", "chr2")
##' chromNames2integers(x=chroms)
##'
##' ## example from apple
##' chroms <- c("Chr15", "Chr01", "Chr02", "Chr00", "Chr02")
##' chromNames2integers(x=chroms)
##'
##' ## example from cherry
##' chroms <- c("Super-Scaffold_14", "Super-Scaffold_4374", "Super-Scaffold_27")
##' chromNames2integers(x=chroms, prefix="Super-Scaffold_",
##'                     thresh.max.chr.int=500)
##'
##' ## example from apricot
##' chroms <- c("Pp08", "scaffold_51", "Pp02", "Pp01", "scaffold_23")
##' chromNames2integers(x=chroms, prefix="Pp", toreplace2=NULL,
##'                     prefix2="scaffold_")
##' }
##' @export
chromNames2integers <- function(x, basic=FALSE,
                                prefix="chr", toreplace2="_random",
                                prefix2=NULL, thresh.max.chr.int=+Inf){
  names.of.x <- names(x)
  x <- as.character(x)
  names(x) <- names.of.x
  stopifnot(is.vector(x),
            is.logical(basic))
  if(basic){
    stopifnot(requireNamespace("gtools"))
  } else{
    stopifnot(is.character(prefix),
              is.numeric(thresh.max.chr.int),
              length(unique(x)) <= thresh.max.chr.int)
    if(! is.null(toreplace2))
      stopifnot(is.null(prefix2))
    if(! is.null(prefix2))
      stopifnot(is.null(toreplace2))
  }

  output <- data.frame(original=x,
                       renamed=NA,
                       stringsAsFactors=FALSE)
  if(! is.null(names(x)))
    rownames(output) <- names(x)

  if(basic){

    sort.uniq.chr.ids <- gtools::mixedsort(unique(x))
    for(i in seq_along(sort.uniq.chr.ids)){
      idx <- which(x == sort.uniq.chr.ids[i])
      output$renamed[idx] <- i
    }

  } else{

    output$renamed <- suppressWarnings(as.integer(gsub(prefix, "", x,
                                                       ignore.case=TRUE)))

    max.chr.int <- max(output$renamed, na.rm=TRUE)

    if(any(is.na(output$renamed))){
      tmp <- data.frame(orig=x[is.na(output$renamed)],
                        idx=which(is.na(output$renamed)),
                        as.chr=NA,
                        as.int=NA,
                        stringsAsFactors=FALSE)
      if(! is.null(toreplace2)){
        tmp$as.chr <- gsub(paste0(prefix, "|", toreplace2), "", tmp$orig)
      } else if(! is.null(prefix2)){ # e.g. for PRPER-Lovell-v2
        tmp$as.chr <- gsub(prefix2, "", tmp$orig)
      }
      for(i in 1:nrow(tmp)){
        suppressWarnings(tmp$as.int[i] <- as.integer(tmp$as.chr[i]))
        if(! is.na(tmp$as.int[i])){ # e.g. chr1_random or 23 (from scaffold_23)
          stopifnot(tmp$as.int[i] > 0) # trigger error if scaffold_0
          output$renamed[tmp$idx[i]] <- max.chr.int + tmp$as.int[i]
        } else{ # e.g. chrUn
          output$renamed[tmp$idx[i]] <- 2 * max.chr.int + 1
        }
      }
    }

    if(any(output$renamed == 0)){ # e.g. for MADOM-GoldenDelicious-dh-v1
      output$renamed[output$renamed == 0] <- 2 * max.chr.int + 1
    }

    if(max.chr.int > thresh.max.chr.int){ # e.g. for PRUAV-Regina-v1
      uniq.chr.names <- sort(unique(x))
      if(requireNamespace("gtools")){
        uniq.chr.names <- gtools::mixedsort(uniq.chr.names)
      } else
        uniq.chr.names <- sort(uniq.chr.names)
      uniq.chr.ints <- 1:length(uniq.chr.names)
      for(i in seq_along(uniq.chr.ints)){
        idx <- which(output$original == uniq.chr.names[i])
        output$renamed[idx] <- uniq.chr.ints[i]
      }
    }
  }

  stopifnot(length(unique(output$original)) ==
            length(unique(output$renamed)))
  stopifnot(sort(as.vector(table(output$original))) ==
            sort(as.vector(table(output$renamed))))

  return(output)
}

##' Fasta file
##'
##' Return key aspects of a fasta file as summary (number of records, length of each record, letter frequency of each record, ...).
##' @param fa.file path to the fasta file (can be gzip-compressed)
##' @param letters vector of characters to get their frequency, or NULL to skip
##' @param algo name of the algorithm used to apply a cryptographical hash function on the sequence of each record (md5/sha1/NULL)
##' @param verbose verbosity level (0/1)
##' @return list
##' @author Timothee Flutre
##' @export
summaryFasta <- function(fa.file, letters=c("A","T","G","C","N"), algo=NULL,
                         verbose=0){
  requireNamespaces(c("Biostrings", "BiocGenerics"))
  stopifnot(file.exists(fa.file))
  if(! is.null(algo)){
    requireNamespace("digest")
    stopifnot(algo %in% c("md5", "sha1"))
  }

  out <- list()

  if(verbose > 0)
    write(paste0("read ", fa.file, " ..."), stdout()); flush(stdout())
  records <- Biostrings::readDNAStringSet(filepath=fa.file, format="fasta")
  out$nb.records <- length(records)
  out$rec.lengths <- stats::setNames(BiocGenerics::width(records),
                                     names(records))
  if(verbose > 0){
    write(paste0("nb of records: ", length(records)), stdout())
    flush(stdout())
  }

  out$letter.freqs <- NULL
  if(! is.null(letters)){
    if(verbose > 0)
      write("calculate letter frequencies ...", stdout()); flush(stdout())
    lf <- Biostrings::letterFrequency(x=records, letters=letters)
    rownames(lf) <- names(records)
    out$letter.freqs <- lf
  }

  if(! is.null(algo)){
    if(verbose > 0){
      write(paste0("apply ", algo, " per record ..."), stdout())
      flush(stdout())
    }
    out[[algo]] <- sapply(records, digest::digest, algo=algo)
  }

  return(out)
}

##' Fasta file
##'
##' Extract a subset of (possibly several) sequences present in fasta file.
##' @param in.fa path to the input fasta file (can be gzipped)
##' @param sub.info information on the subset to extract, as a data.frame with 3 columns, "seq", "start", "end", where the coordinates are 1-based (a 4th column, "name", can specify the names of the resulting subsequences; by default, these are \code{seq:start-end})
##' @param out.fa path to the output fasta file (can be gzipped); if NULL, the subset of sequences isn't written to a file but simply returned
##' @param split.names if not NULL, specify the character at which the sequence names in \code{in.fa} will be split to only keep the first element
##' @param verbose verbosity level (0/1)
##' @return subset of sequences (invisible if \code{out.fa} isn't NULL)
##' @author Timothee Flutre
##' @export
extractFasta <- function(in.fa, sub.info, out.fa, split.names=" ", verbose=1){
  requireNamespaces(c("Biostrings", "GenomicRanges", "S4Vectors", "IRanges",
                      "GenomeInfoDb", "BSgenome"))
  stopifnot(file.exists(in.fa),
            is.data.frame(sub.info),
            ncol(sub.info) >= 3,
            all(c("seq", "start", "end") %in% colnames(sub.info)),
            all(sub.info$start > 0),
            all(sub.info$end > 0),
            all(sub.info$start <= sub.info$end))
  if(! is.null(out.fa))
    stopifnot(! file.exists(out.fa))
  if(is.factor(sub.info$seq))
    sub.info$seq <- as.character(sub.info$seq)
  if(! is.numeric(sub.info$start))
    sub.info$start <- as.numeric(sub.info$start)
  if(! is.numeric(sub.info$end))
    sub.info$end <- as.numeric(sub.info$end)

  ## read the input fasta file
  if(verbose > 0){
    msg <- paste0("read '", in.fa, "' ...")
    write(msg, stdout())
    utils::flush.console()
  }
  records <- Biostrings::readDNAStringSet(filepath=in.fa, format="fasta")
  if(verbose > 0){
    msg <- paste0("nb of records: ", length(records))
    write(msg, stdout())
    utils::flush.console()
  }
  if(! is.null(split.names))
    names(records) <- sapply(strsplit(names(records), split.names),
                             function(x){x[1]})

  ## if(verbose > 0){
  ##   msg <- "convert subsequence coordinates into a GRanges ..."
  ##   write(msg, stdout())
  ## }
  ## idx <- which(sub.info$seq %in% names(records))
  ## if(length(idx) > 0){
  ##   sub.info.gr <- GenomicRanges::GRanges(
  ##       seqnames=S4Vectors::Rle(sub.info$seq[idx]),
  ##       ranges=IRanges::IRanges(
  ##           start=sub.info$start[idx],
  ##           end=sub.info$end[idx]))
  ##   if("name" %in% colnames(sub.info)){
  ##     names(sub.info.gr) <- sub.info$name[idx]
  ##   } else
  ##     names(sub.info.gr) <- paste0(sub.info$seq[idx], ":", sub.info$start[idx],
  ##                                  "-", sub.info$end[idx])
  ## }
  if(verbose > 0){
    msg <- "convert subsequence coordinates into a RangesList/IRangesList..."
    write(msg, stdout())
    utils::flush.console()
  }
  if(unlist(getRversion())[2] < 5){ # assume major R version is 3
    sub.info.rl <- IRanges::RangesList()
  } else
    sub.info.rl <- IRanges::IRangesList()
  if(any(sub.info$seq %in% names(records))){
    seq.names <- unique(sub.info$seq)
    for(i in 1:length(seq.names)){
      seq.name <- seq.names[i]
      idx <- which(sub.info$seq == seq.name)
      r <- IRanges::IRanges(start=sub.info$start[idx],
                            end=sub.info$end[idx])
      if("name" %in% colnames(sub.info)){
        names(r) <- sub.info$name[idx]
      } else
        names(r) <- paste0(sub.info$seq[idx], ":", sub.info$start[idx],
                           "-", sub.info$end[idx])
      sub.info.rl[[seq.name]] <- r
    }
  }

  ## if(length(sub.info.gr) > 0){
  ##   if(verbose > 0){
  ##     msg <- paste0("nb of sequences to extract from: ",
  ##                   nlevels(GenomeInfoDb::seqnames(sub.info.gr)),
  ##                   "\nnb of subsequences to extract: ", length(sub.info.gr))
  ##     write(msg, stdout())
  ##   }
  if(length(sub.info.rl) > 0){
    if(verbose > 0){
      msg <- paste0("nb of sequences to extract from: ", length(sub.info.rl),
                    "\nnb of subsequences to extract: ",
                    sum(sapply(sub.info.rl, length)))
      write(msg, stdout())
      utils::flush.console()
    }

    ## perform the subsequence extractions
    ## ## sub.records <- Biostrings::getSeq(x=records, sub.info.gr)
    ## sub.records <- Biostrings::DNAStringSet()
    ## seq.names <- unique(as.character(GenomeInfoDb::seqnames(sub.info.gr)))
    ## for(i in 1:length(seq.names)){
    ##   seq.name <- seq.names[i]
    ##   sub.rec <- Biostrings::getSeq(x=records[seq.name],
    ##                                 sub.info.gr[seqnames(sub.info.gr) == seq.name])
    ##   sub.records[[seq.name]] <- sub.rec
    ## }
    sub.records <- Biostrings::DNAStringSet()
    seq.names <- names(sub.info.rl)
    for(i in 1:length(seq.names)){
      seq.name <- seq.names[i]
      sub.rec <- Biostrings::extractAt(x=records[seq.name],
                                       at=sub.info.rl[[seq.name]])[[1]]
      sub.records <- c(sub.records, sub.rec)
    }

    if(! is.null(out.fa)){
      if(verbose > 0){
        msg <- paste0("write '", out.fa, "' ...")
        write(msg, stdout())
        utils::flush.console()
      }
      comp <- ifelse(tail(strsplit(out.fa, "\\.")[[1]], 1) == "gz", TRUE, FALSE)
      Biostrings::writeXStringSet(x=sub.records, filepath=out.fa,
                                  compress=comp)
      invisible(sub.records)
    } else
      return(sub.records)
  }
}

##' Calculate the GC content of a set of sequences
##'
##' Requires the Biostrings package.
##' @param x vector of sequences (e.g. "AGGT"), possibly with names
##' @return vector
##' @author Timothee Flutre
##' @export
gcContent <- function(x){
  requireNamespaces("Biostrings")
  stopifnot(is.character(x))

  sapply(x, function(xi){
    sum(Biostrings::alphabetFrequency(Biostrings::DNAString(xi),
                                      baseOnly=TRUE, as.prob=TRUE)[c("C","G")])
  })
}

##' Align each sequence against each other (and itself)
##'
##' Requires the Biostrings package.
##' @title All pairwise alignments
##' @param x vector of sequences (e.g. "AGGT"), possibly with names
##' @param type type of alignment ("global", i.e. Needleman-Wunsch)
##' @param ... arguments to be passed to Biostrings::pairwiseAlignment()
##' @return list of instances of class PairwiseAlignments
##' @author Timothee Flutre
##' @export
allPairAligns <- function(x, type="global", ...){
  requireNamespaces("Biostrings")
	stopifnot(is.character(x))

	aligns <- list()

	for(i in 1:(length(x)-1)){
		for(j in i:length(x)){
			aligns[[paste0(i,"-",j)]] <-
				Biostrings::pairwiseAlignment(
            pattern=x[j],
            subject=x[i],
            type=type,
            substitutionMatrix=
              Biostrings::nucleotideSubstitutionMatrix(match=1,
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
##' @param aligns list of instances of class PairwiseAlignments (see \code{\link{allPairAligns}})
##' @param nb.sequences number of sequences
##' @return list of matrices
##' @author Timothee Flutre
##' @export
statsAllPairAligns <- function(aligns, nb.sequences){
  requireNamespaces("Biostrings")
	stopifnot(is.list(aligns))

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
			scores[i,j] <- Biostrings::score(aligns[[paste0(i,"-",j)]])
			dists[i,j] <- Biostrings::nedit(aligns[[paste0(i,"-",j)]])
			pids[i,j] <- Biostrings::pid(aligns[[paste0(i,"-",j)]])
			nmatchs[i,j] <- Biostrings::nmatch(aligns[[paste0(i,"-",j)]])
			nmismatchs[i,j] <- Biostrings::nmismatch(aligns[[paste0(i,"-",j)]])
			ninss[i,j] <- sum(Biostrings::nindel(aligns[[paste0(i,"-",j)]])@insertion[,"Length"] != 0) # insertions in pattern wrt subject
			ndels[i,j] <- sum(Biostrings::nindel(aligns[[paste0(i,"-",j)]])@deletion[,"Length"] != 0) # idem
		}
	}

	return(list(scores=scores, dists=dists, pids=pids, nmatchs=nmatchs,
							nmismatchs=nmismatchs, ninss=ninss, ndels=ndels))
}

##' Convert data frame to GRanges
##'
##' Convert a data frame of genomic coordinates to GRanges.
##' @param x data frame of genomic coordinates
##' @param seq name of the column containing the sequence of the intervals
##' @param start name of the column containing the start of the intervals
##' @param end name of the column containing the end of the intervals
##' @param strand name of the column containing the strand of the intervals
##' @param si output from the "seqinfo" function from the GenomeInfoDb package
##' @param keep.all.seqlevels if TRUE, all sequence levels from \code{si} are kept, otherwise, only those present in \code{x}
##' @param names2use if NULL, the output will have no names; if "row", it will have the row names of the input (if any); else, it will have the content of the specified column (if it exists)
##' @return GRanges
##' @author Timothee Flutre
##' @export
df2gr <- function(x, seq="chr", start="start", end="end", strand=NULL,
                  si=NULL, keep.all.seqlevels=FALSE, names2use="row"){
  requireNamespace("GenomicRanges")
  requireNamespace("S4Vectors")
  requireNamespace("IRanges")
  requireNamespace("GenomeInfoDb")
  stopifnot(all(c(seq, start, end, strand) %in% colnames(x)),
            is.logical(keep.all.seqlevels))
  if(! is.null(si))
    stopifnot(class(si) == "Seqinfo")

  if(is.null(strand)){
    out.gr <-
      GenomicRanges::GRanges(seqnames=S4Vectors::Rle(x[[seq]]),
                             ranges=IRanges::IRanges(start=x[[start]],
                                                     end=x[[end]]))
  } else
    out.gr <-
      GenomicRanges::GRanges(seqnames=S4Vectors::Rle(x[[seq]]),
                             ranges=IRanges::IRanges(start=x[[start]],
                                                     end=x[[end]]),
                             strand=S4Vectors::Rle(x[[strand]]))
  if(! is.null(names2use)){
    if(names2use == "row"){
      names(out.gr) <- rownames(x)
    } else if(names2use %in% colnames(x))
      names(out.gr) <- x[[names2use]]
  }
  out.gr <- GenomeInfoDb::sortSeqlevels(out.gr)

  if(any(! colnames(x) %in% c(seq, start, end, strand))){
    for(coln in colnames(x))
      if(! coln %in% c(seq, start, end, strand))
        S4Vectors::mcols(out.gr)[[coln]] <- x[[coln]]
  }

  if(! is.null(si)){
    GenomeInfoDb::seqlevels(si) <-
      GenomeInfoDb::sortSeqlevels(GenomeInfoDb::seqlevels(si))
    if(length(GenomeInfoDb::seqlevels(out.gr)) ==
       length(GenomeInfoDb::seqlevels(si))){
      stopifnot(all(GenomeInfoDb::seqlevels(out.gr) ==
                    GenomeInfoDb::seqlevels(si)))
    } else{
      stopifnot(all(GenomeInfoDb::seqlevels(out.gr) %in%
                    GenomeInfoDb::seqlevels(si)))
      if(keep.all.seqlevels){
        GenomeInfoDb::seqlevels(out.gr, pruning.mode="coarse") <-
          GenomeInfoDb::seqlevels(si)
      } else
        si <- GenomeInfoDb::keepSeqlevels(x=si,
                                          value=GenomeInfoDb::seqlevels(out.gr))
    }
    GenomeInfoDb::seqinfo(out.gr) <- si
  }

  return(out.gr)
}

##' Merge overlapping genomic intervals
##'
##' Merges overlapping genomic intervals (as proposed \href{https://support.bioconductor.org/p/68021/}{here}).
##' @param x input \code{GRanges}
##' @param thresh threshold of relative overlap above which \code{GRanges} are merged
##' @return \code{GRanges}
##' @author Herve Pages [aut], Timothee Flutre [ctb]
##' @export
mergeOverlaps <- function(x, thresh=0.0){
  requireNamespace("methods")
  requireNamespace("S4Vectors")
  requireNamespace("BiocGenerics")
  requireNamespace("IRanges")
  requireNamespace("GenomicRanges")

  hits <- GenomicRanges::findOverlaps(x)
  if(thresh >= 0){
    gr.qry <- x[S4Vectors::queryHits(hits)]
    gr.sbj <- x[S4Vectors::subjectHits(hits)]
    relative_overlap <-
      BiocGenerics::width(GenomicRanges::pintersect(gr.qry, gr.sbj)) /
      BiocGenerics::pmin(BiocGenerics::width(gr.qry),
                         BiocGenerics::width(gr.sbj))
    hits <- hits[relative_overlap >= thresh]
  }

  extractClustersFromSelfHits <- function(hits){
    stopifnot(methods::is(hits, "Hits"))
    stopifnot(S4Vectors::queryLength(hits) == S4Vectors::subjectLength(hits))
    hits <- BiocGenerics::union(hits, t(hits))
    qh <- S4Vectors::queryHits(hits)
    sh <- S4Vectors::subjectHits(hits)
    cid <- seq_len(S4Vectors::queryLength(hits))  # cluster IDs
    while(TRUE){
      h <- S4Vectors::Hits(qh, cid[sh],
                           S4Vectors::queryLength(hits),
                           S4Vectors::subjectLength(hits))
      cid2 <- BiocGenerics::pmin(cid, S4Vectors::selectHits(h, "first"))
      if (identical(cid2, cid))
        break
      cid <- cid2
    }
    unname(IRanges::splitAsList(seq_len(S4Vectors::queryLength(hits)), cid))
  }
  clusters <- extractClustersFromSelfHits(hits)

  gr.out <- range(IRanges::extractList(x, clusters))
  if(any(S4Vectors::elementNROWS(gr.out) != 1L))
    stop("some connected ranges are not on the same ",
         "chromosome and strand, and thus cannot be ",
         "merged")
  gr.out <- BiocGenerics::unlist(gr.out)
  S4Vectors::mcols(gr.out)$revmap <- clusters
  if(! is.null(names(x)))
    S4Vectors::mcols(gr.out)$revnames <-
                             IRanges::CharacterList(lapply(clusters, function(clu){
                               paste(names(x)[clu], collapse=", ")
                             }))

  return(gr.out)
}

##' Genomic bins
##'
##' Summarize a given metric per \code{GRanges} over genomic bins.
##' @param gr \code{GRanges} object
##' @param colname character specifying a given column in \code{mcols(gr)}; if NULL, the number of ranges per bin will be counted
##' @param binwidth fixed width of each bin
##' @param which.summary names of the summaryzing function (one or several among \code{c("sum", "mean", "min", "max")})
##' @param plot.it if TRUE, a karyogram will be plotted
##' @param ... other arguments passed to \code{autoplot}, such as \code{main}
##' @return invisible \code{GRanges} containing the bins with the summary
##' @author Timothee Flutre
##' @examples
##' \dontrun{## dummy data
##' set.seed(1)
##' C <- 2
##' P <- C * 10^3
##' df <- data.frame(CHR=c(rep("chr1", 10^3), rep("chr2", 10^3)),
##'                  SNP=paste0("snp", 1:P),
##'                  POS=NA,
##'                  N=abs(rnorm(P)))
##' df$POS[df$CHR == "chr1"] <- sort(sample.int(10^5, 10^3))
##' df$POS[df$CHR == "chr2"] <- sort(sample.int(10^5, 10^3))
##' head(df)
##' str(df)
##'
##' library(GenomicRanges)
##' df.gr <- GRanges(seqnames=df$CHR, ranges=IRanges(start=df$POS, width=1),
##'                  metric=df$N,
##'                  seqlengths=c(chr1=10^5, chr2=10^5))
##' df.gr
##'
##' (bins <- grSummaryPerBin(gr=df.gr, colname=NULL, binwidth=200))
##' (bins <- grSummaryPerBin(gr=df.gr, colname="metric", binwidth=200,
##'                          which.summary="sum"))
##' (bins <- grSummaryPerBin(gr=df.gr, colname="metric", binwidth=200,
##'                          which.summary="sum", which.plot="sum", plot.it=TRUE))
##' }
##' @export
grSummaryPerBin <- function(gr, colname=NULL, binwidth=200,
                            which.summary="sum",
                            plot.it=FALSE, ...){
  requireNamespaces(c("S4Vectors", "GenomicRanges", "IRanges",
                      "GenomeInfoDb", "BiocGenerics"))
  if(! is.null(colname))
    stopifnot(colname %in% colnames(S4Vectors::mcols(gr)))
  stopifnot(all(which.summary %in% c("sum", "mean", "min", "max")))
  if(plot.it)
    requireNamespaces("ggbio")

  bins <- GenomicRanges::tileGenome(seqlengths=GenomeInfoDb::seqlengths(gr),
                                    tilewidth=binwidth,
                                    cut.last.tile.in.chrom=TRUE)
  if(is.null(colname)){
    weight <- 1
    colname <- "count"
  } else
    weight <- colname
  cvg <- GenomicRanges::coverage(gr, weight=weight)
  list.views <- IRanges::RleViewsList(
    lapply(names(cvg), function(seqname){
      if(utils::packageVersion("GenomeInfoDb") < "1.14"){
        tmp <- IRanges::ranges(GenomeInfoDb::keepSeqlevels(bins,
                                                           seqname))
      } else
        tmp <- IRanges::ranges(GenomeInfoDb::keepSeqlevels(bins,
                                                           seqname,
                                                           "coarse"))
      IRanges::Views(cvg[[seqname]], tmp)
    }))
  if("sum" %in% which.summary)
    S4Vectors::mcols(bins)[[paste0(colname, ".sum")]] <-
      BiocGenerics::unlist(IRanges::viewSums(list.views))
  if("mean" %in% which.summary)
    S4Vectors::mcols(bins)[[paste0(colname, ".mean")]] <-
      BiocGenerics::unlist(IRanges::viewMeans(list.views))
  if("min" %in% which.summary)
    S4Vectors::mcols(bins)[[paste0(colname, ".min")]] <-
      BiocGenerics::unlist(IRanges::viewMins(list.views))
  if("max" %in% which.summary)
    S4Vectors::mcols(bins)[[paste0(colname, ".max")]] <-
      BiocGenerics::unlist(IRanges::viewMeans(list.views))

  if(plot.it)
    print(
        ggbio::autoplot(bins, layout="karyogram",
                        ggplot2::aes_string(color=paste0(colname, ".",
                                                         which.summary),
                                            fill=paste0(colname, ".",
                                                        which.summary)),
                        ...))

  invisible(bins)
}

##' Load read counts
##'
##' Load read counts per individual and lane
##' @param lanes.dir vector of paths to "lane" directories, each containing a "demultiplex" directory in which \href{https://github.com/timflutre/quantgen/blob/master/demultiplex.py}{demultiplex.py} was executed
##' @return data frame with (individual,lane) in rows and counts in columns
##' @author Timothee Flutre
##' @export
loadReadCountsPerIndAndLane <- function(lanes.dir){
  stopifnot(is.vector(lanes.dir),
            length(lanes.dir) > 0,
            is.character(lanes.dir),
            all(file.exists(lanes.dir)))

  reads <- list()

  for(lane.dir in lanes.dir){
    lane <- gsub("lane_", "", basename(lane.dir))
    lane.file <- paste0(lane.dir, "/demultiplex/",
                        lane, "_stats-demultiplex.txt.gz")
    if(! file.exists(lane.file))
      stop(paste0("file '", lane.file, "' doesn't exist"))
    reads[[lane]] <- utils::read.table(lane.file, header=TRUE)
    if(! all(colnames(reads[[lane]] %in% c("ind","barcode","assigned","lane"))))
      warning(paste0("look at header line of file '", lane.file, "'"))
  }

  reads <- do.call(rbind, reads)
  reads$flowcell <- as.factor(sapply(strsplit(rownames(reads), "_"),
                                     function(x){x[1]}))
  reads$lane <- as.factor(sapply(strsplit(rownames(reads), "\\."),
                                 function(x){x[1]}))

  return(reads)
}

##' Reformat read counts per lane
##'
##' Reformat as a matrix of counts a data frame formatted, for instance, as the output of \code{\link{loadReadCountsPerIndAndLane}}
##' @param x data frame with columns "lane", "ind" and "assigned"
##' @return matrix with lanes in rows and individuals in columns
##' @author Timothee Flutre
##' @export
formatReadCountsPerLane <- function(x){
  stopifnot(is.data.frame(x),
            "lane" %in% colnames(x),
            "ind" %in% colnames(x),
            "assigned" %in% colnames(x))

  counts <- matrix(data=0, nrow=length(unique(x$lane)),
                   ncol=length(unique(x$ind)),
                   dimnames=list(unique(x$lane), unique(x$ind)))

  for(j in 1:ncol(counts)){
    ind <- colnames(counts)[j]
    for(lane in rownames(counts)){
      if(ind %in% x$ind[x$lane == lane]){
        stopifnot(length(x$assigned[x$ind == ind & x$lane == lane]) == 1)
        counts[lane, ind] <- x$assigned[x$ind == ind &
                                        x$lane == lane]
      }
    }
  }
  counts <- counts[, order(apply(counts, 2, sum))]

  return(counts)
}

##' Stacked barplot of read counts
##'
##' Make a stacked barplot, for instance from the output of \code{\link{formatReadCountsPerLane}}. Default graphical parameters are optimized for few categories (e.g. lanes) but many individuals.
##' @param counts matrix, with categories (e.g. lanes) in rows and individuals in columns; if colnames is not NULL, they will be used as labels of the x-axis
##' @param xlab a title for the x axis
##' @param ylab a title for the y axis
##' @param main an overall title for the plot
##' @param legend.text see \code{\link[graphics]{barplot}}
##' @param col see \code{\link[graphics]{barplot}}
##' @param lines.h the y-value(s) for horizontal line(s)
##' @param mar see \code{\link[graphics]{par}}
##' @param xlab.line see \code{line} from \code{\link[graphics]{mtext}}
##' @param font.lab see \code{\link[graphics]{par}}
##' @param cex.xtxt see \code{cex} from \code{\link[graphics]{text}}
##' @param perc if TRUE, plot data from "counts" as percentages (using \code{\link[base]{prop.table}}), with legend on the topright outside the x-axis range, and columns (of "counts") sorted in the increasing order corresponding of the first row
##' @param inset.x inset distance along the x axis from the margins as a fraction of the plot region
##' @return invisible output of \code{\link[graphics]{barplot}}
##' @author Timothee Flutre
##' @examples
##' \dontrun{## make fake data
##' (x <- matrix(data=c(10,3, 20,0, 21,17, 35,19), nrow=2, ncol=4,
##'              dimnames=list(c("lane1", "lane2"),
##'              c("ind1", "ind2", "ind3", "ind4"))))
##'
##' ## plot them
##' barplotReadCounts(counts=x)
##' barplotReadCounts(counts=x, perc=TRUE, ylab="Percentage of read pairs",
##'                   mar=c(5, 4.5, 4, 6))
##' }
##' @export
barplotReadCounts <- function(counts,
                              xlab=paste0("Individuals (", ncol(counts), ")"),
                              ylab="Number of read pairs",
                              main="Read pairs per individual",
                              legend.text=rownames(counts),
                              col=1:nrow(counts),
                              lines.h=c(10^5, 10^6, 5*10^6),
                              mar=c(5, 4.5, 4, 0),
                              xlab.line=3,
                              font.lab=2,
                              cex.xtxt=0.7,
                              perc=FALSE,
                              inset.x=-0.2){
  stopifnot(is.matrix(counts),
            is.logical(perc))

  if(! is.null(mar)){
    stopifnot(is.numeric(mar),
              length(mar) == 4)
    oldpar <- graphics::par(mar=mar)
  }

  if(! perc){
    height <- counts
    idx <- 1:ncol(height)
    bp <- graphics::barplot(height=height,
                            width=1,
                            col=col,
                            border=NA,
                            xaxt="n",
                            xlab="",
                            ylab=ylab,
                            main=main,
                            font.lab=font.lab,
                            legend.text=legend.text,
                            args.legend=list(x="topleft", bty="n",
                                             fill=rev(col),
                                             border=rev(col)))
  } else{
    height <- 100 * prop.table(counts, 2)
    idx <- order(height[1,])
    bp <- graphics::barplot(height=height[, idx],
                            width=1,
                            col=col,
                            border=NA,
                            las=1,
                            xaxt="n",
                            xlab="",
                            ylab=ylab,
                            main=main,
                            font.lab=font.lab,
                            legend.text=legend.text,
                            args.legend=list(x="topright", bty="n",
                                             xpd=TRUE, inset=c(inset.x, 0),
                                             fill=rev(col),
                                             border=rev(col)))
  }

  graphics::axis(1, at=bp, labels=FALSE)
  if(! is.null(colnames(counts)))
    graphics::text(x=bp, y=graphics::par("usr")[3], srt=45, adj=c(1.2,2),
                   labels=colnames(counts)[idx], xpd=TRUE, cex=cex.xtxt)
  graphics::mtext(text=xlab, side=1, line=xlab.line, font=font.lab)

  if(! is.null(lines.h)){
    for(h in lines.h)
      graphics::abline(h=h, lty=2)
  }

  graphics::par(oldpar)

  invisible(bp)
}

##' Coverage
##'
##' Calculate the breadth (fraction of the genome that is covered by reads) and depth (average number of times each position in the genome has been sequenced) of coverage from a set of BAM file(s) (see \href{http://dx.doi.org/10.1093/bib/bbu029}{Molnar & Ilie (2015)}). Results are identical to those from \code{samtools depth}.
##' @param bamFiles vector of paths to BAM file(s), each of them having an index file with the same name but finishing by ".bai"; if there are several BAM files, all of them should contain reads aligned on the same set of sequences (i.e. the same reference genome)
##' @param yieldSize number of records to yield each time the file is read from (see \code{\link[Rsamtools]{BamFile}})
##' @param seq.ids sequence identifier(s) to focus on, e.g. "chr2" or c("chr2","chr5"); if NULL, all of them
##' @param out.file if not NULL, the output will be transformed into a data.frame and written into the given file
##' @param max.depth see \code{\link[Rsamtools]{PileupParam}}
##' @param min.base.quality see \code{\link[Rsamtools]{PileupParam}}
##' @param min.mapq see \code{\link[Rsamtools]{PileupParam}}
##' @param min.nucleotide.depth see \code{\link[Rsamtools]{PileupParam}}
##' @param distinguish.strands see \code{\link[Rsamtools]{PileupParam}}
##' @param distinguish.nucleotides see \code{\link[Rsamtools]{PileupParam}}
##' @param ignore.query.Ns see \code{\link[Rsamtools]{PileupParam}}
##' @param include.deletions see \code{\link[Rsamtools]{PileupParam}}
##' @param include.insertions see \code{\link[Rsamtools]{PileupParam}}
##' @param is.paired see \code{\link[Rsamtools]{scanBamFlag}}
##' @param is.proper.pair see \code{\link[Rsamtools]{scanBamFlag}}
##' @param is.unmapped.query see \code{\link[Rsamtools]{scanBamFlag}}
##' @param has.unmapped.mate see \code{\link[Rsamtools]{scanBamFlag}}
##' @param is.minus.strand see \code{\link[Rsamtools]{scanBamFlag}}
##' @param is.mate.minus.strand see \code{\link[Rsamtools]{scanBamFlag}}
##' @param is.first.mate.read see \code{\link[Rsamtools]{scanBamFlag}}
##' @param is.second.mate.read see \code{\link[Rsamtools]{scanBamFlag}}
##' @param is.secondary.align see \code{\link[Rsamtools]{scanBamFlag}}
##' @param is.not.passing.qc see \code{\link[Rsamtools]{scanBamFlag}}
##' @param is.dupl see \code{\link[Rsamtools]{scanBamFlag}}
##' @param verbose verbosity level (0/1/2)
##' @param nb.cores the number of cores to use to parallelize over bamFiles, i.e. at most how many child processes will be run simultaneously (not on Windows)
##' @return list with one component per BAM file, each being a matrix with columns "nb.bases", "mapped.bases", "count.sum", "breadth", "depth" and "depth.map" (corresponding to the depth at positions with at least one mapped read)
##' @seealso \code{\link{freadBedtoolsCoverageHist}}
##' @author Timothee Flutre (with the help of Martin Morgan)
##' @export
coverageBams <- function(bamFiles, yieldSize=10^4, seq.ids=NULL, out.file=NULL,
                         max.depth=10^4, min.base.quality=1,
                         min.mapq=5, min.nucleotide.depth=1,
                         distinguish.strands=FALSE,
                         distinguish.nucleotides=FALSE, ignore.query.Ns=FALSE,
                         include.deletions=FALSE, include.insertions=FALSE,
                         is.paired=NA, is.proper.pair=NA,
                         is.unmapped.query=FALSE, has.unmapped.mate=NA,
                         is.minus.strand=NA, is.mate.minus.strand=NA,
                         is.first.mate.read=NA, is.second.mate.read=NA,
                         is.secondary.align=FALSE, is.not.passing.qc=FALSE,
                         is.dupl=FALSE,
                         verbose=1, nb.cores=1){
  requireNamespaces(c("parallel", "IRanges", "Rsamtools"))
  stopifnot(is.character(bamFiles),
            all(sapply(bamFiles, file.exists)))
  if(! is.null(out.file))
    stopifnot(! file.exists(out.file))

  if(verbose > 0)
    write("extract sequence lengths ...", stdout()); flush(stdout())
  tmp <- Rsamtools::scanBamHeader(files=bamFiles, what="targets")
  tmp <- lapply(tmp, function(x){x$targets})
  seq.lengths <- tmp[[1]]
  if(length(tmp) > 1){
    for(i in 2:length(tmp))
      if(! identical(tmp[[i]], seq.lengths))
        stop(paste0(bamFiles[i], " has different target sequences than ",
                    bamFiles[1]))
  }
  if(! is.null(seq.ids)){
    for(seq.id in seq.ids)
      if(! seq.id %in% names(seq.lengths))
        stop(paste0(seq.id, " not in ", bamFiles[1]))
    seq.lengths <- seq.lengths[seq.ids]
  }
  if(verbose > 1)
    print(seq.lengths)
  rl <- IRanges::RangesList()
  for(seq.id in names(seq.lengths))
    rl[[seq.id]] <- IRanges::IRanges(start=1, end=seq.lengths[seq.id],
                                     names=seq.id)

  if(verbose > 0){
    msg <- paste0("compute coverage from ", length(bamFiles), " BAM files ...")
    write(msg, stdout()); flush(stdout())
  }
  fl <- Rsamtools::scanBamFlag(isPaired=is.paired,
                               isProperPair=is.proper.pair,
                               isUnmappedQuery=is.unmapped.query,
                               hasUnmappedMate=has.unmapped.mate,
                               isMinusStrand=is.minus.strand,
                               isMateMinusStrand=is.mate.minus.strand,
                               isFirstMateRead=is.first.mate.read,
                               isSecondMateRead=is.second.mate.read,
                               isSecondaryAlignment=is.secondary.align,
                               isNotPassingQualityControls=is.not.passing.qc,
                               isDuplicate=is.dupl)
  sbp <- Rsamtools::ScanBamParam(which=rl, flag=fl)
  pp <- Rsamtools::PileupParam(max_depth=max.depth,
                               min_base_quality=min.base.quality,
                               min_mapq=min.mapq,
                               min_nucleotide_depth=min.nucleotide.depth,
                               distinguish_strands=distinguish.strands,
                               distinguish_nucleotides=distinguish.nucleotides,
                               ignore_query_Ns=ignore.query.Ns,
                               include_deletions=include.deletions,
                               include_insertions=include.insertions)
  out <- parallel::mclapply(bamFiles, function(bamFile){
    if(verbose > 0)
      write(bamFile, stdout()); flush(stdout())
    ## cmd <- paste0("samtools depth -Q ", min.mapq, " ", bamFile)
    ## cvg <- data.table::fread(input=cmd,
    ##                          col.names=c("chr", "pos", "count"))
    ## tmp <- cbind(nb.bases=seq.lengths,
    ##              mapped.bases=cvg[, length(count), chr][, V1],
    ##              count.sum=cvg[, sum(count), chr][, V1])
    p <- Rsamtools::pileup(file=bamFile, yieldSize=yieldSize,
                           scanBamParam=sbp, pileupParam=pp)
    p$seqnames <- p$seqnames[, drop=TRUE]
    tmp <- cbind(nb.bases=seq.lengths,
                 mapped.bases=tapply(p$count, p$seqnames, length),
                 count.sum=tapply(p$count, p$seqnames, sum))
    cbind(tmp,
          breadth=tmp[,"mapped.bases"] / tmp[,"nb.bases"],
          depth=tmp[,"count.sum"] / tmp[,"nb.bases"],
          depth.map=tmp[,"count.sum"] / tmp[,"mapped.bases"])
  }, mc.cores=nb.cores, mc.silent=ifelse(nb.cores == 1, FALSE, TRUE))
  names(out) <- basename(bamFiles)

  if(! is.null(out.file)){
    if(verbose > 0)
      write("write results into file ...", stdout())
    tmp <- do.call(rbind, lapply(names(out), function(bam.file){
      cbind(bam.file, seq=rownames(out[[bam.file]]), out[[bam.file]])
    }))
    utils::write.table(x=tmp, file=gzfile(out.file), quote=FALSE, sep="\t",
                       row.names=FALSE, col.names=TRUE)
  }

  return(out)
}

##' Bar plot of insert sizes
##'
##' Creates a bar plot with vertical bars of insert sizes from histogram data as calculated by the CollectInsertSizeMetrics program from the \href{https://broadinstitute.github.io/picard/}{Picard} software.
##' @param file path to the output file from \code{java -jar picard.jar CollectInsertSizeMetrics}
##' @param main overall title for the plot
##' @param add.text add total count, as well as Q25, median and Q75 for insert sizes, to the topright of the plot
##' @return invisible data frame of the content of the file
##' @author Timothee Flutre
##' @export
barplotInsertSizes <- function(file, main=NULL, add.text=FALSE){
  stopifnot(file.exists(file))

  dat <- utils::read.table(file=file, skip=10, header=TRUE)
  if(! ncol(dat) == 2)
    stop(paste0("file ", file, " doesn't seem to come from",
                " Picard CollectInsertSizeMetrics"))
  colnames(dat) <- c("insert.size", "count")

  tot.count <- sum(dat$count)
  q25.insert.size <- rev(dat$insert.size[cumsum(dat$count) <=
                                         0.25 * tot.count])[1]
  med.insert.size <- rev(dat$insert.size[cumsum(dat$count) <=
                                         0.5 * tot.count])[1]
  q75.insert.size <- rev(dat$insert.size[cumsum(dat$count) <=
                                         0.75 * tot.count])[1]

  bp <- graphics::barplot(height=dat$count, width=1,
                          xlab="Insert size (in bp)",
                          ylab="Counts",
                          main=main)
  graphics::axis(side=1, at=c(0, bp[seq(100, max(dat$insert.size), 100)]),
                 labels=c(0, seq(100, max(dat$insert.size), 100)))

  graphics::abline(v=bp[q25.insert.size], lty=2)
  graphics::abline(v=bp[med.insert.size], lty=2)
  graphics::abline(v=bp[q75.insert.size], lty=2)

  if(add.text){
    graphics::text(x=bp[floor(0.6*max(dat$insert.size))],
                   y=0.6*max(dat$count), adj=0,
                   labels=paste0("Q25 = ", format(q25.insert.size, digits=2),
                                 " bp"))
    graphics::text(x=bp[floor(0.6*max(dat$insert.size))],
                   y=0.7*max(dat$count), adj=0,
                   labels=paste0("median = ", format(med.insert.size, digits=2),
                                 " bp"))
    graphics::text(x=bp[floor(0.6*max(dat$insert.size))],
                   y=0.8*max(dat$count), adj=0,
                   labels=paste0("Q75 = ", format(q75.insert.size, digits=2),
                                 " bp"))
    graphics::text(x=bp[floor(0.6*max(dat$insert.size))],
                   y=0.9*max(dat$count), adj=0,
                   labels=paste0("total count = ", format(tot.count, digits=2)))
  }

  invisible(dat)
}

##' Read bedtools-coverage-hist as data.table
##'
##' Read output files from bedtools coverage with -hist into a list of data.tables. All lines starting by "all" must have been discarded beforehand.
##' @param files character vector of relative or absolute filepaths (e.g. from Sys.glob).
##' @param verbose verbosity level (0/1)
##' @return list of data.tables
##' @seealso \code{\link{depthsPerSample}}, \code{\link{depthsPerRegion}}
##' @author Timothee Flutre
##' @export
freadBedtoolsCoverageHist <- function(files, verbose=1){
  colClasses <- sapply(utils::read.table(files[1], nrows=5), class)
  ldat <-
    lapply(1:length(files), function (i, verbose){
             if(verbose > 0){
               txt <- paste0(i, "/", length(files))
               write(txt, stdout()); flush(stdout())
             }
             f <- files[i]
             nrows <- as.numeric(system(paste("zcat", f, "| wc -l"),
                                        intern=TRUE))
             dat <- fread(paste("zcat", f),
                          sep="\t",
                          nrows=nrows,
                          colClasses=colClasses,
                          data.table=TRUE)
             setnames(x=dat, old=1:ncol(dat),
                      new=c("chrom", "start", "end", "name", "depth", "nb.bases",
                          "length", "fraction"))
             dat
           }, verbose)

  return(ldat)
}

##' Depth per sample across regions
##'
##' Summarize the depth per sample across regions of a given length range.
##' @param dat data.table (see freadBedtoolsCoverageHist)
##' @param min.reg.len minimum length of a region to be considered
##' @param max.reg.len maximum length of a region to be considered
##' @param min.reg.dep minimum depth of a region when reporting the number of interesting regions
##' @param max.reg.dep maximum depth of a region when reporting the number of interesting regions
##' @param min.reg.frac minimum fraction of a region when reporting the number of interesting regions
##' @return data.table
##' @author Timothee Flutre
##' @export
depthsPerSample <- function(dat, min.reg.len=30, max.reg.len=500,
                            min.reg.dep=10, max.reg.dep=200,
                            min.reg.frac=0.25){
  stopifnot(is.data.table(dat))
  for(col in c("ind", "flowcell", "lane", "chrom", "start", "end", "name",
               "depth", "fraction"))
    stopifnot(col %in% colnames(dat))

  ## http://stackoverflow.com/a/8096882/597069
  ind=flowcell=lane=chrom=start=end=name=depth=fraction=NULL

  setkey(dat, NULL)
  setkeyv(x=dat, cols=c("ind", "flowcell", "lane"))

  depths.sample <- dat[end - start >= min.reg.len &
                         end - start <= max.reg.len,
                       .(depth.n=.N,
                         depth.min=min(.SD[,depth]),
                         depth.med=as.double(stats::median(.SD[,depth])),
                         depth.mean=mean(.SD[,depth]),
                         depth.max=max(.SD[,depth]),
                         depth.q65=stats::quantile(.SD[,depth], 0.65),
                         depth.q70=stats::quantile(.SD[,depth], 0.70),
                         depth.q75=stats::quantile(.SD[,depth], 0.75),
                         depth.q80=stats::quantile(.SD[,depth], 0.80),
                         regions.ok=nrow(unique(.SD[depth >= min.reg.dep &
                                                      depth <= max.reg.dep &
                                                        fraction >= min.reg.frac,
                             .(chrom,start,end,name)]))),
                       by=.(ind,flowcell,lane)]

  setkey(dat, NULL)

  return(depths.sample)
}


##' Depth per region across samples
##'
##' Summarize the depth per region of a given length range across samples.
##' @param dat data.table (see freadBedtoolsCoverageHist)
##' @param min.reg.len minimum length of a region to be considered
##' @param max.reg.len maximum length of a region to be considered
##' @param min.reg.dep minimum depth of a region when reporting the number of interesting regions
##' @param max.reg.dep maximum depth of a region when reporting the number of interesting regions
##' @param min.reg.frac minimum fraction of a region when reporting the number of interesting regions
##' @return data.table
##' @author Timothee Flutre
##' @export
depthsPerRegion <- function(dat, min.reg.len=30, max.reg.len=500,
                            min.reg.dep=10, max.reg.dep=200,
                            min.reg.frac=0.25){
  stopifnot(is.data.table(dat))
  for(col in c("ind", "flowcell", "lane", "chrom", "start", "end", "name",
               "depth", "fraction"))
    stopifnot(col %in% colnames(dat))

  ## http://stackoverflow.com/a/8096882/597069
  ind=flowcell=lane=chrom=start=end=name=depth=fraction=NULL

  setkey(dat, NULL)

  depths.region <- dat[end - start >= min.reg.len &
                         end - start <= max.reg.len,
                       .(depth.n=.N,
                         depth.min=min(.SD[,depth]),
                         depth.med=as.double(stats::median(.SD[,depth])),
                         depth.max=max(.SD[,depth]),
                         samples.ok=nrow(unique(.SD[depth >= min.reg.dep &
                                                      depth <= max.reg.dep &
                                                        fraction >= min.reg.frac,
                             .(ind,flowcell,lane)]))),
                       by=.(chrom,start,end,name)]

  setkey(dat, NULL)

  return(depths.region)
}

##' Plot the covered fraction of regions as a function of depth
##'
##' Need to first run bedtools coverage as in http://www.gettinggeneticsdone.com/2014/03/visualize-coverage-exome-targeted-ngs-bedtools.html by Stephen Turner
##' @param path string
##' @param pattern character string containing a regular expression passed to gsub to be removed to simplify the file basenames used in the legend (e.g. "bedtools-coverage|\\.txt\\.gz")
##' @param covrg list of data frames
##' @param plot.it logical
##' @param xlim vector of limits for the x-axis
##' @param ylim vector of limits for the y-axis
##' @param points.type character indicating the type of plotting ("p"/"l"/"b")
##' @param plot.legend logical
##' @param verbose verbosity level (0/1)
##' @return invisible list of covrg and cumcovrg
##' @author Timothee Flutre
##' @export
coverageRegions <- function(path=NULL, pattern=NULL, covrg=NULL, plot.it=TRUE,
                            xlim=NULL, ylim=c(0,1), points.type="l",
                            plot.legend=TRUE, verbose=1){
  stopifnot(xor(is.null(path), is.null(covrg)))

  if(is.null(covrg)){
    files <- Sys.glob(path)
    if(verbose > 0){
      txt <- paste0("nb of files: ", length(files))
      write(txt, stdout())
    }
    covrg <- list()
    for(i in seq_along(files))
      covrg[[i]] <- utils::read.table(files[i], sep="\t",
                               col.names=c("chrom","depth", "nb.bases",
                                   "size","fraction"))
    if(is.null(pattern)){
      names(covrg) <- basename(files)
    } else
      names(covrg) <- gsub(pattern=pattern, replacement="", x=basename(files))
  }

  cumcovrg <- list()
  for(i in seq_along(covrg))
    cumcovrg[[i]] <- 1 - cumsum(covrg[[i]][, "fraction"])

  if(plot.it){
    if(is.null(xlim))
      xlim <- c(min(sapply(cumcovrg, min)), max(sapply(cumcovrg, max)))
    labs <- names(covrg)
    cols <- 1:length(covrg)
    graphics::plot(x=covrg[[1]][0:500,2], y=cumcovrg[[1]][0:500], type="n",
         xlab="Depth", ylab="Fraction of bases >= depth",
         main="Region coverage",
         xlim=xlim, ylim=ylim, las=1)
    for(i in seq_along(covrg))
      graphics::points(covrg[[i]][0:500, 2], cumcovrg[[i]][0:500], type=points.type, lwd=1,
             col=cols[i])
    if(plot.legend)
      graphics::legend("topright", legend=labs, col=cols, lty=1, lwd=1, bty="n")
    graphics::abline(h=0.10, lty=2, col="grey60")
    graphics::abline(h=0.50, lty=2, col="grey60")
    graphics::abline(h=0.90, lty=2, col="grey60")
  }

  invisible(list(covrg=covrg, cumcovrg=cumcovrg))
}

##' Read SAM dict
##'
##' Reads a file in the "dict" format. See \url{https://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary}.
##' @param file the name of the ".dict" file which the data are to be read from
##' @return data.frame
##' @author Timothee Flutre
##' @export
readSamDict <- function(file){
  stopifnot(file.exists(file))

  lines <- readLines(con=file)
  list.lines <- strsplit(x=lines, split="\t")

  if(length(unique(sapply(list.lines[-1], length)))  > 1){
    msg <- paste0("all @SQ in '", file,
                  "' should have the same number of tags")
    stop(msg)
  }
  df <- do.call(rbind,
                lapply(list.lines[-1], function(x){
                         return(x[-1])
                       }))

  out <- as.data.frame(
      apply(df[,-1], 2, function(x){
              sapply(strsplit(x, ":"), function(y){
                       paste(y[-1], collapse=":")
                     })
            }), stringsAsFactors=FALSE)
  colnames(out) <- sapply(strsplit(df[1,-1], ":"), function(x){x[1]})
  rownames(out) <- sapply(strsplit(df[,1], ":"), function(x){x[2]})
  out$LN <- as.numeric(out$LN)

  return(out)
}

##' Read output from bcftools counts
##'
##' Reads the output from "bcftools plugin counts"
##' @param files vector of file names to read
##' @return data frame
##' @author Timothee Flutre
##' @export
readBcftoolsCounts <- function(files){
  stopifnot(all(file.exists(files)),
            ! anyDuplicated(files))

  out <- lapply(files, function(f){
    lines <- readLines(con=f)
    idx <- grep("WARNING", lines)
    skip <- ifelse(length(idx) == 0, 0, idx)
    tmp <- utils::read.table(file=f, header=FALSE, sep=":", skip=skip,
                             stringsAsFactors=FALSE)
    stats::setNames(object=tmp[,2],
                    nm=sapply(strsplit(tmp[,1], " "), function(x){
                      x[length(x)]
                    }))
  })
  names(out) <- files
  out <- do.call(rbind, out)

  return(out)
}

##' Convert VCF file to dosage file
##'
##' Convert a VCF file to a genotype dosage file with \code{bcftools +dosage} and \code{datamash}.
##' Caution, by default, if PL is 0,3,16 and GT=./., only PL will be used; if tag=GT is specified, then the output will contain -1, which will be converted into NA by \code{\link{readGenoDoseFileFromBcftools}}.
##' @param in.vcf.file path to the input VCF file (gzip-compressed)
##' @param out.txt.file path to the output file (gzip-compressed)
##' @param tag VCF tags to determine the dosage (PL/GL/GT); if NULL, chosen in that order
##' @param save.timestamp if TRUE, the time stamp (and original name) will be saved, making the output file not bit-by-bit reproducible
##' @param nb.cores number of threads for bcftools
##' @return nothing
##' @author Timothee Flutre
##' @seealso \code{\link{readGenoDoseFileFromBcftools}}
##' @export
convertVcfToGenoDoseWithBcftools <- function(in.vcf.file, out.txt.file,
                                             tag=NULL,
                                             save.timestamp=FALSE,
                                             nb.cores=1){
  stopifnot(file.exists(Sys.which("bcftools")),
            file.exists(Sys.which("datamash")),
            file.exists(Sys.which("tr")),
            file.exists(Sys.which("gzip")),
            file.exists(in.vcf.file),
            out.txt.file != in.vcf.file)
  if(! is.null(tag))
    stopifnot(tag %in% c("PL", "GL", "GT"))

  cmd <- paste0("bcftools +dosage",
                " --threads ", nb.cores,
                " ", in.vcf.file,
                ifelse(is.null(tag), "", paste0(" -- -t ", tag)),
                " | tr ' ' '\t'",
                " | datamash transpose",
                " | gzip",
                ifelse(save.timestamp, "", " -n"),
                " > ", out.txt.file)
  system(cmd)
}

##' Read SNP genotypes as dosage from bcftools
##'
##' Read an input file with genotypes as allele dosage as written by \code{bcftools +dosage genos.vcf.gz | tr ' ' '\t' | datamash transpose | gzip > genos_dosage.tsv.gz}.
##' Caution, depending on the file size, this may require the allocation of too much RAM for your system.
##' @param genos.file path to the input file containing the SNP genotype data
##' @param get.coords if TRUE, SNP coordinates will also be extracted
##' @param get.alleles if TRUE, SNP alleles will also be extracted
##' @param use.fread if TRUE, \code{\link[data.table]{fread}} will be used to speed-up
##' @param rm.dup if TRUE, duplicated SNPs (i.e., with the exact same coordinate) will be removed
##' @param verbose verbosity level (0/1)
##' @return list with a matrix of SNP genotypes and optional data frames of SNP coordinates and alleles
##' @author Timothee Flutre
##' @seealso \code{\link{convertVcfToGenoDoseWithBcftools}}
##' @export
readGenoDoseFileFromBcftools <- function(genos.file, get.coords=TRUE,
                                         get.alleles=TRUE, use.fread=TRUE,
                                         rm.dup=FALSE, verbose=0){
  stopifnot(file.exists(genos.file),
            is.logical(get.coords),
            is.logical(get.alleles),
            is.logical(use.fread))
  if(use.fread){
    requireNamespace("data.table")
    requireNamespace("tools")
  }

  if(use.fread){
    if(tools::file_ext(genos.file) == "gz"){
      cmd <- paste0("zcat ", genos.file)
      genos <- data.table::fread(input=cmd, sep="\t", header=FALSE)
    } else
      genos <- data.table::fread(input=genos.file, sep="\t", header=FALSE)
    genos <- as.data.frame(genos)
  } else
    genos <- utils::read.table(file=genos.file, sep="\t", comment.char="",
                               stringsAsFactors=FALSE)
  stopifnot(nrow(genos) >= 5,
            genos[1,1] == "#[1]CHROM",
            genos[2,1] == "[2]POS",
            genos[3,1] == "[3]REF",
            genos[4,1] == "[4]ALT")

  colnames(genos)[-1] <- paste0(t(genos[1,-1]), ":", t(genos[2,-1]), "_",
                                t(genos[3,-1]), "/", t(genos[4,-1]))

  ind.names <- gsub("\\[[0-9]+\\]", "", genos[-c(1:4), 1])
  stopifnot(all(! duplicated(ind.names)))
  genos <- as.matrix(genos[-c(1:4), -1])
  mode(genos) <- "numeric"
  rownames(genos) <- ind.names
  is.neg <- (genos < 0)
  if(any(is.neg)){
    if(verbose > 0){
      perc <- 100 * sum(is.neg) / length(is.neg)
      txt <- paste0("negative entries converted to NA: ",
                    sum(is.neg), " / ", length(is.neg),
                    " (", round(perc, 2), "%)")
      write(txt, stdout())
    }
    genos[which(is.neg)] <- NA
  }

  coords <- NULL
  if(get.coords){
    coords <- utils::strcapture("^([0-9a-zA-Z]+):([0-9]+)_",
                                colnames(genos),
                                list(chr=character(), coord=double()))
    rownames(coords) <- colnames(genos)
  }

  alleles <- NULL
  if(get.alleles){
    alleles <- utils::strcapture("_([A-Z])/([A-Z])$",
                                 colnames(genos),
                                 list(ref=character(), alt=character()))
    rownames(alleles) <- colnames(genos)
  }

  if(all(rm.dup, get.coords)){
    tmp <- paste0(coords$chr, "_", coords$coord)
    is.dup <- duplicated(tmp)
    if(verbose > 0){
      perc <- 100 * sum(is.dup) / length(is.dup)
      txt <- paste0("duplicates: ", sum(is.dup), " / ", length(is.dup),
                    " (", round(perc, 2), "%)")
      write(txt, stdout())
    }
    if(sum(is.dup) > 0){
      coords <- coords[! is.dup,]
      genos <- genos[, ! is.dup]
      if(get.alleles)
        alleles <- alleles[! is.dup,]
    }
  }

  return(list(genos=genos, coords=coords, alleles=alleles))
}

##' VCF dimensions
##'
##' Return the dimensions of a VCF file, i.e. the number of sites (rows) and samples (columns)
##' @param vcf.file path to the VCF file (if the bgzip index doesn't exist in the same directory, it will be created)
##' @param genome genome identifier (e.g. "VITVI_12x2")
##' @param yieldSize number of records to yield each time the file is read from
##' @return vector
##' @author Timothee Flutre
##' @export
dimVcf <- function(vcf.file, genome="", yieldSize=10000){
  requireNamespaces(c("Rsamtools", "VariantAnnotation"))
  stopifnot(file.exists(vcf.file))

  out <- stats::setNames(object=c(NA, NA), nm=c("sites", "samples"))

  if(! file.exists(paste0(vcf.file, ".tbi")))
    Rsamtools::indexTabix(file=vcf.file, format="vcf")
  tabix.file <- Rsamtools::TabixFile(file=vcf.file,
                                     yieldSize=yieldSize)
  open(tabix.file)
  hdr <- VariantAnnotation::scanVcfHeader(file=tabix.file)
  out["samples"] <- length(VariantAnnotation::samples(hdr))
  out["sites"] <- Rsamtools::countTabix(tabix.file)[[1]]
  close(tabix.file)

  return(out)
}

##' Read VCF
##'
##' Read a subset of a VCF file
##' @param vcf.file path to the VCF file (bgzip index should exist in same directory)
##' @param genome genome identifier (e.g. "VITVI_12x2")
##' @param seq.id sequence identifier to work on (e.g. "chr2")
##' @param seq.start start of the sequence to work on (1-based coordinate)
##' @param seq.end end of the sequence to work on (1-based coordinate)
##' @return CollapsedVCF (see pkg \href{http://bioconductor.org/packages/VariantAnnotation/}{VariantAnnotation})
##' @author Timothee Flutre
##' @export
readVcfSubset <- function(vcf.file, genome="", seq.id, seq.start, seq.end){
  requireNamespaces(c("IRanges", "GenomicRanges", "VariantAnnotation",
                      "Rsamtools", "S4Vectors"))
  stopifnot(file.exists(vcf.file))

  tabix.file <- Rsamtools::TabixFile(file=vcf.file,
                                     yieldSize=NA_integer_)
  rngs <- GenomicRanges::GRanges(seqnames=c(seq.id),
                                 ranges=IRanges::IRanges(start=c(seq.start),
                                                         end=c(seq.end)))
  names(rngs) <- c(seq.id)
  vcf.params <- VariantAnnotation::ScanVcfParam(which=rngs)
  vcf <- VariantAnnotation::readVcf(file=tabix.file, genome=genome,
                                    param=vcf.params)
  return(vcf)
}

##' Read VCF
##'
##' Read a VCF file and print the number of alternate alleles over all records (via \code{\link[base]{table}}).
##' @param vcf.file path to the VCF file
##' @param genome genome identifier (e.g. "VITVI_12x2")
##' @param verbose verbosity level (0/1)
##' @return invisible \code{DNAStringSetList} (from pkg Biostrings)
##' @author Timothee Flutre
##' @export
tableVcfAlt <- function(vcf.file, genome="", verbose=1){
  requireNamespaces(c("VariantAnnotation", "S4Vectors"))
  stopifnot(file.exists(vcf.file))

  if(verbose > 0){
    msg <- "read VCF to count alternate alleles ..."
    write(msg, stdout()); flush(stdout())
  }

  svp <- VariantAnnotation::ScanVcfParam(fixed="ALT", info=NA, geno=NA)
  vcf <- VariantAnnotation::readVcf(file=vcf.file, genome=genome, param=svp)
  if(verbose > 0){
    msg <- paste0("nb of records: ", nrow(vcf))
    write(msg, stdout())
  }

  alts <- VariantAnnotation::alt(vcf)
  names(alts) <- names(vcf)

  print(table(S4Vectors::elementNROWS(alts)))

  invisible(alts)
}

##' Fasta file
##'
##' Simulate a reference genome sequence compatible with alleles from genotype data in the VCF format.
##' Useful with the GATK software.
##' @param vcf.file path to the file in the VCF format (will be used only if \code{vcf.rng.all=NULL}
##' @param vcf.rng.all data frame output from \code{\link{rngVcf2df}}
##' @param chr.lengths vector of chromosome lengths (if a single number, will be repeated, i.e. all chromosomes will have the ame length)
##' @param fa.file path to the output file in the fasta format (will be compressed via \code{gzip} if the name ends with \code{gz}); if NULL, the reference sequence will be returned as an object but not written into a file
##' @param seed seed for the pseudo-random number generator (will be set only if non-NULL)
##' @param verbose verbosity level (0/1)
##' @return invisible DNAStringSet
##' @author Timothee Flutre
##' @export
simulRefseqCompatibleWithVcf <- function(vcf.file=NULL, vcf.rng.all=NULL,
                                         chr.lengths=10^5, fa.file=NULL,
                                         seed=NULL, verbose=1){
  requireNamespaces(c("VariantAnnotation", "Biostrings"))
  stopifnot(! all(is.null(vcf.file), is.null(vcf.rng.all)))
  if(is.null(vcf.rng.all))
    stopifnot(file.exists(vcf.file))

  if(is.null(vcf.rng.all)){
    if(verbose > 0){
      msg <- paste0("retrieve ranges and alleles from file '",
                    vcf.file, "' ...")
      write(msg, stdout())
    }
    vcf <- VariantAnnotation::readVcf(file=vcf.file, genome="")
    vcf.rng.all <- rngVcf2df(vcf=vcf, with.coords=TRUE, with.alleles=TRUE,
                             single.ref=TRUE, single.alt=TRUE)
  }

  chr.names <- unique(vcf.rng.all$chr)
  nb.chrs <- length(chr.names)
  if(verbose > 0){
    msg <- paste0("nb of chromosomes: ", nb.chrs)
    write(msg, stdout())
  }
  if(! is.null(names(chr.lengths)))
    stopifnot(all(names(chr.lengths) %in% chr.names),
              all(chr.names %in% chr.lengths))
  if(all(length(chr.names) > 1, length(chr.lengths) == 1)){
    chr.lengths <- rep(chr.lengths, nb.chrs)
  }
  if(is.null(names(chr.lengths)))
    names(chr.lengths) <- chr.names
  chr.lengths <- chr.lengths[chr.names] # same order
  max.pos <- tapply(vcf.rng.all$pos, as.factor(vcf.rng.all$chr), max)
  for(i in 1:nb.chrs){
    if(max.pos[i] > chr.lengths[i]){
      msg <- paste0(chr.names[i], " has length smaller than",
                    " maximum SNP position")
      stop(msg)
    }
  }
  if(! is.null(seed))
    set.seed(seed)

  if(verbose > 0){
    msg <- paste0("simulate chromosomes for the reference sequence ...")
    write(msg, stdout())
  }
  refseq <- lapply(1:nb.chrs, function(i){
    sample(x=c("A","T","G","C"), size=chr.lengths[i], replace=TRUE)
  })
  refseq <- lapply(1:nb.chrs, function(i){
    tmp <- refseq[[i]]
    idx <- which(vcf.rng.all$chr == chr.names[i])
    tmp[vcf.rng.all$pos[idx]] <- vcf.rng.all$allele.ref[idx]
    paste(tmp, collapse="")
  })
  names(refseq) <- chr.names
  refseq <- Biostrings::DNAStringSet(unlist(refseq))

  if(! is.null(fa.file)){
    if(verbose > 0){
      msg <- paste0("write as fasta in file '", fa.file, "' ...")
      write(msg, stdout())
    }
    Biostrings::writeXStringSet(x=refseq, filepath=fa.file,
                                compress=tools::file_ext(fa.file) == "gz",
                                format="fasta")
  }

  invisible(refseq)
}

##' Information on variant-level calls
##'
##' Return some information to help in hard-filtering variant-level calls.
##' See GATK's Best Practices tutorial (\url{https://www.broadinstitute.org/gatk/guide/topic?name=tutorials#tutorials2806}).
##' @param x data.frame, e.g. from "bcftools query --format '\%CHROM\\t\%POS\\t\%TYPE\\t...'"
##' @param type variant type
##' @param thresh.qual exclude variant if QUAL < threshold
##' @param thresh.qd exclude variant if QD < threshold
##' @param thresh.fs exclude variant if FS > threshold
##' @param thresh.mq exclude variant if MQ < threshold
##' @param thresh.mqrs exclude variant if MQRankSum < threshold (must be negative)
##' @param thresh.rprs exclude variant if ReadPosRankSum < threshold (must be negative)
##' @param thresh.hs exclude variant if HaplotypeScore > threshold (must be positive, \url{http://dx.doi.org/10.1038/sdata.2015.11})
##' @param thresh.dp if NULL, exclude variant if DP > 3 * mean(DP) (\url{https://www.biostars.org/p/110670/#110671})
##' @param thresh.crs exclude variant if ClippingRankSum < threshold (must be negative)
##' @param thresh.pnb exclude variant if PercentNBase < threshold (must be in [0,100])
##' @param thresh.bqrs exclude variant if BaseQRankSum < threshold (must be negative)
##' @param thresh.gc exclude variant if GC > threshold (must be in [0,1])
##' @param thresh.hw exclude variant if HW > threshold (must be positive)
##' @param thresh.hr exclude variant if HRun > threshold (must be positive, \url{http://dx.doi.org/10.1038/sdata.2015.11})
##' @param thresh.ic exclude variant if InbreedingCoeff > threshold
##' @param thresh.lrs exclude variant if LikelihoodRankSum < threshold (must be negative)
##' @param thresh.sor exclude variant if SOR > threshold (must be positive)
##' @return list
##' @author Timothee Flutre
##' @export
infoVariantCalls <- function(x, type="SNP", thresh.qual=20, thresh.qd=2,
                             thresh.fs=60, thresh.mq=40, thresh.mqrs=-12.5,
                             thresh.rprs=-8, thresh.hs=13,
                             thresh.dp=NULL, thresh.crs=NULL,
                             thresh.pnb=NULL, thresh.bqrs=NULL,
                             thresh.gc=NULL, thresh.hw=7, thresh.hr=6,
                             thresh.ic=NULL, thresh.lrs=NULL, thresh.sor=NULL){
  stopifnot(is.data.frame(x),
            "type" %in% colnames(x))

  ## http://stackoverflow.com/a/32012103/597069
  my_summary <- function(v){
    if(! any(is.na(v))){
      res <- c(summary(v), "NA's"=0)
    } else{
      res <- summary(v)
    }
    return(res)
  }

  info <- list()

  idx <- which(x$type == type)
  stopifnot(length(idx) != 0)
  info[[paste0("nb.", type)]] <- length(idx)

  if("qual" %in% colnames(x)){
    info[["summary.qual"]] <- my_summary(x$qual[idx])
    info[["thresh.qual"]] <- thresh.qual
    info[["perc.excl.qual"]] <- 100 * sum(x$qual[idx] < info[["thresh.qual"]],
                                          na.rm=TRUE) / sum(! is.na(x$qual[idx]))
  }

  if("qd" %in% colnames(x)){
    info[["summary.qd"]] <- my_summary(x$qd[idx])
    info[["thresh.qd"]] <- thresh.qd
    info[["perc.excl.qd"]] <- 100 * sum(x$qd[idx] < info[["thresh.qd"]],
                                        na.rm=TRUE) / sum(! is.na(x$qd[idx]))
  }

  if("fs" %in% colnames(x)){
    info[["summary.fs"]] <- my_summary(x$fs[idx])
    info[["thresh.fs"]] <- thresh.fs
    info[["perc.excl.fs"]] <- 100 * sum(x$fs[idx] > info[["thresh.fs"]],
                                        na.rm=TRUE) / sum(! is.na(x$fs[idx]))
  }

  if("mq" %in% colnames(x)){
    info[["summary.mq"]] <- my_summary(x$mq[idx])
    info[["thresh.mq"]] <- thresh.mq
    info[["perc.excl.mq"]] <- 100 * sum(x$mq[idx] < info[["thresh.mq"]],
                                        na.rm=TRUE) / sum(! is.na(x$mq[idx]))
  }

  if("mqrs" %in% colnames(x)){
    info[["summary.mqrs"]] <- my_summary(x$mqrs[idx])
    stopifnot(thresh.mqrs < 0)
    info[["thresh.mqrs"]] <- thresh.mqrs
    info[["perc.excl.mqrs"]] <- 100 * sum(x$mqrs[idx] < info[["thresh.mqrs"]],
                                          na.rm=TRUE) / sum(! is.na(x$mqrs[idx]))
  }

  if("rprs" %in% colnames(x)){
    info[["summary.rprs"]] <- my_summary(x$rprs[idx])
    stopifnot(thresh.rprs < 0)
    info[["thresh.rprs"]] <- thresh.rprs
    info[["perc.excl.rprs"]] <- 100 * sum(x$rprs[idx] < info[["thresh.rprs"]],
                                          na.rm=TRUE) / sum(! is.na(x$rprs[idx]))
  }

  if("hs" %in% colnames(x)){
    if(! all(is.na(x$hs[idx]))){
      info[["summary.hs"]] <- my_summary(x$hs[idx])
      stopifnot(thresh.hs > 0)
      info[["thresh.hs"]] <- thresh.hs
      info[["perc.excl.hs"]] <- 100 * sum(x$hs[idx] < info[["thresh.hs"]],
                                          na.rm=TRUE) / sum(! is.na(x$hs[idx]))
    }
  }

  if("dp" %in% colnames(x)){
    info[["summary.dp"]] <- my_summary(x$dp[idx])
    if(is.null(thresh.dp))
      info[["threshold.dp"]] <- 3 * mean(x$dp[idx])
    info[["perc.excl.dp"]] <- 100 * sum(x$dp[idx] > info[["threshold.dp"]],
                                        na.rm=TRUE) / sum(! is.na(x$dp[idx]))
  }

  if("an" %in% colnames(x)){
    info[["summary.an"]] <- my_summary(x$an[idx])
    ## if(is.null(thresh.an)){
    ##   stopifnot(! is.null(nb.samples))
    ##   info[["thresh.an"]] <- 0.8 * 2 * nb.samples # assume samples are diploids
    ## } else
    ##   info[["thresh.an"]] <- thresh.an
    ## info[["perc.excl.an"]] <- 100 * sum(x$an[idx] > info[["thresh.an"]],
    ##                                     na.rm=TRUE) / sum(! is.na(x$an[idx]))
  }

  if("crs" %in% colnames(x)){
    info[["summary.crs"]] <- my_summary(x$crs[idx])
    if(! is.null(thresh.crs)){
      stopifnot(thresh.crs < 0)
      info[["threshold.crs"]] <- thresh.crs
      info[["perc.excl.crs"]] <- 100 * sum(x$crs[idx] < info[["threshold.crs"]],
                                           na.rm=TRUE) / sum(! is.na(x$crs[idx]))
    }
  }

  if("pnb" %in% colnames(x)){
    info[["summary.pnb"]] <- my_summary(x$pnb[idx])
    if(! is.null(thresh.pnb)){
      stopifnot(thresh.pnb >= 0, thresh.pnb <= 100)
      info[["threshold.pnb"]] <- thresh.pnb
      info[["perc.excl.pnb"]] <- 100 * sum(x$pnb[idx] < info[["threshold.pnb"]],
                                           na.rm=TRUE) / sum(! is.na(x$pnb[idx]))
    }
  }

  if("bqrs" %in% colnames(x)){
    info[["summary.bqrs"]] <- my_summary(x$bqrs[idx])
    if(! is.null(thresh.bqrs)){
      stopifnot(thresh.bqrs < 0)
      info[["threshold.bqrs"]] <- thresh.bqrs
      info[["perc.excl.bqrs"]] <- 100 * sum(x$bqrs[idx] < info[["threshold.bqrs"]],
                                            na.rm=TRUE) / sum(! is.na(x$bqrs[idx]))
    }
  }

  if("gc" %in% colnames(x)){
    info[["summary.gc"]] <- my_summary(x$gc[idx])
    if(! is.null(thresh.gc)){
      stopifnot(thresh.gc >= 0, thresh.gc <= 1)
      info[["threshold.gc"]] <- thresh.gc
      info[["perc.excl.gc"]] <- 100 * sum(x$gc[idx] > info[["threshold.gc"]],
                                          na.rm=TRUE) / sum(! is.na(x$gc[idx]))
    }
  }

  if("hw" %in% colnames(x)){
    info[["summary.hw"]] <- my_summary(x$hw[idx])
    if(! is.null(thresh.hw)){
      stopifnot(thresh.hw > 0)
      info[["threshold.hw"]] <- thresh.hw
      info[["perc.excl.hw"]] <- 100 * sum(x$hw[idx] > info[["threshold.hw"]],
                                          na.rm=TRUE) / sum(! is.na(x$hw[idx]))
    }
  }

  if("hr" %in% colnames(x)){
    info[["summary.hr"]] <- my_summary(x$hr[idx])
    if(! is.null(thresh.hr)){
      stopifnot(thresh.hr > 0)
      info[["threshold.hr"]] <- thresh.hr
      info[["perc.excl.hr"]] <- 100 * sum(x$hr[idx] > info[["threshold.hr"]],
                                          na.rm=TRUE) / sum(! is.na(x$hr[idx]))
    }
  }

  if("ic" %in% colnames(x)){
    info[["summary.ic"]] <- my_summary(x$ic[idx])
    if(! is.null(thresh.ic)){
      info[["threshold.ic"]] <- thresh.ic
      info[["perc.excl.ic"]] <- 100 * sum(x$ic[idx] > info[["threshold.ic"]],
                                          na.rm=TRUE) / sum(! is.na(x$ic[idx]))
    }
  }

  if("lrs" %in% colnames(x)){
    info[["summary.lrs"]] <- my_summary(x$lrs[idx])
    if(! is.null(thresh.lrs)){
      stopifnot(thresh.lrs < 0)
      info[["threshold.lrs"]] <- thresh.lrs
      info[["perc.excl.lrs"]] <- 100 * sum(x$lrs[idx] < info[["threshold.lrs"]],
                                           na.rm=TRUE) / sum(! is.na(x$lrs[idx]))
    }
  }

  if("sor" %in% colnames(x)){
    info[["summary.sor"]] <- my_summary(x$sor[idx])
    if(! is.null(thresh.sor)){
      stopifnot(thresh.sor > 0)
      info[["threshold.sor"]] <- thresh.sor
      info[["perc.excl.sor"]] <- 100 * sum(x$sor[idx] < info[["threshold.sor"]],
                                           na.rm=TRUE) / sum(! is.na(x$sor[idx]))
    }
  }

  tmp <- list()
  tmp[["summary"]] <- do.call(rbind, info[grepl("summary", names(info))])
  tmp[["filters"]] <- cbind(do.call(rbind, info[grepl("thresh", names(info))]),
                            do.call(rbind, info[grepl("perc", names(info))]))
  colnames(tmp[["filters"]]) <- c("thresh", "perc.excl")
  rownames(tmp[["filters"]]) <- sapply(strsplit(x=rownames(tmp[["filters"]]),
                                                split="\\."),
                                       function(x){x[2]})

  return(tmp)
}

##' Plot information on variant-level calls
##'
##' Plot histograms of QUAL, DP, BaseQRankSum, MQRankSum, ReadPosRankSum and ClippingRankSum over 2 rows and 3 columns.
##' @param x data frame, e.g. from "bcftools query --format '\%CHROM\\t\%POS\\t\%TYPE\\t...'"
##' @param main main title
##' @return nothing
##' @author Timothee Flutre
##' @seealso \code{\link{infoVariantCalls}}
##' @export
plotInfoVariantCalls <- function(x, main=""){
  stopifnot(is.data.frame(x),
            all(c("type", "qual", "dp", "bqrs", "mqrs", "rprs", "crs") %in%
                colnames(x)))

  def.par <- graphics::par(mfrow=c(2,3), oma=c(0,0,2,0))

  graphics::hist(log10(x[x$type == "SNP", "qual"]), breaks="Sturges",
                 main="log10(QUAL)",
                 xlab=paste(sum(! is.na(x[x$type == "SNP", "qual"])), "SNPs"))
  graphics::hist(log10(x[x$type == "SNP", "dp"]), breaks="Sturges",
                 main="log10(DP)",
                 xlab=paste(sum(! is.na(x[x$type == "SNP", "dp"])), "SNPs"))
  graphics::hist(x[x$type == "SNP", "bqrs"], breaks="Sturges",
                 main="BaseQRankSum",
                 xlab=paste(sum(! is.na(x[x$type == "SNP", "bqrs"])), "SNPs"))
  graphics::abline(v=c(-2,2), col="red")
  graphics::hist(x[x$type == "SNP", "mqrs"], breaks="Sturges",
                 main="MQRankSum",
                 xlab=paste(sum(! is.na(x[x$type == "SNP", "mqrs"])), "SNPs"))
  graphics::abline(v=c(-2,2), col="red")
  graphics::hist(x[x$type == "SNP", "rprs"], breaks="Sturges",
                 main="ReadPosRankSum",
                 xlab=paste(sum(! is.na(x[x$type == "SNP", "rprs"])), "SNPs"))
  graphics::abline(v=c(-2,2), col="red")
  graphics::hist(x[x$type == "SNP", "crs"], breaks="Sturges",
                 main="ClippingRankSum",
                 xlab=paste(sum(! is.na(x[x$type == "SNP", "crs"])), "SNPs"))
  graphics::abline(v=c(-2,2), col="red")
  graphics::title(main=main, outer=TRUE)

  on.exit(graphics::par(def.par))
}

##' Confidence in one variant's genotypes
##'
##' Provide measure of confidence in the genotypes at a given variant.
##' @param x CollapsedVCF corresponding to a SNV (see pkg \href{http://bioconductor.org/packages/VariantAnnotation/}{VariantAnnotation})
##' @param plot.it plot GQ=f(DP) if TRUE
##' @return invisible data.frame
##' @author Timothee Flutre
##' @export
confidenceGenoOneVar <- function(x, plot.it=FALSE){
  requireNamespaces("VariantAnnotation")
  stopifnot(class(x) == "CollapsedVCF",
            nrow(x) == 1,
            VariantAnnotation::isSNV(x))

  genos <- VariantAnnotation::geno(x)
  stopifnot("GQ" %in% names(genos))

  out <- data.frame(GT=genos[["GT"]][1,],
                    GQ=genos[["GQ"]][1,],
                    stringsAsFactors=FALSE)
  if("DP" %in% names(genos))
    out$DP <- genos[["DP"]][1,]
  if("AD" %in% names(genos)){
    out$AD.0 <- do.call(rbind, genos[["AD"]][1,])[,1]
    out$AD.1 <- do.call(rbind, genos[["AD"]][1,])[,2]
  }
  if("PL" %in% names(genos)){
    out$PL.00 <- do.call(rbind, genos[["PL"]])[,1]
    out$PL.01 <- do.call(rbind, genos[["PL"]])[,2]
    out$PL.11=do.call(rbind, genos[["PL"]])[,3]
  }

  if("DP" %in% names(genos) && plot.it)
    graphics::plot(x=genos[["DP"]][1,], y=genos[["GQ"]][1,],
         xlab="DP", ylab="GQ", main=names(x))

  invisible(out)
}

##' Set GT to NA
##'
##' Set genotypes (GT field) to missing if the genotype quality (GQ field) isn't high enough.
##' @param vcf.file path to the VCF file (if the bgzip index doesn't exist in the same directory, it will be created)
##' @param genome genome identifier (e.g. "VITVI_12x2")
##' @param out.file path to the output VCF file (a bgzip index will be created in the same directory)
##' @param yieldSize number of records to yield each time the file is read from (see ?TabixFile) if seq.id is NULL
##' @param dict.file path to the SAM dict file (see \url{https://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary}) if seq.id is specified with no start/end
##' @param seq.id sequence identifier to work on (e.g. "chr2")
##' @param seq.start start of the sequence to work on (if NULL, whole seq)
##' @param seq.end end of the sequence to work on (if NULL, whole seq)
##' @param min.gq minimum GQ below which GT is set to missing
##' @param verbose verbosity level (0/1)
##' @return the destination file path as an invisible character(1)
##' @author Timothee Flutre
##' @seealso \code{\link{filterVariantCalls}}
##' @export
setGt2Na <- function(vcf.file, genome="", out.file,
                     yieldSize=NA_integer_, dict.file=NULL,
                     seq.id=NULL, seq.start=NULL, seq.end=NULL,
                     min.gq=90, verbose=1){
  requireNamespaces(c("IRanges", "GenomicRanges", "VariantAnnotation",
                      "Rsamtools", "S4Vectors"))
  stopifnot(file.exists(vcf.file),
            xor(is.na(yieldSize), is.null(seq.id)))
  if(! is.null(seq.id) & is.null(seq.start) & is.null(seq.end))
    stopifnot(! is.null(dict.file),
              file.exists(dict.file))
  out.file <- sub("\\.gz$", "", out.file)
  out.file <- sub("\\.bgz$", "", out.file)

  dest <- NULL

  if(verbose > 0){
    msg <- paste0("read VCF to set GT to NA if GQ < ", min.gq, " ...")
    write(msg, stdout()); flush(stdout())
  }

  if(! file.exists(paste0(vcf.file, ".tbi")))
    Rsamtools::indexTabix(file=vcf.file, format="vcf")
  tabix.file <- Rsamtools::TabixFile(file=vcf.file,
                                     yieldSize=yieldSize)
  if(! is.null(seq.id)){
    rngs <- seqIdStartEnd2GRanges(seq.id=seq.id, seq.start=seq.start,
                                  seq.end=seq.end, dict.file=dict.file)
    vcf.params <- VariantAnnotation::ScanVcfParam(which=rngs)
    vcf <- VariantAnnotation::readVcf(file=tabix.file, genome=genome,
                                      param=vcf.params)
    nb.variants <- nrow(vcf)
    idx <- (VariantAnnotation::geno(vcf)[["GQ"]] < min.gq)
    VariantAnnotation::geno(vcf)[["GT"]][idx] <- "./."
    dest <- VariantAnnotation::writeVcf(obj=vcf, filename=out.file, index=TRUE)
  } else{
    open(tabix.file)
    if(file.exists(out.file))
      file.remove(out.file)
    con <- file(out.file, open="a")
    nb.variants <- 0
    while(nrow(vcf <- VariantAnnotation::readVcf(file=tabix.file,
                                                 genome=genome))){
      nb.variants <- nb.variants + nrow(vcf)
      idx <- (VariantAnnotation::geno(vcf)[["GQ"]] < min.gq)
      VariantAnnotation::geno(vcf)[["GT"]][idx] <- "./."
      VariantAnnotation::writeVcf(obj=vcf, filename=con)
    }
    close(con)
    close(tabix.file)
    dest <- Rsamtools::bgzip(file=out.file, overwrite=TRUE)
    Rsamtools::indexTabix(file=dest, format="vcf")
    file.remove(out.file)
  }

  if(verbose > 0){
    msg <- paste0("nb of variants: ", nb.variants)
    write(msg, stdout())
  }

  invisible(dest)
}

##' Filter variant calls
##'
##' Filter out variant calls from VCF file according to several criteria (bi-allelic, single nucleotide variant, proper amount of missing genotypes, overall depth and allele frequency).
##' @param vcf.file path to the VCF file (if the bgzip index doesn't exist in the same directory, it will be created)
##' @param genome genome identifier (e.g. "VITVI_12x2")
##' @param out.file path to the output VCF file (a bgzip index will be created in the same directory)
##' @param yieldSize number of records to yield each time the file is read from (see ?TabixFile) if seq.id is NULL
##' @param dict.file path to the SAM dict file (see \url{https://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary}) if seq.id is specified with no start/end
##' @param seq.id sequence identifier to work on (e.g. \code{"chr2"})
##' @param seq.start start of the sequence to work on (if NULL, whole seq)
##' @param seq.end end of the sequence to work on (if NULL, whole seq)
##' @param variants.tokeep character vector of variant names to keep (e.g. \code{c("chr1:35718_C/A","chr1:61125_A/G")})
##' @param is.snv if not NULL but TRUE, filter out the variants which are not SNVs
##' @param is.biall if not NULL but TRUE, filter out the variants with more than one alternative allele
##' @param min.var.dp minimum variant-level DP below which variants are filtered out
##' @param max.var.dp maximum variant-level DP above which variants are filtered out
##' @param min.alt.af minimum variant-level AF below which variants are filtered out
##' @param max.alt.af maximum variant-level AF above which variants are filtered out
##' @param min.spl.dp minimum sample-level DP
##' @param min.perc.spl.dp minimum percentage of samples with DP above threshold
##' @param min.spl.gq minimum sample-level GQ
##' @param min.perc.spl.gq minimum percentage of samples with GQ above threshold
##' @param max.var.nb.gt.na maximum number of samples with missing GT
##' @param max.var.perc.gt.na maximum percentage of samples with missing GT
##' @param verbose verbosity level (0/1)
##' @return the destination file path as a character(1)
##' @author Timothee Flutre
##' @export
filterVariantCalls <- function(vcf.file, genome="", out.file,
                               yieldSize=NA_integer_, dict.file=NULL,
                               seq.id=NULL, seq.start=NULL, seq.end=NULL,
                               variants.tokeep=NULL,
                               is.snv=NULL, is.biall=NULL,
                               min.var.dp=NULL, max.var.dp=NULL,
                               min.alt.af=NULL, max.alt.af=NULL,
                               min.spl.dp=NULL, min.perc.spl.dp=NULL,
                               min.spl.gq=NULL, min.perc.spl.gq=NULL,
                               max.var.nb.gt.na=NULL, max.var.perc.gt.na=NULL,
                               verbose=1){
  requireNamespaces(c("IRanges", "GenomicRanges", "VariantAnnotation",
                      "Rsamtools", "S4Vectors", "BiocInstaller",
                      "SummarizedExperiment"))
  stopifnot(file.exists(vcf.file))
  if(! is.null(seq.id) & is.null(seq.start) & is.null(seq.end))
    stopifnot(! is.null(dict.file),
              file.exists(dict.file))
  if(! is.null(variants.tokeep))
    stopifnot(is.vector(variants.tokeep),
              is.character(variants.tokeep))
  if(! is.null(is.snv))
    stopifnot(is.logical(is.snv))
  if(! is.null(is.biall))
    stopifnot(is.logical(is.biall))
  if(! is.null(min.perc.spl.dp))
    stopifnot(min.perc.spl.dp >= 0, min.perc.spl.dp <= 100)
  if(! is.null(min.perc.spl.gq))
    stopifnot(min.perc.spl.gq >= 0, min.perc.spl.gq <= 100)
  if(! is.null(max.var.nb.gt.na))
    stopifnot(max.var.nb.gt.na >= 0)
  if(! is.null(max.var.perc.gt.na))
    stopifnot(max.var.perc.gt.na >= 0, max.var.perc.gt.na <= 100)

  dest <- NULL

  if(verbose > 0){
    msg <- "read VCF to filter variants ..."
    write(msg, stdout()); flush(stdout())
  }

  out.file <- sub("\\.gz$", "", out.file)
  out.file <- sub("\\.bgz$", "", out.file)

  ## return TRUE if in subset of variants to keep
  filterVariantsToKeep <- function(x, vtk=variants.tokeep){
    (names(SummarizedExperiment::rowRanges(x)) %in% vtk)
  }

  ## return TRUE if SNV (single nucleotide variant)
  filterSnv <- function(x){
    (VariantAnnotation::isSNV(x))
  }

  ## return TRUE if at most one alternate allele
  filterBiall <- function(x){
    if(utils::compareVersion(as.character(BiocInstaller::biocVersion()),
                             "3.4") < 0){
      (S4Vectors::elementLengths(VariantAnnotation::alt(x)) <= 1)
    } else
      (S4Vectors::elementNROWS(VariantAnnotation::alt(x)) <= 1)
  }

  ## return TRUE if variant-level DP inside of given range
  filterVariantDp <- function(x, min.dp=min.var.dp, max.dp=max.var.dp){
    if(nrow(x) == 0){
      logical(0)
    } else{
      dp <- VariantAnnotation::info(x)$DP
      (dp >= min.dp) || (dp <= max.dp)
    }
  }

  ## return TRUE if allele frequency inside of given range
  filterAf <- function(x, min.af=min.alt.af, max.af=max.alt.af){
    if(nrow(x) == 0){
      logical(0)
    } else{
      af <- unlist(VariantAnnotation::info(x)$AF)
      (af >= min.af) & (af <= max.af)
    }
  }

  ## return TRUE if high-enough percentage of samples with DP above threshold
  filterSampleDp <- function(x, min.dp=min.spl.dp,
                             min.perc.dp=min.perc.spl.dp){
    if(nrow(x) == 0){
      logical(0)
    } else{
      dp <- VariantAnnotation::geno(x)$DP
      (100 * rowSums(dp >= min.dp) / ncol(dp) >= min.perc.dp)
    }
  }

  ## return TRUE if high-enough percentage of samples with GQ above threshold
  filterSampleGq <- function(x, min.gq=min.spl.gq,
                             min.perc.gq=min.perc.spl.gq){
    if(nrow(x) == 0){
      logical(0)
    } else{
      gq <- VariantAnnotation::geno(x)$GQ
      (100 * rowSums(gq >= min.gq) / ncol(gq) >= min.perc.gq)
    }
  }

  ## return TRUE if high-enough number of samples without missing genotypes
  filterGtNb <- function(x, max.nb.gt.na=max.var.nb.gt.na){
    if(nrow(x) == 0){
      logical(0)
    } else{
      gt <- VariantAnnotation::geno(x)$GT
      is.missing <- matrix(gt %in% c(".","./."), nrow(gt), ncol(gt))
      (rowSums(is.missing) <= max.nb.gt.na)
    }
  }

  ## return TRUE if high-enough percentage of samples without missing genotypes
  filterGtPerc <- function(x, max.perc.gt.na=max.var.perc.gt.na){
    if(nrow(x) == 0){
      logical(0)
    } else{
      gt <- VariantAnnotation::geno(x)$GT
      is.missing <- matrix(gt %in% c(".","./."), nrow(gt), ncol(gt))
      (100 * rowSums(is.missing) / ncol(gt) <= max.perc.gt.na)
    }
  }

  ## set the filters
  tmp <- list()
  if(! is.null(variants.tokeep)){
    tmp[["filterVariantsToKeep"]] <- filterVariantsToKeep
  }
  if(! is.null(is.snv)){
    if(is.snv)
      tmp[["filterSnv"]] <- filterSnv
  }
  if(! is.null(is.biall)){
    if(is.biall)
      tmp[["filterBiall"]] <- filterBiall
  }
  if(! is.null(min.var.dp) & ! is.null(max.var.dp)){
    tmp[["filterVariantDp"]] <- filterVariantDp
  }
  if(! is.null(min.alt.af) & ! is.null(max.alt.af)){
    tmp[["filterAf"]] <- filterAf
  }
  if(! is.null(min.spl.dp) & ! is.null(min.perc.spl.dp)){
    tmp[["filterSampleDp"]] <- filterSampleDp
  }
  if(! is.null(min.spl.gq) & ! is.null(min.perc.spl.gq)){
    tmp[["filterSampleGq"]] <- filterSampleGq
  }
  if(! is.null(max.var.nb.gt.na)){
    tmp[["filterGtNb"]] <- filterGtNb
  }
  if(! is.null(max.var.perc.gt.na)){
    tmp[["filterGtPerc"]] <- filterGtPerc
  }

  if(length(tmp) > 0){
    filters <- S4Vectors::FilterRules(tmp)
    if(verbose > 0)
      print(filters)

    ## filter the VCF file
    if(! file.exists(paste0(vcf.file, ".tbi")))
      Rsamtools::indexTabix(file=vcf.file, format="vcf")
    tabix.file <- Rsamtools::TabixFile(file=vcf.file,
                                       yieldSize=yieldSize)
    if(! is.null(seq.id)){
      rngs <- seqIdStartEnd2GRanges(seq.id=seq.id, seq.start=seq.start,
                                    seq.end=seq.end, dict.file=dict.file)
      vcf.params <- VariantAnnotation::ScanVcfParam(which=rngs)
      dest <- VariantAnnotation::filterVcf(file=tabix.file, genome=genome,
                                           destination=out.file, index=TRUE,
                                           verbose=(verbose > 0),
                                           filters=filters,
                                           param=vcf.params)
    } else{
      dest <- VariantAnnotation::filterVcf(file=tabix.file, genome=genome,
                                           destination=out.file, index=TRUE,
                                           verbose=(verbose > 0),
                                           filters=filters)
    }
  }

  invisible(sub(".tbi", "", dest))
}

##' Summary per variant
##'
##' Compute the mean, sd, min, Q1, med, mean, Q3, max of the genotype qualities per variant, also reporting the number of samples and the number of missing data.
##' @param vcf CollapsedVCF (see pkg \href{http://bioconductor.org/packages/VariantAnnotation/}{VariantAnnotation})
##' @param fields genotype field(s) of the VCF to parse (\code{"DP"}/\code{"GQ"}/\code{c("DP","GQ")})
##' @return list of matrices (one per field) with one row per variant and 9 columns (n, na, mean, sd, min, q1, med, q3, max)
##' @author Timothee Flutre
##' @seealso \code{\link{summaryVariant}}
##' @export
varqual2summary <- function(vcf, fields="GQ"){
  requireNamespace("VariantAnnotation")
  stopifnot(is.character(fields))

  output <- list()

  for(field in fields){
    mat <- VariantAnnotation::geno(vcf)[[field]]
    if(class(mat[1,1]) != "list"){ # DP, GQ
      output[[field]] <-
        cbind(n=rep(ncol(mat), nrow(mat)),
              na=rowSums(t(apply(mat, 1, is.na))),
              mean=suppressWarnings(apply(mat, 1, mean, na.rm=TRUE)),
              sd=suppressWarnings(apply(mat, 1, stats::sd, na.rm=TRUE)),
              min=suppressWarnings(apply(mat, 1, min, na.rm=TRUE)),
              q1=suppressWarnings(apply(mat, 1, stats::quantile,
                                        probs=0.25, na.rm=TRUE)),
              med=suppressWarnings(apply(mat, 1, stats::median, na.rm=TRUE)),
              q3=suppressWarnings(apply(mat, 1, stats::quantile,
                                        probs=0.75, na.rm=TRUE)),
              max=suppressWarnings(apply(mat, 1, max, na.rm=TRUE)))
    } else{ # AD
      ## TODO
    }
    rownames(output[[field]]) <- rownames(mat)
  }

  return(output)
}

##' Summary per variant
##'
##' Compute the mean, sd, min, Q1, med, mean, Q3, max of the genotype qualities per variant, also reporting the number of samples and the number of missing data.
##' @param vcf.file path to the VCF file (if the bgzip index doesn't exist in the same directory, it will be created)
##' @param genome genome identifier (e.g. "VITVI_12x2")
##' @param yieldSize number of records to yield each time the file is read from (see ?TabixFile) if seq.id is NULL
##' @param dict.file path to the SAM dict file (see \url{https://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary}) if seq.id is specified with no start/end
##' @param seq.id sequence identifier to work on (e.g. "chr2")
##' @param seq.start start of the sequence to work on
##' @param seq.end end of the sequence to work on
##' @param fields genotype field(s) of the VCF to parse (\code{"DP"}/\code{"GQ"}/\code{c("DP","GQ")})
##' @param verbose verbosity level (0/1)
##' @return list of matrices (one per field) with one row per variant and 9 columns (n, na, mean, sd, min, q1, med, q3, max)
##' @author Timothee Flutre
##' @seealso \code{\link{varqual2summary}}
##' @export
summaryVariant <- function(vcf.file, genome="", yieldSize=NA_integer_,
                           dict.file=NULL, seq.id=NULL, seq.start=NULL,
                           seq.end=NULL, fields="GQ", verbose=1){
  requireNamespaces(c("IRanges", "GenomicRanges", "VariantAnnotation",
                      "Rsamtools"))
  stopifnot(file.exists(vcf.file),
            xor(is.na(yieldSize), is.null(seq.id)),
            is.character(fields),
            all(fields %in% c("GQ", "DP")))
  if(! is.null(seq.id) & is.null(seq.start) & is.null(seq.end))
    stopifnot(! is.null(dict.file),
              file.exists(dict.file))

  output <- NULL

  if(verbose > 0){
    msg <- paste0("read VCF to summarize ", paste(fields, collapse="-"),
                  " from variants ...")
    write(msg, stdout()); flush(stdout())
  }

  if(! file.exists(paste0(vcf.file, ".tbi")))
    Rsamtools::indexTabix(file=vcf.file, format="vcf")
  tabix.file <- Rsamtools::TabixFile(file=vcf.file,
                                     yieldSize=yieldSize)
  if(! is.null(seq.id)){
    rngs <- seqIdStartEnd2GRanges(seq.id=seq.id, seq.start=seq.start,
                                  seq.end=seq.end, dict.file=dict.file)
    vcf.params <- VariantAnnotation::ScanVcfParam(which=rngs)
    vcf <- VariantAnnotation::readVcf(file=tabix.file, genome=genome,
                                      param=vcf.params)
    output <- varqual2summary(vcf, fields)
  } else{
    open(tabix.file)
    while(nrow(vcf <- VariantAnnotation::readVcf(file=tabix.file,
                                                 genome=genome))){
      if(is.null(output)){
        output <- varqual2summary(vcf, fields)
      } else{
        tmp <- varqual2summary(vcf, fields)
        for(field in fields)
          output[[field]] <- rbind(output[[field]], tmp[[field]])
      }
    }
    close(tabix.file)
  }

  if(verbose > 0){
    msg <- paste0("nb of variants: ", nrow(output[[1]]))
    write(msg, stdout())
  }

  return(output)
}

##' Make GRanges
##'
##' Return a \code{GRanges} object from sequence identifiers, starts and ends.
##' @param seq.id sequence identifier(s) to work on (e.g. \code{"chr2"} or \code{c("chr2","chr7")})
##' @param seq.start start(s) of the sequence(s) to work on (if NULL, whole seq; see \code{dict.file})
##' @param seq.end end(s) of the sequence(s) to work on (if NULL, whole seq; see \code{dict.file})
##' @param dict.file path to the SAM dict file (see \url{https://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary}) if \code{seq.id} is specified with no start/end
##' @param subseq.name names of the subsequence(s) (optional; e.g. \code{"gene34"} or \code{"snp005"})
##' @return \code{GRanges} object from the \code{GenomicRanges} package
##' @author Timothee Flutre
##' @seealso \code{\link{vcf2dosage}}, \code{\link{vcf2genoClasses}}
##' @export
seqIdStartEnd2GRanges <- function(seq.id, seq.start=NULL, seq.end=NULL,
                                  dict.file=NULL, subseq.name=NULL){
  requireNamespaces(c("IRanges", "GenomicRanges"))
  if(all(is.null(seq.start), is.null(seq.end)))
    stopifnot(! is.null(dict.file),
              file.exists(dict.file))
  if(all(! is.null(seq.start), ! is.null(seq.end)))
    stopifnot(all(length(seq.start) == length(seq.id),
                  length(seq.end) == length(seq.id)))
  if(! is.null(subseq.name))
    stopifnot(is.character(subseq.name),
              length(subseq.name) == length(seq.id))

  for(i in seq_along(seq.id)){
    if(any(is.null(seq.start), is.null(seq.end))){
      dict <- readSamDict(file=dict.file)
      if(! seq.id[i] %in% rownames(dict)){
        msg <- paste0("seq '", seq.id[i], "' not in '", dict.file, "'")
        stop(msg)
      }
      if(is.null(seq.start))
        seq.start <- c(seq.start, 1)
      if(is.null(seq.end))
        seq.end <- c(seq.end, dict$LN[rownames(dict) == seq.id])
    }
  }

  rngs <- GenomicRanges::GRanges(seqnames=c(seq.id),
                                 ranges=IRanges::IRanges(start=c(seq.start),
                                                         end=c(seq.end)))
  if(! is.null(subseq.name))
    names(rngs) <- subseq.name

  return(rngs)
}

##' Parse VCF
##'
##' Returns a VCF object subsetted for allelicity.
##' @param vcf CollapsedVCF (see pkg \href{http://bioconductor.org/packages/VariantAnnotation/}{VariantAnnotation})
##' @param single.ref if TRUE, only records with a single 'ref' are kept
##' @param single.alt if TRUE, only records with a single 'alt' are kept
##' @return CollapsedVCF
##' @author Gautier Sarah, Timothee Flutre
##' @export
subsetVcfOnAllelicity <- function(vcf, single.ref=TRUE, single.alt=TRUE){
  stopifnot(all(is.logical(single.ref), is.logical(single.alt)))

  idx <- 1:nrow(vcf)

  if(any(single.ref, single.alt)){
    requireNamespaces(c("S4Vectors", "VariantAnnotation", "BiocInstaller"))

    if(utils::compareVersion(as.character(BiocInstaller::biocVersion()),
                             "3.4") < 0){
      idxRef <- S4Vectors::elementLengths(VariantAnnotation::ref(vcf)) == 1L
      idxAlt <- S4Vectors::elementLengths(VariantAnnotation::alt(vcf)) == 1L

    } else{
      idxRef <- S4Vectors::elementNROWS(VariantAnnotation::ref(vcf)) == 1L
      idxAlt <- S4Vectors::elementNROWS(VariantAnnotation::alt(vcf)) == 1L
    }

    if(single.ref){
      idx <- idxRef
      if(single.alt){
        idx <- idx & idxAlt
        if(! all(idx))
          warning("keep records with single 'ref' and single 'alt'")
      } else if(! all(idx))
        warning("keep records with single 'ref' and possibly multiple 'alt'")
    } else if(single.alt){
      idx <- idxAlt
      if(! all(idx))
        warning("keep records with single 'alt' and possibly multiple 'ref'")
    }
  }

  return(vcf[idx])
}

##' Convert GT to dosage
##'
##' Non-bi-allelic variants are discarded.
##' From Martin Morgan (see http://grokbase.com/t/r/bioconductor/135b460s2b/bioc-how-to-convert-genotype-snp-matrix-to-nucleotide-genotypes).
##' @param vcf CollapsedVCF (see pkg \href{http://bioconductor.org/packages/VariantAnnotation/}{VariantAnnotation})
##' @return matrix (in integer \code{\link{mode}}) with variants in rows and samples in columns
##' @author Timothee Flutre
##' @seealso \code{\link{vcf2dosage}}, \code{\link{gtVcf2genoClasses}}, \code{\link{dsVcf2dose}}
##' @export
gtVcf2dose <- function(vcf){
  requireNamespace("VariantAnnotation")

  vcf <- subsetVcfOnAllelicity(vcf=vcf, single.ref=TRUE, single.alt=TRUE)

  gt <- sub("\\|", "/", VariantAnnotation::geno(vcf)$GT) # ignore phasing
  gt[which(gt == "0/0")] <- 0
  gt[which(gt == "0/1")] <- 1
  gt[which(gt == "1/1")] <- 2
  gt[which(gt %in% c(".","./."))] <- NA
  mode(gt) <- "integer"

  return(gt)
}

##' Convert DS to dosage
##'
##' Non-bi-allelic variants are discarded.
##' From Martin Morgan (see http://grokbase.com/t/r/bioconductor/135b460s2b/bioc-how-to-convert-genotype-snp-matrix-to-nucleotide-genotypes).
##' @param vcf CollapsedVCF (see pkg \href{http://bioconductor.org/packages/VariantAnnotation/}{VariantAnnotation})
##' @return matrix (in numeric \code{\link{mode}}) with variants in rows and samples in columns
##' @author Timothee Flutre
##' @seealso \code{\link{vcf2dosage}}, \code{\link{gtVcf2dose}}
##' @export
dsVcf2dose <- function(vcf){
  requireNamespace("VariantAnnotation")

  vcf <- subsetVcfOnAllelicity(vcf=vcf, single.ref=TRUE, single.alt=TRUE)

  ds <- VariantAnnotation::geno(vcf)$DS
  mode(ds) <- "numeric"

  return(ds)
}

##' Convert ranges to data.frame
##'
##' Non-bi-allelic variants are discarded.
##' @param vcf CollapsedVCF (see pkg \href{http://bioconductor.org/packages/VariantAnnotation/}{VariantAnnotation})
##' @param with.coords if TRUE, the output will contain variant coordinates
##' @param with.alleles if TRUE, the output will contain variant alleles
##' @param single.ref if TRUE, only records with a single 'ref' are kept
##' @param single.alt if TRUE, only records with a single 'alt' are kept
##' @return data.frame with variants in rows
##' @author Timothee Flutre
##' @seealso \code{\link{subsetVcfOnAllelicity}}
##' @export
rngVcf2df <- function(vcf, with.coords=TRUE, with.alleles=TRUE,
                      single.ref=FALSE, single.alt=FALSE){
  requireNamespaces(c("VariantAnnotation", "GenomeInfoDb",
                      "BiocGenerics", "SummarizedExperiment"))
  stopifnot(all(is.logical(with.coords), is.logical(with.alleles)),
            any(with.coords, with.alleles))

  vcf <- subsetVcfOnAllelicity(vcf=vcf, single.ref=single.ref,
                               single.alt=single.alt)

  tmp <- SummarizedExperiment::rowRanges(vcf)

  if(all(with.coords, with.alleles)){
    out <- data.frame(chr=as.character(GenomeInfoDb::seqnames(tmp)),
                      pos=BiocGenerics::end(tmp),
                      allele.ref=as.character(tmp$REF),
                      allele.alt=as.character(BiocGenerics::unlist(tmp$ALT)),
                      row.names=names(tmp),
                      stringsAsFactors=FALSE)
  } else if(all(with.coords, ! with.alleles)){
    out <- data.frame(chr=as.character(GenomeInfoDb::seqnames(tmp)),
                      pos=BiocGenerics::end(tmp),
                      row.names=names(tmp),
                      stringsAsFactors=FALSE)
  } else if(all(! with.coords, with.alleles)){
    out <- data.frame(allele.ref=as.character(tmp$REF),
                      allele.alt=as.character(BiocGenerics::unlist(tmp$ALT)),
                      row.names=names(tmp),
                      stringsAsFactors=FALSE)
  }

  return(out)
}

##' Convert VCF to dose
##'
##' Convert genotypes at bi-allelic variants from a VCF file into allele doses.
##' This is also doable with \code{bcftools +dosage}.
##' @param vcf.file path to the VCF file (if the bgzip index doesn't exist in the same directory, it will be created)
##' @param genome genome identifier (e.g. "VITVI_12x2")
##' @param gdose.file path to the output file to record genotypes as allele doses (will be gzipped); variants will be in rows and samples in columns
##' @param ca.file path to the output file to record SNP 1-based coordinates and alleles (will be gzipped)
##' @param yieldSize number of records to yield each time the file is read from (see ?TabixFile) if seq.id is NULL
##' @param dict.file path to the SAM dict file (see \url{https://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary}) if seq.id is specified with no start/end
##' @param seq.id see \code{\link{seqIdStartEnd2GRanges}}
##' @param seq.start see \code{\link{seqIdStartEnd2GRanges}}
##' @param seq.end see \code{\link{seqIdStartEnd2GRanges}}
##' @param field the genotypes to convert should come from the "GT" or "DS" fields
##' @param verbose verbosity level (0/1)
##' @return an invisible list with both output file paths
##' @author Timothee Flutre
##' @seealso \code{\link{gtVcf2dose}}, \code{\link{dsVcf2dose}}, \code{\link{filterVariantCalls}}
##' @export
vcf2dosage <- function(vcf.file, genome="", gdose.file, ca.file,
                       yieldSize=NA_integer_, dict.file=NULL,
                       seq.id=NULL, seq.start=NULL, seq.end=NULL,
                       field="GT", verbose=1){
  requireNamespaces(c("IRanges", "GenomicRanges", "VariantAnnotation",
                      "Rsamtools"))
  stopifnot(file.exists(vcf.file),
            xor(is.na(yieldSize), is.null(seq.id)),
            field %in% c("GT", "DS"))
  if(! is.null(seq.id) & is.null(seq.start) & is.null(seq.end))
    stopifnot(! is.null(dict.file),
              file.exists(dict.file))

  if(verbose > 0){
    msg <- paste0("read VCF to convert genotypes into allele doses",
                  " (field=", field, ")...")
    write(msg, stdout()); flush(stdout())
  }

  for(out.file in c(gdose.file, ca.file))
    if(file.exists(out.file))
      file.remove(out.file)
  gdose.file <- sub("\\.gz$", "", gdose.file)
  ca.file <- sub("\\.gz$", "", ca.file)
  gdose.con <- file(gdose.file, open="a")
  ca.con <- file(ca.file, open="a")
  cat("chr\tpos\tallele.ref\tallele.alt\n", file=ca.con, append=TRUE)

  if(! file.exists(paste0(vcf.file, ".tbi")))
    Rsamtools::indexTabix(file=vcf.file, format="vcf")
  tabix.file <- Rsamtools::TabixFile(file=vcf.file,
                                     yieldSize=yieldSize)
  if(! is.null(seq.id)){
    rngs <- seqIdStartEnd2GRanges(seq.id=seq.id, seq.start=seq.start,
                                  seq.end=seq.end, dict.file=dict.file)
    vcf.params <- VariantAnnotation::ScanVcfParam(which=rngs)
    vcf <- VariantAnnotation::readVcf(file=tabix.file, genome=genome,
                                      param=vcf.params)
    nb.variants <- nrow(vcf)
    if(field == "GT"){
      gtmp <- gtVcf2dose(vcf=vcf)
    } else if(field == "DS")
      gtmp <- dsVcf2dose(vcf=vcf)
    catmp <- rngVcf2df(vcf=vcf, with.coords=TRUE, with.alleles=TRUE,
                       single.ref=TRUE, single.alt=TRUE)
    cat(paste(colnames(gtmp), collapse="\t"), file=gdose.con,
        append=TRUE, sep="\n")
    utils::write.table(x=gtmp,
                       file=gdose.con, append=TRUE,
                       quote=FALSE, sep="\t", row.names=TRUE,
                       col.names=FALSE)
    utils::write.table(x=catmp, file=ca.con, append=TRUE,
                       quote=FALSE, sep="\t", row.names=TRUE,
                       col.names=FALSE)
  } else{
    open(tabix.file)
    nb.variants <- 0
    while(nrow(vcf <- VariantAnnotation::readVcf(file=tabix.file,
                                                 genome=genome))){
      if(field == "GT"){
        gtmp <- gtVcf2dose(vcf=vcf)
      } else if(field == "DS")
        gtmp <- dsVcf2dose(vcf=vcf)
      catmp <- rngVcf2df(vcf=vcf, with.coords=TRUE, with.alleles=TRUE,
                         single.ref=TRUE, single.alt=TRUE)
      if(nb.variants == 0)
        cat(paste(colnames(gtmp), collapse="\t"), file=gdose.con,
            append=TRUE, sep="\n")
      nb.variants <- nb.variants + nrow(vcf)
      utils::write.table(x=gtmp, file=gdose.con, append=TRUE,
                         quote=FALSE, sep="\t", row.names=TRUE,
                         col.names=FALSE)
      utils::write.table(x=catmp, file=ca.con, append=TRUE,
                         quote=FALSE, sep="\t", row.names=TRUE,
                         col.names=FALSE)
    }
    close(tabix.file)
  }

  if(verbose > 0){
    msg <- paste0("nb of variants: ", nb.variants)
    write(msg, stdout())
  }

  close(gdose.con)
  close(ca.con)
  for(gz.out.file in c(paste0(gdose.file, ".gz"), paste0(ca.file, ".gz")))
    if(file.exists(gz.out.file))
      file.remove(gz.out.file)
  if(file.exists(Sys.which("gzip"))){ # should work on Linux and Mac OS
    system(command=paste("gzip -n", gdose.file))
    system(command=paste("gzip -n", ca.file))
  } else{
    warning("'gzip' not available on this computer")
  }

  invisible(list(gdose.file=paste0(gdose.file, ".gz"),
                 ca.file=paste0(ca.file, ".gz")))
}

##' Samples from VCF file
##'
##' Return the samples from a VCF file using \code{\link{pipe}} and \code{cat} or \code{zcat} (depending on the file extension).
##' @param vcf.file path to the VCF file
##' @return vector of samples
##' @author Timothee Flutre
##' @export
getSamplesFromVcfFile <- function(vcf.file){
  stopifnot(file.exists(vcf.file))

  cmd <- ""
  if(grepl("gz", tools::file_ext(vcf.file)))
    cmd <- "z"
  cmd <- paste0(cmd, "cat")
  stopifnot(file.exists(Sys.which(cmd)))
  cmd <- paste0(cmd, " ", vcf.file, " 2>/dev/null",
                " | grep -m 1 '^#CHROM'")
  con <- pipe(cmd)
  samples <- scan(con, what="character", quiet=TRUE)
  close(con)
  samples <- samples[-c(1:9)]

  return(samples)
}

##' Sort a VCF file
##'
##' Sort a VCF file with \code{bcftools sort} (requires bcftools >= 1.6)
##' @param in.vcf.file path to the input VCF file
##' @param out.vcf.file path to the output VCF file (if it is the same as \code{in.vcf.file}, the input VCF file will be sorted into a temporary file which will then be copied back to it)
##' @return invisible path to the output file
##' @author Timothee Flutre
##' @seealso \code{\link{indexVcfFile}}
##' @export
sortVcfFile <- function(in.vcf.file, out.vcf.file){
  exe.name <- "bcftools"
  stopifnot(file.exists(Sys.which(exe.name)))
  stopifnot(file.exists(in.vcf.file),
            tools::file_ext(in.vcf.file) == "gz",
            tools::file_ext(out.vcf.file) == "gz")

  tmp.file <- out.vcf.file
  if(out.vcf.file == in.vcf.file)
    tmp.file <- tempfile(pattern="sortVcfFile_", fileext=".vcf.gz")

  cmd <- paste0(exe.name, " view",
                " -O u", # export as BCF
                " ", in.vcf.file,
                " | ", exe.name, " sort",
                " -O z", # export as compressed VCF
                " -o ", tmp.file,
                " -")
  ret <- system(cmd)

  if(out.vcf.file == in.vcf.file){
    file.copy(from=tmp.file, to=out.vcf.file)
    file.remove(tmp.file)
  }

  invisible(out.vcf.file)
}

##' Index a VCF file
##'
##' Index a VCF file if it exists and if its index file doesn't already exist or is older.
##' @param vcf.file path to a VCF file (must be sorted!)
##' @param software name of the sofwatre with which to index the file (Rsamtools/bcftools)
##' @param nb.cores number of threads for "bcftools"
##' @return invisible path to the index file
##' @author Timothee Flutre
##' @seealso \code{\link{sortVcfFile}}
##' @export
indexVcfFile <- function(vcf.file, software="Rsamtools", nb.cores=1){
  idx.file <- paste0(vcf.file, ".tbi")
  need.to.index <- FALSE

  if(file.exists(vcf.file)){
    if(! file.exists(idx.file)){
      need.to.index <- TRUE
    } else{
      mtime.vcf <- file.mtime(vcf.file)
      mtime.idx <- file.mtime(idx.file)
      if(as.POSIXlt(mtime.idx) < as.POSIXlt(mtime.vcf)){
        need.to.index <- TRUE
      }
    }
  } else{
    msg <- paste0("can't index ", vcf.file, " because it doesn't exist")
    warning(msg)
  }

  if(need.to.index){
    stopifnot(software %in% c("Rsamtools", "bcftools"))
    if(software == "Rsamtools"){
      requireNamespace("Rsamtools")
      Rsamtools::indexTabix(file=vcf.file, format="vcf")
    } else if(software == "bcftools"){
      cmd <- paste0("bcftools index",
                    " -f -t",
                    " --threads ", nb.cores,
                    " ", vcf.file)
      ret <- system(cmd)
    }
  }

  invisible(idx.file)
}

##' Rename VCF samples
##'
##' Rename samples in a VCF file with \code{bcftools reheader}.
##' Variants can also be sorted.
##' @param in.vcf.file path to the input VCF file (will be automatically indexed if necessary)
##' @param samples data frame with (at least) two columns, the first corresponding to the old names (same order as in the VCF file) and the second to the new names; spaces inside names are forbidden
##' @param out.vcf.file path to the output VCF file
##' @param skip.check if TRUE, no comparison is done to check that \code{samples} is coherent with the data in the VCF file
##' @param sort.variants if TRUE, variants are also sorted, with \code{bcftools sort} (requires bcftools >= 1.6)
##' @return invisible \code{out.vcf.file}
##' @author Timothee Flutre
##' @export
renameVcfSamples <- function(in.vcf.file, samples, out.vcf.file,
                             skip.check=FALSE, sort.variants=TRUE){
  exe.name <- "bcftools"
  stopifnot(file.exists(Sys.which(exe.name)))
  if(is.matrix(samples))
    samples <- as.data.frame(samples)
  stopifnot(file.exists(in.vcf.file),
            is.data.frame(samples),
            ncol(samples) >= 2,
            all(! grepl(" ", samples[,1])),
            all(! grepl(" ", samples[,2])),
            all(! is.na(samples[,2])),
            out.vcf.file != in.vcf.file,
            is.logical(skip.check),
            is.logical(sort.variants))

  samples <- convertFactorColumnsToCharacter(samples)
  if(! skip.check){
    samples.init <- getSamplesFromVcfFile(in.vcf.file)
    if(nrow(samples) != length(samples.init)){
      msg <- "number of rows in 'samples' different than the number of samples in the VCF file"
      stop(msg)
    }
    if(any(! samples[,1] %in% samples.init)){
      msg <- "some old names in 'samples' are absent from the VCF file"
      stop(msg)
    }
    if(! all(samples[,1] == samples.init)){
      msg <- "some old names different than names in the VCF file"
      stop(msg)
    }
  }

  tmpf <- tempfile(pattern="renameVcfSamples", fileext=".txt")
  tmp.df <- data.frame(old_name=samples[,1],
                       new_name=samples[,2])
  utils::write.table(tmp.df, tmpf, quote=FALSE, sep=" ", row.names=FALSE)

  if(sort.variants){
    ## (1) convert vcf to bcf, (2) sort, (3) rename samples
    cmd <- paste0(exe.name, " view",
                  " -O u", # export as BCF
                  " ", in.vcf.file,
                  " | ", exe.name, " sort",
                  " -O z", # export as compressed VCF
                  " -",
                  " | ", exe.name, " reheader",
                  " -s ", tmpf,
                  " -o ", out.vcf.file, # will be compressed VCF
                  " -")
  } else{
    cmd <- paste0(exe.name, " reheader",
                  " -s ", tmpf,
                  " -o ", out.vcf.file,
                  " ", in.vcf.file)
  }
  if(file.exists(out.vcf.file))
    file.remove(out.vcf.file)
  ret <- system(cmd)
  if(ret != 0){
    msg <- paste0(exe.name, " returned ", ret)
    warning(msg)
  }

  file.remove(tmpf)

  invisible(out.vcf.file)
}

##' Parse VCF
##'
##' Return a matrix of genotypes for SNPs (with possibly multiple alternative alleles).
##' Phasing information is ignored.
##' With \href{http://grokbase.com/t/r/bioconductor/135b460s2b/bioc-how-to-convert-genotype-snp-matrix-to-nucleotide-genotypes}{help} from Martin Morgan.
##' @param vcf CollapsedVCF (see pkg \href{http://bioconductor.org/packages/VariantAnnotation/}{VariantAnnotation})
##' @param na.string a symbol to indicate missing genotypes (e.g. NA, "NN", "--", etc)
##' @param single.alt if TRUE, only records with a single 'alt' are kept
##' @return matrix with variants in rows and samples in columns
##' @author Gautier Sarah, Timothee Flutre
##' @seealso \code{\link{vcf2genoClasses}}, \code{\link{gtVcf2dose}}
##' @export
gtVcf2genoClasses <- function(vcf, na.string=NA, single.alt=TRUE){
  requireNamespaces(c("VariantAnnotation", "GenomicRanges"))
  stopifnot(is.logical(single.alt))

  if(single.alt)
    vcf <- subsetVcfOnAllelicity(vcf=vcf, single.ref=TRUE,
                                 single.alt=TRUE)

  ref <- VariantAnnotation::ref(vcf)
  ref <- as.character(ref)
  alt <- VariantAnnotation::alt(vcf)
  alt <- GenomicRanges::as.data.frame(alt)
  gt <- sub("\\|", "/", VariantAnnotation::geno(vcf)$GT) # ignore phasing

  idx.var <- (which(gt == "0/0")-1) %% nrow(gt) + 1
  gt[which(gt == "0/0")] <- paste0(ref[idx.var], ref[idx.var])

  idx.var <- (which(gt == "0/1")-1) %% nrow(gt) + 1
  gt[which(gt == "0/1")] <- paste0(ref[idx.var], alt[idx.var, 3])

  idx.var <- (which(gt == "1/1")-1) %% nrow(gt) + 1
  gt[which(gt == "1/1")] <- paste0(alt[idx.var, 3], alt[idx.var, 3])

  gt[which(gt %in% c(".","./."))] <- na.string

  if(! single.alt){
    gclasses <- names(table(gt))
    gclasses <- gclasses[grepl("/", gclasses)] # keep only the unmodified
    gclasses <- gclasses[! grepl("\\.", gclasses)] # discard the missing
    nb.alts <- tapply(alt$group, factor(alt$group), length)
    tmp <- data.frame(record=rep(NA, sum(nb.alts + 1)),
                      allele=NA,
                      stringsAsFactors=FALSE)
    tmp$record <- rep(as.integer(names(nb.alts)), nb.alts + 1)
    is.dupl <- duplicated(tmp$record)
    tmp$allele[! is.dupl] <- ref
    tmp$allele[is.dupl] <- alt$value
    for(gclass in gclasses){
      message(gclass)
      idx.geno <- which(gt == gclass)
      idx.rec <- (idx.geno - 1) %% nrow(gt) + 1
      gt[idx.geno] <- sapply(1:length(idx.geno), function(i){
        idx.alleles <- as.numeric(strsplit(gt[idx.geno[i]], "/")[[1]])
        paste0(tmp$allele[tmp$record == idx.rec[i]][idx.alleles + 1], collapse="")
      })
    }
  }

  return(gt)
}

##' Convert VCF to genotypic classes
##'
##' Convert genotypes at SNPs from a VCF file into genotypic classes.
##' @param vcf.file path to the VCF file (if the bgzip index doesn't exist in the same directory, it will be created)
##' @param genome genome identifier (e.g. "VITVI_12x2")
##' @param gclasses.file path to the output file to record genotypes into genotypic classes (will be gzipped if suffix is ".gz")
##' @param ca.file path to the output file to record SNP 1-based coordinates and alleles (will be gzipped if suffix is ".gz")
##' @param yieldSize number of records to yield each time the file is read from (see \code{?TabixFile}) if seq.id is NULL
##' @param dict.file path to the SAM dict file (see \url{https://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary}) if seq.id is specified with no start/end
##' @param seq.id see \code{\link{seqIdStartEnd2GRanges}}
##' @param seq.start see \code{\link{seqIdStartEnd2GRanges}}
##' @param seq.end see \code{\link{seqIdStartEnd2GRanges}}
##' @param na.string a symbol to indicate missing genotypes (e.g. NA, "NN", "--", etc)
##' @param single.alt if TRUE, only records with a single 'alt' are kept
##' @param verbose verbosity level (0/1)
##' @return invisible vector with the path to the output file
##' @author Gautier Sarah, Timothee Flutre
##' @seealso \code{\link{gtVcf2genoClasses}}, \code{\link{filterVariantCalls}}
##' @export
vcf2genoClasses <- function(vcf.file, genome="", gclasses.file, ca.file,
                            yieldSize=NA_integer_, dict.file=NULL,
                            seq.id=NULL, seq.start=NULL, seq.end=NULL,
                            na.string=NA, single.alt=TRUE, verbose=1){
  requireNamespaces(c("IRanges", "GenomicRanges", "VariantAnnotation",
                      "Rsamtools"))
  stopifnot(file.exists(vcf.file),
            xor(is.na(yieldSize), is.null(seq.id)))
  if(! is.null(seq.id) & is.null(seq.start) & is.null(seq.end))
    stopifnot(! is.null(dict.file),
              file.exists(dict.file))

  if(verbose > 0){
    msg <- "read VCF to convert genotypes into genotypic classes ..."
    write(msg, stdout()); flush(stdout())
  }

  for(out.file in c(gclasses.file, ca.file))
    if(file.exists(out.file))
      file.remove(out.file)
  gzip.gclasses <- ifelse(grepl("\\.gz$", gclasses.file), TRUE, FALSE)
  gclasses.file <- sub("\\.gz$", "", gclasses.file)
  gzip.ca <- ifelse(grepl("\\.gz$", ca.file), TRUE, FALSE)
  ca.file <- sub("\\.gz$", "", ca.file)
  gclasses.con <- file(gclasses.file, open="a")
  ca.con <- file(ca.file, open="a")
  cat("chr\tpos\tallele.ref\tallele.alt\n", file=ca.con, append=TRUE)

  if(! file.exists(paste0(vcf.file, ".tbi")))
    Rsamtools::indexTabix(file=vcf.file, format="vcf")
  tabix.file <- Rsamtools::TabixFile(file=vcf.file,
                                     yieldSize=yieldSize)
  if(! is.null(seq.id)){
    rngs <- seqIdStartEnd2GRanges(seq.id=seq.id, seq.start=seq.start,
                                  seq.end=seq.end, dict.file=dict.file)
    vcf.params <- VariantAnnotation::ScanVcfParam(which=rngs)
    vcf <- VariantAnnotation::readVcf(file=tabix.file, genome=genome,
                                      param=vcf.params)
    nb.variants <- nrow(vcf)
    gtmp <- gtVcf2genoClasses(vcf=vcf, na.string=na.string,
                              single.alt=single.alt)
    catmp <- rngVcf2df(vcf=vcf, with.coords=TRUE, with.alleles=TRUE,
                       single.ref=TRUE, single.alt=single.alt)
    cat(paste(colnames(gtmp), collapse="\t"), file=gclasses.con,
        append=TRUE, sep="\n")
    utils::write.table(x=gtmp,
                       file=gclasses.con, append=TRUE,
                       quote=FALSE, sep="\t", row.names=TRUE,
                       col.names=FALSE)
    utils::write.table(x=catmp, file=ca.con, append=TRUE,
                       quote=FALSE, sep="\t", row.names=TRUE,
                       col.names=FALSE)
  } else{
    open(tabix.file)
    nb.variants <- 0
    while(nrow(vcf <- VariantAnnotation::readVcf(file=tabix.file,
                                                 genome=genome))){
      gtmp <- gtVcf2genoClasses(vcf=vcf, na.string=na.string,
                                single.alt=single.alt)
      catmp <- rngVcf2df(vcf=vcf, with.coords=TRUE, with.alleles=TRUE,
                         single.ref=TRUE, single.alt=single.alt)
      if(nb.variants == 0)
        cat(paste(colnames(gtmp), collapse="\t"), file=gclasses.con,
            append=TRUE, sep="\n")
      nb.variants <- nb.variants + nrow(vcf)
      utils::write.table(x=gtmp, file=gclasses.con, append=TRUE,
                         quote=FALSE, sep="\t", row.names=TRUE,
                         col.names=FALSE)
      utils::write.table(x=catmp, file=ca.con, append=TRUE,
                         quote=FALSE, sep="\t", row.names=TRUE,
                         col.names=FALSE)
    }
    close(tabix.file)
  }

  if(verbose > 0){
    msg <- paste0("nb of variants: ", nb.variants)
    write(msg, stdout())
  }

  close(gclasses.con)
  close(ca.con)
  for(gz.out.file in c(paste0(gclasses.file, ".gz"),
                       paste0(ca.file, ".gz")))
    if(file.exists(gz.out.file))
      file.remove(gz.out.file)
  if(file.exists(Sys.which("gzip"))){ # should work on Linux and Mac OS
    if(gzip.gclasses)
      system(command=paste("gzip -n", gclasses.file))
    if(gzip.ca)
      system(command=paste("gzip -n", ca.file))
  } else{
    warning("'gzip' not available on this computer")
  }

  invisible(list(gclasses.file=paste0(gclasses.file, ".gz"),
                 ca.file=paste0(ca.file, ".gz")))
}

##' Missing genotypes in VCF
##'
##' Calculate the frequency of missing genotypes for each marker.
##' @param vcf CollapsedVCF (see pkg \href{http://bioconductor.org/packages/VariantAnnotation/}{VariantAnnotation})
##' @param with.coords if TRUE, a GRanges object is returned
##' @return numeric vector or GRanges
##' @author Timothee Flutre
##' @export
calcFreqNaVcf <- function(vcf, with.coords=FALSE){
  requireNamespaces(c("VariantAnnotation", "SummarizedExperiment", "S4Vectors"))

  gt <- VariantAnnotation::geno(vcf)$GT
  nb.samples <- ncol(gt)
  isMissing <- matrix(data=gt %in% c("./.", "."), nrow=nrow(gt),
                      ncol=ncol(gt), dimnames=dimnames(gt))
  freqNa <-  rowSums(isMissing) / nb.samples

  if(with.coords){
    gr <- SummarizedExperiment::rowRanges(vcf)
    S4Vectors::mcols(gr)$freqNa <- freqNa
    return(gr)
  } else
    return(freqNa)
}

##' Plot VCF
##'
##' Plot markers from a VCF file along chromosomes, colored according to their percentage of missing genotypes.
##' @param vcf.file path to the VCF file (if the bgzip index doesn't exist in the same directory, it will be created)
##' @param yieldSize number of records to yield each time the file is read from (see \code{?TabixFile}) if seq.id is NULL
##' @param dict.file path to the SAM dict file (see \url{https://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary}) if seq.id is specified with no start/end
##' @param seq.id sequence identifier to work on (e.g. "chr2")
##' @param seq.start start of the sequence to work on (if NULL, whole seq)
##' @param seq.end end of the sequence to work on (if NULL, whole seq)
##' @param main main title of the plot
##' @param plot.it if TRUE, \code{autoplot} from the ggbio package will be used
##' @return GenomicRange
##' @author Timothee Flutre
##' @export
plotVcfPercNa <- function(vcf.file, yieldSize=NA_integer_, dict.file=NULL,
                          seq.id=NULL, seq.start=NULL, seq.end=NULL,
                          main="Density of missing genotypes",
                          plot.it=TRUE){
  requireNamespaces(c("Rsamtools", "VariantAnnotation", "SummarizedExperiment",
                      "S4Vectors"))
  if(plot.it)
    requireNamespaces(c("ggplot2", "ggbio"))
  stopifnot(file.exists(vcf.file),
            xor(is.na(yieldSize), is.null(seq.id)))
  if(! is.null(seq.id) & is.null(seq.start) & is.null(seq.end))
    stopifnot(! is.null(dict.file),
              file.exists(dict.file))

  if(! file.exists(paste0(vcf.file, ".tbi")))
    Rsamtools::indexTabix(file=vcf.file, format="vcf")
  tabix.file <- Rsamtools::TabixFile(file=vcf.file,
                                     yieldSize=yieldSize)
  open(tabix.file)

  if(! is.null(seq.id)){
    rngs <- seqIdStartEnd2GRanges(seq.id=seq.id, seq.start=seq.start,
                                  seq.end=seq.end, dict.file=dict.file)
    vcf.params <- VariantAnnotation::ScanVcfParam(which=rngs)
    vcf <- VariantAnnotation::readVcf(file=tabix.file, genome="",
                                      param=vcf.params)
    gr <- calcFreqNaVcf(vcf=vcf, with.coords=TRUE)
  } else{
    vcf <- VariantAnnotation::readVcf(file=tabix.file, genome="")
    gr <- calcFreqNaVcf(vcf=vcf, with.coords=TRUE)
    vcf <- VariantAnnotation::readVcf(file=tabix.file, genome="")
    while(nrow(vcf)){
      gr <- c(gr, calcFreqNaVcf(vcf=vcf, with.coords=TRUE))
      vcf <- VariantAnnotation::readVcf(file=tabix.file, genome="")
    }
  }

  close(tabix.file)

  if(plot.it){
    percNa <- 100 * S4Vectors::mcols(gr)$freqNa # avoid NOTE in "R CMD check"
    S4Vectors::mcols(gr)$percNa <- percNa
    p <- ggbio::autoplot(gr, layout="karyogram",
                         ggplot2::aes(color=percNa, fill=percNa),
                         main=main)
    print(p)
  }

  invisible(gr)
}

##' Read PLINK
##'
##' Read output files from the \href{https://www.cog-genomics.org/plink2}{PLINK} software used with the \code{--mendel} option.
##' For potentially large files ("lmendel" and "mendel"), the \code{fread} function from the \href{https://cran.r-project.org/package=data.table}{data.table} package is used as it is faster than \code{\link[utils]{read.table}}.
##' @param prefix prefix of the file to read (e.g. \code{"~/work/output_plink"})
##' @param suffix suffix of the file to read
##' \itemize{
##' \item fmendel: one row per parental pair
##' \item imendel: one row per individual per parental pair
##' \item lmendel: one row per variant (SNP)
##' \item mendel: one row per error (Mendelian violation)
##' }
##' @param verbose verbosity level (0/1)
##' @return data frame
##' @author Timothee Flutre
##' @examples
##' \dontrun{## assuming properly-formatted data, launch PLINK via a system call
##' cmd <- "plink --mendel --bed input.bed --bim input.bim --fam input.fam --out output"
##' system(cmd)
##'
##' pl.mend <- readPlinkMendel(prefix="output", suffix="mendel")
##' str(pl.mend)
##'
##' pl.imend <- readPlinkMendel(prefix="output", suffix="imendel")
##' head(pl.imend[order(pl.imend$N, decreasing=TRUE),])
##' }
##' @export
readPlinkMendel <- function(prefix, suffix="", verbose=1){
  stopifnot(suffix %in% c("mendel", "imendel", "fmendel", "lmendel"))

  file <- paste0(prefix, ".", suffix)
  if(verbose > 0){
    msg <- paste0("read PLINK's output ", file, " ...")
    write(msg, stdout())
  }

  if(suffix == "mendel"){
    requireNamespace("data.table")
    line1 <- readLines(con=file, n=1)
    cn <- strsplit(line1, "\\s+")[[1]][-1]
    ## tmp <- utils::read.table(file, sep="", skip=1, stringsAsFactors=TRUE)
    tmp <- as.data.frame(
        data.table::fread(file, skip=1, stringsAsFactors=TRUE))
    data <- cbind(tmp[, 1:5],
                  t(t(apply(tmp[,6:ncol(tmp)], 1, paste, collapse=" "))),
                  stringsAsFactors=TRUE)
    colnames(data) <- cn
    data$CHR <- as.factor(data$CHR)
    data$SNP <- as.factor(data$SNP)
    data$CODE <- as.factor(data$CODE)
  } else if(suffix == "lmendel"){
    ## data <- utils::read.table(file, header=TRUE, sep="")
    data <- as.data.frame(data.table::fread(file))
    data$CHR <- as.factor(data$CHR)
    data$SNP <- as.factor(data$SNP)
  } else if(suffix == "imendel"){
    data <- utils::read.table(file, header=TRUE, sep="", stringsAsFactors=TRUE)
  } else if(suffix == "fmendel"){
    data <- utils::read.table(file, header=TRUE, sep="", stringsAsFactors=TRUE)
  }

  return(data)
}

##' BLAST
##'
##' Convert a data.frame containing alignments coordinates from BLAST into a GRanges object.
##' @param coords data.frame with 14 columns (see \code{\link{loadBlast}})
##' @return GRanges
##' @author Timothee Flutre
##' @seealso \code{\link{loadBlast}}
##' @examples
##' \dontrun{## BLAST should be run beforehand, see the example of `loadBlast`:
##' coords <- loadBlast("megablast_chroms_seqs.txt.gz", asGRanges=FALSE)
##' coords.gr <- blast2granges(coords)
##' }
##' @export
blast2granges <- function(coords){
  requireNamespaces(c("GenomicRanges", "S4Vectors", "IRanges"))
  stopifnot(is.data.frame(coords),
            ncol(coords) == 14,
            all(colnames(coords) == c("qseqid", "sseqid", "pident",
                                      "length", "mismatch", "gapopen",
                                      "qstart", "qend", "sstart", "send",
                                      "evalue", "bitscore", "qlen",
                                      "sstrand")))

  if(any(c("plus" %in% coords$sstrand, "minus" %in% coords$sstrand))){
    coords$sstrand <- gsub("plus", "+", coords$sstrand)
    coords$sstrand <- gsub("minus", "-", coords$sstrand)
  }

  gr <-
    GenomicRanges::GRanges(
        seqnames=S4Vectors::Rle(coords$sseqid),
        ranges=IRanges::IRanges(start=ifelse(coords$sstart <= coords$send,
                                             coords$sstart,
                                             coords$send),
                                end=ifelse(coords$sstart <= coords$send,
                                           coords$send,
                                           coords$sstart)),
        strand=ifelse(coords$sstart <= coords$send, "+", "-"))
  names(gr) <- coords$qseqid
  S4Vectors::mcols(gr) <- coords[, c("pident", "length", "mismatch", "gapopen",
                                     "qstart", "qend", "evalue", "bitscore",
                                     "qlen")]

  return(gr)
}

##' BLAST
##'
##' Load alignment coordinates from the \href{http://blast.ncbi.nlm.nih.gov/Blast.cgi}{NCBI BLAST}.
##' @param file.coords path to the file with alignment coordinates obtained with \code{-outfmt "6 std qlen sstrand"}
##' @param reformat.strand if TRUE, "plus" is replaced by "+", and "minus" by "-"
##' @param apply.sort if TRUE, return the alignments sorted in increasing order for qseqid, then increasing for evalue, then decreasing for bitscore, then decreasing for length, then decreasing for pident, then increasing for sseqid
##' @param asGRanges if TRUE, returns a \code{\link[GenomicRanges]{GRanges}} object
##' @param verbose verbosity level (0/1)
##' @return GRanges or data.frame with 14 columns
##' @author Timothee Flutre
##' @seealso \code{\link{blast2granges}}
##' @examples
##' \dontrun{## BLAST should be run beforehand, for instance:
##' ## zcat chroms.fa.gz | makeblastdb -dbtype nucl -in - -title chroms -out chroms -logfile makeblastdb_chroms.log
##' ## echo "zcat seqs.fa.gz | blastn -query - -task megablast -db chroms -out /dev/stdout -outfmt \"6 std qlen sstrand\" | gzip > megablast_chroms_seqs.txt.gz"  | qsub ...
##' coords <- loadBlast("megablast_chroms_seqs.txt.gz")
##' }
##' @export
loadBlast <- function(file.coords, reformat.strand=TRUE, apply.sort=TRUE,
                      asGRanges=TRUE, verbose=1){
  stopifnot(file.exists(file.coords),
            is.logical(reformat.strand),
            is.logical(apply.sort))

  coords <- utils::read.table(file.coords, sep="\t", stringsAsFactors=FALSE,
                              col.names=c("qseqid", "sseqid", "pident",
                                          "length", "mismatch", "gapopen",
                                          "qstart", "qend", "sstart", "send",
                                          "evalue", "bitscore", "qlen", "sstrand"))
  if(verbose > 0){
    msg <- paste0("nb of alignments: ", nrow(coords),
                  "\nnb of queries: ", length(unique(coords[,"qseqid"])),
                  "\nnb of subjects: ", length(unique(coords[,"sseqid"])))
    write(msg, stdout())
  }

  if(reformat.strand){
    coords$sstrand <- gsub("plus", "+", coords$sstrand)
    coords$sstrand <- gsub("minus", "-", coords$sstrand)
  }

  if(apply.sort)
    coords <- coords[order(coords$qseqid,
                           coords$evalue,
                           -coords$bitscore,
                           -coords$length,
                           -coords$pident,
                           coords$sseqid), ]

  if(asGRanges)
    coords <- blast2granges(coords)

  return(coords)
}

##' MUMmer
##'
##' Convert a data.frame containing alignments coordinates from MUMmer into a GRanges object.
##' @param coords data.frame with 13 columns (see \code{\link{loadMummer}})
##' @return GRanges
##' @author Timothee Flutre
##' @seealso \code{\link{loadMummer}}
##' @examples
##' \dontrun{## MUMmer should be run beforehand, see the example of `loadMummer`:
##' coords <- loadMummer("out-nucmer_filter_coords.txt.gz", asGRanges=FALSE)
##' coords.gr <- mummer2granges(coords)
##' }
##' @export
mummer2granges <- function(coords){
  requireNamespaces(c("GenomicRanges", "S4Vectors", "IRanges"))
  stopifnot(is.data.frame(coords),
            ncol(coords) == 13,
            all(colnames(coords) == c("ref.start", "ref.end",
                                      "qry.start", "qry.end",
                                      "ref.aln.len", "qry.aln.len",
                                      "perc.id",
                                      "ref.len", "qry.len",
                                      "ref.cov", "qry.cov",
                                      "ref", "qry")))

  gr <-
    GenomicRanges::GRanges(
        seqnames=S4Vectors::Rle(coords$ref),
        ranges=IRanges::IRanges(start=ifelse(coords$ref.start <= coords$ref.end,
                                             coords$ref.start,
                                             coords$ref.end),
                                end=ifelse(coords$ref.start <= coords$ref.end,
                                           coords$ref.end,
                                           coords$ref.start)),
        strand=ifelse(coords$ref.start <= coords$ref.end, "+", "-"))
  names(gr) <- coords$qry
  S4Vectors::mcols(gr) <- coords[, c("qry.start", "qry.end", "perc.id",
                                     "ref.len", "ref.cov", "qry.cov")]

  return(gr)
}

##' MUMmer
##'
##' Load alignment coordinates from \href{http://mummer.sourceforge.net/}{MUMmer}.
##' @param file.coords path to the file with alignment coordinates obtained via \code{show-coords} (see the example)
##' @param algo nucmer or promer (the latter isn't available yet)
##' @param asGRanges if TRUE, returns a \code{\link[GenomicRanges]{GRanges}} object
##' @param keep.nlong.qry keep the N queries having the longest cumulated alignments on the reference (keep all of them by default)
##' @param ref.lim limits on reference coordinates outside which alignments are discarded (keep all of them by default)
##' @param verbose verbosity level (0/1)
##' @return GRanges or data.frame
##' @author Timothee Flutre
##' @seealso \code{\link{mummer2granges}}
##' @examples
##' \dontrun{## MUMmer should be run beforehand, for instance:
##' ## nucmer --maxmatch -p out-nucmer chroms.fa seqs.fa
##' ## delta-filter -l 1000 -q out-nucmer.delta > out-nucmer_filter.delta
##' ## show-coords -c -l -L 1000 -r -T out-nucmer_filter.delta | gzip > out-nucmer_filter_coords.txt.gz
##' coords <- loadMummer("out-nucmer_filter_coords.txt.gz")
##' library(GenomicRanges)
##' loc <- GRanges(seqnames=Rle("chr1"),
##'                ranges=IRanges(start=5, end=17))
##' idx <- subjectHits(findOverlaps(query=loc, subject=coords.gr))
##' plotAligns(coords=as.data.frame(ranges(coords.gr[idx,])))
##' }
##' @export
loadMummer <- function(file.coords, algo="nucmer", asGRanges=TRUE,
                       keep.nlong.qry=NULL, ref.lim=NULL, verbose=1){
  if(asGRanges)
    requireNamespaces("GenomicRanges")
  stopifnot(file.exists(file.coords),
            algo %in% c("nucmer"))
  if(! is.null(keep.nlong.qry))
    stopifnot(is.numeric(keep.nlong.qry),
              length(keep.nlong.qry) == 1)
  if(! is.null(ref.lim))
    stopifnot(is.numeric(ref.lim),
              length(ref.lim) == 2)

  coords <- utils::read.table(file.coords, skip=4, stringsAsFactors=FALSE)
  stopifnot(ncol(coords) == 13)
  colnames(coords) <- c("ref.start", "ref.end",
                        "qry.start", "qry.end",
                        "ref.aln.len", "qry.aln.len",
                        "perc.id",
                        "ref.len", "qry.len",
                        "ref.cov", "qry.cov",
                        "ref", "qry")
  if(verbose > 0){
    msg <- paste0("nb of references: ", length(unique(coords[,"ref"])),
                  "\nnb of queries: ", length(unique(coords[,"qry"])),
                  "\nnb of alignments: ", nrow(coords))
    write(msg, stdout())
  }

  if(! is.null(keep.nlong.qry) &&
     keep.nlong.qry < length(unique(coords[,"qry"]))){
    cumlen <- tapply(coords[,"qry.aln.len"], factor(coords[,"qry"]), sum)
    cumlen <- sort(cumlen, decreasing=TRUE)
    coords <- coords[coords[,"qry"] %in% names(cumlen)[1:keep.nlong.qry],]
    if(verbose > 0){
      msg <- paste0("after query filtering:",
                    "\nnb of references: ", length(unique(coords[,"ref"])),
                    "\nnb of queries: ", length(unique(coords[,"qry"])),
                    "\nnb of alignments: ", nrow(coords))
      write(msg, stdout())
    }
  }

  if(! is.null(ref.lim)){
    ref.strand.plus <- coords[,"ref.start"] < coords[,"ref.end"]
    idx.plus <- which(coords[ref.strand.plus, "ref.start"] >= ref.lim[1] &
                      coords[ref.strand.plus, "ref.end"] <= ref.lim[2])
    idx.minus <- which(coords[!ref.strand.plus, "ref.end"] >= ref.lim[1] &
                       coords[!ref.strand.plus, "ref.start"] <= ref.lim[2])
    coords <- coords[c(idx.plus, idx.minus),]
    if(verbose > 0){
      msg <- paste0("after alignment filtering:",
                    "\nnb of references: ", length(unique(coords[,"ref"])),
                    "\nnb of queries: ", length(unique(coords[,"qry"])),
                    "\nnb of alignments: ", nrow(coords))
      write(msg, stdout())
    }
  }

  if(asGRanges)
    coords <- mummer2granges(coords)

  return(coords)
}

##' Plot alignments
##'
##' Plot alignments of several queries on the same reference.
##' @param coords data.frame with at least three columns, "start", "end" and "names"
##' @param main main title
##' @param xlab label of the x-axis
##' @param xlim x-axis limits
##' @param col segment color(s)
##' @return nothing
##' @author Timothee Flutre
##' @seealso \code{\link{plotGRanges}}
##' @examples
##' \dontrun{coords <- data.frame(start=c(2, 21, 29, 50),
##'                      end=c(10, 25, 45, 53),
##'                      names=c("seq1", "seq2", "seq2", "seq3"),
##'                      stringsAsFactors=FALSE)
##' par(mar=c(4, 7, 3, 1))
##' plotAligns(coords, xlab="chr1", col=c(1, 1, 2, 1))
##' }
##' @export
plotAligns <- function(coords, main="Alignments", xlab="reference", xlim=NULL,
                       col="black"){
  stopifnot(is.data.frame(coords),
            all(c("start", "end", "names")
                %in% colnames(coords)))

  ## identify alignment strands
  ref.strand.plus <- coords[,"start"] < coords[,"end"]

  ## determine axes limits
  if(is.null(xlim))
    xlim <- c(min(coords[ref.strand.plus, "start"],
                  coords[!ref.strand.plus, "end"]),
              max(coords[ref.strand.plus, "end"],
                  coords[!ref.strand.plus, "start"]))
  ylim <- c(1, length(unique(coords[,"names"])))

  ## plot empty box
  graphics::plot(x=0, y=0, type="n", xlim=xlim, ylim=ylim, yaxt="n",
                 main=main, xlab=xlab, ylab="")
  graphics::axis(side=2, at=seq(ylim[1], ylim[2], 1),
                 labels=unique(coords[,"names"]), las=1)

  ## add alignments
  y <- match(coords[,"names"], unique(coords[,"names"]))
  graphics::segments(x0=coords[,"start"], y0=y,
                     x1=coords[,"end"], y1=y,
                     col=col)
}

##' Invert GRanges
##'
##' Transform a GRanges by inverting references and queries.
##' @param in.gr GRanges
##' @return GRanges which references are gr.in's queries, and queries are gr.in's references.
##' @author Timothee Flutre
##' @export
invertGRanges <- function(in.gr){
  requireNamespaces(c("S4Vectors", "BiocGenerics", "IRanges", "GenomicRanges",
                      "GenomeInfoDb"))
  stopifnot(all(c("qry.start", "qry.end") %in%
                colnames(S4Vectors::mcols(in.gr))))

  ## get input query starts, ends and strands
  in.qry.start <- S4Vectors::mcols(in.gr)[, "qry.start"]
  in.qry.end <- S4Vectors::mcols(in.gr)[, "qry.end"]
  in.qry.strand <- rep("*", length(in.qry.start))
  in.qry.strand[in.qry.start < in.qry.end] <- "+"
  in.qry.strand[in.qry.start > in.qry.end] <- "-"

  ## make output reference starts, ends and strands
  out.ref.start <- in.qry.start
  out.ref.start[in.qry.strand == "-"] <- in.qry.end[in.qry.strand == "-"]
  out.ref.end <- in.qry.end
  out.ref.end[in.qry.strand == "-"] <- in.qry.start[in.qry.strand == "-"]

  ## make output GRanges
  out.gr <- GenomicRanges::GRanges(
      seqnames=S4Vectors::Rle(names(in.gr)),
      ranges=IRanges::IRanges(start=out.ref.start,
                              end=out.ref.end),
      strand=in.qry.strand)
  names(out.gr) <- as.character(GenomeInfoDb::seqnames(in.gr))

  ## get input reference starts, ends and strands
  in.ref.start <- BiocGenerics::start(in.gr)
  in.ref.end <- BiocGenerics::end(in.gr)
  in.ref.strand <- as.character(BiocGenerics::strand(in.gr))

  ## make output query starts and ends
  out.qry.start <- in.ref.start
  out.qry.start[in.ref.strand == "-"] <- in.ref.end[in.ref.strand == "-"]
  out.qry.end <- in.ref.end
  out.qry.end[in.ref.strand == "-"] <- in.ref.start[in.ref.strand == "-"]

  ## make mcols of output GRanges
  S4Vectors::mcols(out.gr) <- data.frame("qry.start"=out.qry.start,
                                         "qry.end"=out.qry.end)

  return(out.gr)
}

##' Plot GRanges
##'
##' Plot GRanges of several queries on the same reference, with one y-axis line per query. Each alignment can have a different color, specified via the "col" column of \code{mcols(gr)}. If this column doesn't exist, it will be black by default.
##' @param gr GRanges
##' @param main main title
##' @param xlab label of the x-axis (by default, will be the first level of \code{seqnames(gr)})
##' @param xlim x-axis limits
##' @param col.qry.lab color for the query labels used along the y-axis
##' @param shape shape used to represent the alignments (segments/arrows)
##' @return nothing
##' @author Timothee Flutre
##' @seealso \code{\link{plotAligns}}
##' @examples
##' \dontrun{## make fake coordinates to be plotted with a single color
##' library(GenomicRanges)
##' gr <- GRanges(seqnames=c("chr1"),
##'               ranges=IRanges(start=c(1, 15, 6), end=c(10, 20, 12)),
##'               strand=c("+", "-", "+"))
##' names(gr) <- c("qry1", "qry1", "qry2")
##' plotGRanges(gr)
##' plotGRanges(gr, shape="arrows")
##'
##' ## use color to distinguish different alignments of the same query
##' mcols(gr)["col"] <- c("black", "red", "black")
##' plotGRanges(gr, shape="arrows")
##'
##' ## outside of this function, one can always add other things to the plot
##' legend("right", legend=c("a", "b"), col=c("black", "red"), lwd=2, bty="n")
##' }
##' @export
plotGRanges <- function(gr, main="Alignments", xlab=NULL, xlim=NULL,
                        col.qry.lab="black", shape="segments"){
  requireNamespaces(c("BiocGenerics", "GenomeInfoDb", "S4Vectors"))
  GenomeInfoDb::seqlevels(gr) <- GenomeInfoDb::seqlevelsInUse(gr)
  stopifnot(nlevels(GenomeInfoDb::seqnames(gr)) == 1,
            shape %in% c("segments", "arrows"))

  ## determine axes
  if(is.null(xlab))
    xlab <- levels(GenomeInfoDb::seqnames(gr))
  if(is.null(xlim))
    xlim <- range(BiocGenerics::start(gr), BiocGenerics::end(gr))
  ylim <- c(1, length(unique(names(gr))))

  ## plot empty box
  graphics::plot(x=0, y=0, type="n", xlim=xlim, ylim=ylim, yaxt="n",
                 main=main, xlab=xlab, ylab="")
  if(length(col.qry.lab) == 1){
    graphics::axis(side=2, at=seq(ylim[1], ylim[2], 1),
                   labels=unique(names(gr)), las=1)
  } else{
    for(i in seq_along(unique(names(gr)))){
      graphics::axis(side=2, at=i, labels=unique(names(gr))[i], las=1,
                     col=col.qry.lab[i], col.axis=col.qry.lab[i])
    }
  }

  ## add alignments
  if(shape == "segments"){
    y <- match(names(gr), unique(names(gr)))
    if("col" %in% colnames(S4Vectors::mcols(gr))){
      col <- S4Vectors::mcols(gr)$col
    } else
      col <- "black"
    graphics::segments(x0=BiocGenerics::start(gr), y0=y,
                       x1=BiocGenerics::end(gr), y1=y,
                       col=col)
  } else if(shape == "arrows"){
    is.s <- BiocGenerics::strand(gr) == "+"
    if(any(is.s)){
      y <- match(names(gr[is.s]), unique(names(gr)))
      if("col" %in% colnames(S4Vectors::mcols(gr))){
        col <- S4Vectors::mcols(gr[is.s])$col
      } else
        col <- "black"
      graphics::arrows(x0=BiocGenerics::start(gr[is.s]), y0=y,
                       x1=BiocGenerics::end(gr[is.s]), y1=y,
                       col=col)
    }
    is.s <- BiocGenerics::strand(gr) == "-"
    if(any(is.s)){
      y <- match(names(gr[is.s]), unique(names(gr)))
      if("col" %in% colnames(S4Vectors::mcols(gr))){
        col <- S4Vectors::mcols(gr[is.s])$col
      } else
        col <- "black"
      graphics::arrows(x0=BiocGenerics::end(gr[is.s]), y0=y,
                       x1=BiocGenerics::start(gr[is.s]), y1=y,
                       col=col)
    }
    is.s <- BiocGenerics::strand(gr) == "*"
    if(any(is.s)){
      y <- match(names(gr[is.s]), unique(names(gr)))
      if("col" %in% colnames(S4Vectors::mcols(gr))){
        col <- S4Vectors::mcols(gr[is.s])$col
      } else
        col <- "black"
      graphics::segments(x0=BiocGenerics::start(gr[is.s]), y0=y,
                         x1=BiocGenerics::end(gr[is.s]), y1=y,
                         col=col)
    }
  }
}
