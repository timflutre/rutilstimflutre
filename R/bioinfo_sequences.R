## Contains functions handling sequences often used in bioinformatics.

##' Rename chromosomes
##'
##' Rename chromosomes into integers.
##' @param x vector of chromosome names
##' @return data.frame with original and new names
##' @author Timothee Flutre
##' @examples
##' \dontrun{chroms <- c("chr1", "chr1_random", "chr10", "chr10_random", "chrUn", "chr2")
##' chromNames2integers(x=chroms)
##' }
##' @export
chromNames2integers <- function(x){
  stopifnot(is.character(x),
            is.vector(x))

  output <- data.frame(original=x,
                       renamed=NA,
                       stringsAsFactors=FALSE)

  output$renamed <- suppressWarnings(as.integer(gsub("chr", "", x)))

  if(any(is.na(output$renamed))){
    max.chr.int <- max(output$renamed, na.rm=TRUE)

    tmp <- data.frame(orig=x[is.na(output$renamed)],
                      idx=which(is.na(output$renamed)),
                      as.chr=NA,
                      as.int=NA,
                      stringsAsFactors=FALSE)
    tmp$as.chr <- gsub("chr|_random", "", tmp$orig)
    tmp <- tmp[order(tmp$as.chr),]

    for(i in 1:nrow(tmp)){
      suppressWarnings(tmp$as.int[i] <- as.integer(tmp$as.chr[i]))
      if(! is.na(tmp$as.int[i])){ # e.g. chr1_random
        output$renamed[tmp$idx[i]] <- max.chr.int + tmp$as.int[i]
      } else{ # e.g. chrUn
        output$renamed[tmp$idx[i]] <- 2 * max.chr.int + 1
      }
    }
  }

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
    write(paste0("load ", fa.file, " ..."), stdout()); flush(stdout())
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
    msg <- "convert subsequence coordinates into a RangesList ..."
    write(msg, stdout())
    utils::flush.console()
  }
  sub.info.rl <- IRanges::RangesList()
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

##' Load read counts
##'
##' Load read counts per individual and lane
##' @param lanes.dir vector of paths to "lane" directories, each containing a "demultiplex" directory in which \href{https://github.com/timflutre/quantgen/blob/master/demultiplex.py}{demultiplex.py} was executed
##' @return data.frame with (individual,lane) in rows and counts in columns
##' @author Timothee Flutre
##' @export
loadReadCountsPerIndAndLane <- function(lanes.dir){
  stopifnot(is.vector(lanes.dir),
            length(lanes.dir) > 0,
            is.character(lanes.dir))

  reads <- list()

  for(lane.dir in lanes.dir){
    lane <- gsub("lane_", "", basename(lane.dir))
    lane.file <- paste0(lane.dir, "/demultiplex/",
                        lane, "_stats-demultiplex.txt.gz")
    if(! file.exists(lane.file))
      stop(paste0("file ", lane.file, " doesn't exist"))
    reads[[lane]] <- utils::read.table(lane.file, header=TRUE)
    if(! all(colnames(reads[[lane]] %in% c("ind","barcode","assigned","lane"))))
      warning(paste0("look at header line of file '", lane.file, "'"))
  }

  reads <- do.call(rbind, reads)
  reads$lane <- as.factor(sapply(strsplit(rownames(reads), "\\."),
                                 function(x){x[1]}))

  return(reads)
}

##' Reformat read counts per lane
##'
##' Reformat the data.frame from \code{\link{loadReadCountsPerIndAndLane}} in a matrix of counts
##' @param x data.frame
##' @return matrix
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
      if(ind %in% x$ind[x$lane == lane])
        counts[lane, ind] <- x$assigned[x$ind == ind &
                                          x$lane == lane]
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

##' Read VCF
##'
##' Read a subset of a VCF file
##' @param vcf.file path to the VCF file (bgzip index should exist in same directory)
##' @param genome genome identifier (e.g. "VITVI_12x2")
##' @param seq.id sequence identifier to work on (e.g. "chr2")
##' @param seq.start start of the sequence to work on (1-based coordinate)
##' @param seq.end end of the sequence to work on (1-based coordinate)
##' @return CollapsedVCF (see pkg VariantAnnotation)
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
##' Reads a VCF file and returns the number of alternate alleles over all records (via \code{\link[base]{table}}).
##' @param vcf.file path to the VCF file
##' @param genome genome identifier (e.g. "VITVI_12x2")
##' @param verbose verbosity level (0/1)
##' @return table
##' @author Timothee Flutre
##' @export
tableVcfAlt <- function(vcf.file, genome="", verbose=1){
  requireNamespaces(c("VariantAnnotation", "S4Vectors"))
  stopifnot(file.exists(vcf.file))

  if(verbose > 0){
    msg <- "read VCF..."
    write(msg, stdout())
  }
  svp <- VariantAnnotation::ScanVcfParam(fixed="ALT", info=NA, geno=NA)
  vcf <- VariantAnnotation::readVcf(file=vcf.file, genome=genome, param=svp)
  if(verbose > 0){
    msg <- paste0("nb of records: ", nrow(vcf))
    write(msg, stdout())
  }

  return(table(S4Vectors::elementNROWS(VariantAnnotation::alt(vcf))))
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

##' Confidence in one variant's genotypes
##'
##' Provide measure of confidence in the genotypes at a given variant.
##' @param x CollapsedVCF corresponding to a SNV (see pkg VariantAnnotation)
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
##' Set genotypes (GT field) to missing (".") if the genotype quality (GQ field) isn't high enough.
##' @param vcf.file path to the VCF file (bgzip index should exist in same directory; should be already filtered for INFO tags such as QD, FS, MQ, etc)
##' @param genome genome identifier (e.g. "VITVI_12x2")
##' @param out.file path to the output VCF file (a bgzip index will be created in the same directory)
##' @param yieldSize number of records to yield each time the file is read from (see ?TabixFile) if seq.id is NULL
##' @param dict.file path to the SAM dict file (see \url{https://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary}) if seq.id is specified with no start/end
##' @param seq.id sequence identifier to work on (e.g. "chr2")
##' @param seq.start start of the sequence to work on (if NULL, whole seq)
##' @param seq.end end of the sequence to work on (if NULL, whole seq)
##' @param min.gq minimum GQ below which GT is set to "."
##' @param verbose verbosity level (0/1)
##' @return the destination file path as an invisible character(1)
##' @author Timothee Flutre
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
    stopifnot(! is.null(dict),
              file.exists(dict.file))
  out.file <- sub("\\.gz$", "", out.file)
  out.file <- sub("\\.bgz$", "", out.file)

  dest <- NULL

  tabix.file <- Rsamtools::TabixFile(file=vcf.file,
                                     yieldSize=yieldSize)
  if(! is.null(seq.id)){
    if(is.null(seq.start) & is.null(seq.end)){
      dict <- readSamDict(file=dict.file)
      if(! seq.id %in% rownames(dict)){
        msg <- paste0("seq '", seq.id, "' not in '", dict.file, "'")
        stop(msg)
      }
      rngs <- GenomicRanges::GRanges(seqnames=c(seq.id),
                                     ranges=IRanges::IRanges(start=c(1),
                                         end=c(dict$LN[rownames(dict) == seq.id])))
    } else
      rngs <- GenomicRanges::GRanges(seqnames=c(seq.id),
                                     ranges=IRanges::IRanges(start=c(seq.start),
                                         end=c(seq.end)))
    names(rngs) <- c(seq.id)
    vcf.params <- VariantAnnotation::ScanVcfParam(which=rngs)
    vcf <- VariantAnnotation::readVcf(file=tabix.file, genome=genome,
                                      param=vcf.params)
    nb.variants <- nrow(vcf)

    idx <- (VariantAnnotation::geno(vcf)[["GQ"]] < min.gq)
    VariantAnnotation::geno(vcf)[["GT"]][idx] <- "."

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
      VariantAnnotation::geno(vcf)[["GT"]][idx] <- "."
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
##' @param vcf.file path to the VCF file (bgzip index should exist in same directory; should be already filtered for INFO tags such as QD, FS, MQ, etc)
##' @param genome genome identifier (e.g. "VITVI_12x2")
##' @param out.file path to the output VCF file (a bgzip index will be created in the same directory)
##' @param yieldSize number of records to yield each time the file is read from (see ?TabixFile) if seq.id is NULL
##' @param dict.file path to the SAM dict file (see \url{https://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary}) if seq.id is specified with no start/end
##' @param seq.id sequence identifier to work on (e.g. "chr2")
##' @param seq.start start of the sequence to work on (if NULL, whole seq)
##' @param seq.end end of the sequence to work on (if NULL, whole seq)
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
                               is.snv=NULL, is.biall=NULL,
                               min.var.dp=NULL, max.var.dp=NULL,
                               min.alt.af=NULL, max.alt.af=NULL,
                               min.spl.dp=NULL, min.perc.spl.dp=NULL,
                               min.spl.gq=NULL, min.perc.spl.gq=NULL,
                               max.var.nb.gt.na=NULL, max.var.perc.gt.na=NULL,
                               verbose=1){
  requireNamespaces(c("IRanges", "GenomicRanges", "VariantAnnotation",
                      "Rsamtools", "S4Vectors", "BiocInstaller"))
  stopifnot(file.exists(vcf.file))
  if(! is.null(seq.id) & is.null(seq.start) & is.null(seq.end))
    stopifnot(! is.null(dict),
              file.exists(dict.file))
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

  out.file <- sub("\\.gz$", "", out.file)
  out.file <- sub("\\.bgz$", "", out.file)

  ##' @return TRUE if SNV (single nucleotide variant)
  filterSnv <- function(x){
    (VariantAnnotation::isSNV(x))
  }

  ##' @return TRUE if at most one alternate allele
  filterBiall <- function(x){
    if(utils::compareVersion(as.character(BiocInstaller::biocVersion()),
                             "3.4") < 0){
      (S4Vectors::elementLengths(VariantAnnotation::alt(x)) <= 1)
    } else
      (S4Vectors::elementNROWS(VariantAnnotation::alt(x)) <= 1)
  }

  ##' @return TRUE if variant-level DP inside of given range
  filterVariantDp <- function(x, min.dp=min.var.dp, max.dp=max.var.dp){
    if(nrow(x) == 0){
      logical(0)
    } else{
      dp <- VariantAnnotation::info(x)$DP
      (dp >= min.dp) || (dp <= max.dp)
    }
  }

  ##' @return TRUE if allele frequency inside of given range
  filterAf <- function(x, min.af=min.alt.af, max.af=max.alt.af){
    if(nrow(x) == 0){
      logical(0)
    } else{
      af <- unlist(VariantAnnotation::info(x)$AF)
      (af >= min.af) & (af <= max.af)
    }
  }

  ##' @return TRUE if high-enough percentage of samples with DP above threshold
  filterSampleDp <- function(x, min.dp=min.spl.dp,
                             min.perc.dp=min.perc.spl.dp){
    if(nrow(x) == 0){
      logical(0)
    } else{
      dp <- VariantAnnotation::geno(x)$DP
      (100 * rowSums(dp >= min.dp) / ncol(dp) >= min.perc.dp)
    }
  }

  ##' @return TRUE if high-enough percentage of samples with GQ above threshold
  filterSampleGq <- function(x, min.gq=min.spl.gq,
                             min.perc.gq=min.perc.spl.gq){
    if(nrow(x) == 0){
      logical(0)
    } else{
      gq <- VariantAnnotation::geno(x)$GQ
      (100 * rowSums(gq >= min.gq) / ncol(gq) >= min.perc.gq)
    }
  }

  ##' @return TRUE if high-enough number of samples without missing genotypes
  filterGtNb <- function(x, max.nb.gt.na=max.var.nb.gt.na){
    if(nrow(x) == 0){
      logical(0)
    } else{
      gt <- VariantAnnotation::geno(x)$GT
      (rowSums(gt == ".") <= max.nb.gt.na)
    }
  }

  ##' @return TRUE if high-enough percentage of samples without missing genotypes
  filterGtPerc <- function(x, max.perc.gt.na=max.var.perc.gt.na){
    if(nrow(x) == 0){
      logical(0)
    } else{
      gt <- VariantAnnotation::geno(x)$GT
      (100 * rowSums(gt == ".") / ncol(gt) <= max.perc.gt.na)
    }
  }

  ## set the filters
  tmp <- list()
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
    tabix.file <- Rsamtools::TabixFile(file=vcf.file,
                                       yieldSize=yieldSize)
    if(! is.null(seq.id)){
      if(is.null(seq.start) & is.null(seq.end)){
        dict <- readSamDict(file=dict.file)
        if(! seq.id %in% rownames(dict)){
          msg <- paste0("seq '", seq.id, "' not in '", dict.file, "'")
          stop(msg)
        }
        rngs <- GenomicRanges::GRanges(
            seqnames=c(seq.id),
            ranges=IRanges::IRanges(start=c(1),
                                    end=c(dict$LN[rownames(dict) == seq.id])))
      } else
        rngs <- GenomicRanges::GRanges(
            seqnames=c(seq.id),
            ranges=IRanges::IRanges(start=c(seq.start),
                                    end=c(seq.end)))
      names(rngs) <- c(seq.id)
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

  invisible(dest)
}

##' Summary per variant
##'
##' Compute the mean, sd, min, Q1, med, mean, Q3, max of the genotype qualities per variant, also reporting the number of samples and the number of missing data.
##' @param vcf CollapsedVCF (see pkg VariantAnnotation)
##' @param field genotype field of the VCF to parse (GQ/DP)
##' @return matrix with one row per variant
##' @author Timothee Flutre
##' @export
varqual2summary <- function(vcf, field="GQ"){
  requireNamespace("VariantAnnotation")

  output <- c()

  mat <- VariantAnnotation::geno(vcf)[[field]]

  output <- cbind(n=rep(ncol(mat), nrow(mat)),
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
  rownames(output) <- rownames(mat)

  return(output)
}

##' Summary per variant
##'
##' Compute the mean, sd, min, Q1, med, mean, Q3, max of the genotype qualities per variant, also reporting the number of samples and the number of missing data.
##' @param vcf.file path to the VCF file
##' @param genome genome identifier (e.g. "VITVI_12x2")
##' @param yieldSize number of records to yield each time the file is read from (see ?TabixFile) if seq.id is NULL
##' @param seq.id sequence identifier to work on (e.g. "chr2")
##' @param seq.start start of the sequence to work on
##' @param seq.end end of the sequence to work on
##' @param field genotype field of the VCF to parse (GQ/DP)
##' @param verbose verbosity level (0/1)
##' @return matrix with one row per variant and 9 columns (n, na, mean, sd, min, q1, med, q3, max)
##' @author Timothee Flutre
##' @export
summaryVariant <- function(vcf.file, genome, yieldSize=NA_integer_,
                           seq.id=NULL, seq.start=NULL, seq.end=NULL,
                           field="GQ", verbose=1){
  requireNamespaces(c("IRanges", "GenomicRanges", "VariantAnnotation",
                      "Rsamtools"))
  stopifnot(file.exists(vcf.file),
            field %in% c("GQ", "DP"))
  if(! is.null(seq.id))
    stopifnot(all(! is.null(seq.start), ! is.null(seq.end)))

  output <- NULL

  tabix.file <- Rsamtools::TabixFile(file=vcf.file,
                                     yieldSize=yieldSize)
  if(! is.null(seq.id)){
    rngs <- GenomicRanges::GRanges(seqnames=c(seq.id),
                                   ranges=IRanges::IRanges(start=c(seq.start),
                                                           end=c(seq.end)))
    names(rngs) <- c(seq.id)
    vcf.params <- VariantAnnotation::ScanVcfParam(which=rngs)
    vcf <- VariantAnnotation::readVcf(file=tabix.file, genome=genome,
                                      param=vcf.params)
    nb.variants <- nrow(vcf)
    output <- varqual2summary(vcf, field)
  } else{
    open(tabix.file)
    while(nrow(vcf <- VariantAnnotation::readVcf(file=tabix.file,
                                                 genome=genome))){
      if(is.null(output)){
        output <- varqual2summary(vcf, field)
      } else
        output <- rbind(output, varqual2summary(vcf, field))
    }
    close(tabix.file)
  }

  if(verbose > 0){
    msg <- paste0("nb of variants: ", nrow(output))
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
##' @param dict.file path to the SAM dict file (see \url{https://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary}) if seq.id is specified with no start/end
##' @return \code{GRanges} object from the \code{GenomicRanges} package
##' @author Timothee Flutre
##' @seealso \code{\link{vcf2dosage}}, \code{\link{vcf2genoClasses}}
##' @export
seqIdStartEnd2GRanges <- function(seq.id, seq.start=NULL, seq.end=NULL,
                                  dict.file=NULL){
  requireNamespace("GenomicRanges")
  if(is.null(seq.start) & is.null(seq.end))
    stopifnot(! is.null(dict.file),
              file.exists(dict.file))

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
  names(rngs) <- c(seq.id)

  return(rngs)
}

##' Keep single-element VCF records
##'
##' Returns a vcf object in which all multi-element 'ref' or 'alt' records and indels were discarded.
##' @param vcf CollapsedVCF (see pkg VariantAnnotation)
##' @return CollapsedVCF
##' @author Gautier Sarah, Timothee Flutre
##' @export
getVcfOnlySingleElem <- function(vcf){
  requireNamespaces(c("S4Vectors", "VariantAnnotation", "BiocInstaller"))

  if(utils::compareVersion(as.character(BiocInstaller::biocVersion()),
                           "3.4") < 0){
    idxRef <- S4Vectors::elementLengths(VariantAnnotation::ref(vcf)) == 1L
    idxAlt <- S4Vectors::elementLengths(VariantAnnotation::alt(vcf)) == 1L

  } else{
    idxRef <- S4Vectors::elementNROWS(VariantAnnotation::ref(vcf)) == 1L
    idxAlt <- S4Vectors::elementNROWS(VariantAnnotation::alt(vcf)) == 1L
  }

  idx <- idxRef & idxAlt
  if(! all(idx))
    warning("only coercing single-element records")

  return(vcf[idx])
}

##' Convert GT to dosage
##'
##' Non-bi-allelic variants are discarded.
##' From Martin Morgan (see http://grokbase.com/t/r/bioconductor/135b460s2b/bioc-how-to-convert-genotype-snp-matrix-to-nucleotide-genotypes).
##' @param vcf CollapsedVCF (see pkg VariantAnnotation)
##' @return matrix with variants in rows and samples in columns
##' @author Timothee Flutre
##' @export
gtVcf2dose <- function(vcf){
  requireNamespace("VariantAnnotation")
  vcf <- getVcfOnlySingleElem(vcf)
  gt <- sub("\\|", "/", VariantAnnotation::geno(vcf)$GT) # ignore phasing
  gt[which(gt == "0/0")] <- 0
  gt[which(gt == "0/1")] <- 1
  gt[which(gt == "1/1")] <- 2
  gt[which(gt == ".")] <- NA
  mode(gt) <- "integer"
  return(gt)
}

##' Convert ranges to data.frame
##'
##' Non-bi-allelic variants are discarded.
##' @param vcf CollapsedVCF (see pkg VariantAnnotation)
##' @param with.coords if TRUE, the output will contain variant coordinates
##' @param with.alleles if TRUE, the output will contain variant alleles
##' @return data.frame with variants in rows
##' @author Timothee Flutre
##' @export
rngVcf2df <- function(vcf, with.coords=TRUE, with.alleles=TRUE){
  requireNamespaces(c("VariantAnnotation", "GenomeInfoDb",
                      "BiocGenerics", "SummarizedExperiment"))
  stopifnot(all(is.logical(with.coords), is.logical(with.alleles)),
            any(with.coords, with.alleles))

  vcf <- getVcfOnlySingleElem(vcf)

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
##' @param vcf.file path to the VCF file (bgzip index should exist in same directory; should only contain SNPs and be already filtered for QD, FS, MQ, etc)
##' @param genome genome identifier (e.g. "VITVI_12x2")
##' @param gdose.file path to the output file to record genotypes as allele doses (will be gzipped)
##' @param ca.file path to the output file to record SNP 1-based coordinates and alleles (will be gzipped)
##' @param yieldSize number of records to yield each time the file is read from (see ?TabixFile) if seq.id is NULL
##' @param dict.file path to the SAM dict file (see \url{https://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary}) if seq.id is specified with no start/end
##' @param seq.id see \code{\link{seqIdStartEnd2GRanges}}
##' @param seq.start see \code{\link{seqIdStartEnd2GRanges}}
##' @param seq.end see \code{\link{seqIdStartEnd2GRanges}}
##' @param uncertain logical indicating whether the genotypes to convert should come from the "GT" field (uncertain=FALSE) or the "GP" or "GL" field (uncertain=TRUE)
##' @param verbose verbosity level (0/1)
##' @return an invisible list with both output file paths
##' @author Timothee Flutre
##' @export
vcf2dosage <- function(vcf.file, genome, gdose.file, ca.file,
                       yieldSize=NA_integer_, dict.file=NULL,
                       seq.id=NULL, seq.start=NULL, seq.end=NULL,
                       uncertain=FALSE, verbose=1){
  requireNamespaces(c("IRanges", "GenomicRanges", "VariantAnnotation",
                      "Rsamtools"))
  stopifnot(file.exists(vcf.file),
            xor(is.na(yieldSize), is.null(seq.id)))
  if(! is.null(seq.id) & is.null(seq.start) & is.null(seq.end))
    stopifnot(! is.null(dict.file),
              file.exists(dict.file))

  for(out.file in c(gdose.file, ca.file))
    if(file.exists(out.file))
      file.remove(out.file)
  gdose.file <- sub("\\.gz$", "", gdose.file)
  ca.file <- sub("\\.gz$", "", ca.file)
  gdose.con <- file(gdose.file, open="a")
  ca.con <- file(ca.file, open="a")
  cat("chr\tpos\tallele.ref\tallele.alt\n", file=ca.con, append=TRUE)

  tabix.file <- Rsamtools::TabixFile(file=vcf.file,
                                     yieldSize=yieldSize)
  if(! is.null(seq.id)){
    rngs <- seqIdStartEnd2GRanges(seq.id=seq.id, seq.start=seq.start,
                                  seq.end=seq.end, dict.file=dict.file)
    vcf.params <- VariantAnnotation::ScanVcfParam(which=rngs)
    vcf <- VariantAnnotation::readVcf(file=tabix.file, genome=genome,
                                      param=vcf.params)
    nb.variants <- nrow(vcf)
    gtmp <- gtVcf2dose(vcf=vcf)
    catmp <- rngVcf2df(vcf=vcf, with.coords=TRUE, with.alleles=TRUE)
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
      gtmp <- gtVcf2dose(vcf=vcf)
      catmp <- rngVcf2df(vcf=vcf, with.coords=TRUE, with.alleles=TRUE)
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
  system(command=paste("gzip", gdose.file))
  system(command=paste("gzip", ca.file))

  invisible(list(gdose.file=paste0(gdose.file, ".gz"),
                 ca.file=paste0(ca.file, ".gz")))
}


##' Convert GT to genotypic classes
##'
##' Non-bi-allelic variants are discarded.
##' From Martin Morgan (see http://grokbase.com/t/r/bioconductor/135b460s2b/bioc-how-to-convert-genotype-snp-matrix-to-nucleotide-genotypes).
##' @param vcf CollapsedVCF (see pkg VariantAnnotation)
##' @param na.string a symbol to indicate missing genotypes (e.g. NA, "NN", "--", etc)
##' @return matrix with variants in rows and samples in columns
##' @author Gautier Sarah
##' @export
gtVcf2genoClasses <- function(vcf, na.string=NA){
  requireNamespaces(c("VariantAnnotation", "GenomicRanges"))

  vcf <- getVcfOnlySingleElem(vcf)
  alt <- VariantAnnotation::alt(vcf)
  ref <- VariantAnnotation::ref(vcf)

  gt <- sub("\\|", "/", VariantAnnotation::geno(vcf)$GT) # ignore phasing

  rowId <- (which(gt == "0/0")-1)  %% nrow(gt) + 1
  gt[which(gt == "0/0")] <- paste0(as.character(ref[rowId]),
                                   as.character(ref[rowId]))

  rowId <- (which(gt == "0/1")-1)  %% nrow(gt) + 1
  gt[which(gt == "0/1")] <- paste0(as.character(ref[rowId]),
                                   GenomicRanges::as.data.frame(alt)[rowId,3])

  rowId <- (which(gt == "1/1")-1)  %% nrow(gt) + 1
  gt[which(gt == "1/1")] <- paste0(GenomicRanges::as.data.frame(alt)[rowId,3],
                                   GenomicRanges::as.data.frame(alt)[rowId,3])

  gt[which(gt == ".")] <- na.string

  return(gt)
}

##' Convert VCF to genotypic classes
##'
##'Convert genotypes at bi-allelic SNPs from a VCF file into genotypic classes.
##' @param vcf.file path to the VCF file (bgzip index should exist in same directory; should only contain SNPs and be already filtered for QD, FS, MQ, etc)
##' @param genome genome identifier (e.g. "VITVI_12x2")
##' @param gclasses.file path to the output file to record genotypes into genotypic classes (will be gzipped)
##' @param ca.file path to the output file to record SNP 1-based coordinates and alleles (will be gzipped)
##' @param yieldSize number of records to yield each time the file is read from (see \code{?TabixFile}) if seq.id is NULL
##' @param dict.file path to the SAM dict file (see \url{https://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary}) if seq.id is specified with no start/end
##' @param seq.id see \code{\link{seqIdStartEnd2GRanges}}
##' @param seq.start see \code{\link{seqIdStartEnd2GRanges}}
##' @param seq.end see \code{\link{seqIdStartEnd2GRanges}}
##' @param na.string a symbol to indicate missing genotypes (e.g. NA, "NN", "--", etc)
##' @param verbose verbosity level (0/1)
##' @return invisible vector with the path to the output file
##' @author Gautier Sarah, Timothee Flutre
##' @export
vcf2genoClasses <- function(vcf.file, genome, gclasses.file, ca.file,
                            yieldSize=NA_integer_, dict.file=NULL,
                            seq.id=NULL, seq.start=NULL, seq.end=NULL,
                            na.string=NA, verbose=1){
  requireNamespaces(c("IRanges", "GenomicRanges", "VariantAnnotation",
                      "Rsamtools"))
  stopifnot(file.exists(vcf.file),
            xor(is.na(yieldSize), is.null(seq.id)))
  if(! is.null(seq.id) & is.null(seq.start) & is.null(seq.end))
    stopifnot(! is.null(dict.file),
              file.exists(dict.file))

  for(out.file in c(gclasses.file, ca.file))
    if(file.exists(out.file))
      file.remove(out.file)
  gclasses.file <- sub("\\.gz$", "", gclasses.file)
  ca.file <- sub("\\.gz$", "", ca.file)
  gclasses.con <- file(gclasses.file, open="a")
  ca.con <- file(ca.file, open="a")
  cat("chr\tpos\tallele.ref\tallele.alt\n", file=ca.con, append=TRUE)

  tabix.file <- Rsamtools::TabixFile(file=vcf.file,
                                     yieldSize=yieldSize)
  if(! is.null(seq.id)){
    rngs <- seqIdStartEnd2GRanges(seq.id=seq.id, seq.start=seq.start,
                                  seq.end=seq.end, dict.file=dict.file)
    vcf.params <- VariantAnnotation::ScanVcfParam(which=rngs)
    vcf <- VariantAnnotation::readVcf(file=tabix.file, genome=genome,
                                      param=vcf.params)
    nb.variants <- nrow(vcf)
    gtmp <- gtVcf2genoClasses(vcf=vcf, na.string=na.string)
    ctmp <- rngVcf2df(vcf=vcf, with.coords=TRUE, with.alleles=TRUE)
    cat(paste(colnames(gtmp), collapse="\t"), file=gclasses.con,
        append=TRUE, sep="\n")
    utils::write.table(x=gtmp,
                       file=gclasses.con, append=TRUE,
                       quote=FALSE, sep="\t", row.names=TRUE,
                       col.names=FALSE)
    utils::write.table(x=ctmp, file=ca.con, append=TRUE,
                       quote=FALSE, sep="\t", row.names=TRUE,
                       col.names=FALSE)
  } else{
    open(tabix.file)
    nb.variants <- 0
    while(nrow(vcf <- VariantAnnotation::readVcf(file=tabix.file,
                                                 genome=genome))){
      gtmp <- gtVcf2genoClasses(vcf, na.string)
      ctmp <- rngVcf2df(vcf=vcf, with.coords=TRUE, with.alleles=TRUE)
      if(nb.variants == 0)
        cat(paste(colnames(gtmp), collapse="\t"), file=gclasses.con,
            append=TRUE, sep="\n")
      nb.variants <- nb.variants + nrow(vcf)
      utils::write.table(x=gtmp, file=gclasses.con, append=TRUE,
                         quote=FALSE, sep="\t", row.names=TRUE,
                         col.names=FALSE)
      utils::write.table(x=ctmp, file=ca.con, append=TRUE,
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
  system(command=paste("gzip", gclasses.file))
  system(command=paste("gzip", ca.file))

  invisible(list(gclasses.file=paste0(gclasses.file, ".gz"),
                 ca.file=paste0(ca.file, ".gz")))
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
