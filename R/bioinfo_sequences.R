## Contains functions handling sequences often used in bioinformatics.

##' Calculate the GC content of a set of sequences.
##'
##' Requires the Biostrings package.
##' @param x vector of sequences (e.g. "AGGT"), possibly with names
##' @return vector
##' @author Timothee Flutre
gc.content <- function(x){
  if(! requireNamespace("Biostrings", quietly=TRUE))
    stop("Pkg needed for this function to work. Please install it.",
         call.=FALSE)
  stopifnot(is.character(x))

  sapply(x, function(xi){
    sum(Biostrings::alphabetFrequency(Biostrings::DNAString(xi),
                                      baseOnly=TRUE, as.prob=TRUE)[c("C","G")])
  })
}

##' Align each sequence against each other (and itself).
##'
##' Requires the Biostrings package.
##' @title All pairwise alignments
##' @param x vector of sequences (e.g. "AGGT"), possibly with names
##' @param type type of alignment (default="global", i.e. Needleman-Wunsch)
##' @param ... arguments to be passed to Biostrings::pairwiseAlignment()
##' @return list of instances of class PairwiseAlignments
##' @author Timothee Flutre
all.pair.aligns <- function(x, type="global", ...){
  if(! requireNamespace("Biostrings", quietly=TRUE))
    stop("Pkg needed for this function to work. Please install it.",
         call.=FALSE)
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
##' @param aligns list of instances of class PairwiseAlignments (see all.pair.aligns())
##' @param nb.sequences number of sequences
##' @return list of matrices
##' @author Timothee Flutre
stats.all.pair.aligns <- function(aligns, nb.sequences){
  if(! requireNamespace("Biostrings", quietly=TRUE))
    stop("Pkg needed for this function to work. Please install it.",
         call.=FALSE)
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

##' Bar plot of insert sizes
##'
##' Creates a bar plot with vertical bars of insert sizes from histogram data as calculated by Picard CollectInsertSizeMetrics.
##' @param file path to the output file from Picard CollectInsertSizeMetrics
##' @param main overall title for the plot
##' @param add.text add total count, as well as Q25, median and Q75 for insert sizes, to the topright of the plot
##' @return invisible data frame of the content of the file
##' @author Timothee Flutre
barplot.insert.sizes <- function(file, main=NULL, add.text=FALSE){
  stopifnot(file.exists(file))

  dat <- read.table(file=file, skip=10, header=TRUE)
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

  bp <- barplot(height=dat$count, width=1,
                xlab="Insert size (in bp)",
                ylab="Counts",
                main=main)
  axis(side=1, at=c(0, bp[seq(100, max(dat$insert.size), 100)]),
       labels=c(0, seq(100, max(dat$insert.size), 100)))

  abline(v=bp[q25.insert.size], lty=2)
  abline(v=bp[med.insert.size], lty=2)
  abline(v=bp[q75.insert.size], lty=2)

  if(add.text){
    text(x=bp[floor(0.6*max(dat$insert.size))],
         y=0.6*max(dat$count), adj=0,
         labels=paste0("Q25 = ", format(q25.insert.size, digits=2), " bp"))
    text(x=bp[floor(0.6*max(dat$insert.size))],
         y=0.7*max(dat$count), adj=0,
         labels=paste0("median = ", format(med.insert.size, digits=2), " bp"))
    text(x=bp[floor(0.6*max(dat$insert.size))],
         y=0.8*max(dat$count), adj=0,
         labels=paste0("Q75 = ", format(q75.insert.size, digits=2), " bp"))
    text(x=bp[floor(0.6*max(dat$insert.size))],
         y=0.9*max(dat$count), adj=0,
         labels=paste0("total count = ", format(tot.count, digits=2)))
  }

  invisible(dat)
}

##' Read bedtools-coverage-hist as data.table.
##'
##' Read output files from bedtools coverage with -hist into a list of data.tables. All lines starting by "all" must have been discarded beforehand.
##' @param files character vector of relative or absolute filepaths (e.g. from Sys.glob).
##' @param verbose verbosity level (0/default=1)
##' @return list of data.tables
##' @author Timothee Flutre
fread.bedtools.coverage.hist <- function(files, verbose=1){
  colClasses <- sapply(read.table(files[1], nrows=5), class)
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
##' @param dat data.table (see fread.bedtools.coverage.hist)
##' @param min.reg.len minimum length of a region to be considered
##' @param max.reg.len maximum length of a region to be considered
##' @param min.reg.dep minimum depth of a region when reporting the number of interesting regions
##' @param max.reg.dep maximum depth of a region when reporting the number of interesting regions
##' @param min.reg.frac minimum fraction of a region when reporting the number of interesting regions
##' @return data.table
##' @author Timothee Flutre
depths.per.sample <- function(dat, min.reg.len=30, max.reg.len=500,
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
                         depth.med=as.double(median(.SD[,depth])),
                         depth.mean=mean(.SD[,depth]),
                         depth.max=max(.SD[,depth]),
                         depth.q65=quantile(.SD[,depth], 0.65),
                         depth.q70=quantile(.SD[,depth], 0.70),
                         depth.q75=quantile(.SD[,depth], 0.75),
                         depth.q80=quantile(.SD[,depth], 0.80),
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
##' @param dat data.table (see fread.bedtools.coverage.hist)
##' @param min.reg.len minimum length of a region to be considered
##' @param max.reg.len maximum length of a region to be considered
##' @param min.reg.dep minimum depth of a region when reporting the number of interesting regions
##' @param max.reg.dep maximum depth of a region when reporting the number of interesting regions
##' @param min.reg.frac minimum fraction of a region when reporting the number of interesting regions
##' @return data.table
##' @author Timothee Flutre
depths.per.region <- function(dat, min.reg.len=30, max.reg.len=500,
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
                         depth.med=as.double(median(.SD[,depth])),
                         depth.max=max(.SD[,depth]),
                         samples.ok=nrow(unique(.SD[depth >= min.reg.dep &
                                                      depth <= max.reg.dep &
                                                        fraction >= min.reg.frac,
                             .(ind,flowcell,lane)]))),
                       by=.(chrom,start,end,name)]

  setkey(dat, NULL)

  return(depths.region)
}

##' Plot the covered fraction of regions as a function of depth.
##'
##' Need to first run bedtools coverage as in http://www.gettinggeneticsdone.com/2014/03/visualize-coverage-exome-targeted-ngs-bedtools.html by Stephen Turner
##' @param path string
##' @param pattern character string containing a regular expression passed to gsub to be removed to simplify the file basenames used in the legend (e.g. "bedtools-coverage|\\.txt\\.gz")
##' @param covrg list of data frames
##' @param plot.it boolean (default=TRUE)
##' @param xlim vector of limits for the x-axis
##' @param ylim vector of limits for the y-axis
##' @param points.type character indicating the type of plotting ("p"/default="l"/"b")
##' @param plot.legend boolean (default=TRUE)
##' @param verbose verbosity level (0/default=1)
##' @return invisible list of covrg and cumcovrg
##' @author Timothee Flutre
coverage.regions <- function(path=NULL, pattern=NULL, covrg=NULL, plot.it=TRUE,
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
      covrg[[i]] <- read.table(files[i], sep="\t",
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
    plot(x=covrg[[1]][0:500,2], y=cumcovrg[[1]][0:500], type="n",
         xlab="Depth", ylab="Fraction of bases >= depth",
         main="Region coverage",
         xlim=xlim, ylim=ylim, las=1)
    for(i in seq_along(covrg))
      points(covrg[[i]][0:500, 2], cumcovrg[[i]][0:500], type=points.type, lwd=1,
             col=cols[i])
    if(plot.legend)
      legend("topright", legend=labs, col=cols, lty=1, lwd=1, bty="n")
    abline(h=0.10, lty=2, col="grey60")
    abline(h=0.50, lty=2, col="grey60")
    abline(h=0.90, lty=2, col="grey60")
  }

  invisible(list(covrg=covrg, cumcovrg=cumcovrg))
}
