## Contains functions handling sequences often used in bioinformatics.

##' Calculate the GC content of a set of sequences
##'
##' Requires the Biostrings package.
##' @param x vector of sequences (e.g. "AGGT"), possibly with names
##' @return vector
##' @author Timothee Flutre
gc.content <- function(x){
  if(! requireNamespace("Biostrings", quietly=TRUE))
    stop("Pkg Biostrings needed for this function to work. Please install it.",
         call.=FALSE)
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
##' @param type type of alignment (default="global", i.e. Needleman-Wunsch)
##' @param ... arguments to be passed to Biostrings::pairwiseAlignment()
##' @return list of instances of class PairwiseAlignments
##' @author Timothee Flutre
all.pair.aligns <- function(x, type="global", ...){
  if(! requireNamespace("Biostrings", quietly=TRUE))
    stop("Pkg Biostrings needed for this function to work. Please install it.",
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
    stop("Pkg Biostrings needed for this function to work. Please install it.",
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
barplotInsertSizes <- function(file, main=NULL, add.text=FALSE){
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

##' Read bedtools-coverage-hist as data.table
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

##' Plot the covered fraction of regions as a function of depth
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

##' Read SAM dict
##'
##' Reads a file in the "dict" format. See \url{https://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary}.
##' @param file the name of the ".dict" file which the data are to be read from
##' @return data.frame
##' @author Timothee Flutre
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

##' Information on variant-level calls
##'
##' Return some information to help in hard-filtering variant-level calls. See GATK's Best Practices tutorial (\url{https://www.broadinstitute.org/gatk/guide/topic?name=tutorials#tutorials2806}).
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
confidenceGenoOneVar <- function(x, plot.it=FALSE){
  if(! requireNamespace("VariantAnnotation", quietly=TRUE))
    stop("Pkg VariantAnnotation needed for this function to work. Please install it.",
         call.=FALSE)
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
    plot(x=genos[["DP"]][1,], y=genos[["GQ"]][1,],
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
##' @param verbose verbosity level (0/default=1)
##' @return the destination file path as a character(1)
##' @author Timothee Flutre
setGt2Na <- function(vcf.file, genome, out.file,
                     yieldSize=NA_integer_, dict.file=NULL,
                     seq.id=NULL, seq.start=NULL, seq.end=NULL,
                     min.gq=90, verbose=1){
  if(! requireNamespace("IRanges", quietly=TRUE))
    stop("Pkg IRanges needed for this function to work. Please install it.",
         call.=FALSE)
  if(! requireNamespace("GenomicRanges", quietly=TRUE))
    stop("Pkg GenomicRanges needed for this function to work. Please install it.",
         call.=FALSE)
  if(! requireNamespace("VariantAnnotation", quietly=TRUE))
    stop("Pkg VariantAnnotation needed for this function to work. Please install it.",
         call.=FALSE)
  if(! requireNamespace("Rsamtools", quietly=TRUE))
    stop("Pkg Rsamtools needed for this function to work. Please install it.",
         call.=FALSE)
  if(! requireNamespace("S4Vectors", quietly=TRUE))
    stop("Pkg S4Vectors needed for this function to work. Please install it.",
         call.=FALSE)
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
    nb.records <- nrow(vcf)

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
##' @param min.prop.spl.dp minimum proportion of samples with DP above threshold
##' @param min.spl.gq minimum sample-level GQ
##' @param min.prop.spl.gq minimum proportion of samples with GQ above threshold
##' @param max.var.nb.gt.na maximum number of samples with missing GT
##' @param max.var.prop.gt.na maximum proportion of samples with missing GT
##' @param verbose verbosity level (0/default=1)
##' @return the destination file path as a character(1)
##' @author Timothee Flutre
filterVariantCalls <- function(vcf.file, genome, out.file,
                               yieldSize=NA_integer_, dict.file=NULL,
                               seq.id=NULL, seq.start=NULL, seq.end=NULL,
                               is.snv=NULL, is.biall=NULL,
                               min.var.dp=NULL, max.var.dp=NULL,
                               min.alt.af=NULL, max.alt.af=NULL,
                               min.spl.dp=NULL, min.prop.spl.dp=NULL,
                               min.spl.gq=NULL, min.prop.spl.gq=NULL,
                               max.var.nb.gt.na=NULL, max.var.prop.gt.na=NULL,
                               verbose=1){
  if(! requireNamespace("IRanges", quietly=TRUE))
    stop("Pkg IRanges needed for this function to work. Please install it.",
         call.=FALSE)
  if(! requireNamespace("GenomicRanges", quietly=TRUE))
    stop("Pkg GenomicRanges needed for this function to work. Please install it.",
         call.=FALSE)
  if(! requireNamespace("VariantAnnotation", quietly=TRUE))
    stop("Pkg VariantAnnotation needed for this function to work. Please install it.",
         call.=FALSE)
  if(! requireNamespace("Rsamtools", quietly=TRUE))
    stop("Pkg Rsamtools needed for this function to work. Please install it.",
         call.=FALSE)
  if(! requireNamespace("S4Vectors", quietly=TRUE))
    stop("Pkg S4Vectors needed for this function to work. Please install it.",
         call.=FALSE)
  stopifnot(file.exists(vcf.file))
  if(! is.null(seq.id) & is.null(seq.start) & is.null(seq.end))
    stopifnot(! is.null(dict),
              file.exists(dict.file))
  if(! is.null(min.prop.spl.dp))
    stopifnot(min.prop.spl.dp >= 0, min.prop.spl.dp <= 1)
  if(! is.null(min.prop.spl.gq))
    stopifnot(min.prop.spl.gq >= 0, min.prop.spl.gq <= 1)
  if(! is.null(max.var.nb.gt.na))
    stopifnot(max.var.nb.gt.na >= 0)
  if(! is.null(max.var.prop.gt.na))
    stopifnot(max.var.prop.gt.na >= 0, max.var.prop.gt.na <= 1)

  dest <- NULL

  out.file <- sub("\\.gz$", "", out.file)
  out.file <- sub("\\.bgz$", "", out.file)

  ##' @return TRUE if SNV (single nucleotide variant)
  filterSnv <- function(x){
    (VariantAnnotation::isSNV(x))
  }

  ##' @return TRUE if at most one alternate allele
  filterBiall <- function(x){
    (S4Vectors::elementLengths(VariantAnnotation::alt(x)) <= 1)
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

  ##' @return TRUE if high-enough proportion of samples with DP above threshold
  filterSampleDp <- function(x, min.dp=min.spl.dp,
                             min.prop.dp=min.prop.spl.dp){
    if(nrow(x) == 0){
      logical(0)
    } else{
      dp <- VariantAnnotation::geno(x)$DP
      (rowSums(dp >= min.dp) / ncol(dp) >= min.prop.dp)
    }
  }

  ##' @return TRUE if high-enough proportion of samples with GQ above threshold
  filterSampleGq <- function(x, min.gq=min.spl.gq,
                             min.prop.gq=min.prop.spl.gq){
    if(nrow(x) == 0){
      logical(0)
    } else{
      gq <- VariantAnnotation::geno(x)$GQ
      (rowSums(gq >= min.gq) / ncol(gq) >= min.prop.gq)
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

  ##' @return TRUE if high-enough proportion of samples without missing genotypes
  filterGtProp <- function(x, max.prop.gt.na=max.var.prop.gt.na){
    if(nrow(x) == 0){
      logical(0)
    } else{
      gt <- VariantAnnotation::geno(x)$GT
      (rowSums(gt == ".") / ncol(gt) <= max.prop.gt.na)
    }
  }

  ## set the filters
  tmp <- list()
  if(! is.null(is.snv)){
    stopifnot(is.logical(is.snv))
    if(is.snv)
      tmp[["filterSnv"]] <- filterSnv
  }
  if(! is.null(is.biall)){
    stopifnot(is.logical(is.biall))
    if(is.biall)
      tmp[["filterBiall"]] <- filterBiall
  }
  if(! is.null(min.var.dp) & ! is.null(max.var.dp)){
    tmp[["filterVariantDp"]] <- filterVariantDp
  }
  if(! is.null(min.alt.af) & ! is.null(max.alt.af)){
    tmp[["filterAf"]] <- filterAf
  }
  if(! is.null(min.spl.dp) & ! is.null(min.prop.spl.dp)){
    tmp[["filterSampleDp"]] <- filterSampleDp
  }
  if(! is.null(min.spl.gq) & ! is.null(min.prop.spl.gq)){
    tmp[["filterSampleGq"]] <- filterSampleGq
  }
  if(! is.null(max.var.nb.gt.na)){
    tmp[["filterGtNb"]] <- filterGtNb
  }
  if(! is.null(max.var.prop.gt.na)){
    tmp[["filterGtProp"]] <- filterGtProp
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
        rngs <- GenomicRanges::GRanges(seqnames=c(seq.id),
                                       ranges=IRanges::IRanges(start=c(1),
                                           end=c(dict$LN[rownames(dict) == seq.id])))
      } else
        rngs <- GenomicRanges::GRanges(seqnames=c(seq.id),
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

##' Summary of GQ per variant
##'
##' Compute the min, Q1, med, mean, Q3, max of the genotype qualities per variant.
##' @param vcf.file path to the VCF file
##' @param genome genome identifier (e.g. "VITVI_12x2")
##' @param yieldSize number of records to yield each time the file is read from (see ?TabixFile)
##' @param verbose verbosity level (0/default=1)
##' @return matrix with one row per variant and 7 columns (min, q1, med, mean, q3, max, na)
##' @author Timothee Flutre
summaryGq <- function(vcf.file, genome, yieldSize=10^4, verbose=1){
  stopifnot(file.exists(vcf.file))

  all.smry.gq <- list(min=c(),
                      q1=c(),
                      med=c(),
                      mean=c(),
                      q3=c(),
                      max=c(),
                      na=c())
  var.names <- c()
  tabix.file <- Rsamtools::TabixFile(file=vcf.file,
                                     yieldSize=yieldSize)
  open(tabix.file)
  while(nrow(vcf <- VariantAnnotation::readVcf(file=tabix.file,
                                               genome=genome))){
    gq <- VariantAnnotation::geno(vcf)[["GQ"]]
    all.smry.gq$min <- append(all.smry.gq$min,
                              suppressWarnings(apply(gq, 1, min, na.rm=TRUE)))
    all.smry.gq$q1 <- append(all.smry.gq$q1,
                             suppressWarnings(apply(gq, 1, quantile,
                                                    probs=0.25, na.rm=TRUE)))
    all.smry.gq$med <- append(all.smry.gq$med,
                              suppressWarnings(apply(gq, 1, median, na.rm=TRUE)))
    all.smry.gq$mean <- append(all.smry.gq$mean,
                               suppressWarnings(apply(gq, 1, mean, na.rm=TRUE)))
    all.smry.gq$q3 <- append(all.smry.gq$q3,
                             suppressWarnings(apply(gq, 1, quantile,
                                                    probs=0.75, na.rm=TRUE)))
    all.smry.gq$max <- append(all.smry.gq$max,
                              suppressWarnings(apply(gq, 1, max, na.rm=TRUE)))
    all.smry.gq$na <- append(all.smry.gq$na, rowSums(t(apply(gq, 1, is.na))))
    var.names <- append(var.names, rownames(gq))
  }
  close(tabix.file)

  if(verbose > 0){
    msg <- paste0("nb of variants: ", length(var.names))
    write(msg, stdout())
  }

  all.smry.gq <- matrix(data=do.call(cbind, all.smry.gq),
                        nrow=length(var.names),
                        ncol=length(all.smry.gq),
                        dimnames=list(var.names,
                            names(all.smry.gq)))
  return(all.smry.gq)
}

##' Convert VCF to dose
##'
##' Convert genotypes from a VCF file into allele doses.
##' @param vcf.file path to the VCF file (bgzip index should exist in same directory; should only contain SNPs and be already filtered for QD, FS, MQ, etc)
##' @param dict.file path to the SAM dict file (see \url{https://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary})
##' @param genome genome identifier (e.g. "VITVI_12x2")
##' @param seqs sequence identifier(s) to work on (e.g. "chr2")
##' @param gdose.file path to the output file to record genotypes as allele doses (will be gzipped)
##' @param amap.file path to the output file to record SNP positions and alleles (will be gzipped)
##' @param uncertain boolean indicating whether the genotypes to convert should come from the "GT" field (uncertain=FALSE) or the "GP" or "GL" field (uncertain=TRUE)
##' @param verbose verbosity level (0/default=1)
##' @return nothing
##' @author Timothee Flutre
vcf2dosage <- function(vcf.file, dict.file, genome, seqs=NULL, gdose.file,
                       amap.file, uncertain=FALSE, verbose=1){
  if(! requireNamespace("IRanges", quietly=TRUE))
    stop("Pkg IRanges needed for this function to work. Please install it.",
         call.=FALSE)
  if(! requireNamespace("GenomicRanges", quietly=TRUE))
    stop("Pkg GenomicRanges needed for this function to work. Please install it.",
         call.=FALSE)
  if(! requireNamespace("VariantAnnotation", quietly=TRUE))
    stop("Pkg VariantAnnotation needed for this function to work. Please install it.",
         call.=FALSE)
  if(! requireNamespace("Rsamtools", quietly=TRUE))
    stop("Pkg Rsamtools needed for this function to work. Please install it.",
         call.=FALSE)
  if(! requireNamespace("snpStats", quietly=TRUE))
    stop("Pkg snpStats needed for this function to work. Please install it.",
         call.=FALSE)
  stopifnot(file.exists(vcf.file),
            file.exists(dict.file),
            ! file.exists(gdose.file),
            ! file.exists(amap.file))

  dict <- readSamDict(file=dict.file)

  if(! file.exists(paste(vcf.file, "tbi", sep=".")))
    Rsamtools::indexTabix(file=vcf.file, format="vcf")
  tabix.file <- Rsamtools::TabixFile(vcf.file)
  if(is.null(seqs)){
    seqs <- Rsamtools::seqnamesTabix(tabix.file)
  } else{
    for(seq in seqs){
      if(! seq %in% Rsamtools::seqnamesTabix(tabix.file)){
        msg <- paste0("seq '", seq, "' not in '", vcf.file, "'")
        stop(msg)
      }
    }
  }

  gdose.file <- gsub(".gz", "", gdose.file)
  amap.file <- gsub(".gz", "", amap.file)
  gdose.con <- file(gdose.file, open="a")
  amap.con <- file(amap.file, open="a")

  for(seq in seqs){
    if(! seq %in% rownames(dict)){
      msg <- paste0("seq '", seq, "' not in '", dict.file, "'")
      stop(msg)
    }

    rngs <- GenomicRanges::GRanges(seqnames=c(seq),
                                   ranges=IRanges::IRanges(start=c(1),
                                       end=c(dict$LN[rownames(dict) == seq])))
    names(rngs) <- c(seq)
    vcf.params <- VariantAnnotation::ScanVcfParam(which=rngs)
    if(verbose > 0){
      msg <- paste0("read seq '", seq, "'...")
      write(msg, stdout())
    }
    vcf <- VariantAnnotation::readVcf(file=tabix.file, genome=genome,
                                      param=vcf.params)
    if(verbose > 0){
      msg <- paste0("seq '", seq, "': ", nrow(vcf), " variants x ",
                    ncol(vcf), " samples")
      write(msg, stdout())
    }

    res <- VariantAnnotation::genotypeToSnpMatrix(vcf)

    snpStats::write.SnpMatrix(x=res$genotypes, file=gdose.con, append=TRUE,
                              quote=FALSE, sep="\t", row.names=TRUE,
                              col.names=TRUE)
    write.table(x=res$map, file=amap.con, append=TRUE,
                quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
  }

  close(gdose.con)
  close(amap.con)
  system(command=paste("gzip", gdose.file))
  system(command=paste("gzip", amap.file))
}
