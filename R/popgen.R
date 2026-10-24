## Contains functions useful for population genetics.

##' SNP heterozygosity
##'
##' Compute the observed and expected heterozygosity (H_o and H_e), averaged over all SNPs.
##' @param X matrix of bi-allelic SNP genotypes encoded in allele doses in \{0,1,2\}, with genotypes in rows and SNPs in columns; missing values should be encoded as NA; SNPs with NAs will be discarded
##' @param verbose verbosity level (0/1)
##' @return named vector of two numeric values
##' @seealso \code{\link{estimSnpMaf}}
##' @author Timothee Flutre
##' @examples
##' \dontrun{set.seed(1859)
##' library(scrm)
##' nb.genos <- 100
##' Ne <- 10^4
##' chrom.len <- 10^5
##' mu <- 10^(-8)
##' c <- 10^(-8)
##'
##' ## under HWE:
##' genomes <- simulCoalescent(nb.inds=nb.genos,
##'                            pop.mut.rate=4 * Ne * mu * chrom.len,
##'                            pop.recomb.rate=4 * Ne * c * chrom.len,
##'                            chrom.len=chrom.len)
##' avgSnpHeterozygosity(genomes$genos)
##' ## Ho ~= He as expected under HWE
##'
##' ## with medium pop struct:
##' mig_rates <- c("high"=10^2, "med"=10, "low"=0.5)
##' nb_pops <- 10
##' level_mig_rate <- "low"
##' genomes <- simulCoalescent(nb.inds=nb.genos,
##'                            pop.mut.rate=4 * Ne * mu * chrom.len,
##'                            pop.recomb.rate=4 * Ne * c * chrom.len,
##'                            chrom.len=chrom.len,
##'                            nb.pops=nb_pops,
##'                            mig.rate=mig_rates[level_mig_rate])
##' avgSnpHeterozygosity(genomes$genos)
##' ## Ho < He: deficit of heterozygotes (Wahlund effect)
##' }
##' @export
avgSnpHeterozygosity <- function(X, verbose = 1) {
  out <- setNames(rep(NA, 2), c("exp", "obs"))

  if (any(is.na(X))) {
    X <- discardMarkersMissGenos(X, verbose)
  }
  nbGenos <- nrow(X)

  afs <- estimSnpAf(X)
  He_perSNP <- 2 * afs * (1 - afs)
  out["exp"] <- mean(He_perSNP)

  Ho_perSNP <- apply(X, 2, function(x) {
    sum(x == 1) / nbGenos
  })
  out["obs"] <- mean(Ho_perSNP)

  return(out)
}

##' Genetic diversity
##'
##' Compute Watterson's theta (Watterson, 1975; \doi{10.1016/0040-5809(75)90020-9}).
##' @param X matrix of bi-allelic SNP genotypes encoded in allele doses in \{0,1,2\}, with genotypes in rows and SNPs in columns; missing values should be encoded as NA; these NAs will be ignored
##' @param verbose verbosity level (0/1)
##' @return numeric
##' @author Timothee Flutre
##' @export
thetaWatterson <- function(X, verbose = 1) {
  out <- NA

  n <- nrow(X)
  stopifnot(n >= 2)

  ## nb of segregating sites
  S <- sum(apply(X, 2, var, na.rm = TRUE) > 0, na.rm = TRUE)

  ## (n-1)-th harmonic number
  a_n <- sum(1 / (1:(n - 1)))

  out <- S / a_n

  attr(out, "n") <- n
  attr(out, "S") <- S
  attr(out, "a_n") <- a_n

  return(out)
}

##' Clusteredness
##'
##' Computes the statistic which "measures the extent to which a randomly chosen individual is inferred to have ancestry in only one cluster (clusteredness = 1), with the other extreme being equal membership in all clusters (clusteredness = 0)", from \href{http://dx.plos.org/10.1371/journal.pgen.0010070}{Rosenberg et al (2005)}.
##' @param q.ik matrix of assignment probabilities with individuals in rows and clusters in columns
##' @return numeric
##' @author Timothee Flutre
##' @export
clusteredness <- function(q.ik) {
  I <- nrow(q.ik)
  K <- ncol(q.ik)
  (1 / I) * sum(sqrt((K / (K - 1)) * rowSums((q.ik - (1 / K))^2)))
}

##' Genetic clustering
##'
##' Describes the strong assignment of individuals per cluster.
##' @param assign.probs vector of assignment probabilities
##' @param assign.grps vector of assigned clusters
##' @param thresh.strong.assign threshold
##' @return list
##' @author Timothee Flutre
##' @export
statsStrongAssign <- function(assign.probs, assign.grps, thresh.strong.assign = 0.8) {
  stopifnot(
    nrow(assign.probs) == length(assign.grps),
    all(unique(assign.grps) %in% colnames(assign.probs))
  )
  out <- list()
  assign.grps <- as.character(assign.grps)
  id.grps <- sort(unique(assign.grps))
  nb.grps <- length(id.grps)
  out[["nb.samples.per.grp"]] <- stats::setNames(rep(NA, nb.grps), id.grps)
  out[["nb.samples.strong.assign.per.grp"]] <- stats::setNames(rep(NA, nb.grps), id.grps)
  out[["perc.samples.strong.assign.per.grp"]] <- stats::setNames(rep(NA, nb.grps), id.grps)
  for (id in id.grps) {
    idx <- which(assign.grps == id)
    out[["nb.samples.per.grp"]][id] <- length(idx)
    out[["nb.samples.strong.assign.per.grp"]][id] <-
      sum(assign.probs[idx, id] >= thresh.strong.assign)
  }
  out[["perc.samples.strong.assign.per.grp"]] <- 100 * out[["nb.samples.strong.assign.per.grp"]] /
    out[["nb.samples.per.grp"]]
  return(out)
}

##' STRUCTURE
##'
##' Read the output file from STRUCTURE.
##' @param out.file path to the output file
##' @return list
##' @author Timothee Flutre
##' @export
readOutSTRUCTURE <- function(out.file) {
  stopifnot(file.exists(out.file))

  out <- list()

  tmp <- readLines(out.file, warn = FALSE)

  idx.start <- grep("^Inferred ancestry of individuals:$", tmp)
  idx.end <- grep("^Estimated Allele Frequencies in each cluster$", tmp)
  ancestries <- tmp[(idx.start + 2):(idx.end - 3)]
  ancestries <- as.data.frame(
    t(sapply(
      strsplit(ancestries, " "),
      function(x) {
        x[!x == ""]
      }
    )),
    stringsAsFactors = FALSE
  )
  stopifnot(all(ancestries[, 4] == ":"))
  nb.clusters <- ncol(ancestries) - 4
  colnames(ancestries) <- c(
    "index", "label", "perc.miss", "sep",
    paste0("cluster", 1:nb.clusters)
  )
  rownames(ancestries) <- ancestries$label
  ancestries$index <- as.integer(ancestries$index)
  for (i in 1:nb.clusters) {
    col.idx <- ncol(ancestries) - i + 1
    ancestries[, col.idx] <- as.numeric(ancestries[, col.idx])
  }
  out$ancestries <- ancestries

  return(out)
}

##' STRUCTURE to CLUMPP
##'
##' Gather output files from STRUCTURE (several runs at the same number of clusters, K) into a single file for input to CLUMPP.
##' @param struct.files vector of file names
##' @param clumpp.file file name
##' @param verbose verbosity level (0/1/2)
##' @return invisible list of lists, one per STRUCTURE file
##' @author Timothee Flutre
##' @seealso \code{\link{readOutSTRUCTURE}}
##' @export
structure2clumpp <- function(struct.files, clumpp.file, verbose = 1) {
  stopifnot(
    all(file.exists(struct.files)),
    !clumpp.file %in% struct.files
  )
  if (file.exists(clumpp.file)) {
    file.remove(clumpp.file)
  }

  out <- list()

  if (verbose > 0) {
    msg <- paste0(
      "gather ", length(struct.files), " output files",
      " from STRUCTURE",
      " to use with CLUMPP"
    )
    write(msg, stdout())
  }

  for (i in seq_along(struct.files)) {
    tmp <- readOutSTRUCTURE(struct.files[i])
    out[[struct.files[i]]] <- tmp
    tmp <- tmp$ancestries
    tmp <- cbind(
      tmp$index, # col 1) ignored by CLUMPP
      tmp$index, # col 2) integer identifying the individuals
      tmp$perc.miss, # col 3) ignored by CLUMPP
      rep(1, nrow(tmp)), # col 4) ignored by CLUMPP
      tmp$sep, # col 5) ignored by CLUMPP
      tmp[, 5:ncol(tmp)]
    ) # cols 6...) membership coefficients
    utils::write.table(
      x = tmp, file = clumpp.file, append = TRUE, quote = FALSE,
      sep = " ", row.names = FALSE, col.names = FALSE
    )
  }

  invisible(out)
}
