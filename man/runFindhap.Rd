% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findhap.R
\name{runFindhap}
\alias{runFindhap}
\title{Genotype imputation via findhap}
\usage{
runFindhap(
  X = NULL,
  chips = NULL,
  snp.coords = NULL,
  vcf.file = NULL,
  yieldSize = 10000,
  ped = NULL,
  work.dir = getwd(),
  task.id = "findhap",
  params.findhap = NULL,
  nb.threads = 1,
  clean = "none",
  verbose = 1
)
}
\arguments{
\item{X}{matrix of bi-allelic SNP genotypes encoded in number of copies of the 2nd allele, i.e. as allele doses in \{0,1,2\}, with genotypes in rows and SNPs in columns; missing values should be encoded as NA; if not NULL, will be used in priority even if \code{vcf.file} is not NULL}

\item{chips}{if several microarrays were used, provide a vector with chip numbers which names are genotype identifiers (same as in \code{X})}

\item{snp.coords}{data.frame with SNP identifiers as row names, and two compulsory columns in that order, chromosome identifiers and coordinate/position identifiers; the maximum length of SNP identifiers is 50 characters; chromosome identifiers should be numeric; other optional column(s) contain the SNP order on the microarray(s) (maximum 10); compulsory if \code{X} is specified}

\item{vcf.file}{path to the VCF file (if the bgzip index doesn't exist in the same directory, it will be created); used only if \code{X=NULL}}

\item{yieldSize}{number of records to yield each time the file is read from  (see \code{?TabixFile})}

\item{ped}{data frame of pedigree with four columns in that order, genotype identifiers (same as in \code{X}), parent1 identifiers (considered as father/sire), parent2 identifiers (considered as mother/dam), and sex (as M or F)}

\item{work.dir}{directory in which the input and output files will be saved}

\item{task.id}{identifier of the task (used in temporary and output file names)}

\item{params.findhap}{list of additional parameters to pass to findhap}

\item{nb.threads}{number of threads to pass to findhap}

\item{clean}{remove files: none, some (temporary only), all (temporary and results)}

\item{verbose}{verbosity level (0/1)}
}
\value{
list
}
\description{
Impute SNP genotypes via findhap (\href{http://dx.doi.org/10.3168/jds.2012-5702}{VanRaden et al, 2013}).
}
\examples{
\dontrun{## simulate haplotypes and genotypes in a single population
set.seed(1859)
nb.inds <- 5*10^2
nb.chrs <- 1
Ne <- 10^4
chrom.len <- 1*10^6
mu <- 10^(-8)
c.rec <- 10^(-7)
genomes <- simulCoalescent(nb.inds=nb.inds,
                           nb.reps=nb.chrs,
                           pop.mut.rate=4 * Ne * mu * chrom.len,
                           pop.recomb.rate=4 * Ne * c.rec * chrom.len,
                           chrom.len=chrom.len,
                           get.alleles=TRUE)
nb.snps <- nrow(genomes$snp.coords)
plotHaplosMatrix(genomes$haplos[[1]]) # quick view of the amount of LD

## discard some genotypes according to a "microarray" design:
## some inds with high density of genotyped SNPs, and the others with
## low density of SNPs, these being on both microarrays
ind.names <- rownames(genomes$genos)
inds.high <- sample(x=ind.names, size=floor(0.4 * nb.inds))
inds.low <- setdiff(ind.names, inds.high)
snp.names <- colnames(genomes$genos)
mafs <- estimSnpMaf(X=genomes$genos)
length(idx <- which(mafs >= 0.05))
snps.high <- names(idx)
snps.high.only <- sample(x=snps.high, size=floor(0.7*length(snps.high)))
dim(X <- genomes$genos[ind.names, snps.high])
X.na <- X
X.na[inds.low, snps.high.only] <- NA
sum(is.na(X.na)) / length(X.na)
plotGridMissGenos(X=X.na)

## perform imputation
snp.coords <- genomes$snp.coords[colnames(X.na),]
snp.coords$chr <- as.numeric(sub("chr", "", snp.coords$chr))
out.fh <- runFindhap(X=X.na, snp.coords=snp.coords, clean="all")
X.impfh <- outfh$genos.imp

## assess imputation accuracy
genomes$haplos[[1]][1:4, 1:7]
genomes$genos[1:2, 1:7]
X.na[1:2, 1:6]
X.impfh[1:2, 1:6]
100 * sum(X.impfh != X) / sum(is.na(X.na))
}
}
\seealso{
\code{\link{writeInputsFindhap}}, \code{\link{readOutputsFindhap}}, \code{\link{runFastphase}}
}
\author{
Timothee Flutre
}
