## Contains functions useful to perform QTL detection by interval mapping
## with R/qtl.

##' QTL detection by SIM
##'
##' Perform QTL detection by simple interval mapping
##' @param cross object
##' @param response.in.cross logical
##' @param response character
##' @param numeric.chr.format logical
##' @param geno.joinmap genotypes at markers in the JoinMap format
##' @param phase phase
##' @param pheno.col character
##' @param nperm number of permutations to be done in \code{qtl::scanone}
##' @param threshold.alpha vector of length 1 or 2 (optional) with thresholds for (1) the significance of QTL presence (based on permutations) and (2) the significance of QTL effects
##' @param QTL_position useful for simulated data
##' @param method method to detect QTL in \code{qtl::scanone}
##' @param plot logical
##' @param verbose verbosity level (0/1/2)
##' @return list
##' @author Charlotte Brault [aut], Timothee Flutre [ctb]
##' @export
SIMQTL <- function (cross, response.in.cross=TRUE, response=NULL, numeric.chr.format=TRUE,
                    geno.joinmap=NULL, phase=NULL,
                    pheno.col="y", nperm=100, threshold.alpha=0.05, QTL_position=NULL,
                    method="em", plot=FALSE, verbose=0){

  requireNamespace(c("qtl", "rutilstimflutre", "caret"))

  ## Reformat chromosome names
  if(! numeric.chr.format){
    nb_chr <- length(lapply(cross$geno, attributes))
    chr_renamed <- seq(1:nb_chr)
    names(cross$geno) <- chr_renamed
  }

  ## Verifications
  stopifnot(all(threshold.alpha <= 0.2),
            threshold.alpha[1] %in% c(0.05,0.1),
            xor(all(!response.in.cross, !is.null(response)),
                all(response.in.cross & is.null(response))))
  if(! is.null(QTL_position))
    stopifnot(c("linkage.group", "genetic.distance") %in% colnames(QTL_position))
  ## Check the conformity between geno.joinmap and cross
  if(!is.null(geno.joinmap)){
    jm <- data.frame(seg=rutilstimflutre::getJoinMapSegregs(geno.joinmap), phase=phase)
    jm <- cbind(jm, geno.joinmap)
    jm[jm == "--"] <- NA
    converted_qtl <- rutilstimflutre::phasedJoinMapCP2qtl(jm, verbose = verbose)

    tmp1 <- t(cross$geno$`1`$data) ; geno_qtl <- tmp1
    for(i in 2:length(cross$geno)){
      tmp <- t(cross$geno[[i]]$data)
      geno_qtl <- rbind(geno_qtl, tmp)
    }
    geno_qtl <- geno_qtl[order(rownames(geno_qtl)),]
    converted_qtl <- converted_qtl[order(rownames(converted_qtl)),]
    colnames(geno_qtl) <- NULL ; colnames(converted_qtl) <- NULL
    stopifnot(identical(as.numeric(converted_qtl), as.numeric(geno_qtl)))
  }

  ## if no 2nd threshold provided, take the same as the first one
  if(length(threshold.alpha) == 1)
    threshold.alpha[2] = threshold.alpha[1]

  cross <- qtl::calc.genoprob(cross, step=1, map.function="kosambi")
  if(! response.in.cross){
    cross$pheno[[pheno.col]] <- NA
    stopifnot(!is.null(rownames(response)))
    for(i in 1:nrow(cross$pheno)){
      ind <- as.character(cross$pheno$indiv[i])
      if(ind %in% rownames(response)){
        cross$pheno[[pheno.col]][i] <- response[ind,]
      }
    }
  }

  ## Find LOD threshold to apply with permutations
  qtl.em.perm <- qtl::scanone(cross, pheno.col=pheno.col, method=method, n.perm = nperm, verbose=verbose)
  threshold.LOD <- ifelse(threshold.alpha[1] == 0.05,
                          threshold.LOD <- summary(qtl.em.perm)[1], # 5% alpha error
                          threshold.LOD <- summary(qtl.em.perm)[2]) # 10% alpha error

  qtl.em <- qtl::scanone(cross, pheno.col=pheno.col, method=method, verbose=verbose)
  qtl.em$chr <- as.numeric(qtl.em$chr)
  LOD.max <- 0 ; pos <- 0
  qtl.em$is.qtl <- FALSE ; qtl.em$nb_interval <- NA
  qtl.em$is.qtl[which(qtl.em$lod> threshold.LOD)] <- TRUE
  nb_interval <- 1

  for(i in 2:nrow(qtl.em)){
    if(!qtl.em$is.qtl[i-1] & qtl.em$is.qtl[i]){  # if beginning new qtl (F -> T)
      qtl.em$nb_interval[i] <- nb_interval
      nb_interval <- nb_interval + 1
    }
    if(qtl.em$is.qtl[i-1] & qtl.em$is.qtl[i]){ # same qtl (T -> T)
      qtl.em$nb_interval[i] <- qtl.em$nb_interval[i-1]
      if(i==2){ # QTL at the beginning
        qtl.em$nb_interval[c(i-1,i)] <- 1
        nb_interval <- nb_interval +1
      }
    }
  }
  ifelse(!all(is.na(qtl.em$nb_interval)),
         nb_interval <- max(qtl.em$nb_interval, na.rm=TRUE),
         nb_interval <- 0)
  stopifnot(nrow(qtl.em[qtl.em$is.qtl,]) == nrow(qtl.em[! is.na(qtl.em$nb_interval),]))
  if(nb_interval != 0){
    qtl.df <- data.frame(linkage.group=rep(0, nb_interval),
                         position.peak=rep(0, nb_interval), nearest.mrk = rep(0, nb_interval),
                         LOD=rep(0, nb_interval), interval.inf=rep(0, nb_interval),
                         interval.sup=rep(0,nb_interval))
    for(i in 1:nb_interval){
      qtl.df$linkage.group[i] <- unique(qtl.em$chr[which(qtl.em$nb_interval == i)])
      qtl.df$LOD[i] <- max(qtl.em[which(qtl.em$nb_interval == i),"lod"])
      qtl.df$position.peak[i] <- qtl.em$pos[qtl.em$nb_interval == i & qtl.em$lod == qtl.df$LOD[i]]
      qtl.df$nearest.mrk[i] <- qtl::find.marker(cross, qtl.df$linkage.group[i],  qtl.df$position.peak[i])
      qtl.df$interval.inf[i] <- min(qtl.em$pos[qtl.em$nb_interval == i & qtl.em$lod > qtl.df$LOD[i] - 1], na.rm=TRUE)
      qtl.df$interval.sup[i] <- max(qtl.em$pos[qtl.em$nb_interval == i & qtl.em$lod > qtl.df$LOD[i] - 1], na.rm=TRUE)
    }
    ## Plot LOD profile
    if(plot){
      for(link in unique(qtl.df$linkage.group)){
        plot(qtl.em, chr=link,
             main=paste0("Profile LOD score \n linkage group ",link),
             ylim=c(0, max(qtl.df$LOD[qtl.df$linkage.group == link])+2))
        panel.first=graphics::grid()
        graphics::abline(h=threshold.LOD, col="red")
        graphics::abline(v=qtl.df$interval.inf[qtl.df$linkage.group == link], col="blue", lty=2)
        graphics::abline(v=qtl.df$interval.sup[qtl.df$linkage.group == link], col="blue", lty=2)
        if(! is.null(QTL_position))
          graphics::abline(v=QTL_position$genetic.distance[QTL_position$linkage.group == link], col="green")
        graphics::text(x=qtl.df$position[qtl.df$linkage.group == link],
                       y=qtl.df$LOD[qtl.df$linkage.group == link]+1,
                       labels = paste0("Max LOD = ", round(qtl.df$LOD[qtl.df$linkage.group == link], 2),
                                       " \n at ", qtl.df$position[qtl.df$linkage.group == link], " cM"))
      }
    }

    ## Provide list of markers in the interval
    selected.markers <- NA
    for(i in 1:nb_interval){
      tmp <- colnames(cross$geno[[qtl.df$linkage.group[i]]]$map)[cross[["geno"]][[qtl.df$linkage.group[i]]]$map[1,] >=
                                                                 qtl.df$interval.inf[i] &
                                                                 cross$geno[[qtl.df$linkage.group[i]]]$map[1,] <= qtl.df$interval.sup[i]]
      selected.markers <- append(selected.markers, tmp)
    }
    selected.markers <- as.character(selected.markers)
    selected.markers <- subset(selected.markers, subset= !is.na(selected.markers))

  } else { # no QTL found
    qtl.df <- data.frame(chromosome=0,
                         position.peak=0, nearest.mrk = 0,
                         LOD=0, interval.inf=0, interval.sup=0)
    selected.markers <- 0 ; length(selected.markers) <- 0
  }

  ## Fit a linear model to estimate allelic effects
  if(! is.null(geno.joinmap)){
    if(length(selected.markers) != 0) {
      ## Convert cross object to a design matrix for selected markers
      tmp <- data.frame(locus=rownames(geno.joinmap), seg=rutilstimflutre::getJoinMapSegregs(geno.joinmap),
                        phase=phase, clas=0)
      jm <- cbind(tmp, geno.joinmap)
      jm[is.na(jm)] <- "--"
      jm$seg <- as.character(jm$seg)
      jm2 <- jm[selected.markers,]
      colnames(jm2)[5:ncol(jm2)] <- cross$pheno$indiv
      jm2[,c(5:ncol(jm2))] <- lapply(jm2[,c(5:ncol(jm2))], as.character)
      X <- rutilstimflutre::joinMap2designMatrix(jm=jm2, verbose=verbose)
      ## Find linear combination between columns and remove them
      lin.comb <- caret::findLinearCombos(X)
      X <- X[,-lin.comb$remove]
      X <- X[,-1] # remove intercept

      ## Estimate allele effects by multiple linear regression model
      fit.lm <- stats::lm(cross[["pheno"]][[pheno.col]]  ~ X)
      summary.fit <- summary(fit.lm)
      signif.coeff <- as.matrix(summary.fit$coefficients[,1][summary.fit$coefficients[,4] < threshold.alpha[2]])
      rownames(signif.coeff) <- substr(rownames(signif.coeff), start = 2, stop = nchar(signif.coeff))

      allelic.effects <- data.frame(predictor=colnames(rutilstimflutre::joinMap2designMatrix(jm=jm, verbose=0))[-1], effect=0)

      for(i in 1:nrow(signif.coeff)){
        allelic.effects[allelic.effects$predictor %in% rownames(signif.coeff)[i], "effect"] <- signif.coeff[i]
      }

    } else { # There is no QTL found
      tmp <- data.frame(locus=rownames(geno.joinmap),
                        seg=rutilstimflutre::getJoinMapSegregs(geno.joinmap),
                        phase=phase, clas=0)
      jm <- cbind(tmp, geno.joinmap)
      jm[is.na(jm)] <- "--"
      jm$seg <- as.character(jm$seg)
      allelic.effects <- data.frame(predictor=colnames(rutilstimflutre::joinMap2designMatrix(jm=jm, verbose=0))[-1], effect=0)
    }
  } else { # No genetic map entered
    allelic.effects <- NULL
  }
  out <- list(qtl.df=qtl.df, selected.markers=selected.markers, allelic.effects=allelic.effects)
  return(out)
}
