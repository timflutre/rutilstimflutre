## Contains functions useful to perform QTL detection by interval mapping
## with R/qtl.

##' QTL detection by SIM
##'
##' Perform QTL detection by simple interval mapping
##' @param cross object
##' @param numeric.chr.format logical to indicate if chromosome names are numeric
##' @param response.in.cross logical to indicate if response studied is in \code{cross$pheno}, default is TRUE
##' @param pheno.col character indicating column to study in \code{cross$pheno}
##' @param response named numeric or vector for response if not in \code{cross$pheno}, default is NULL
##' @param method method to detect QTL in \code{qtl::scanone}, default is "em"
##' @param geno.joinmap genotypes at markers in the JoinMap format, if NULL (default), no estimation of allelic effects is given
##' @param phase marker phases
##' @param threshold genomewide significance LOD threshold, if NULL (default), is found by permutations (with nperm parameter).
##' @param nperm number of permutations to be done in \code{qtl::scanone}, default is 100
##' @param alpha vector of length 1 or 2 (optional) with thresholds for (1) the significance of QTL presence (based on permutations) and (2) the significance of linear regression on QTL effect.
##' @param QTL_position matrix with genetic.distance and linkage.group as columns indicating QTL positions for plotting
##' @param plot logical, default is FALSE.
##' @param verbose verbosity level (0/1/2)
##' @return list of 3 elements: qtl.df is a data frame with QTL informations (linkage.group, position, LOD, interval.inf and interval.sup) / 
##' selected markers is a character vector for markers inside confidence interval /allelic effetcs is a data frame with a column predictor and a column effect with estimated allelic effects.
##' @author Charlotte Brault [aut], Timothee Flutre [ctb]
##' @seealso \code{\link{MIMQTL}}
##' @export
SIMQTL <- function (cross, numeric.chr.format=TRUE, method="em", 
                    response.in.cross=TRUE, pheno.col="y", response=NULL, 
                    geno.joinmap=NULL, phase, 
                    threshold=NULL,  nperm=100, alpha=c(0.05,0.2),
                    plot=FALSE, QTL_position=NULL, verbose=0){
  requireNamespace(c("qtl", "caret"))
  
  ## Reformat chromosome names
  if(! numeric.chr.format){
    nb_chr <- length(lapply(cross$geno, attributes))
    chr_renamed <- seq(1:nb_chr)
    names(cross$geno) <- chr_renamed
  }
  
  ## Verifications
  stopifnot(all(is.numeric(alpha)),
            xor(all(!response.in.cross, !is.null(response)),
                all(response.in.cross & is.null(response))))
  if(! is.null(QTL_position))
    stopifnot(c("linkage.group", "genetic.distance") %in% colnames(QTL_position))
  
  ## Check the conformity between geno.joinmap and cross
  if(!is.null(geno.joinmap)){
    jm <- data.frame(seg=getJoinMapSegregs(geno.joinmap), phase=phase)
    jm <- cbind(jm, geno.joinmap)
    jm[jm == "--"] <- NA
    converted_qtl <- phasedJoinMapCP2qtl(jm, verbose=verbose)
    
    tmp1 <- t(cross$geno$`1`$data) ; geno_qtl <- tmp1
    for(i in 2:length(cross$geno)){
      tmp <- t(cross$geno[[i]]$data)
      geno_qtl <- rbind(geno_qtl, tmp)
    }
    colnames(geno_qtl) <-cross$pheno$indiv 
    geno_qtl <- geno_qtl[order(rownames(geno_qtl)),]
    converted_qtl <- converted_qtl[order(rownames(converted_qtl)),]
    
    ### Keep same genotypes 
    converted_qtl <- converted_qtl[, colnames(converted_qtl) %in% colnames(geno_qtl)]
    geno_qtl <- geno_qtl[, colnames(geno_qtl) %in% colnames(converted_qtl)]
    colnames(geno_qtl) <- NULL ; colnames(converted_qtl) <- NULL
    converted_qtl <- apply(converted_qtl, 2, as.numeric)
    geno_qtl <- apply(geno_qtl, 2, as.numeric)
    stopifnot(identical(converted_qtl, geno_qtl))
  }
  
  ## if no 2nd threshold provided, take the same as the first one
  if(length(alpha) == 1)
    alpha[2] <- alpha[1]
  
  ## Apply calc.genoprob to get genotype probabilities
  cross <- qtl::calc.genoprob(cross, step=1, map.function="kosambi")
  
  ## If response studied is not in object cross, add to it
  if(! response.in.cross){
    pheno.col <- ifelse(is.null(colnames(response)), pheno.col, colnames(response))
    cross$pheno[[pheno.col]] <- NA
    stopifnot(!is.null(rownames(response)))
    for(i in 1:nrow(cross$pheno)){
      ind <- as.character(cross$pheno$indiv[i])
      if(ind %in% rownames(response)){
        cross$pheno[[pheno.col]][i] <- response[ind,]
      }
    }
  }
  
  ## If threshold is not given, apply permutations to find it with error rate alpha
  if(is.null(threshold)){
    qtl.em.perm <- qtl::scanone(cross, pheno.col=pheno.col, method=method,
                                n.perm=nperm, verbose=verbose)
    threshold <- summary(qtl.em.perm, alpha=alpha[1])
  } 
  
  ## Apply scanone to get LOD profile
  qtl.em <- qtl::scanone(cross, pheno.col=pheno.col, method=method, verbose=verbose)
  qtl.em$chr <- as.numeric(qtl.em$chr)
  LOD.max <- 0 ; pos <- 0
  qtl.em$is.qtl <- FALSE ; qtl.em$nb_interval <- NA
  qtl.em$is.qtl[which(qtl.em$lod > as.numeric(threshold))] <- TRUE
  nb_interval <- 1
  
  ## Find QTL intervals
  for(i in 2:nrow(qtl.em)){
    if(!qtl.em$is.qtl[i-1] & qtl.em$is.qtl[i] | # if beginning new qtl (F -> T)
       qtl.em$is.qtl[i-1] & qtl.em$is.qtl[i] & qtl.em$chr[i-1] != qtl.em$chr[i]){ # T -> T new chr = new qtl
      qtl.em$nb_interval[i] <- nb_interval
      nb_interval <- nb_interval + 1
    }
    if(qtl.em$is.qtl[i-1] & qtl.em$is.qtl[i] & qtl.em$chr[i-1] == qtl.em$chr[i]){ # same qtl (T -> T)
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
                         position=0, nearest.mrk=0,
                         LOD=0, interval.inf=0,
                         interval.sup=0)
    for(i in 1:nb_interval){
      qtl.df$linkage.group[i] <- unique(qtl.em$chr[which(qtl.em$nb_interval == i)])
      qtl.df$LOD[i] <- max(qtl.em[which(qtl.em$nb_interval == i),"lod"])
      qtl.df$position[i] <- qtl.em$pos[qtl.em$nb_interval == i & qtl.em$lod == qtl.df$LOD[i]]
      qtl.df$nearest.mrk[i] <- qtl::find.marker(cross, qtl.df$linkage.group[i],  qtl.df$position[i])
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
        graphics::abline(h=threshold, col="red")
        graphics::abline(v=qtl.df$interval.inf[qtl.df$linkage.group == link], col="blue", lty=2)
        graphics::abline(v=qtl.df$interval.sup[qtl.df$linkage.group == link], col="blue", lty=2)
        if(! is.null(QTL_position))
          graphics::abline(v=QTL_position$genetic.distance[QTL_position$linkage.group == link], col="green")
        graphics::text(x=qtl.df$position[qtl.df$linkage.group == link],
                       y=qtl.df$LOD[qtl.df$linkage.group == link]+1,
                       labels=paste0("Max LOD = ", round(qtl.df$LOD[qtl.df$linkage.group == link], 2),
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
    qtl.df <- data.frame(linkage.group=NA,
                         position=NA, nearest.mrk=NA,
                         LOD=NA, interval.inf=NA, interval.sup=NA)
    selected.markers <- 0 ; length(selected.markers) <- 0
  }
  
  ## Fit a linear model to estimate allelic effects
  if(! is.null(geno.joinmap)){
    if(length(selected.markers) != 0) {
      ## Convert cross object to a design matrix for selected markers
      tmp <- data.frame(locus=rownames(geno.joinmap),
                        seg=getJoinMapSegregs(geno.joinmap),
                        phase=phase, clas=0)
      jm <- cbind(tmp, geno.joinmap)
      jm[is.na(jm)] <- "--"
      jm$seg <- as.character(jm$seg)
      jm2 <- jm[selected.markers,]
      colnames(jm2)[5:ncol(jm2)] <- cross$pheno$indiv
      jm2[,c(5:ncol(jm2))] <- lapply(jm2[,c(5:ncol(jm2))], as.character)
      X <- joinMap2designMatrix(jm=jm2, verbose=verbose)
      ## Find linear combination between columns and remove them
      lin.comb <- caret::findLinearCombos(X)
      X <- X[,-lin.comb$remove]
      X <- X[,-1] # remove intercept
      
      ## Estimate allele effects by multiple linear regression model
      fit.lm <- stats::lm(cross[["pheno"]][[pheno.col]]  ~ X)
      coeff <- fit.lm$coefficients[-1]
      names(coeff) <- substr(names(coeff), start=2, stop=nchar(names(coeff)))
      if(!is.na(alpha[2])){
        summary.fit <- summary(fit.lm)$coefficient
        rownames(summary.fit)[-1] <- substr(rownames(summary.fit)[-1], start=2, stop=nchar(rownames(summary.fit)[-1]))
        coeff <- as.matrix(summary.fit[,1][summary.fit[,4] < alpha[2]])
      } 
      
      allelic.effects <- data.frame(predictor=colnames(joinMap2designMatrix(jm=jm, verbose=0))[-1], effect=0)
      for(i in 1:length(coeff)){
        allelic.effects[allelic.effects$predictor %in% rownames(coeff)[i], "effect"] <- coeff[i]
      }
      
    } else { # There is no QTL found
      tmp <- data.frame(locus=rownames(geno.joinmap),
                        seg=getJoinMapSegregs(geno.joinmap),
                        phase=phase, clas=0)
      jm <- cbind(tmp, geno.joinmap)
      jm[is.na(jm)] <- "--"
      jm$seg <- as.character(jm$seg)
      allelic.effects <- data.frame(predictor=colnames(joinMap2designMatrix(jm=jm, verbose=0))[-1], effect=0)
    }
  } else { # No genetic map entered
    allelic.effects <- NULL
  }
  
  out <- list(qtl.df=qtl.df,
              selected.markers=selected.markers,
              allelic.effects=allelic.effects)
  return(out)
}

##' make.formula
##'
##' Useful function for MIMQTL.
##' @param cross object
##' @param big_list list
##' @param range_qtl vector
##' @param nb_run integer
##' @param pLOD logical
##' @seealso \code{\link{MIMQTL}}
##' @return list
make.formula <- function(cross, big_list, range_qtl, nb_run, pLOD=FALSE){
  requireNamespace("qtl")

  qtltemp <- qtl::makeqtl(cross,
                          chr=big_list[[range_qtl]][[nb_run]]$chr,
                          pos=big_list[[range_qtl]][[nb_run]]$pos, what="prob")
  if(pLOD)
    pLOD <- attr(big_list[[range_qtl]][[nb_run]], "pLOD")

  form.tmp <- attr(big_list[[range_qtl]][[nb_run]], "formula")
  temp <- strsplit(form.tmp," + ",fixed=TRUE)

  # If there are interaction terms, retrieve the formula without interactions
  if (length(temp[[1]]) > length(qtltemp$chr)) {
    form.int.tmp <- temp[[1]][1]

    for (k in 2:length(big_list[[range_qtl]][[1]]$chr)) {
      form.int.tmp <- paste0(form.int.tmp," + ", temp[[1]][k])
    }
  } else {
    form.int.tmp <- form.tmp
  }
  out <- list(qtltemp=qtltemp, pLOD=pLOD,
              form.tmp=form.tmp, temp=temp, form.int.tmp=form.int.tmp)
  return(out)
}

##' QTL detection by MIM
##'
##' Perform QTL detection by multiple interval mapping
##' @param cross object
##' @param response.in.cross logical
##' @param response character
##' @param numeric.chr.format logical
##' @param geno.joinmap genotypes at markers in the JoinMap format
##' @param phase phase
##' @param pheno.col character
##' @param nperm number of permutations to be done in \code{qtl::scantwo}
##' @param threshold.alpha vector of length 1 or 2 (optional) with thresholds for (1) the significance of QTL presence (based on permutations) and (2) the significance of QTL effects
##' @param QTL_position useful for simulated data
##' @param method method to detect QTL in \code{qtl::scantwo}
##' @param plot logical
##' @param range_nb_qtl_max vector
##' @param nrun integer
##' @param additive.only logical
##' @param p2d character
##' @param scan2file optional path to file
##' @param type.CI type of confidence interval
##' @param nb.cores number of cores
##' @param verbose verbosity level (0/1/2)
##' @return list
##' @author Agnes Doligez [aut], Charlotte Brault [ctbt], Timothee Flutre [ctb]
##' @seealso \code{\link{SIMQTL}}
##' @export
MIMQTL <- function(cross, response.in.cross=TRUE, response=NULL,
                   numeric.chr.format=FALSE, geno.joinmap=NULL,
                   phase=NULL, pheno.col="y", nperm=100,
                   threshold.alpha =0.05, QTL_position=NULL,
                   method="hk", plot=c(FALSE, FALSE),
                   range_nb_qtl_max=seq(1:5),
                   nrun=10, additive.only=TRUE, p2d="",
                   scan2file="", type.CI="LOD-1", nb.cores=1,
                   verbose=0){
  requireNamespace(c("qtl", "caret"))

  ## Reformat chromosome names
  if(! numeric.chr.format){
    nb_chr <- length(lapply(cross$geno, attributes))
    chr_renamed <- seq(1:nb_chr)
    names(cross$geno) <- chr_renamed
  }

  ## Verifications
  stopifnot(max(range_nb_qtl_max) < 15,
            method %in% c("hk", "em"),
            type.CI %in% c("LOD-1", "LOD-2"),
            xor(all(!response.in.cross, !is.null(response)),
                all(response.in.cross & is.null(response))))
  if(! is.null(QTL_position))
    stopifnot(c("linkage.group", "genetic.distance") %in% colnames(QTL_position))
  ## Check the conformity between geno.joinmap and cross
  if(!is.null(geno.joinmap)){
    jm <- data.frame(seg=getJoinMapSegregs(geno.joinmap), phase=phase)
    jm <- cbind(jm, geno.joinmap)
    jm[jm == "--"] <- NA
    converted_qtl <- phasedJoinMapCP2qtl(jm, verbose=verbose)

    tmp1 <- t(cross$geno$`1`$data) ; geno_qtl <- tmp1
    for(i in 2:length(cross$geno)){
      tmp <- t(cross$geno[[i]]$data)
      geno_qtl <- rbind(geno_qtl, tmp)
    }
    geno_qtl <- geno_qtl[order(rownames(geno_qtl)),]
    converted_qtl <- converted_qtl[order(rownames(converted_qtl)),]
    colnames(geno_qtl) <- NULL ; colnames(converted_qtl) <- NULL

    stopifnot(identical(converted_qtl, geno_qtl))
  }

  ## if no 2nd threshold provided, take the same as the first one
  if(length(threshold.alpha) == 1)
    threshold.alpha[2] <- threshold.alpha[1]

  ## Add phenotypes to cross object
  if( ! response.in.cross){
    cross$pheno[[pheno.col]] <- NA
    for(i in 1:nrow(cross$pheno)){
      ind <- as.character(cross$pheno$indiv[i])
      if(ind %in% rownames(response)){
        cross$pheno[[pheno.col]][i] <- response[ind, pheno.col]
      }
    }
  }
  stopifnot(!all(is.na(cross$pheno[[pheno.col]])))
  cross <- qtl::calc.genoprob(cross, step=1, map.function="kosambi")

  ## Run scantwo
  if(file.exists(paste0(p2d, "/", scan2file)) & scan2file != ""){
    load(paste0(p2d, "/", scan2file))
  } else {
    system.time(sc2 <- qtl::scantwo(cross=cross, pheno.col=pheno.col,
                                    model="normal", n.perm=nperm, method=method,
                                    verbose=FALSE,
                                        #maxit=100000, tol=1e-2,
                                    n.cluster=nb.cores))
    if(scan2file != "" && p2d != "")
      save(sc2, file=paste0(p2d, "/", scan2file))
  }
  stopifnot("sc2" %in% ls())

  ## Run stepwiseqtl
  calc_pen <- qtl::calc.penalties(perms=sc2, alpha=threshold.alpha[1])
  big_list <- list()
  small_list <- list()
  tmp <- list()
  small_list[[1]] <- tmp
  step <- 1
  for (nb_qtl in range_nb_qtl_max){
    for(nb_run in 1:nrun){
      small_list[[nb_run]] <- qtl::stepwiseqtl(cross, pheno.col=pheno.col, max.qtl=nb_qtl,
                                               method=method,  penalties= calc_pen,
                                               model="normal", incl.markers=TRUE,
                                               refine.locations=TRUE, additive.only=additive.only,
                                               scan.pairs=FALSE, keeplodprofile=FALSE, keeptrace=T,
                                               verbose=verbose, tol=1e-4, maxit=1000)
    }
    big_list[[step]] <- small_list
    if(verbose > 0)
      write(paste("number of max qtl tested: ", step, "out of", length(range_nb_qtl_max)), stderr())
    step <- step + 1
  }

  ## Plot model Comparison
  if(plot[1]){
    for (nb_run in 1:nrun) {
      toprint <- ""
      for (range_qtl in 1:length(range_nb_qtl_max)){
        if (attr(big_list[[range_qtl]][[nb_run]],"pLOD") == 0) {    # if Null model
          toprint <- paste(toprint,"0")
        } else {
          toprint <- paste(toprint,length(big_list[[range_qtl]][[nb_run]]$name))
        }
      }
    }
    ## vizualize model selection steps for the nrun with a given value of nb max qtl
    m <- ceiling(mean(range_nb_qtl_max)) # select a value for nb max qtl
    for (nb_run in 1:1) { # why ?
      if (attr(big_list[[m]][[nb_run]],"pLOD") == 0) {    # if Null model
        writeLines("Null QTL model")
      } else {
        thetrace<-attr(big_list[[m]][[nb_run]],"trace")
        for(k in seq(along=thetrace))
          print(qtl::plotModel(thetrace[[k]], chronly=T,
                               main=paste(k,": pLOD = ", round(attr(thetrace[[k]], "pLOD"), 2),
                                          "\n for ", m, " number of max qtl")))
      }
    }
  }

  ## Model selection
  prec <- list()  #  nb of QTLs for the preceding value of max.qtl
  qtl <- list()  # selected qtl objects
  form <- list()  # all formulas without interaction
  form.int <- list() #  all formulas with interaction

  ## 0: For each value of max.qtl parameter, create a list containing all qtl objects for fitqtl +
  ## a list containing all formulas for fitqtl + a second list containing all formulas for fitqtl with no interaction terms
  for (range_qtl in range_nb_qtl_max) {
    qtltemp <- list()  #  selected qtl objects
    pLOD <- list() #  compare pLOD among outcomes from different runs with the same nb of QTLs
    form.temp <- list()  #  all formulas with no interaction
    form.int.temp <- list()  # all formulas with interaction terms
    nbid1 <- 0  # nb of runs with an outcome identical to the 1st outcome
    nbnull <- 0 # nb of runs with a Null model as outcome
    temp <- list()
    pLOD <- 0

    ## 1:If the outcome of 1st run is Null model
    if (attr(big_list[[range_qtl]][[1]],"pLOD") == 0) {
      nbnull <- nbnull + 1
    }

    for (nb_run in 2:nrun) {
      ## 2.1: If other outcome is  Null model
      if (attr(big_list[[range_qtl]][[nb_run]],"pLOD") == 0) {
        nbnull <- nbnull + 1
        ## 2.2: Else if all is identical to the 1st model
      } else if(identical(big_list[[range_qtl]][[1]]$pos, big_list[[range_qtl]][[nb_run]]$pos) &
                identical(big_list[[range_qtl]][[1]]$chr, big_list[[range_qtl]][[nb_run]]$chr) &
                identical(attr(big_list[[range_qtl]][[1]],"pLOD"), attr(big_list[[range_qtl]][[nb_run]],"pLOD")) &
                identical(attr(big_list[[range_qtl]][[1]],"formula"),
                          attr(big_list[[range_qtl]][[nb_run]],"formula"))) {
        nbid1 <- nbid1 + 1 # outcome identical to the 1st
      }
    }

    ## 3.1: All or some outcomes are Null model
    if (nbnull == nrun | nbnull > 0) {
      qtltemp <- "NO QTL"
      form.tmp <- "NO QTL"
      qtl <- qtltemp #### Agnes : nb_run'ai rajoute cette ligne pour ne pas avoir de liste nulle pour qtl quand il n'y a pas de QTL, est-ce que c'est correct ?

      ## 3.2: No Null model, several possible cases
    } else if (nbnull == 0) {
      ##  3.2.1: All nrun outcomes are identical, keep the 1st one
      if (nbid1 == nrun - 1) {
        out <- make.formula(cross, big_list, range_qtl, nb_run=1, pLOD=FALSE)
        ## 3.2.2: Some outcomes are different from the 1st one
                                        #TODO no stability
      }else {
        out <- make.formula(cross, big_list, range_qtl, nb_run=1, pLOD=TRUE)

        for (nb_run in 2:nrun) {
          ## 3.2.2.1 If nb of QTLs < for the currently kept outcome, keep this outcome
          if (length(big_list[[range_qtl]][[nb_run]]$chr) < length(qtltemp$chr)) {
            out <- make.formula(cross, big_list, range_qtl, nb_run=nb_run, pLOD=FALSE)

            ## 3.2.2.2: If same nb of QTLs but pLOD is > for the currently outcome, keep this outcome
          } else if ((length(big_list[[range_qtl]][[nb_run]]$chr) == length(qtltemp$chr)) &
                     (attr(big_list[[range_qtl]][[nb_run]],"pLOD") > pLOD)) {
            out <- make.formula(cross, big_list, range_qtl, nb_run=nb_run, pLOD=TRUE)

          } # end else if pLOD is > fir this outcome
        } # end for nb_run in 2:nrun
      } # end of else: some outcomes are different from the first one

      qtltemp=out$qtltemp ; pLOD=out$pLOD
      form.tmp=out$form.tmp ; temp=out$temp ; form.int.tmp=out$form.int.tmp

      ## Keep the results of the current max.qtl value only if the nb of QTLs has increased by one compared to the preceding value of max.qtl
      ## (assumes that qtl nb cannot decrease when max.qtl value increases)

      if (range_qtl == range_nb_qtl_max[1]) {
        ## keep the 1st max.qtl value as a default (includes any Null model)
        qtl <- qtltemp
        form <- form.tmp
        form.int <- form.int.tmp
        if (all(qtl == "NO QTL")) {
          prec <- 0
        } else {
          prec <- length(qtltemp$chr)
        }
      } else if (all(qtl == "NO QTL")) {
        prec <- 0

        ## If the nb of QTLs is one more than with the preceding max.qtl value
      } else if (length(qtltemp$chr) == prec + 1) {
        qtl <- qtltemp
        form <- form.tmp
        form.int <- form.int.tmp
        prec <- length(qtltemp$chr)
      }
    } # end else if no null model
  }  # end for range max qtl

  ## Fit the best model
  fit <- list()  # list of fit qtl results without interaction
  fit.int <- list() # list of fit qtl results with intarction
  if (all(qtl == "NO QTL")) {   # if Null model
    fit <- "NO QTL"
    fit.int <- "NO QTL"
  } else {
    fit <- suppressWarnings(qtl::fitqtl(cross, pheno.col=pheno.col, qtl=qtl,
                                        formula=form, method=method, get.ests=T))
    fit.int <- suppressWarnings(qtl::fitqtl(cross, pheno.col=pheno.col, qtl=qtl,
                                            formula=form.int, method=method, get.ests=T))
  }

  ## Put all model information together
  if (all(fit[1] == "NO QTL")) {
    temp <- as.data.frame(t(rep("NO QTL", 6)))
    names(temp) <-
      c("df", "SS", "MS", "LOD", "perc.var", "Pvalue(Chi2)")# add pvalue ??
  } else {
    colnames(fit$result.full)[5] <- "perc.var"
    if (!is.null(fit$result.drop))
      colnames(fit$result.drop)[4] <- "perc.var"
    temp <- as.data.frame(fit$result.full[, 1:6])
  }
  temp$Trait <- pheno.col  # adds a column with Trait name
  ## converts the content of the 6 first columns to "factor" type for rbind
  temp[,c(1:6)] <- lapply(temp[, c(1:6)], as.factor)
  mod <- temp
  mod$ord <- seq(from=1, by=1, to=nrow(mod))

  ## Summary
  if (length(attr(fit, "names")) < 3) {
    ## Null Model
    temp <- as.data.frame(t(rep("NO QTL", 5)))            ## why t ?
    names(temp) <- c("df", "Type III SS", "LOD", "perc.var")
  } else if (length(attr(fit, "names")) == 3) {
                                        # model with only one QTL
    temp <- as.data.frame(t(fit$result.full[1, c(1,3:6)]));
    names(temp)[2] <- "Type III SS"
  } else if (length(attr(fit, "names")) > 3) {
    ## model with more than one QTL
    temp <- as.data.frame(fit$result.drop[,1:5])
  }
  temp$Trait <- pheno.col #  add a variable Trait because rownames are not explicit for traits with no QTLs

  if (length(attr(fit, "names")) == 0 ) {
    ## Null Model
    temp$Terms <- "NO QTL"
    temp$est.BC <- NA
    temp$est.AD <- NA
    temp$est.BD <- NA

  }else if (form == form.int) {
    nbQTLSeg <- length(names(fit$est$ests))
    namefit <- names(fit$est$ests)

    ## If there is no interaction term in the model
    temp$Terms <- substr(namefit[seq(from=2, by=3, length.out=((nbQTLSeg-1)/3))], 1,
                         nchar(names((fit$est$ests)[seq(from=2, by=3,
                                                        length.out=((nbQTLSeg-1)/3))]))-3)
    if(verbose > 0)
      print(temp$Terms)
    temp$est.BC <- t(fit$est$ests)[seq(from=2, by=3, length.out=((nbQTLSeg - 1) /3))]
    temp$est.AD <- t(fit$est$ests)[seq(from=3, by=3, length.out=((nbQTLSeg - 1) /3))]
    temp$est.BD <- t(fit$est$ests)[seq(from=4, by=3, length.out=((nbQTLSeg - 1) /3))]

  }else {
    ## If there is an interaction term in the model
    temp$Terms <- rownames(temp);
    temp$est.BC <- c(t(fit$est$ests)[seq(from=2,by=3,
                                         length.out=((length(names(fit$est$ests))-1)/3))],
                     rep(NA,nrow(temp)-nrow(fit$result.drop)))
    temp$est.AD <- c(t(fit$est$ests)[seq(from=3,by=3,
                                         length.out=((length(names(fit$est$ests))-1)/3))],
                     rep(NA,nrow(temp) - nrow(fit$result.drop)))
    temp$est.BD <- c(t(fit$est$ests)[seq(from=4,by=3,
                                         length.out=((length(names(fit$est$ests))-1)/3))],
                     rep(NA,nrow(temp) - nrow(fit$result.drop)))
  }

  QTLest <- temp
  tmp <- strsplit(QTLest$Terms, split="@")
  QTLest$linkage.group <- suppressWarnings(as.numeric(unlist(lapply(tmp, `[`, 1))))
  QTLest$position <- suppressWarnings(as.numeric(unlist(lapply(tmp, `[`, 2))))
  ## Add a column for further order purpose
  QTLest$ord <- seq(from=1, to=nrow(QTLest), by=1)
  if(all(QTLest$df == "NO QTL")){
    QTLest[,c(1:ncol(QTLest))] <- NA
  }

  ## Refine QTL location
  refine <- list()
  if (all(qtl == "NO QTL")) {   # if Null model
    refine <- "NO QTL"
  } else {
    refine <- qtl::refineqtl(cross, pheno.col=pheno.col, qtl=qtl,
                             formula=form, method=method, verbose=verbose,
                             keeplodprofile=TRUE)
  }

  ## Confidence intervals
  if(length(attr(refine, "names")) < 9) {    # if Null model
    temp4 <- as.data.frame(cbind(Trait=pheno.col, linkage.group="NO QTL", position="NO QTL",
                                 interval.inf="NO QTL", interval.sup="NO QTL", CI.int="NO QTL"))
  } else {
    interval.inf <- 0
    interval.sup <- 0
    for (j in 1:length(refine$name)) {       # for each QTL
      Trait <- pheno.col
      Terms <- names(attr(refine,"lodprofile")[j])

      ## Type of Confidence interval is LOD - 1
      if(type.CI == "LOD-1"){
        temp <- qtl::lodint(refine, qtl.index=j, drop=1)

        ## Type of Confidence interval is LOD - 2
      } else if(type.CI == "LOD-2"){
        temp <- qtl::lodint(refine, qtl.index=j, drop=2)
      }
      interval.inf <- as.numeric(min(temp$pos))
      interval.sup <- as.numeric(max(temp$pos))
      position <- as.numeric(temp$pos[2])
      CI.int <- interval.sup - interval.inf
      linkage.group <- as.numeric(strsplit(Terms, split="@")[[1]][1])
      temp3 <- as.data.frame(cbind(Trait, linkage.group, position,
                                   interval.inf, interval.sup, CI.int))

      if (j == 1) {
        temp4 <- temp3
      } else {
        temp4 <- as.data.frame(rbind(temp4, temp3))
      }
    }
  }
  CI <- temp4
  CI[,c(2:6)] <- suppressWarnings(lapply(CI[,c(2:6)], as.character))
  CI[,c(2:6)] <- suppressWarnings(lapply(CI[,c(2:6)], as.numeric))

  ## QTL table
  if(all(apply(QTLest, 1, is.na)) & all(apply(CI[-1], 1, is.na))){
    qtl.df <- data.frame("Trait"=NA, "linkage.group"=NA, "position"=NA, "LOD"=NA,
                         "perc.var"=NA, "nearest.marker"=NA, "interval.inf"=NA, "interval.sup"=NA,
                         "est.BC"=NA, "est.AD"=NA, "est.BD"=NA,
                         "effect.Af"=NA, "effect.Am"=NA, "effect.D"=NA)
  } else {
    qtl.df <- merge(QTLest, CI, by=c("Trait", "linkage.group", "position"), all=TRUE)
    if(all(is.numeric(qtl.df[,c("LOD", "perc.var")])))
      qtl.df[,c("LOD", "perc.var")] <- round(qtl.df[,c("LOD", "perc.var")], 2)
    qtl.df <- qtl.df[order(qtl.df$ord),]
  }

  ## Formulas from V.Segura
  qtl.df$Af <- 0.25 * (qtl.df$est.AD - qtl.df$est.BD - qtl.df$est.BC) #female
  qtl.df$Am <- 0.25 * (qtl.df$est.BC - qtl.df$est.AD - qtl.df$est.BD) #male
  qtl.df$D <- 0.25 * (qtl.df$est.BD - qtl.df$est.BC - qtl.df$est.AD) #dominance

  qtl.df$tot <- abs(qtl.df$Af) + abs(qtl.df$Am) + abs(qtl.df$D)

  qtl.df[is.na(qtl.df$Af) == F &
         abs(qtl.df$Af) / qtl.df$tot > 0.3, "effect.Af"] <- "Af"
  qtl.df[is.na(qtl.df$Am) == F &
         abs(qtl.df$Am) / qtl.df$tot > 0.3, "effect.Am"] <- "Am"
  qtl.df[is.na(qtl.df$D) == F &
         abs(qtl.df$D) / qtl.df$tot > 0.3, "effect.D"] <- "D"

  qtl.df$effect <- paste0(qtl.df$effect.Af, ", ",
                          qtl.df$effect.Am, ", ",
                          qtl.df$effect.D)


  ## Marker in Confidence Interval
  if(all(! is.na(qtl.df$linkage.group))){
    selected.markers <- NA
    for(i in 1:nrow(qtl.df)){
      chr <- qtl.df$linkage.group[i]
      tmp <- colnames(cross[["geno"]][[chr]]$map)[cross[["geno"]][[chr]]$map[1,] >= qtl.df$interval.inf[i] &
                                                  cross[["geno"]][[chr]]$map[1,] <= qtl.df$interval.sup[i]]
      selected.markers <- as.character(c(selected.markers, tmp))
    }
    selected.markers <- subset(selected.markers, subset= !is.na(selected.markers))

  } else { # no QTL found
    selected.markers <- 0 ; length(selected.markers) <- 0
  }

  ## Allelic effect estimation
  if(all(! is.na(qtl.df$linkage.group))){
    qtl.df$nearest.marker <- qtl::find.marker(cross, qtl.df$linkage.group, qtl.df$position)
  } else {
    qtl.df[c("interval.inf", "interval.sup", "position", "nearest.marker")] <- NA
  }

  qtl.df <- subset(qtl.df, select=c("Trait", "linkage.group", "position", "LOD",
                                    "perc.var", "nearest.marker", "interval.inf", "interval.sup",
                                    "est.BC", "est.AD", "est.BD",
                                    "effect.Af", "effect.Am", "effect.D"))
  if(p2d != "")
    utils::write.table(qtl.df, file=paste0(p2d, "/summary_qtl_table-", pheno.col,".tsv"),
                       sep="\t", na="NA", dec=".", col.names=TRUE, row.names=FALSE)

  ## Plot LOD Profile
  if(plot[2] & all(!is.na(qtl.df$linkage.group))){
    for(link in unique(levels(as.factor(qtl.df$linkage.group)))){
      qtl::plotLodProfile(qtl=refine, qtl.labels=FALSE, chr=link,
                          main=paste0("Profile LOD score \n linkage group ",link),
                          ylim=c(0, max(qtl.df$LOD[qtl.df$linkage.group == link])+2))
      panel.first=graphics::grid()
      graphics::abline(h=calc_pen[1], col="red") # display main penalty
      graphics::abline(v=qtl.df$interval.inf[qtl.df$linkage.group == link], col="blue", lty=2)
      graphics::abline(v=qtl.df$interval.sup[qtl.df$linkage.group == link], col="blue", lty=2)
      if(! is.null(QTL_position))
        graphics::abline(v=QTL_position$genetic.distance[QTL_position$linkage.group == link], col="green")
      graphics::text(x=qtl.df$position[qtl.df$linkage.group == link],
                     y=qtl.df$LOD[qtl.df$linkage.group == link]+1,
                     labels=paste0("Max LOD = ", round(qtl.df$LOD[qtl.df$linkage.group == link],2),
                                     " \n at ", qtl.df$position[qtl.df$linkage.group == link], " cM"))
    }
  }

  ## Marker in Confidence Interval
  if(all(!is.na(qtl.df$linkage.group))){
    selected.markers <- NA
    for(i in 1:nrow(qtl.df)){
      chr <- qtl.df$linkage.group[i]
      tmp <- colnames(cross[["geno"]][[chr]]$map)[cross[["geno"]][[chr]]$map[1,] >= qtl.df$interval.inf[i] &
                                                  cross[["geno"]][[chr]]$map[1,] <= qtl.df$interval.sup[i]]
      selected.markers <- as.character(c(selected.markers, tmp))
    }
    selected.markers <- subset(selected.markers, subset= !is.na(selected.markers))

  } else { # no QTL found
    selected.markers <- 0 ; length(selected.markers) <- 0
  }

  ## Estimated allelic effect
  ## Fit a linear model to estimate allelic effects
  if(! is.null(geno.joinmap)){
    tmp <- data.frame(locus=rownames(geno.joinmap),
                      seg=getJoinMapSegregs(geno.joinmap),
                      phase=phase, clas=0)
    jm <- cbind(tmp, geno.joinmap)
    jm[is.na(jm)] <- "--"
    jm$seg <- as.character(jm$seg)
    if(length(selected.markers) != 0) {
      ## Convert cross object to a design matrix for selected markers
      jm2 <- subset(jm, subset=rownames(geno.joinmap) %in% selected.markers)
      X <- joinMap2designMatrix(jm=jm2, verbose=0)
      ## Find linear combination between columns and remove them
      lin.comb <- caret::findLinearCombos(X)
      X <- X[,-lin.comb$remove]
      X <- X[,-1] # remove intercept

      ## Estimate allele effects by multiple linear regression model
      fit.lm <- stats::lm(cross[["pheno"]][[pheno.col]]  ~ X)
      summary.fit <- summary(fit.lm)
      signif.coeff <- as.matrix(summary.fit$coefficients[,1][summary.fit$coefficients[,4] < threshold.alpha[2]])
      rownames(signif.coeff) <- substr(rownames(signif.coeff), start=2, stop=nchar(signif.coeff))

      allelic.effects <- data.frame(predictor=colnames(joinMap2designMatrix(jm=jm, verbose=0))[-1], effect=0)

      for(i in 1:nrow(signif.coeff)){
        allelic.effects[allelic.effects$predictor %in% rownames(signif.coeff)[i], "effect"] <- signif.coeff[i]
      }

    }else { # There is no QTL found
      allelic.effects <- data.frame(predictor=colnames(joinMap2designMatrix(jm=jm,
                                                                            verbose=0))[-1],
                                    effect=0)
    }
  } else { # No genetic map entered
    allelic.effects <- NULL
  }
  if(verbose > 0)
    if(all(! is.na(qtl.df$linkage.group))){
      write(paste0("QTL is found in linkage group ", qtl.df$linkage.group), stderr())
    } else if (is.na(qtl.df$linkage.group))
      write(paste0("No QTL found"), stderr())

  out <- list(qtl.df=qtl.df,
              selected.markers=selected.markers,
              allelic.effects=allelic.effects)

  return(out)
}
