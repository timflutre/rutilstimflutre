## Contains functions useful for grapevine projects.

##' Plot the map of SxG
##'
##' reciprocal cross Syrah x Grenache at the Domaine du Chapitre (near Montpellier, France)
##' @param map data.frame with columns rank, location, x and y
##' @param main main title (there is one by default in both languages)
##' @param lang language for the title and legend (default=FR/EN)
##' @param col.block vector with block colors (set to NULL to ignore)
##' @param show.plant.types show the plants corresponding to F1, parents and controls
##' @param show.pu.type show the "parcelles unitaires" corresponding to F1, parents and controls
##' @param x.off offset on the x-axis to let some space for the legend
##' @param mar margins (see the "mar" option of "par()")
##' @return the "map" data.frame as an invisible object
##' @author Timothee Flutre
##' @export
plotMapSxGDomaineChapitre <- function(map,
                                      main="default",
                                      lang="FR",
                                      col.block=stats::setNames(
                                          c("grey81", "bisque1"),
                                          c("A","B")),
                                      show.plant.types=TRUE,
                                      show.pu.type=FALSE,
                                      x.off=4,
                                      mar=c(5,4,4,1)){
  stopifnot(is.data.frame(map),
            "x" %in% colnames(map),
            "y" %in% colnames(map))
  if("rang" %in% colnames(map) & ! "rank" %in% colnames(map))
    colnames(map)[colnames(map) == "rang"] <- "rank"
  if("placette" %in% colnames(map) & ! "location" %in% colnames(map))
    colnames(map)[colnames(map) == "placette"] <- "location"
  stopifnot("rank" %in% colnames(map),
            "location" %in% colnames(map))
  if(show.plant.types)
    stopifnot(! show.pu.type)
  if(show.pu.type)
    stopifnot(! show.plant.types)
  if(show.plant.types || show.pu.type){
    if("t\u00E9moin" %in% colnames(map) & ! "control" %in% colnames(map))
      colnames(map)[colnames(map) == "t\u00E9moin"] <- "control"
    stopifnot("control" %in% colnames(map),
              "parent" %in% colnames(map),
              "F1" %in% colnames(map))
  }
  if(! is.null(col.block)){
    if("bloc" %in% colnames(map) & ! "block" %in% colnames(map))
      colnames(map)[colnames(map) == "bloc"] <- "block"
    stopifnot("block" %in% colnames(map),
              is.vector(col.block),
              "A" %in% names(col.block),
              "B" %in% names(col.block))
  }

  if(! is.null(mar))
    graphics::par(mar=mar)

  range.ranks <- c(min(map$rank), max(map$rank))
  range.locations <- c(min(map$location), max(map$location))
  range.x <- c(min(map$x), max(map$x))
  range.y <- c(min(map$y), max(map$y))
  if(lang == "FR"){
    if(! is.null(main))
      if(main == "default")
        main <- paste0("Carte du croisement r\u00E9ciproque Syrah x Grenache (INRA)",
                       "\nDomaine du Chapitre, Villeneuve-l\u00E8s-Maguelone, France")
    xlab <- "Rangs"
    ylab <- "Placettes"
  } else if(lang == "EN"){
    if(! is.null(main))
      if(main == "default")
        main <- paste0("Map of the reciprocal cross Syrah x Grenache (INRA)",
                       "\nDomaine du Chapitre, Villeneuve-l\u00E8s-Maguelone, France")
    xlab <- "Ranks"
    ylab <- "Locations"
  }
  xlab <- paste0(xlab, " (", length(unique(map$rank)), ")")
  ylab <- paste0(ylab, " (", length(unique(map$location)), ")")
  graphics::plot(0, type="n", axes=FALSE, xlim=c(range.x[1], range.x[2]+x.off), ylim=range.y,
       main=NULL, xlab=xlab, ylab=ylab)
  if(! is.null(main))
    graphics::title(main=main)
  graphics::axis(side=1, at=range.x[1]:range.x[2], labels=sort(unique(map$rank)))
  graphics::axis(side=2, at=seq(range.y[1], range.y[2]+1, by=5), labels=FALSE, tick=TRUE)
  graphics::axis(side=2, at=seq(range.y[1], range.y[2], by=5)+2.5, las=1,
       labels=sort(unique(map$location)), tick=FALSE)

  if(! is.null(col.block)){
    off <- 0.5
    graphics::rect(xleft=min(map$x[map$block == "A"]) - off,
         ybottom=min(map$y[map$block == "A"]) - off,
         xright=max(map$x[map$block == "A"]) + off,
         ytop=max(map$y[map$block == "A"]) + off,
         col=col.block["A"], border=NA)
    graphics::rect(xleft=min(map$x[map$block == "B"]) - off,
         ybottom=min(map$y[map$block == "B"]) - off,
         xright=max(map$x[map$block == "B"]) + off,
         ytop=max(map$y[map$block == "B"]) + off,
         col=col.block["B"], border=NA)
    if(lang == "FR"){
      graphics::legend("right", legend=paste("bloc", names(col.block)),
             col=col.block, fill=col.block, border=col.block,
             pch=NA, bty="n")
    } else if(lang == "EN")
      graphics::legend("right", legend=paste("block", names(col.block)),
               col=col.block, fill=col.block, border=col.block,
               pch=NA, bty="n")
  }

  symbols <- data.frame(pch=c(4, 2, 1),
                        col=c("blue", "darkgreen", "black"),
                        stringsAsFactors=FALSE)
  rownames(symbols) <- c("control", "parent", "F1")
  if(lang == "FR"){
    symbols$leg <- c("t\u00E9moin", "parent", "F1")
  } else if(lang == "EN"){
    symbols$leg <- rownames(symbols)
  }
  which.y <- "y"
  if(show.pu.type){
    if(! "pu" %in% colnames(map))
      map$pu <- paste(map$id, map$rank, map$location, sep="-")
    if(! "y.pu" %in% colnames(map)){
      map$y.pu <- NA
      for(pu in unique(map$pu)){
        idx <- which(map$pu == pu)
        y <- map$y[idx]
        map$y.pu[idx] <- rep(stats::median(y-0.5), length(y))
      }
    }
    which.y <- "y.pu"
  }
  if(show.plant.types || show.pu.type){
    graphics::points(x=map$x[map$control],
           y=map[[which.y]][map$control],
           col=symbols["control","col"],
           pch=symbols["control","pch"], cex=0.8)
    graphics::points(x=map$x[map$parent],
           y=map[[which.y]][map$parent],
           col=symbols["parent","col"],
           pch=symbols["parent","pch"], cex=0.8)
    graphics::points(x=map$x[!map$control & !map$parent],
           y=map[[which.y]][!map$control & !map$parent],
           col=symbols["F1","col"],
           pch=symbols["F1","pch"], cex=0.2)
    graphics::legend("topright", legend=symbols$leg, col=symbols$col,
           pch=symbols$pch, bty="n", pt.cex=c(0.8,0.8,0.2))
  }

  return(invisible(map))
}

##' SSR genotypes from PlantGrape
##'
##' Extract SSR genotypes from a PDF produced by \href{http://plantgrape.plantnet-project.org/fr/}{PlantGrape} for a given grape variety.
##' @param file path to the PDF file
##' @param var.name name of the grape variety
##' @param sep character to use to separate both alleles of a given SSR
##' @return matrix with one row and as many columns as SSRs, usable by "df2genind" from the "adegenet" package
##' @author Timothee Flutre
##' @export
getSsrGenosFromPlantGrapePdf <- function(file, var.name, sep="/"){
  requireNamespace("pdftools")
  stopifnot(file.exists(file))

  text <- pdftools::pdf_text(file)

  text <- do.call(c, strsplit(text, "\n"))
  idx <- (grep("Profil G\u00E9n\u00E9tique", text)+1):(grep("Ph\u00E9nologie", text)-1)
  stopifnot(length(idx) == 3)
  text <- text[idx]

  ssr.names <- strsplit(text[1], split=" ")[[1]]
  ssr.names <- ssr.names[! ssr.names == ""][-1]
  stopifnot(length(ssr.names) == 9)

  all1 <- strsplit(text[2], split=" ")[[1]]
  all1 <- as.numeric(all1[! all1 == ""][-c(1,2)])
  stopifnot(length(all1) == 9)

  all2 <- strsplit(text[3], split=" ")[[1]]
  all2 <- as.numeric(all2[! all2 == ""][-c(1,2)])
  stopifnot(length(all2) == 9)

  out <- matrix(data=paste0(all1, sep, all2),
                nrow=1, ncol=9,
                dimnames=list(var.name, ssr.names))
  return(out)
}
