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
plotMapSxGDomaineChapitre <- function(map,
                                      main="default",
                                      lang="FR",
                                      col.block=setNames(
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
    par(mar=mar)

  range.ranks <- c(min(map$rank), max(map$rank))
  range.locations <- c(min(map$location), max(map$location))
  range.x <- c(min(map$x), max(map$x))
  range.y <- c(min(map$y), max(map$y))
  if(lang == "FR"){
    if(! is.null(main))
      if(main == "default")
        main <- paste0("Carte du croisement r\u00E9ciproque Syrah x Grenache (INRA)",
                       "\nDomaine du Chapitre, Villeneuve-l\u00E8s-Maguelone, France")
    xlab <- "rangs"
    ylab <- "placette"
  } else if(lang == "EN"){
    if(! is.null(main))
      if(main == "default")
        main <- paste0("Map of the reciprocal cross Syrah x Grenache (INRA)",
                       "\nDomaine du Chapitre, Villeneuve-l\u00E8s-Maguelone, France")
    xlab <- "ranks"
    ylab <- "locations"
  }
  plot(0, type="n", axes=FALSE, xlim=c(range.x[1], range.x[2]+x.off), ylim=range.y,
       main=NULL, xlab=xlab, ylab=ylab)
  if(! is.null(main))
    title(main=main)
  axis(side=1, at=range.x[1]:range.x[2], labels=sort(unique(map$rank)))
  axis(side=2, at=seq(range.y[1], range.y[2]+1, by=5), labels=FALSE, tick=TRUE)
  axis(side=2, at=seq(range.y[1], range.y[2], by=5)+2.5, las=1,
       labels=sort(unique(map$location)), tick=FALSE)

  if(! is.null(col.block)){
    off <- 0.5
    rect(xleft=min(map$x[map$block == "A"]) - off,
         ybottom=min(map$y[map$block == "A"]) - off,
         xright=max(map$x[map$block == "A"]) + off,
         ytop=max(map$y[map$block == "A"]) + off,
         col=col.block["A"], border=NA)
    rect(xleft=min(map$x[map$block == "B"]) - off,
         ybottom=min(map$y[map$block == "B"]) - off,
         xright=max(map$x[map$block == "B"]) + off,
         ytop=max(map$y[map$block == "B"]) + off,
         col=col.block["B"], border=NA)
    if(lang == "FR"){
      legend("right", legend=paste("bloc", names(col.block)),
             col=col.block, fill=col.block, border=col.block,
             pch=NA, bty="n")
    } else if(lang == "EN")
        legend("right", legend=paste("block", names(col.block)),
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
        map$y.pu[idx] <- rep(median(y-0.5), length(y))
      }
    }
    which.y <- "y.pu"
  }
  if(show.plant.types || show.pu.type){
    points(x=map$x[map$control],
           y=map[[which.y]][map$control],
           col=symbols["control","col"],
           pch=symbols["control","pch"], cex=0.8)
    points(x=map$x[map$parent],
           y=map[[which.y]][map$parent],
           col=symbols["parent","col"],
           pch=symbols["parent","pch"], cex=0.8)
    points(x=map$x[!map$control & !map$parent],
           y=map[[which.y]][!map$control & !map$parent],
           col=symbols["F1","col"],
           pch=symbols["F1","pch"], cex=0.2)
    legend("topright", legend=symbols$leg, col=symbols$col,
           pch=symbols$pch, bty="n", pt.cex=c(0.8,0.8,0.2))
  }

  return(invisible(map))
}
