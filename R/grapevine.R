## Contains functions useful for grapevine projects.

##' Plot the map of SxG
##'
##' reciprocal cross Syrah x Grenache at the Domaine du Chapitre (near Montpellier, France)
##' @param map data.frame with columns rank, location, x and y
##' @param lang language for the title and legend (default=FR/EN)
##' @param show.plant.types show the plants corresponding to F1, parents and controls
##' @param col.block vector with block colors (set to NULL to ignore)
##' @return nothing
##' @author Timothee Flutre
plotMapSxGDomaineChapitre <- function(map,
                                      lang="FR",
                                      show.plant.types=TRUE,
                                      col.block=setNames(
                                          c("grey81", "bisque1"),
                                          c("A","B"))){
  stopifnot(is.data.frame(map),
            "x" %in% colnames(map),
            "y" %in% colnames(map))
  if("rang" %in% colnames(map) & ! "rank" %in% colnames(map))
    colnames(map)[colnames(map) == "rang"] <- "rank"
  if("placette" %in% colnames(map) & ! "location" %in% colnames(map))
    colnames(map)[colnames(map) == "placette"] <- "location"
  stopifnot("rank" %in% colnames(map),
            "location" %in% colnames(map))
  if(show.plant.types){
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

  par(mar=c(5,4,4,1))

  range.ranks <- c(min(map$rank), max(map$rank))
  range.locations <- c(min(map$location), max(map$location))
  range.x <- c(min(map$x), max(map$x))
  range.y <- c(min(map$y), max(map$y))
  if(lang == "FR"){
    main <- paste0("Carte du croisement r\u00E9ciproque Syrah x Grenache (INRA)",
                   "\nDomaine du Chapitre, Villeneuve-l\u00E8s-Maguelone, France")
    xlab <- "rangs"
    ylab <- "placette"
  } else if(lang == "EN"){
      main <- paste0("Map of the reciprocal cross Syrah x Grenache (INRA)",
                     "\nDomaine du Chapitre, Villeneuve-l\u00E8s-Maguelone, France")
      xlab <- "ranks"
      ylab <- "locations"
    }
  plot(0, type="n", axes=FALSE, xlim=c(range.x[1], range.x[2]+4), ylim=range.y,
       main=main, xlab=xlab, ylab=ylab)
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

  if(show.plant.types){
    symbols <- data.frame(pch=c(4, 2, 1),
                          col=c("blue", "darkgreen", "black"),
                          stringsAsFactors=FALSE)
    rownames(symbols) <- c("control", "parent", "F1")
    if(lang == "FR"){
      symbols$leg <- c("t\u00E9moin", "parent", "F1")
    } else if(lang == "EN"){
      symbols$leg <- rownames(symbols)
    }
    points(x=map$x[map$control], y=map$y[map$control],
           col=symbols["control","col"], pch=symbols["control","pch"], cex=0.8)
    points(x=map$x[map$parent], y=map$y[map$parent], col=symbols["parent","col"],
           pch=symbols["parent","pch"], cex=0.8)
    points(x=map$x[!map$control & !map$parent], y=map$y[!map$control & !map$parent],
           col=symbols["F1","col"], pch=symbols["F1","pch"], cex=0.2)
    legend("topright", legend=symbols$leg, col=symbols$col,
           pch=symbols$pch, bty="n", pt.cex=c(0.8,0.8,0.2))
  }
}
