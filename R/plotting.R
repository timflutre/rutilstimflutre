## Contains functions useful for plotting.

##' Make a scatter plot of y as a function of x, along with the regression
##' line from lm() as well as both confidence and both prediction lines from
##' predict().
##'
##'
##' @param x vector of points
##' @param y vector of points
##' @param ... arguments to be passed to plot()
##' @return nothing
##' @author Timothee Flutre
##' @export
regplot <- function(x, y, ...){
  x <- as.numeric(x)
  y <- as.numeric(y)
  graphics::plot(x, y, ...)
  fit <- stats::lm(y ~ x)
  graphics::abline(fit, col="red")
  newx <- seq(min(x), max(x), length.out=length(x))
  pred.ci <- stats::predict(fit, newdata=data.frame(x=newx), interval="confidence")
  graphics::lines(newx, pred.ci[,"lwr"], lty=2)
  graphics::lines(newx, pred.ci[,"upr"], lty=2)
  pred.pi <- stats::predict(fit, newdata=data.frame(x=newx), interval="prediction")
  graphics::lines(newx, pred.pi[,"lwr"], lty=3)
  graphics::lines(newx, pred.pi[,"upr"], lty=3)
}

##' Plot a Hinton diagram
##'
##' Modified from http://www.cs.princeton.edu/~mimno/factor-analysis.R
##' @title Hinton diagram
##' @param m matrix
##' @param main main title
##' @param max.sqrt.m to play with the scaling
##' @author Timothee Flutre
##' @export
hinton <- function(m, main="", max.sqrt.m=NULL){
  rows <- dim(m)[1]
  cols <- dim(m)[2]

  left <- rep(0, rows * cols)
  right <- rep(0, rows * cols)
  bottom <- rep(0, rows * cols)
  top <- rep(0, rows * cols)

  box.colors <- rep("white", rows * cols)

  if(is.null(max.sqrt.m))
    max.sqrt.m <- max(sqrt(abs(m)))
  scale <- 0.9 / (2 * max.sqrt.m)

  position <- 1

  for(row in 1:rows){
    for(col in 1:cols){
      if(m[row,col] < 0)
        box.colors[position] <- "black"
      x <- sqrt(abs(m[row,col]))
      left[position] <- col - (x * scale)
      right[position] <- col + (x * scale)
      top[position] <- -(row - (x * scale))
      bottom[position] <- -(row + (x * scale))
      position <- position + 1
    }
  }

  xlab <- ""
  ylab <- ""
  if(! is.null(names(dimnames(m)))){
    xlab <- names(dimnames(m))[2]
    ylab <- names(dimnames(m))[1]
  }

  opar <- graphics::par(mar=c(ifelse(xlab == "", 1, 3),
                    ifelse(ylab == "", 1, 3), 5, 1) + 0.1)

  graphics::plot(0, xlim=c(0.25,cols+0.75), ylim=c(-rows-0.75, -0.25),
       type="n", xaxt="n", yaxt="n", xlab="", ylab="")

  graphics::rect(left, bottom, right, top, col=box.colors)

  if(main != "")
    graphics::title(main=main, line=3)

  if(! is.null(colnames(m))){
    graphics::axis(side=3, at=1:cols, labels=FALSE)
    graphics::text(x=1:cols, y=graphics::par("usr")[4] + 0.25, labels=colnames(m), adj=0, srt=45, xpd=TRUE)
  } else
    graphics::axis(side=3, at=1:cols, labels=1:cols)

  if(! is.null(rownames(m))){
    graphics::axis(side=2, at=-(1:rows), labels=rownames(m), las=1)
  } else
    graphics::axis(side=2, at=-(1:rows), labels=1:rows, las=1)

  if(xlab != "")
    graphics::mtext(xlab, side=1, line=1)
  if(ylab != "")
    graphics::mtext(ylab, side=2, line=2)

  on.exit(graphics::par(opar))
}

##' Plot a scale, e.g. to add on the side of image()
##'
##' Takes some time to draw (there is one polygon per break...)
##' http://menugget.blogspot.de/2011/08/adding-scale-to-image-plot.html
##' @param z vector
##' @param zlim lim
##' @param col color
##' @param breaks vector
##' @param horiz boolean
##' @param ylim lim
##' @param xlim lim
##' @param ... arguments to be passed to plot()
##' @author Timothee Flutre
##' @export
plotWithScale <- function(z, zlim, col = grDevices::heat.colors(12),
                          breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
  if(! missing(breaks))
    if(length(breaks) != (length(col)+1))
      stop("must have one more break than colour")

  if(missing(breaks) & ! missing(zlim))
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))

  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2] + c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1] - c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }

  poly <- vector(mode="list", length(col))
  for(i in seq(poly))
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])

  xaxt <- ifelse(horiz, "s", "n")
  yaxt <- ifelse(horiz, "n", "s")
  if(horiz){
    YLIM <- c(0,1)
    XLIM <- range(breaks)
  } else{
    YLIM <- range(breaks)
    XLIM <- c(0,1)
  }
  if(missing(xlim))
    xlim <- XLIM
  if(missing(ylim))
    ylim <- YLIM

  graphics::plot(1, 1, t="n", ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt,
       xaxs="i", yaxs="i", bty="n", ...)

  for(i in seq(poly)){
    if(horiz){
      graphics::polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    } else
      graphics::polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
  }
}

##' Plot a matrix as a heatmap in its natural orientation, with a colored
##' scale on the right side, and optionally using its dimension names for
##' rows and columns
##'
##' To print all row names, choose idx.rownames=1:nrow(z). To print a subset
##' of 10 row names, choose idx.rownames=floor(seq(1, nrow(z), length.out=10)).
##' Similarly for column names.
##' @param z matrix to be plotted
##' @param main title to appear above the heatmap
##' @param idx.rownames vector giving the indices of the row names of z to be added on the left side of the plot
##' @param idx.colnames vector giving the indices of the column names of z to be added on top of the plot
##' @param breaks vector (default=seq(min(z), max(z), length.out=100))
##' @param left.text.at vector which names and values will be used to label the left side of the plot; if not NULL, takes precedence over idx.rownames
##' @author Timothee Flutre
##' @examples
##' \dontrun{set.seed(1859)
##' genomes <- simulCoalescent(nb.inds=200, nb.pops=3, mig.rate=3)
##' X <- genomes$genos
##' A <- estimGenRel(X=X, relationships="additive", method="vanraden1")
##' imageWithScale(z=A, main="Additive genetic relationships", breaks=seq(0,1,length.out=20),
##'                left.text.at=setNames(c(0.83, 0.5, 0.17), c("pop1", "pop2", "pop3")))
##' }
##' @export
imageWithScale <- function(z, main=NULL, idx.rownames=NULL, idx.colnames=NULL,
                           breaks=NULL, left.text.at=NULL){
  stopifnot(is.matrix(z))
  if(! is.null(left.text.at))
    stopifnot(is.null(idx.rownames))
  if(! is.null(idx.rownames) & is.null(rownames(z)))
    stop("non-null idx.rownames requires z to have row names")
  if(! is.null(idx.colnames) & is.null(colnames(z)))
    stop("non-null idx.colnames requires z to have column names")
  if(is.null(breaks))
    breaks <- seq(min(z), max(z), length.out=100)

  graphics::layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(7,1))
  ## layout.show(2) # for debugging purposes

  col.pal <- grDevices::colorRampPalette(c("black", "red", "yellow"), space="rgb")

  ## plot the heatmap
  custom.mar <- c(1, 5, 6, 1)
  if(is.null(idx.rownames) & is.null(left.text.at))
      custom.mar[2] <- 1
  if(is.null(idx.colnames))
      custom.mar[3] <- 3
  opar <- graphics::par(mar=custom.mar)
  graphics::image(t(z)[,nrow(z):1], axes=FALSE, col=col.pal(length(breaks)-1))
  if(! is.null(main))
    graphics::mtext(text=main, side=3, line=ifelse(is.null(idx.colnames), 1, 4),
                    font=2, cex=1.3)
  if(! is.null(idx.colnames))
    graphics::text(x=seq(0,1,length.out=length(idx.colnames)), y=graphics::par("usr")[4]+0.02,
                   srt=45, adj=0, labels=colnames(z)[idx.colnames], xpd=TRUE)
  if(! is.null(idx.rownames))
    graphics::mtext(text=rev(rownames(z)[idx.rownames]), side=2, line=1,
                    at=seq(0,1,length.out=length(idx.rownames)),
                    las=2)
  if(! is.null(left.text.at))
    graphics::mtext(text=names(left.text.at), side=2, line=1,
                    at=left.text.at, las=2)
  on.exit(graphics::par(opar))

  ## plot the scale
  opar <- graphics::par(mar=c(1,0,6,3))
  plotWithScale(z, col=col.pal(length(breaks)-1), breaks=breaks, horiz=FALSE,
                yaxt="n")
  graphics::axis(4, at=format(breaks[seq.int(1, length(breaks), length.out=5)],
                              digits=2),
       las=2, lwd=0, lwd.ticks=1)
  on.exit(graphics::par(opar))
}

##' Principal component analysis
##'
##' Plot the two first principal components from a PCA.
##' @param rotation rotated matrix which columns corresponds to "principal components"
##' @param prop.vars vector with the proportion of variance explained per PC
##' @param plot use "points" to show a plot with \code{\link[graphics]{points}} of PC1 versus PC2, and "text" to use \code{\link[graphics]{text}} with row names of \code{rotation} as labels
##' @param main main title of the plot
##' @param cols N-vector of colors
##' @return nothing
##' @author Timothee Flutre
##' @seealso \code{\link{pca}}
##' @export
plotPca <- function(rotation, prop.vars, plot="points", main="PCA",
                    cols=rep("black", nrow(rotation))){
  stopifnot(is.matrix(rotation),
            is.vector(prop.vars),
            is.numeric(prop.vars),
            all(prop.vars >= 0),
            all(prop.vars <= 1),
            plot %in% c("points", "text"),
            is.vector(cols),
            length(cols) == nrow(rotation))

  graphics::plot(x=rotation[,1], y=rotation[,2], las=1,
                 xlab=paste0("PC1 (", format(100 * prop.vars[1], digits=3), "%)"),
                 ylab=paste0("PC2 (", format(100 * prop.vars[2], digits=3), "%)"),
                 main=main, type="n")
  graphics::abline(h=0, lty=2)
  graphics::abline(v=0, lty=2)

  if(plot == "points"){
    graphics::points(x=rotation[,1], y=rotation[,2], col=cols, pch=20)
  } else if(plot == "text")
    graphics::text(x=rotation[,1], y=rotation[,2], labels=rownames(rotation),
                   col=cols, pch=20)
}
