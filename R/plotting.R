## Contains functions useful for plotting.

##' Scatter plot with regression lines
##'
##' Make a scatter plot of y as a function of x, along with regression line(s).
##' @param x vector or 1-column matrix (missing data coded as NA will be automatically discarded)
##' @param y vector or 1-column matrix (missing data coded as NA will be automatically discarded)
##' @param reg specifies which model(s) to use to add regression lines(s) (lm/loess/rlm; can be \code{c("lm", "rlm")})
##' @param col named vector specifying the color(s) of the regression line(s) specified via \code{reg}
##' @param show.cor if TRUE, Pearson and Spearman correlation coefficients are shown in the top left corner
##' @param ci.int if \code{reg="lm"}, add lines corresponding to confidence intervals
##' @param pred.int if \code{reg="lm"}, add lines corresponding to prediction intervals
##' @param legend.x the x co-ordinate to be used to position the legend (see \code{\link[graphics]{legend}})
##' @param legend.y the y co-ordinate to be used to position the legend (see \code{\link[graphics]{legend}})
##' @param ... arguments to be passed to \code{\link[graphics]{plot}}
##' @return list of object(s) returned by the function(s) specified via \code{reg}
##' @author Timothee Flutre
##' @examples
##' set.seed(1859)
##' n <- 500
##' x <- rnorm(n=n, mean=37, sd=3)
##' y <- 50 + 1.2 * x + rnorm(n=n, mean=0, sd=3)
##' fit <- regplot(x=x, y=y, reg="lm", las=1, main="Linear regression")
##' fit <- regplot(x=x, y=y, reg="loess", las=1, col=c(loess="red"),
##'                main="Locally weighted scatterplot smoothing (loess)")
##' y2 <- y + sample(x=c(rep(0, 0.7*floor(n)),
##'                      rnorm(n=ceiling(0.3*n), mean=c(7,13), sd=20)), size=n)
##' fit <- regplot(x=x, y=y2, reg=c("lm","rlm"), las=1,
##'                col=c(lm="red", rlm="blue"), legend.x="bottomright",
##'                main="(Robust) linear regressions")
##' @export
regplot <- function(x, y, reg="lm", col=c(lm="red"), show.cor=TRUE,
                    ci.int=TRUE, pred.int=TRUE,
                    legend.x="right", legend.y=NULL,
                    ...){
  stopifnot(is.numeric(x) || is.vector(x) || (is.matrix(x) && ncol(x) == 1),
            is.numeric(y) || is.vector(y) || (is.matrix(y) && ncol(y) == 1),
            is.character(reg),
            all(reg %in% c("lm", "loess", "rlm")),
            all(names(col) %in% c("lm", "loess", "rlm")),
            length(col) == length(reg),
            is.logical(ci.int),
            is.logical(pred.int))
  if("rlm" %in% reg){
    requireNamespace("MASS")
  }

  fit <- list()

  x <- as.numeric(x)
  y <- as.numeric(y)
  stopifnot(length(x) == length(y))
  tmp <- data.frame(x=x, y=y)
  tmp <- tmp[stats::complete.cases(tmp),]

  graphics::plot(formula=y ~ x, data=tmp, ...)
  legd <- c()
  col.order <- c()

  if("lm" %in% reg){
    fit$lm <- stats::lm(formula=y ~ x, data=tmp)
    graphics::abline(fit$lm, col=col["lm"])
    legd <- append(legd, "lm")
    col.order <- append(col.order, col["lm"])

    newx <- seq(min(tmp$x), max(tmp$x), length.out=length(tmp$x))
    if(ci.int){
      pred.ci <- stats::predict(fit$lm, newdata=data.frame(x=newx),
                                interval="confidence")
      graphics::lines(newx, pred.ci[,"lwr"], lty=2)
      graphics::lines(newx, pred.ci[,"upr"], lty=2)
    }
    if(pred.int){
      pred.pi <- stats::predict(fit$lm, newdata=data.frame(x=newx),
                                interval="prediction")
      graphics::lines(newx, pred.pi[,"lwr"], lty=3)
      graphics::lines(newx, pred.pi[,"upr"], lty=3)
    }
  }
  if("loess" %in% reg){
    fit$loess <- stats::loess(formula=y ~ x, data=tmp)
    tmp$fitted.loess <- fit$loess$fitted
    graphics::lines(x=tmp$x[order(tmp$x)],
                    y=tmp$fitted.loess[order(tmp$x)],
                    col=col["loess"])
    legd <- append(legd, "loess")
    col.order <- append(col.order, col["loess"])
  }
  if("rlm" %in% reg){
    fit$rlm <- MASS::rlm(formula=y ~ x, data=tmp)
    tmp$fitted.rlm <- stats::fitted(fit$rlm)
    graphics::lines(x=tmp$x[order(tmp$x)],
                    y=tmp$fitted.rlm[order(tmp$x)],
                    col=col["rlm"])
    legd <- append(legd, "rlm")
    col.order <- append(col.order, col["rlm"])
  }

  graphics::legend(x=legend.x, y=legend.y, legend=legd, col=col.order, lty=1,
                   bty="n")

  if(show.cor){
    fit$cor.p <- stats::cor(x=x, y=y, method="pearson")
    fit$cor.s <- stats::cor(x=x, y=y, method="spearman")
    graphics::text(x=graphics::par("usr")[1] +
                     0.03 * abs(graphics::par("usr")[2] -
                                graphics::par("usr")[1]),
                   y=graphics::par("usr")[4] -
                     0.03 * abs(graphics::par("usr")[4] -
                                graphics::par("usr")[3]),
                   adj=c(0, 1),
                   labels=paste0("cor.Pearson=", format(fit$cor.p, digits=2),
                                 "\n",
                                 "cor.Spearman=", format(fit$cor.s, digits=2)))
  }

  invisible(fit)
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
