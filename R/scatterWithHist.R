scatterWithHist = function(x, breaks, barcols, xlab, ylab, ...) {
  stopifnot(is.matrix(x), ncol(x)==2)
  xmax = breaks[length(breaks)]
  xmid = breaks[length(breaks)/2]
  x[x>xmax]=NA
  xhist = hist(x[,1], breaks=breaks, plot=FALSE)
  yhist = hist(x[,2], breaks=breaks, plot=FALSE)
  topx  = max(xhist$counts)
  topy  = max(yhist$counts)
  top   = max(topx, topy)
  xrange = yrange = breaks[c(1, length(breaks))]
  nf = layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)
 
  par(mar=c(3,3,1,1))
  plot(x, xlim=xrange, ylim=yrange, xlab="", ylab="", ...)
  par(mar=c(0,3,1,1))
  barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0, col=barcols[1])
  text(length(xhist$counts)/2, topx, adj=c(0.5, 1), labels=xlab, cex=1)
  par(mar=c(3,0,1,1))
  barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, col=barcols[2], horiz=TRUE)
  text(topy, length(yhist$counts)/2, adj=c(0.5, 1), labels=ylab, cex=1, srt=270)
} 
