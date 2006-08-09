showDens = function(z, breaks, xat, xtickLabels=paste(xat), col, ylab="",  ...) {

  y  = matrix(NA, nrow=length(breaks)-1, ncol=length(z))
  scaleFac = numeric(length(z))
  for(k in seq(along=z)) {
    h = hist(z[[k]], breaks=breaks, plot=FALSE)
    scaleFac[k] = 1/max(h$counts) * 0.92 
    y[,k] = h$counts * scaleFac[k]
    if(k==1) {
      mids = h$mids
    } else {
      stopifnot(identical(mids, h$mids))
    }
  }

  plot(breaks[c(1, length(breaks))], c(0, length(z)), ylab=ylab, type="n",
       yaxt="n", xaxt="n", ...)
  axis(1, at=xat, labels=xtickLabels)
       
  for(k in seq(along=z)) {
    poy = (length(z)-k) + c(0, rep(y[,k], each=2), 0, 0)
    pox = breaks[c(rep(seq(along=breaks), each=2), length(breaks))]
    polygon(pox, poy, col=col[k])
  }

  return(scaleFac)
}

