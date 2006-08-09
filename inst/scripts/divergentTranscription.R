##
## Explore divergent transcription
## (c) W. Huber 2006
##
## the file categorizedSegments.rda is produced by makecs.rda
## and it contains the data objects "cs"
if(!exists("cs"))
  load("categorizedSegments.rda")

options(error=recover)
rt = "seg-polyA-050909"
segScore = cs[[rt]] 

## jl - index of segment to the left, jr - to the right
jl= c(NA, 1:(nrow(segScore)-1))
jr= c(2:nrow(segScore), NA)
jl[ !((segScore$chr==segScore$chr[jl]) & (segScore$strand==segScore$strand[jl])) ] = NA
jr[ !((segScore$chr==segScore$chr[jr]) & (segScore$strand==segScore$strand[jr])) ] = NA

stopifnot(all(segScore$end <= segScore$start[jr]+17, na.rm=TRUE),
          all(segScore$start+17 >= segScore$end[jl], na.rm=TRUE))

## selWS: segments on W strand whose left boundary appears to be a transcription start site
## selWE: segments on W strand whose right boundary appears to be a transcription stop site
## selCS: segments on C strand whose right boundary appears to be a transcription start site
## selCE: segments on C strand whose left boundary appears to be a transcription stop site
minLev = 1
sel   = ((segScore$level>minLev) & (segScore$frac.dup<0.5) & (segScore$chr<=16))
selWS = (sel & (segScore$level[jl]<0) & (segScore$strand=="+"))
selWE = (sel & (segScore$level[jr]<0) & (segScore$strand=="+"))
selCS = (sel & (segScore$level[jr]<0) & (segScore$strand=="-"))
selCE = (sel & (segScore$level[jl]<0) & (segScore$strand=="-"))

## for each segment, find its 5' and 3' end
stopifnot(all(segScore$strand%in%c("+","-")))
p5  = ifelse(segScore$strand=="+", segScore$start, segScore$end) + 1e8*segScore$chr
p3  = ifelse(segScore$strand=="+", segScore$end, segScore$start) + 1e8*segScore$chr

## For each segment in 'wh', find the nearest neighbours to its 5' end,
##   on same and opposite strands, looking in the set of segments indicated
##   by 'whsame' and 'whoppo'. Segments are identified by their integer index
##   in segScore.
## The return value is a matrix whose rows correspond to 'wh',
##   1st column is the distance to the closest 5' end of a segment from 'whoppo'
##   2nd column is that segment's index
##   3rd column is the distance to the closest 3' end of a segment from 'whsame'
##   4th column is that segment's index
##
nearestNeighbourDists = function(wh, whsame, whoppo, offset=0) {
  
  d = matrix(as.integer(NA), nrow=length(wh), ncol=5)
  colnames(d) = c("which", "dminoppo", "whichminoppo", "dminsame", "whichminsame")
  
  for(i in seq(along=wh)){
    ## 5' end of segment on opposite strand (candidate for divergent transcription)
    dop = p5[wh[i]] - (p5[whoppo] + offset)
    ## 3' end of segment on same strand
    dsm = p5[wh[i]] - (p3[whsame] + offset)

    dop[ (dop > 1e7) | (dop < -25) ] = +Inf
    dsm[ (dsm > 1e7) | (dsm < -25) ] = +Inf
    
    wop = which.min(dop)
    wsm = which.min(dsm)

    d[i, "which"] = wh[i]
    
    if(length(wop)==1) {
      d[i, "dminoppo"]     = dop[wop]
      d[i, "whichminoppo"] = whoppo[wop]
    } else {
      stopifnot(length(wop)==0)
    }
    if(length(wsm)==1) {
      d[i, "dminsame"]     = dsm[wsm]
      d[i, "whichminsame"] = whsame[wsm]
    } else {
      stopifnot(length(wsm)==0)
    }
  }
  return(d)
}

## write a table with the best 30
##
dW = nearestNeighbourDists(which(selWS), which(selWE), which(selCS))
dC = nearestNeighbourDists(which(selCS), which(selCE), which(selWS))
if(TRUE){
  source(scriptsDir("writeSegmentTable.R"))
  if(!exists("gff"))load("probeAnno.rda")

  myTable = function(d, ...) {
    ord = order(d[, "dminoppo"])[1:30]
    out = cbind("distance" = d[ord, "dminoppo"],
      segScore[d[ord, "which"], c("overlappingFeature", "chr", "strand", "start", "end", "length", "level","frac.dup")])
    
    writeSegmentTable(out, colOrd=colnames(out), 
                      fn = file.path(rt, "viz", "divergentTranscription"), HTML=TRUE, 
                      ..., interact=TRUE)
  }
  
  myTable(dW, title = paste("div. transcription (", rt, ") Watson strand", sep=""), open="wt")
  myTable(dC, title = paste("div. transcription (", rt, ") Crick strand", sep=""),  open="at")
    
}

## are there many short segments?
##

makeLengths = function(d) {
   ord = order(d[, "dminoppo"])[1:50]
   cbind(segScore$length[ d[ord, "which"] ],
         segScore$length[ d[ord, "whichminoppo"] ])
 }

lengths = makeLengths(dW)
## lengths = makeLengths(dC)

## minl = 200
## tab = table(lengths[,1]>minl, lengths[,2]>minl)
## print(tab)
## print(fisher.test(tab))

pdf(file="divergentTranscription-lengths.pdf", width=12, height=4.5)
layout(cbind(c(1,1),c(2,3),c(4,5)), widths=1, heights=1)
breaks = seq(0, max(segScore$length, na.rm=TRUE), length=21)

plot(lengths, xlab="length (WS)", ylab="length (CS)", pch=21)
hist(lengths[,1], breaks=breaks, col="orange", main="length (divergent WS)")
hist(segScore$length[which(selWS)], breaks=breaks, col="orange", main="length (all WS)")
hist(lengths[,2], breaks=breaks, col="lightblue", main="length (divergent CS)")
hist(segScore$length[which(selCS)], breaks=breaks, col="lightblue", main="length (all CS)")
dev.off()


## is the nearest 5' ORF to a filtered antisense is on the same
## or opposite strand.?
isAS   = (segScore$category %in% c("novel antisense - filtered", "novel antisense - unassigned")[1])
isORF  = (segScore$simpleCatg=="annotated ORF")
dW = nearestNeighbourDists(which(selWS & isAS), which(selWE), which(selCS))
dC = nearestNeighbourDists(which(selCS & isAS), which(selCE), which(selWS))

pdf(file="divergentTranscription-sameOrOppo.pdf", width=8, height=4.5)
par(mfrow=c(1,2))
myPlot = function(d) {
  plot(d[, c("dminoppo", "dminsame")], pch=16, xlab="opposite strand", ylab="same strand", xlim=c(0,1.6e4), ylim=c(0,1.6e4))
  abline(a=0,b=1,col="blue")
}

myPlot(dW)
myPlot(dC)

dev.off()


stop()



## Statistical significance of the minimal distances:
##   by shifting the C strand back and forth
if(!exists("vn")) {
  ## offsets = 150*seq(-8, 6, by=1)
  offsets = 1000*seq(-10, 10, by=1)
  vn = vector(mode="list", length=length(offsets))
  for(ioff in seq(along=offsets)) {
    cat(ioff,"")
    vn[[ioff]] = nearestNeighbourDists(offsets[ioff])
  }
}

if(!TRUE){
  library("geneplotter")
  dists = lapply(vn, function(x) x[, "dmin"])
  each  = unique(listLen(dists))
  stopifnot(length(each)==1)
  dat  = cbind(x=unlist(dists), y=rep(offsets, each=each))
  dat  = dat[!is.na(dat[,1]) & (dat[,1]<1000), ]
  smoothScatter(dat, bandwidth=c(10, 100))
}

breaks  = seq(-400, 1200, by=25)
mind    = breaks[1]
maxd    = breaks[length(breaks)]
dh = sapply(vn, function(x) {
  px = x[, "dmin"]
  px = px[px<maxd & px>mind]
  h = hist(px, plot=FALSE, breaks=breaks)
  rv = h$counts
  names(rv) = paste(h$mids)
  return(rv)
})
colnames(dh) = paste(offsets)
  
## to do: maybe we can use a cool density estimate instead of histogram
## image(x=breaks, y=offsets, z=h, xlab="distance (bp)", col=gray(seq(1,0,by=-0.01)))

library("grid")
histoPlot = function(h, breaks) {
  pushViewport(viewport(layout=grid.layout(ncol(h)+1, 1, height=1)))
  ymax = max(h)*1.2
  xmin = -4
  for(i in 1:ncol(h)){
    pushViewport(dataViewport(xscale=c(xmin, nrow(h)+1),
                              yscale=c(0, ymax), 
                              ##yscale=c(0, max(h[,i])), 
                              layout.pos.col=1, layout.pos.row=i))
    
    grid.lines(c(0, rep(1:(nrow(h)-1), each=2), nrow(h)),
               rep(h[,i], each=2), default.units="native")
    grid.text(colnames(h)[i], x=unit(0, "native"), just=c("right", "centre"))
    popViewport()
  }
  stopifnot(length(breaks)==nrow(h)+1)
  pushViewport(dataViewport(xscale=c(xmin, nrow(h)+1), yscale=c(0,1),
                            layout.pos.col=1, layout.pos.row=ncol(h)+1))
  px  = unit(0:nrow(h), "native")
  grid.segments(x0=px, x1=px, y0=unit(0.7, "npc"), y1=unit(0.9, "npc"))
  sel = (0:nrow(h))%%4==0
  grid.text(paste(breaks[sel]), x=px[sel], y=unit(0, "npc"), just=c("centre", "bottom"))
  popViewport(2)  
}

pdf(file="divergentTranscription.pdf", width=9, height=9)
grid.newpage()
histoPlot(dh, breaks)
dev.off()


