## This script produces Fig. 1 and 2 for the paper
##  (do not be mislead by its name)
## (C) W. Huber 2006
##
doSave = TRUE
interact = !doSave

library("davidTiling")
source("setScriptsDir.R")
source(scriptsDir("categorizeSegments.R"))
graphics.off()
options(error=recover, warn=0)

rnaTypes  = c("seg-polyA-050909")
what = c("fig1", "fig2", "cover")[3]  ## 1=big overview plot, 2=seven detailed plots
## source(functionsDir("plotAlongChrom.R")) 

stopifnot(length(rnaTypes)==1)
so = get(rnaTypes)

if(!exists("myGff"))
  myGff = gff[ gff$Name!="tR(UCU)E", ]

ylim    = c(-4.5, 3.5)
switch(what,
       "fig1"= {
         width   = 16
         height  = 6
       },
       "fig2"= {
         width   = 11 ## width of A4 is 8.3''
         height  = 5.5
       },
       "cover"= {
         width   = 8.46 ## this is what PNAS says: 21.5cm x 22.5cm high
         height  = 8.86
       },
       stop("Zapperlot"))

if(doSave) {
  fn = sprintf("Figures/%s.pdf", what)
  pdf(file=fn, width=width, height=height)
} else {
  x11(width=width, height=height)
}
grid.newpage()

switch(what,
       fig1 = {
         plotAlongChrom(segObj=so, chr=1, coord = c(30, 130)*1e3, ylim=ylim, 
                        probeAnno = probeAnno, gff=myGff,
                        featureNoLabel = c("uORF", "binding_site", "TF_binding_site"),
                        doLegend=FALSE)  
         
       },
       fig2 = {
         ## -------------------------------------------------------------
         ## layout
         dx = 0.20
         dy = 0.05
         pushViewport(viewport(x=0.01, width=0.97, height=0.97, just=c("left", "center"),
                               layout=grid.layout(3, 8,
                                 height=c(1, dy, 1),
                                 width =c(dx, 1, dx, 1, dx, 1, dx, 1))))
         
         myPlot = function(row, col, ...) {
           pushViewport(viewport(layout.pos.row=row, layout.pos.col=col))
           grid.rect(x=-0.1, width=1.15, y=0.0, height=1.02, just=c("left", "bottom"),
                     default.units="npc", gp=gpar(lwd=0.2))
           plotAlongChrom(..., segObj=so, ylim=ylim,  
                          probeAnno = probeAnno, gff=myGff,
                          featureNoLabel = c("CDS", "uORF", "binding_site", "TF_binding_site"),
                          doLegend=FALSE)  
           popViewport()
         }
         
         ##
         ## A) 13:550k splicing RPS16A, RPL13B
         myPlot(1, 2, chr=13, coord = c(550044, 553360), main="a")
         
         ## B) GCN4
         myPlot(1, 4, chr=5, coord = c(138660, 141880), main="b")
         
         ## C) MET7, novel architecture
         myPlot(1, 6, chr=15, coord = c(784700, 790000), main="c")
         
         ## D) overlapping transcripts
         myPlot(1, 8, chr=14, coord = c(342200, 347545), main="d")
         
         ## E) SER3
         myPlot(3, 2, chr=5, coord = c(321900, 326100), main="e")
         
         ## F) 2:360.5-366.5: novel isolated
         myPlot(3, 4, chr=2, coord = c(360500, 365970), main="f")
         
         ## G) 9:221-227: novel antisense SPO22
         myPlot(3, 6, chr=9, coord = c(221000, 226500), main="g")
         
         ## legend
         fc = tilingArray:::featureColors(1)[c("CDS", "CDS_dubious", "uORF", "ncRNA", "TF_binding_site"), ]
         pc = c("Watson strand probe" = "#00441b",
           "Crick strand probe"  = "#081d58",
           "Non-unique probe"    = "grey")
         sc = c("Segment boundary" =  "#777777")
         
         pushViewport(dataViewport(xscale=c(0,1), yscale=c(-7,nrow(fc)+1), layout.pos.col=8, layout.pos.row=3))
         h1  = nrow(fc):1
         h2  = 0:(1-length(pc))
         h3  = -length(pc)
         
         w   = 0.2
         grid.rect(x=0, width=w, y=h1, height = unit(1, "native")- unit(2, "mm"), 
                   just  = c("left", "center"), default.units="native",
                   gp    = do.call("gpar", fc))
         grid.circle(x = w/2, y=h2, r=0.2, default.units="native",
                     gp = gpar(col=pc, fill=pc))
         grid.lines(x=w/2, y=h3+c(-.3,+.3), default.units="native", gp=gpar(col=sc))
         grid.text(label = c(rownames(fc), names(pc), names(sc)), x = w*1.1, y = c(h1,h2,h3),
                   just  = c("left", "center"), default.units="native",
                   gp=gpar(cex=.7))
         popViewport(2)
       },
       cover = {
         w = rbind(c(4,  6e4, 12e4),
                   c(4, 12e4, 18e4),
                   c(4, 18e4, 24e4))
         
         pushViewport(viewport(layout=grid.layout(nrow(w), 1), height=1, width=1))
         for(j in 1:nrow(w)) {
           pushViewport(viewport(layout.pos.row=j, layout.pos.col=1))
           plotAlongChrom(chr=w[j,1], coord=w[j,2:3], segObj=so, ylim=ylim, main="",
                          probeAnno = probeAnno, gff=myGff,
                          featureNoLabel = c("all"),
                          doLegend=FALSE)  
           popViewport()
         }
         popViewport()             
       }
 
       )

if(doSave) {
  dev.off()
  cmd = paste("convert -density 600", fn, "-compress RLE", sub(".pdf", ".tiff", fn))
  cat(cmd, "\n")
  if(!TRUE) {
    cat("PLEASE RUN!\n")
  } else {
    system(cmd)
  }
}
         
