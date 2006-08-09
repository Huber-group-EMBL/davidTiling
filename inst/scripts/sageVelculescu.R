## compare Velculescu et al. (Cell 1997) SAGE tages to our novel categories
##
## From: Victor Velculescu
## To: 'Lior David'
## Sent: Tuesday, January 24, 2006 9:00 AM
## Subject: RE: SAGE tags - Data request
## 
## The yeast SAGE data are available at
## ftp://genome-ftp.stanford.edu/pub/yeast/data_download/systematic_results/SAGE/
##
## and you can query individual SAGE tags directly at
## http://db.yeastgenome.org/cgi-bin/SAGE/querySAGE
##
## Good of luck with your experiments.
##
## Best regards,
## Victor Velculescu, M.D., Ph.D.
## Assistant Professor of Oncology
## Ludwig Center for Cancer Genetics and Therapeutics
## The Sidney Kimmel Comprehensive Cancer Center
## Johns Hopkins University School of Medicine
## 1650 Orleans St., Room 5M05
## Baltimore, MD 21231
## Tel      410-614-4582
## Fax     410-955-0548
##
## From: Lior David [mailto:liord@stanford.edu]
## Sent: Monday, January 23, 2006 9:57 PM
## To: velculescu@jhmi.edu
## Subject: SAGE tags - Data request
##
## Dear Victor,
## I am a post doc at the lab. of Ron Davis at Stanford and we are
## interested in the data of the SAGE tags from your 1997 Cell paper
## on the yeast transcriptome. Since I could not find any other way of
## downloading the data, would you be so kind as to send me the SAGE
## tags data and their expression levels?
##
## thanks in advance,
## Lior David

library("davidTiling")
library("matchprobes")


if(TRUE){
  source("setScriptsDir.R")
  source(scriptsDir("categorizeSegments.R"))
  
  if(!exists("sage")) {
    sage=read.table("sageVelculescu/genome-ftp.stanford.edu/pub/yeast/data_download/systematic_results/SAGE/sage.tab", header=FALSE, sep="\t", as.is=TRUE)
    colnames(sage) = c("sequence", "chromosome", "class", "start", "strand", "featurename", "genomehits", "conditiontested", "expressionvalue")
  }
  
  rt = "seg-polyA-050909"
  if(!exists("segSeq"))
    load(file.path(rt, "segSeq.rda"))
  segScore = cs[[rt]]
  stopifnot(nrow(segSeq)==nrow(segScore))
  
  ## The "class" of a SAGE tag is in column 3
  ## We care about class 3
  ##   table(sage[, 3])
  ##
  ##     1     2     3     4 
  ## 17049  6228  6663  4140 
  
  ## collapse to unique Tags
  sageSeq    = unique(sage[,"sequence"])

  ## sageSeq = sageSeq[1:1000]
  
  cat(length(sageSeq), "unique SAGE tags.\n")

  sageSeqRev = reverseSeq(complementSeq(sageSeq))
  whW = which(segScore$strand=="+")
  whC = which(segScore$strand=="-")
  stopifnot(setequal( union(whW, whC), 1:nrow(segScore)))

  ## sagS: sage tag sequences (possibly reverse-complement)
  ## segS: segment sequences (only those from one strand)
  ## return value: integer vector of same size as query (sagS)
  myFun = function(sagS, segS) {
    hit = rep(as.integer(NA), length(sagS))
    for(i in seq(along=sagS)) {
      if (i%%100==0) { cat(i,"") }
      g = grep(sagS[i], segS, fixed=TRUE, perl=TRUE)
      if(length(g)==1)
        hit[i] = g
    }
    return(hit)
  }

  v1 = myFun( sageSeq,    segSeq[whW])
  v2 = myFun( sageSeqRev, segSeq[whC]) 
  good = xor(is.na(v1), is.na(v2))
  
  hitSeq = rep(as.integer(NA), length(sageSeq))
  hitSeq[ good & !is.na(v1) ] =  whW[ v1[ good & !is.na(v1) ] ]
  hitSeq[ good & !is.na(v2) ] =  whC[ v2[ good & !is.na(v2) ] ]
  
  save(hitSeq, v1, v2, sage, segScore, file="sageVelculescu.rda")
} else {
  load("sageVelculescu.rda")
}


px = table(segScore$category[hitSeq])
py = table(segScore$category)
xl = "no. of SAGE tag hits"
yl = "found by array"
plot(as.vector(px), as.vector(py), log="xy", pch=21, col="orange", xlab=xl, ylab=yl)
text(as.vector(px), as.vector(py), names(px))

dev.copy(pdf, file="sageVelculescu.pdf", width=12, height=12)
dev.off()

sink("sageVelculescu.txt")
out = cbind(px, py, signif(100*px/py,3))
out = out[px>0,]
colnames(out) = c(xl, yl, "percent")
print(out)
sink()
