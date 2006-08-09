## This script does the UTR related analyses
## It produces the following output
## 
## Fig 2a: scatterplot 5' vs 3' UTR lengths and marginal distributions
##    print the basic distribution statistics
## Fig 2b: barplot of UTR lengths of specific GO categories
##    write the table of significant GO categories: utrmap-GOcategories.txt
##
## Optional parts of the script
## explen: expression vs length -> utrmap-expression-vs-length.pdf
## 1vs2:   difference in estimated UTR lengths between hybeset 1 and 2
## wst:    html table of all UTR estimates

library("davidTiling")
library("geneplotter")

interact=(!TRUE)
options(error=recover, warn=0)
graphics.off()
if(!interact)
  sink("utrmap.txt")

source("setScriptsDir.R")
source(scriptsDir("categorizeSegments.R"))


## The variable rnaTypes contains the names of the hybe sets to look at.
## Fig 2a+2b are produced only for the first of these
## The segment tables are produced for both


what = c("wst", "explen", "1vs2")[numeric(0)]

##
## Fig 2a: scatterplot 5' vs 3' UTR lengths and marginal distributions
##         print the basic distribution statistics
##

for(rt in rnaTypes) {
  ul = utr[[rt]]
  cat("\n-------", rt, "-------\n")
  cat("Length distribution summary of", nrow(ul), "5'-UTRs:\n")
  print(summary(ul[, "5' UTR"]))
  cat("Length distribution summary of", nrow(ul), "3'-UTRs:\n")
  print(summary(ul[, "3' UTR"]))

  pdf(file=file.path("Figures", paste("utrmap-", rt, ".pdf", sep="")), width = 4, height = 4)
  scatterWithHist(ul,
      xlab=paste("Length of", colnames(ul)[1]),
                  ylab=paste("Length of", colnames(ul)[2]),
                  breaks = seq(0, 700, by=20), pch=20, barcols=rep("#d0d0d0", 2))
  dev.off()
}

## common:
comUTR = intersect(rownames(utr[[1]]), rownames(utr[[2]]))
allUTR = union(rownames(utr[[1]]), rownames(utr[[2]]))
cat("\n", length(comUTR)," in both ", paste(rnaTypes, collapse=" and "), ", ", 
    length(allUTR), " altogether.\n", sep="")

  
##
## GO analysis of UTR lengths
##
if(!exists("pGO")) {
  myUTR = utr[[1]]
  medAll = apply(myUTR, 2, median)

  goCat = getAllGO(rownames(myUTR), gff)
    allGO = unique(unlist(goCat))
  gm    = matrix(FALSE, nrow=length(allGO), ncol=length(goCat))
  rownames(gm) = allGO
  colnames(gm) = names(goCat)
  
  for(i in seq(along=goCat))
    gm[ goCat[[i]], i] = TRUE
  
  w.all = which(rownames(gm)=="all")
  stopifnot(length(w.all)==1, all(gm[w.all, ]))
  gm = gm[-w.all,]
  
  GOterms = mget(rownames(gm), GOTERM)

  ## Wilcoxon test
  wtfun = function() {
    function(z) {
      sz = sum(z)
      if(sz>=5&&(length(z)-sz)>=5) {
        w5 = wilcox.test(myUTR[, "5' UTR"] ~ z)$p.value
        w3 = wilcox.test(myUTR[, "3' UTR"] ~ z)$p.value
      } else {
        w5 = w3 = as.numeric(NA)
      }
      m5 = median(myUTR[z, "5' UTR"])
      m3 = median(myUTR[z, "3' UTR"])
      c(w5, m5, w3, m3, sz)
    }
  }
  wt = wtfun()

  pGO = t(apply(gm, 1, wt))
  colnames(pGO) = c("p 5' UTR", "median 5' UTR", "p 3' UTR", "median 3' UTR", "nrGenes")
}


## Select those with p<=pthresh
pthresh = 2e-3
pp  = pmin(pGO[, "p 5' UTR"], pGO[, "p 3' UTR"])
ord = order(pp)
ord = ord[ which(pp[ord]<= pthresh) ]
sGOterms = GOterms[ord]

for(j in 1:ncol(myUTR)) {
  pname = paste("p",      colnames(myUTR)[j])
  mname = paste("median", colnames(myUTR)[j])
  cat("--------------------------------------------------\n",
      colnames(myUTR)[j], "\n",
      "--------------------------------------------------\n", sep="")
  for(s in ord) {
    g   = rownames(gm)[s]
    cat(g, " p=", format.pval(pGO[s, pname]), "  median=",
        pGO[s, mname], " (versus ", medAll[j], ")\n", sep="")
    print(GOterms[[s]])
    nr = pGO[s, "nrGenes"]
    cat(nr, " genes", sep="")
    if(nr<=20)
      cat(":", replaceSystematicByCommonName(names(which(gm[s,]))))
    cat("\n\n")
  }
}

## Write tab-delimited table
## GO-ID, Term, Median 3' UTR, Median 5' UTR, No. Genes

utrGO = cbind(data.frame(
  GOID = I(sapply(sGOterms, GOID)),
  Term = I(sapply(sGOterms, Term)), 
  Ontology = I(sapply(sGOterms, Ontology))), 
  pGO[ord, ])
write.table(utrGO, file="utrmap-GOcategories.txt", sep="\t", row.names=FALSE)


## categories to be shown in the plot
gonames = rev(c("GO:0007154" = "cell communication",
       "GO:0007264" = "small GTPase mediated signal transduction",
       "GO:0004672" = "protein kinase activity",
       "GO:0016310" = "phosphorylation",
       "GO:0050790" = "regulation of enzyme activity",
       "GO:0000271" = "polysaccharide biosynthesis",
       "GO:0005886" = "plasma membrane",
       "GO:0005618" = "cell wall",
       "GO:0005215" = "transporter activity",
       "GO:0006811" = "ion transport",
       "GO:0005746" = "mitochondrial electron transport chain",
       "GO:0000082" = "G1/S transition of mitotic cell cycle",
       "GO:0000086" = "G2/M transition of mitotic cell cycle",
       "GO:0051325" = "interphase",
       "GO:0016072" = "rRNA metabolism",
       "GO:0006396" = "RNA processing",
       "GO:0007046" = "ribosome biogenesis",
       "GO:0005830" = "cytosolic ribosome (sensu Eukaryota)",
       "GO:0003735" = "structural constituent of ribosome",
       "GO:0004532" = "exoribonuclease activity",
       "GO:0005730" = "nucleolus",
       "GO:0005681" = "spliceosome complex"))

## select and sort
mt = match(names(gonames), utrGO$GOID)
stopifnot(!any(is.na(mt)), all(utrGO$Term[mt]==gonames))
utrGO = utrGO[mt, ]
          
utrGO$Term = gsub(" (sensu Eukaryota)", "", utrGO$Term, extended=FALSE)
utrGO$Term = gsub("transduction", "transd.", utrGO$Term, extended=FALSE)
utrGO$Term = gsub("mitochondrial", "mitoch.", utrGO$Term, extended=FALSE)
utrGO$Term = gsub("of mit", "mit", utrGO$Term, extended=FALSE)
utrGO$Term = gsub("mediated", "med.", utrGO$Term, extended=FALSE)

totmedian5 = -68
totmedian3 =  91

p.cutoff = 0.05

signif5 <- (-1)^as.numeric(-utrGO[,"median 5' UTR"] > totmedian5) * as.numeric(utrGO[,"p 5' UTR"]< p.cutoff)
signif3 <- (-1)^as.numeric( utrGO[,"median 3' UTR"] < totmedian3) * as.numeric(utrGO[,"p 3' UTR"]< p.cutoff)

upperQuart <-  90
lowerQuart <- -upperQuart

table1stats <- rbind(-utrGO[,"median 5' UTR"]+lowerQuart, lowerQuart, Inf, upperQuart, upperQuart+utrGO[,"median 3' UTR"])

tab1 <- list(stats=table1stats, n=rep(100, nrow(utrGO)))

transformColors <- function(COL, H, S, V){
  colhsv <- as(hex2RGB(COL),"HSV")
  colhsv@coords[,1] <- colhsv@coords[,1]+H
  colhsv@coords[,2] <- S
  colhsv@coords[,3] <- V
  return(hex(colhsv))
}

bpcols <- transformColors(brewer.pal(9,"Greens")[3:5],H=c(-10,0,10),S=0.5, V=1)
cccols <- transformColors(brewer.pal(9,"Oranges")[2:4],H=c(-10,0,10),S=0.5, V=1)
mfcols <- transformColors(brewer.pal(9,"Blues")[2:4],H=c(-10,0,10),S=0.5, V=1)

boxcolors = rep("",nrow(utrGO))
for (i in 1:nrow(utrGO)){
  boxcolors[i] <- switch(utrGO$Ontology[i],
         "BP"=bpcols[3],
         "CC"=cccols[3],
         "MF"=mfcols[3] 
         )
}

mypsize  <- 8
textcex  <- 1.3
textfont <- 1

##
## Fig 2b: boxplot GO categroy UTR lengths
##
pdf(file=file.path("Figures", "utrLengthsGOCateg.pdf"), width=7, height=5.5, pointsize=mypsize)
par(mar=c(3.6,0.1,0,1))

bxp(tab1, horiz=T, col=boxcolors, medlty="blank", yaxt="n", xaxt="n", staplelty=0, lwd=0.5, frame.plot=FALSE)

## lines for median UTR lengths (taken from table 1 header):
lines(x=rep(totmedian5+lowerQuart,2),y=c(0,nrow(utrGO)+0.4), lty=2, col="grey")
lines(x=rep(totmedian3+upperQuart,2),y=c(0,nrow(utrGO)+0.4), lty=2, col="grey")

## x-axis
axis(1, at=c(c(-180,-150,-120, -90, -60,-30, 0)+lowerQuart,upperQuart+c(0, 30,60,90,120,150,180)),
     labels=c("180","150","120","90","60","30","0","0","30","60","90","120","150","180"), line=-1, cex=textcex)

medutr5 <- -utrGO[,"median 5' UTR"]+lowerQuart
medutr3 <-  utrGO[,"median 3' UTR"]+upperQuart
mtext(side=1, at=c(lowerQuart-100, upperQuart+100), text=c("Length of 5' UTR", "Length of 3'UTR"), line=2, cex=textcex)
text(x=0,y=1:nrow(utrGO), utrGO$Term, cex=textcex, font=textfont) 
text(x=medutr5-0.2,y=1:length(medutr3), "|", cex=0.5, col="black") # ends of whiskers
text(x=medutr3+0.5,y=1:length(medutr5), "|", cex=0.5, col="black")

# mark significant differences:
text(x=medutr3[signif3>0], y=which(signif3>0), "*", col="red",cex=2)
text(x=medutr3[signif3<0], y=which(signif3<0), "*", col="blue",cex=2)
text(x=medutr5[signif5>0], y=which(signif5>0), "*", col="red",cex=2)
text(x=medutr5[signif5<0], y=which(signif5<0), "*", col="blue",cex=2)

## finished Fig. 2b
dev.off()
if(!interact)
  sink()


##
## expression vs length
##
if("explen" %in% what){
  cols = brewer.pal(12, "Paired")

  investigateExpressionVersusLength = function(lev, len, main) {
    theCut = cut(len, breaks=quantile(len, probs=c(0, 0.95, 1)))
    e = tapply(lev, theCut, ecdf)
    theCol = cols[1]
    plot(e[[1]], pch=".", xlab="level", main=main)
    for(i in 2:length(e)) {
      theCol = cols[i*2]
      lines(e[[i]], col.hor=theCol, col.points=theCol, col.vert=theCol)
    }
  }

  investigateLengthVersusLength = function(l1, l2, ...) {
    cc = cor(l1, l2, method="kendall")
    plot(l1+1, l2+1, log="xy", pch=".", main=paste("length, cor=", signif(cc, 3)), ...)
  }
  
  if(!interact) {
    pdf(file=paste("utrmap-expression-vs-length.pdf", sep=""), width = 10.5, height = 14)
  } else {
    x11(width = 10.5, height = 14)
  }
  par(mfrow = c(4, 3))

  for(rt in rnaTypes) {
    s = csu[[rt]]
    ## get the CDS length
    mt = match(s[,"featureInSegment"], gff[,"Name"])
    stopifnot(!any(is.na(mt)))
    cdslen = gff[mt, "end"]-gff[mt, "start"]

    investigateExpressionVersusLength(s[,"level"], s[,"utr3"], paste(rt, ": length of 3' UTR", sep=""))
    investigateExpressionVersusLength(s[,"level"], s[,"utr5"], paste(rt, ": length of 5' UTR", sep=""))
    investigateExpressionVersusLength(s[,"level"], cdslen, paste(rt, ": length of CDS", sep=""))
    investigateLengthVersusLength(s[,"utr3"], s[,"utr5"], xlab="3' UTR", ylab="5' UTR")
    investigateLengthVersusLength(cdslen, s[,"utr3"], xlab="CDS", ylab="3' UTR")
    investigateLengthVersusLength(cdslen, s[,"utr5"], xlab="CDS", ylab="5' UTR")
  }
  
  if(!interact)
    dev.off()
}

##
## difference between hybeset 1 and 2
## 
if("1vs2" %in% what){
  if(interact) {
    x11(width=6.6, height=10)
  } else {
    pdf(file=paste("utrmap-scatter.pdf", sep=""), width=6.6, height=10)
  }
  
  par(mfrow=c(3,2))
  for(i in 1:2){
    px = utr[[1]][comUTR,i]
    py = utr[[2]][comUTR,i]
    axlim = c(0, quantile(c(px, py), 0.8))
    plot(px, py,
         main = paste("length of ", colnames(utr[[1]])[i], " (", length(comUTR), " common)", sep=""),
         xlab = longNames[rnaTypes[1]], ylab=longNames[rnaTypes[2]],
         xlim = axlim, ylim = axlim, pch=20)
    abline(a=0, b=1, col="red")
  }

  vec = c("5' UTR", "3' UTR")
  d = utr[[1]][comUTR, vec] - utr[[2]][comUTR, vec]
  colnames(d)=vec
  
  ex1 = csu[[1]][comUTR, "level"]
  ex2 = csu[[2]][comUTR, "level"]

  ex = cbind(difference=ex1-ex2, average=(ex1+ex2)/2)

  col2 = "lightblue"; col3="pink"
  for(j in 1:ncol(ex)) {
    for(i in 1:ncol(d)){
      ec1 = ecdf(ex[d[,i]==0, j])
      ec2 = ecdf(ex[d[,i]>0, j])
      ec3 = ecdf(ex[d[,i]<0, j])
      plot(ec1, pch=".", xlab=colnames(ex)[j], main=colnames(d)[i])
      lines(ec2, pch=".", col.points=col2, col.hor=col2, col.vert=col2)
      lines(ec3, pch=".", col.points=col3, col.hor=col3, col.vert=col3)
    }
  }
  if(!interact)
    dev.off()
}


##
## WRITE THE SEGMENT TABLE
##
if("wst" %in% what){
  for(rt in rnaTypes) {
    fn = file.path(indir[rt], "viz", "utrmap")
    nr = nrow(csu[[rt]])
    cat("Writing", nr, "UTRs to", fn, "\n")
    writeSegmentTable(csu[[rt]],
                      title    = paste(nr, "UTRs from", longNames[rt], "data"),
                      subtitle = paste("<i>Version ",  strsplit(rt, "-")[[1]][3], "</i>", sep=""),
                      fn=fn, sortBy = "goodUTR", sortDecreasing=TRUE, interact=interact)
  }
}

