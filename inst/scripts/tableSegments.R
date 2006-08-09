library("davidTiling")
source("setScriptsDir.R")
source(scriptsDir("categorizeSegments.R"))
graphics.off()
options(error=recover, warn=0)

interact = TRUE
what     = c("fig4", "compare", "overlap", "lvsx", "wst")[1]

outfile = "tableSegments"

if(!interact){
  sink(paste(outfile, ".txt", sep=""))
  cat("Made on", date(), "\n\n")
}


fillColors = c(c(brewer.pal(10, "Paired")[c(1, 2, 6, 8, 5:8, 2, 10)]), "#d0d0d0")
names(fillColors) = c("overlap <50%", "overlap >=50%",
       "novel antisense", "novel isolated",
       "novel isolated - unassigned",  "novel isolated - filtered",
       "novel antisense - unassigned", "novel antisense - filtered",
         "annotated ORF", "ncRNA(all)", "untranscribed")

lineColors = c(brewer.pal(8, "Paired")[c(1,2,6,8)], "grey")
names(lineColors) =c("annotated ORF", "ncRNA(all)", 
       "novel antisense - filtered", "novel isolated - filtered",
       "unexpressed isolated")

##
## piechart
##
if("fig4" %in% what){

  if(interact) {
    x11(width=5*length(rnaTypes), height=4*3.2)
  } else {
    pdf("fig4.pdf", width=2.5*length(rnaTypes), height=2.6*3)
  }

  layout(matrix(1:8, ncol=2, byrow=TRUE), widths=c(1,1), height=c(2,1,1,0.5))
  counts = NULL

  cat("\nSegment overlap with known features (genes):\n",
        "============================================\n\n", sep="")
  mai.old = par(mai=c(0.3,0.1,0.3,0.25))
  for(irt in seq(along=rnaTypes)) {
    rt = rnaTypes[irt]
    s  = cs[[rt]] 
   
    px = table(s[, "pieCat"])
    ##labels = LETTERS[ match(names(px), levels(s[, "pieCat"])) ]
    ##labels = names(px)[ match(names(px), levels(s[, "pieCat"])) ]
    
    stopifnot(all(names(px) %in% names(fillColors)))
    counts = cbind(counts, px)
    pie(px, radius=0.9, main=longNames[rt], col = fillColors[names(px)], labels = paste(px))

    category = s[, "category"]
    levels(category) = sub("ncRNA", "ncRNA(all)", levels(category))
    category[ s[, "simpleCatg"]=="ncRNA(all)" ] = "ncRNA(all)"

    cat(rt, ":\n", sep="")
    tab = table(category, s[, "overlap"])
    tab = tab[rowSums(tab)!=0, ]
    print(tab)
    cat("\n\n")
  } ## for rt
  par(mai.old)
  
  colnames(counts)=rnaTypes
  cat("\nSegment counts (pie charts):\n",
        "============================\n\n", sep="")
  print(counts)
  cat("\n")
  
  ##
  ## LENGTH & LEVEL DISTRIBUTIONS
  ##
  cat("\n\nLength distributions:\n",
          "=====================\n", sep="")
  maxlen=4000
  mai = par("mai")
  mai[2:3] = c(0.5,0.1)
  par(mai=mai)
  br  = seq(1, 4, by=0.2)
  xat = seq(1, 4, by=1)
  xtickLabels = paste(10^xat)
  for(irt in seq(along=rnaTypes)) {
    s   = cs[[rnaTypes[irt]]]
    plotCat = s[, "pieCat"]
    stopifnot(all(levels(plotCat) %in% names(fillColors)))
    len = split(s[, "length"], plotCat)
    ## len = lapply(len, function(z) {z[z>maxlen]=maxlen; z})
    slen = lapply(len, function(z) {z[z>maxlen]=NA; log(z, 10)})
    cols = fillColors[names(len)]

    showDens(slen, breaks=br, xat=xat, xtickLabels=xtickLabels, col=cols, main="",
             xlab=expression(plain(Length)~~plain((Nucleotides))), ylab="")
    text(2*br[1]-br[3], length(slen)/2, "Frequency", adj=c(0.5, 0.5), srt=90, xpd=NA)
    cat("\n", rnaTypes[irt], "\n")
    print(sapply(len, summary))
  }


  br  = seq(-0.2, 6.6, by=0.2)
  xat = seq(0, 6, by=1)
  lvall = lapply(rnaTypes, function(rt) cs[[rt]][, "level"] )
  rg    = quantile(unlist(lvall), probs=c(0.001, 0.999), na.rm=TRUE)
  stopifnot(rg[2]<=br[length(br)])
  
  for(irt in seq(along=rnaTypes)) {
    s       = cs[[rnaTypes[irt]]]
    plotCat = s[, "pieCat"]
    ## levels(plotCat) = c(levels(plotCat), "untranscribed")
    ## plotCat[ s[, "simpleCatg"]=="untranscribed" ] = "untranscribed"
    stopifnot(all(levels(plotCat) %in% names(fillColors)))

    lv = split(lvall[[irt]], plotCat)
    lv = lapply(lv, function(z)
      ifelse(z<=rg[2], ifelse(z>=rg[1], z, rg[1]), rg[2]))
    showDens(lv, breaks=br, xat=xat, col=fillColors[names(lv)], main="",
             xlab=expression(log[2]*~~plain(Level)), ylab="")
    text(2*br[1]-br[5], length(lv)/2, "Frequency", adj=c(0.5, 0.5), srt=90, xpd=NA)
  }

  ## legend
  opar = par(mai=rep(0,4))
  dy = 0.45
  
  for(j in 1:2) {
    w  = list(1:2, 3:6)[[j]]
    x0 = c(0.3, 0)[j]
    plot(c(0,1), c(0,5), type="n", bty="n", xaxt="n", yaxt="n")
    y = 5-seq(along=w)
    rect(rep(x0, length(y)), y-dy, rep(x0+0.22, length(y)), y+dy, col=fillColors[names(px)[w]], border="black")
    text(rep(x0+0.25, length(y)), y, names(px)[w], adj=c(0, 0.5))
  }
  par(opar)
  
  if(!interact)
    dev.off()
}

##
## Compare total to poly-A, the goal is: which transcripts do we find 
## specifically in total RNA?
##
if("compare" %in% what){
  stopifnot(length(rnaTypes)==2)
  s1     = cs[[rnaTypes[1]]]
  s2     = cs[[rnaTypes[2]]]
  start1 = s1[, "start"]
  end1   = s1[, "end"]
  start2 = s2[, "start"]
  end2   = s2[, "end"]

  unTrCatgs = c("excluded", "untranscribed")
    
  isTr1  = !(s1[, "category"] %in% unTrCatgs)
  jstart = jend = 1
  ov     = numeric(nrow(s2))

  for(k in 1:nrow(s2)) {
    ## make sure that jstart points to a segment in s1 whose
    ## start <= ks <= end, where ks=start of current segment in s2.
    ks = start2[k]
    ke = end2[k]
    while(!((start1[jstart]<=ks) && (end1[jstart]>=ks)))
      jstart = jstart+1
    ## Similarly, make sure that jend points to a segment in s1 whose
    ## start <= ke <= end, where ke=end of current segment in s2.
    while(!((start1[jend]<=ke) && (end1[jend]>=ke)))
      jend = jend+1
    stopifnot(jstart<=nrow(s1), jend<=nrow(s1))
    ## cat(k, jstart, jend, "\n")
    
    if(jstart==jend) {
      lne = lns = 0
      lni = (ke-ks+1) * as.numeric(isTr1[jstart])
    } else {
      lns = (end1[jstart] - ks + 1) * as.numeric(isTr1[jstart])
      lne = (ke - start1[jend] + 1) * as.numeric(isTr1[jend])
      if(jend-jstart>1) {
        j = (jstart+1):(jend-1)
        lni = sum( (end1[j]-start1[j]+1) * as.numeric(isTr1[j]) )
      } else {
        lni = 0
      }
    }
    stopifnot(lns>=0, lne>=0, lni>=0)
    ov[k] = (lns+lni+lne) / (ke-ks+1)
  }
  
  cat("\n\nOverlap of segments from total RNA with those from poly-A\n",
          "(not counting: ", paste(unTrCatgs, collapse=", "), ")\n",
          "=========================================================\n", sep="")
  
  tab = table(s2[, "category"], ov > .5)
  tab = tab[ !(rownames(tab)%in%c("excluded", "untranscribed")), 2:1]
  print(tab)
  cat("\n\n")
}



##
##
if("overlap" %in% what){
  selectedCategories = c(
    "(1): only 'overlap <50%' (i.e. in 'overlappingFeature' but not 'mostOfFeatureinSegment')", 
    "(2): only 'overlap >=50%', but not 'complete' (i.e. in 'mostOfFeatureinSegment' but not 'featureInSegment')",
    "(3): 'complete' (i.e. in 'featureInSegment')",
    "(4): (1) AND (2), i.e. 'overlap <50%' in one segment and 'overlap >=50%' in another segment.\n")
  
  nsc = length(selectedCategories)
  nfn = length(featNames)
  tab = matrix(NA, nrow=nsc*nfn, ncol=length(rnaTypes)+1)
  rownames(tab) = paste("(", rep(1:nsc, nfn), ") ", rep(names(featNames), each=nsc), sep="")
  colnames(tab) = c(rnaTypes, "in genome")
  
  multiGenesPerSegment = matrix(as.numeric(NA), nrow=2, ncol=length(rnaTypes))
  colnames(multiGenesPerSegment) = rnaTypes
  rownames(multiGenesPerSegment) = c("any feature", "only annotated ORFs")
        
  for(irt in seq(along=rnaTypes)) {
    isT = !(cs[[irt]][, "category"] %in% c("excluded", "untranscribed"))
    ovf = strsplit(cs[[irt]][isT, "overlappingFeature"],     split=", ")
    mof = strsplit(cs[[irt]][isT, "mostOfFeatureInSegment"], split=", ")
    fis = strsplit(cs[[irt]][isT, "featureInSegment"],       split=", ")
    s1  = mapply(setdiff, ovf, mof)
    s2  = mapply(setdiff, mof, fis)
    s3  = fis
    for(isc in seq(along=selectedCategories)) {
      fIDs = unique(switch(isc,
        unlist(s1),
        unlist(s2),
        unlist(s3),
        intersect(unlist(s1), unlist(s2))))
      for(k in seq(along=featNames)) {
        m =  intersect(fIDs, featNames[[k]])
        tab[ (k-1)*length(selectedCategories)+isc, irt ] = length(m)
        if(irt==1 && isc==4 && k==1)
          writeLines(replaceSystematicByCommonName(m), con="tableSegments-unusual-architecture.txt")
      }
    }
    multiGenesPerSegment[1, irt] = sum(listLen(mof)>=2)
    multiGenesPerSegment[2, irt] = sum(listLen(mof)>=2 & sapply(fis, function(g) all(g %in% featNames$"annotated ORFs")))
  }
  tab[ , "in genome" ] = rep(listLen(featNames), each=length(selectedCategories))

  cat("\nHow many unique known features (SGD Names) do we find that occur with\n",
      paste(selectedCategories, collapse="\n"), "\n",
      "=====================================================================\n", sep="")
  print(tab)  
  cat("\n\n")

  cat("How many segments have more than one verified/uncharaterized gene in featureInSegment:\n",
      "======================================================================================\n", sep="")
  print(multiGenesPerSegment)
}


##
## WRITE THE SEGMENT TABLE
##
if("wst" %in% what){
  for(rt in rnaTypes) {
    s = cs[[rt]]
    s$segID = paste(1:nrow(s))
    drop =  (s[,"category"]=="excluded") | (s[,"category"]=="untranscribed"&(!s[,"isUnIso"]))
    writeSegmentTable(s[!drop, ],
      fn = file.path(indir[rt], "viz", "index"), HTML=TRUE, 
      sortBy = "category-level",
      title    = paste("Segmentation table for", longNames[rt], "sample"),
      subtitle = paste("<i>Version ",  strsplit(rt, "-")[[1]][3], "</i>", sep=""),
      interact=interact)
  }
  cat("\n")
}

##
## LENGTH VERSUS EXPRESSION LEVEL
##
if("lvsx" %in% what){
  if(!interact) {
    pdf(paste(outfile, "lvsx.pdf", sep="-"), width=14, height=length(rnaTypes)*3)
    pch="."
  } else {
    pch=18
  }
  par(mfrow=c(length(rnaTypes), 2))
  maxlen=5000
  br = seq(0, maxlen, by=200)
  selectedCategories = c("annotated ORF", "novel isolated - filtered")
  for(rt in rnaTypes) {
    s = cs[[rt]]
    stopifnot(all(selectedCategories %in% levels(s[,"simpleCatg"])))
    ylim = quantile(s[,"level"][s[,"simpleCatg"] %in% selectedCategories],
      probs=c(0.01, 0.99), na.rm=TRUE)
    for(lev in selectedCategories) {
      len = s[s[,"simpleCatg"] == lev, "length"]
      exl = s[s[,"simpleCatg"] == lev, "level"]
      len[len>maxlen] = maxlen
      plot(len, exl, pch=pch, ylim=ylim,
           main=paste(longNames[rt], ": ", lev, sep=""),
           ylab="expression level", xlab="length")
      lf = loess(exl ~ len)
      slen = sort(len)
      lines(slen, predict(lf, newdata=slen), col="blue")
    }
    rm(list=c("s", "ylim", "len", "lev", "exl", "lf", "slen"))
  }
  if(!interact)
    dev.off()
}

if(!interact)
  sink()
