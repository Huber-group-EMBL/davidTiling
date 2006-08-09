## look for cases where transcription over an active promoter
##
## We use the ends of 5' UTRs as proxy for promoters

library("davidTiling")

source("setScriptsDir.R")
source(scriptsDir("categorizeSegments.R"))

rt = "seg-polyA-050909"
segScore = cs[[rt]] 
wh = which(!is.na(segScore$goodUTR))

segCrit = ((segScore$simpleCatg %in% c("annotated ORF", "ncRNA(all)")) &
           (segScore$overlap %in% c(">=50%", "complete")) &
           (segScore$level >= 1))

levs = matrix(as.numeric(NA), nrow=length(wh), ncol=3)
colnames(levs)= c("level", "same strand adjacent", "opposite")

for(i in seq(along=wh)) {
  if(i%%100==0)cat(i, "")
  js     = wh[i]  
  chr    = segScore$chr[js]
  strand = segScore$strand[js]
  switch(segScore$strand[js],
    "+" = {
      pro = segScore$start[js]
      jsame = js-1
    },
    "-" = {
      pro = segScore$end[js]
      jsame = js+1
    },
    stop("Sapperlot"))

  levs[i, 1] = segScore$level[js]
  if(segCrit[jsame])
    levs[i, 2] = segScore$level[jsame]

  ## joppo = which((segScore$chr==chr) & segCrit & (segScore$start<=pro-25) & (segScore$end>=pro+25) & (segScore$strand==otherStrand(strand)))
  ## if(length(joppo)>0)
  ##   levs[i, 3] = mean(segScore$level[joppo])

  ## if(segScore$featureInSegment[js] == "YDL181W") browser()
  
}

## pairs(levs, panel=function(...){points(...,pch=21);abline(v=0,col="blue");abline(h=0,col="blue")})


selSame = levs[,2]>0
selOppo = levs[,3]>0
whTot   = which(selSame|selOppo)

con = file("transcriptionThruPromoters.txt", open="w")
cat("Segments above background adjacent to mapped UTR:", sum(selSame, na.rm=TRUE), "\n",
    "                      ... opposite ...          :", sum(selOppo, na.rm=TRUE), "\n",
    "Together                                        :", length(whTot), "\n",
    "Out of total                                    :", length(selSame), "\n", file=con)
close(con)

out = segScore[wh[whTot], c("featureInSegment", "chr", "strand", "start", "end", "length", "level","frac.dup")]
ord = order(pmax(levs[whTot,2], levs[whTot,3], na.rm=TRUE), decreasing=TRUE)

writeSegmentTable(out[ord,],
      fn = file.path(rt, "viz", "transcriptionThruPromoters"), HTML=TRUE, 
      title  = paste("Instances of transcription through promoters (", rt, " sample)", sep=""),
      interact=TRUE)
