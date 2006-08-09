library("davidTiling")

source("setScriptsDir.R")
source(scriptsDir("categorizeSegments.R"))

options(error=recover)
rt = "seg-polyA-050909"
segScore = cs[[rt]] 


sp   = strsplit(segScore$featureInSegment, ", ")
whOp = which(listLen(sp)>1 & (segScore$level>0) & (segScore$frac.dup<0.5))

# annotated ORFs, no introns, no non-coding RNA
isAnnoORF = gff$feature=="CDS"

isGood = logical(length(whOp))
for(i in seq(along=whOp)) {
  j = whOp[i]
  good = sapply(sp[[j]], function(nm)
    any((gff$Name==nm) & isAnnoORF & (gff$start>=segScore$start[j]) & (gff$end<=segScore$end[j])))
  isGood[i] = all(good)
  ## if("YPR117W" %in% sp[[j]])browser()
}

writeSegmentTable(segScore[whOp[isGood],
                           c("featureInSegment", "chr", "strand", "start", "end", "length", "level","frac.dup")],
      sortBy = "level",
      fn = file.path(rt, "viz", "multifeatureSegments"), HTML=TRUE, 
      title  = paste("Segments with multiple features (", rt, " sample)", sep=""),
      interact=TRUE)
