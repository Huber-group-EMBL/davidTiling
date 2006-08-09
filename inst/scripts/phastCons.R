library("tilingArray")
source("setScriptsDir.R")

graphics.off()
options(error=recover, warn=2)

rnaTypes = c("seg-polyA-050909")
outfile  = "phastCons"

source(scriptsDir("readSegments.R"))
source(scriptsDir("categorizeSegments.R")) 

interact=TRUE
if(!interact){
  sink(paste(outfile, ".txt", sep=""))
  cat("Made on", date(), "\n\n")
}

source(scriptsDir("calcThreshold.R"))

##
## CATEGORIZE
##
if(!exists("cs")) {
  cs = vector(mode="list", length=length(rnaTypes))
  names(cs)=rnaTypes

  cat("\n\nCategorization of segments:\n",
      "===========================\n", sep="")
  for(rt in rnaTypes) {
    cat(rt, ":\n", sep="")
    s = categorizeSegments(get(rt))
    cs[[rt]] =s
  }
} else {
  cat("\n**************************************************\n",
        "*      NOT REDOING categorizeSegments            *\n",
        "**************************************************\n", sep="")
}

## read phastCons conserved segments
if(!exists("pc"))
  pc = read.table("phastConsElements.txt", sep="\t", as.is=TRUE, header=FALSE)

rt  = rnaTypes[1]
res = numeric(nlevels(cs[[rt]]$simpleCatg))
names(res) = levels(cs[[rt]]$simpleCatg)

for(chr in 1:16){
  p  = pc[ pc[, 2]==paste("chr", chr, sep=""), ]
  mp = max(p[, 4])
  sp = mp + 5000
  vp = logical(sp)
  for(j in 1:nrow(p))
    vp[ p[j,3]:p[j,4] ] = TRUE
  
  for(str in c("+", "-")) {
    s  = cs[[rt]][ cs[[rt]]$chr==chr & cs[[rt]]$strand==str, ]
    ms = max(s$end) 
    cat(sprintf("%2d  %8d  %8d\n", chr, ms, ms-mp))
    stopifnot(sp >= ms)
    vs = matrix(as.logical(FALSE), nrow=ms, ncol=length(res))
    colnames(vs) = names(res)
    for(j in which(!is.na(s$simpleCatg)))
      vs[ s$start[j]:s$end[j], as.character(s$simpleCatg[j]) ] = TRUE
    for(i in 1:ncol(vs))
      res[i] = res[i] + mean(vs[, i]*vp[1:ms])
    stopifnot(!any(is.na(res)))
  }
}


if(!interact)
  sink()
