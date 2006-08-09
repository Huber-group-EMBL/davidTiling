library("davidTiling")

## this will load the data
source("setScriptsDir.R")
source(scriptsDir("categorizeSegments.R"))
       
for(i in seq(along=cs))
  write.table(cs[[i]], file=paste(names(cs)[i], "txt", sep="."),
              col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
