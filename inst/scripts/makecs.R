library("davidTiling")

source("setScriptsDir.R")
source(scriptsDir("categorizeSegments.R"))

save(list=c("cs", "featNames", "functionsDir","scriptsDir"),
     file="categorizedSegments.rda")
