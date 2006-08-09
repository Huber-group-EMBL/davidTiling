##
##  Export a table of antisense transcripts
##
library("davidTiling")

source("setScriptsDir.R")
source(scriptsDir("categorizeSegments.R"))

hybeSet ="seg-polyA-050909"
s = cs[[hybeSet]]

sel = s$category %in% c("novel antisense - filtered", "novel antisense - unassigned")
antiSenseTable = s[ sel, c("category", "chr", "strand", "start", "end", "level","oppositeFeature")]
antiSenseTable = cbind(antiSenseTable,
  oppositeName=I(replaceSystematicByCommonName(antiSenseTable$oppositeFeature)))

save(antiSenseTable, file=file.path(hybeSet, "antiSenseTable.rda"))
