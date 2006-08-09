## write FASTA files with segments so they can be blasted against
## another genome

## options
seqDir = "SGD-0508"
outdir = "fasta"
interact = TRUE
rnaTypes  = "seg-polyA-050909"

library("davidTiling")
source("setScriptsDir.R")
source(scriptsDir("categorizeSegments.R"))

options(error=recover)
if(!exists("gff"))
  load("probeAnno.rda")

if(!exists("fsa")) {
  fsa = new.env()
  fsa.files = paste("chr", c(sapply(1:16, function(n) sprintf("%02d", n)), "mt"),
    ".fsa", sep="")
  for(i in seq(along=fsa.files)) {
    s = readLines(file.path(seqDir, fsa.files[i]))
    s = paste(s[-1], collapse="")
    assign(paste(i), s, envir=fsa)
    cat(fsa.files[i], ": ", nchar(s), "\n", sep="")
  }
}


## Check if the sequence lengths found here coincide with the end of
## the telomere in the GFF table. If yes, all is well!
chrLengths = sapply(fsa, nchar)
chrLengths = chrLengths[order(as.numeric(names(chrLengths)))]

## double-check
sgff = gff[ gff$feature=="telomere", ]
for(i in 1:16) {
  w = which(sgff$chr==i)
  stopifnot(length(w)==2)
  stopifnot(chrLengths[i]==sgff$end[w[2]])
}

## write segment seqs:
for(rt in rnaTypes) {
  ss = cs[[rt]]
  stopifnot(all(ss$strand %in% c("+", "-")))
  
  fout = file.path(rt, "segSeq.fsa")
  con = file(fout, open="wt")

  segSeq = character(nrow(ss))
  for(i in 1:nrow(ss)) {
    sequence = substr(fsa[[paste(ss$chr[i])]],
      start = ss$start[i],
      stop  = ss$end[i])
    cat(">", i, "\n", sequence, "\n", sep="", file=con, append=TRUE)
    segSeq[i] = sequence
  }
  close(con)
  save(segSeq, file=file.path(rt, "segSeq.rda"))
}
