# script to transfer posterior-probabilties downloaded from A.Siepel's site to usable RData object files.

# change paths accordingly to your system layout:
indir <- file.path("~","yeast","POSTPROBS.forPaper")
outdir <- file.path("~","yeast","data")

# postprob files are separated by chrom:
for (chr in c(1:16,"M")){
  cat(".")
  infile  <- file.path(indir,paste("chr",chr,".pp",sep=""))
  outfile <- file.path(outdir,paste("Chr",chr,"Posteriors.RData",sep=""))
  inread  <- read.delim(infile, header=FALSE, as.is=TRUE)
  basePosition   <- as.numeric(inread[,1])
  basePosterior  <- as.numeric(inread[,2])
  chrPost <- numeric(max(basePosition))
  chrPost[basePosition] <- basePosterior
  save(chrPost, file=outfile, compress=TRUE)
  #browser()
}#for



