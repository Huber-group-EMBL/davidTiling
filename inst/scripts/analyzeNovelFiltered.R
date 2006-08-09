# script to analyze novel isolated (filtered) transcripts
load(file.path("~","downloads","polyASegments.RData"))

posteriorDir <- file.path("~","yeast","data") #where to find the save posterior probabilities

## now: for each segment in categorizedSegments, get conservation score:
#   take median posterior prob for being in state "conserved" over the segment
currentChrom <- 0
conserv <- numeric(nrow(polyA))

for (i in 1:nrow(polyA)){
  if (i %% 1000 == 0){
    cat(i,"... ")
    save(conserv, file=file.path("~","yeast","data","conserv.RData"), compress=TRUE)
  }
  thisChrom <- polyA[i,1]
  if (thisChrom != currentChrom){
    cat("next chromosome ... ")
    currentChrom <- thisChrom
    if (thisChrom==17) thisChrom <- "M"
    postFile <- file.path(posteriorDir,paste("Chr",thisChrom,"Posteriors.RData",sep=""))
    load(postFile)
  }# if (thisChrom != currentChrom)
  thesePosts <- chrPost[polyA[i,"start"]:polyA[i,"end"]]
  conserv[i] <- median(thesePosts)
}#for i

# now compare the categories:
# the novel categories:
summary(conserv[which(polyA$category=="novel isolated - filtered")])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.00000 0.00880 0.03425 0.21970 0.31270 0.99550
summary(conserv[which(polyA$category=="novel antisense - filtered")])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.0000  0.1593  0.5667  0.5476  0.9191  0.9947
## to contranst: the verified transcribed genes show much higher conservation:
summary(conserv[which(polyA$category=="verified gene")])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.0000  0.5373  0.8399  0.7135  0.9632  1.0000
## and the group of untranscribed regions without any features:
summary(conserv[which(polyA$category=="untranscribed" & polyA$overlapFeatAll=="")])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.0000  0.0149  0.0433  0.1648  0.1921  0.9989

## shows slightly lower conservation scores than the novel ones
        
idxNIF <- which(polyA$category == "novel isolated - filtered")
idxNIFCons <- idxNIF[which(conserv[idxNIF]>0.8)]
    
NIFcons <- polyA2[idxNIFCons,c("chr","strand","start","end","level","category","conserv")]
write.table(NIFcons, file.path("~","hubersvn","projects","Rpacks","davidTiling","inst","doc","strongConservedNovelIsolated.txt"),row.names=FALSE,  sep="\t", quote=F)
