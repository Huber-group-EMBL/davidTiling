
##
## Categorize segments 
##
allncRNA = c("ncRNA","snoRNA","snRNA","tRNA","rRNA")
  
categorizeSegments = function(env, minNewSegmentLength=48, zThresh=0,
  maxDuplicated=0.5, zThreshUTR=2 ) {

  require(tilingArray)
  data(yeastFeatures)

  s = get("segScore", env)
  stopifnot(is(s, "data.frame"))

  
  feat1 = c("verified gene", "uncharacterized gene", "dubious gene")
  feat2 = rev(rownames(yeastFeatures)[yeastFeatures$isTranscribed]) ## trumping order
  
  ## results data structure: a factor which assigns a category to each segment:
  overlap = factor(rep(NA, nrow(s)),
    levels = c("<50%", ">=50%", "complete"))
  
  catg = factor(rep(NA, nrow(s)), 
    levels = c(feat1, feat2, 
      "novel isolated - filtered", "novel isolated - unassigned",
      "novel antisense - filtered", "novel antisense - unassigned",
      "excluded", "untranscribed"))

  simpleCategories = c("annotated ORF", "ncRNA(all)", 
    "novel isolated - filtered",  "novel isolated - unassigned",
    "novel antisense - filtered", "novel antisense - unassigned",
    "excluded", "untranscribed", "dubious gene")

  ## Step 1: frac.dup; exclude features with too many non-unique probes
  sel = s[,"frac.dup"] >= maxDuplicated
  catg[ sel ] = "excluded"
    
  ## Step 2: identify untranscribed regions
  sel = (is.na(catg) & s[,"level"] < 0)
  catg[ sel ] = "untranscribed"
  s$isUnIso = (sel & (s[, "overlapFeatAll"]==""))
  
  ## step 3: identify already annotated regions (misses some so far!)
  wh  = which(is.na(catg))
  attrName = c("<50%"  = "overlappingFeature",
               ">=50%" = "mostOfFeatureInSegment",
               "complete"  = "featureInSegment")

  categIDs = vector(mode="list", length = length(feat1)+length(feat2))
  names(categIDs) = c(feat1, feat2)
  
  for(f in feat2)
    categIDs[[f]] = gff[ gff[, "feature"]==f, "Name" ]

  sel = gff[, "feature"]=="gene"
  categIDs[["dubious gene"]]         = gff[ sel & gff[, "orf_classification"]=="Dubious", "Name"]
  categIDs[["uncharacterized gene"]] = gff[ sel & gff[, "orf_classification"]=="Uncharacterized", "Name"]
  categIDs[["verified gene"]]        = gff[ sel & gff[, "orf_classification"]=="Verified", "Name"]

  stopifnot(all(listLen(categIDs)>0))
  # categIDs <- categIDs[-which(listLen(categIDs)<1)] 
  #ignore non-found features, seems to apply only to uORFs so far

  ## Loop over <50%, >=50%, complete:
  for(i in seq(along=attrName)) {
    ovF = strsplit(s[wh, attrName[i]],  split=", ")

    ## Loop over three annotation classes
    for(j in rev(seq(along=categIDs))) {
      ## find features that this segment is contained in
      sel = sapply(ovF, function(x) any(x %in% categIDs[[j]]))
      if(any(sel)) {
        catg[wh[sel]]    = names(categIDs)[j]
        overlap[wh[sel]] = names(attrName)[i]
      }
    } ## for j
  } ## i

  ## step 4: novelty filter
  zmin = pmin(s[, "z3"], s[, "z5"])
  xOp  = s[,"oppositeExpression"]
  filt1 = (is.na(zmin) | (zmin <zThresh))
  filt2 = (s[,"length"] < minNewSegmentLength)
  filt3 = (!is.na(xOp) & (xOp > s[,"level"]))
  filt  = (filt1|filt2|filt3)
  stopifnot(!any(is.na(filt)))
  
  ## step 5: novel - isolated or antisense
  iso  = (s[,"oppositeFeature"]=="")
  isna = is.na(catg)

  cat("Novelty filter: Considering ", sum(isna), " segments.\n",
    "1. z-scores < ", zThresh, ": ", sum(isna&filt1), "\n",
    "2. length < ", minNewSegmentLength, ": ", sum(isna&filt2), "\n",
    "3. oppositeExpression > segment level: ", sum(isna&filt3), "\n",
    "Rejected by (1 or 2 or 3): ", sum(isna&filt), ".\n\n", sep="")

  catg[isna &  iso & !filt ] = "novel isolated - filtered"
  catg[isna &  iso &  filt ] = "novel isolated - unassigned"
  catg[isna & !iso & !filt ] = "novel antisense - filtered"
  catg[isna & !iso &  filt ] = "novel antisense - unassigned"

  stopifnot(!any(is.na(catg)))
  s$category = catg
  s$overlap  = overlap

  ## simpleCategory
  simc = factor(rep(NA, nrow(s)), levels=simpleCategories)
  simc[ catg %in% c("uncharacterized gene", "verified gene")] = "annotated ORF"
  simc[ catg %in% c(allncRNA)]  = "ncRNA(all)"
  for(lev in simpleCategories[-(1:2)])
    simc[ s[,"category"]==lev] = lev
  s$simpleCatg = simc

  ## piechart category
  pieNames = c(A="overlap >=50%", B="overlap <50%",
  C="novel isolated - filtered", D="novel isolated - unassigned",
  E="novel antisense - filtered", F="novel antisense - unassigned")
  
  pc = factor(rep(NA, nrow(s)), levels=pieNames)
  pc[overlap  %in% c(">=50%", "complete")] = "overlap >=50%"
  pc[overlap  %in% c("<50%")] = "overlap <50%"
  stopifnot(all(levels(pc)[3:6] %in% levels(catg)))
  for(k in levels(pc)[3:6])
    pc[catg == k] = k

  s$pieCat = pc


  ##
  ## for the UTR mapping
  ##
  zmin = pmin(s[, "z3"], s[, "z5"])
  hasGoodFlanks = (!is.na(zmin) & (zmin >= zThreshUTR))

  ## annotated ORF, z-score criterion, nuclear
  stopifnot("annotated ORF" %in% levels(simc))
  sel = ((simc=="annotated ORF") & hasGoodFlanks &
         !is.na(s[, "utr3"]) & !is.na(s[, "utr5"]) & (s[,"chr"] <= 16))
  s$goodUTR = ifelse(sel, zmin, as.numeric(NA))

  
  return(s)
}




##### Start of "main" ####


## these are the hybe sets that we care about in the paper:
rnaTypes = c("seg-polyA-050909", "seg-tot-050909")[1:2]

source(scriptsDir("readSegments.R"))
source(scriptsDir("writeSegmentTable.R"))
source(scriptsDir("calcThreshold.R"))

##
## CATEGORIZE
##
## cs: tables of all segments
## csu: tables only of those segments that have UTR mapping
##      (is redundant with cs, but we have it for convenience) 
## utr: matrix with the UTR length info
##      (is redundant with csu, but we have it for convenience)
##
if(!exists("cs")) {
  cs = csu = vector(mode="list", length=length(rnaTypes))
  names(cs) = names(csu) = rnaTypes
  utr = vector(mode="list", length=length(rnaTypes)+1)
  names(utr) = c(rnaTypes, "combined")

  cat("\n\nCategorization of segments:\n",
      "===========================\n", sep="")
  for(rt in rnaTypes) {
    cat(rt, ":\n", sep="")
    s = categorizeSegments(get(rt))
    cs[[rt]] =s

    s = s[!is.na(s[,"goodUTR"]), ]
    z = as.matrix(s[, c("utr5", "utr3")])
    rownames(z) = rownames(s) = s[, "featureInSegment"]
    colnames(z) = c("5' UTR", "3' UTR")
    csu[[rt]] = s
    utr[[rt]] = z
  }
  ##
  ## UTR lengths combined from poly-A and total RNA samples
  ##
  utr[["combined"]] = rbind(utr[[1]], utr[[2]][ !(rownames(utr[[2]])%in%rownames(utr[[1]])), ])
  save(utr, file=paste("utr-", date(), ".rda", sep=""))
} else {
  cat("\n**************************************************\n",
        "*      NOT REDOING categorizeSegments            *\n",
        "**************************************************\n", sep="")
}

##
##  define 'featNames' list
##

isGene = ((gff[,"feature"]=="gene") & (gff[, "orf_classification"] %in% c("Verified", "Uncharacterized")))
featNames = list("annotated ORFs" = unique(gff[ isGene, "Name"]),
  "ncRNA(all)" = unique(gff[ gff[, "feature"] %in% allncRNA, "Name"]))





