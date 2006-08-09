
getAllGO = function(x, gff) {
  require("GO")

  ## create environment of ancestors
  e = new.env(hash=TRUE)

  for(j in ls(GOMFANCESTOR))
    assign(j, get(j, GOMFANCESTOR), envir=e)
  for(j in ls(GOBPANCESTOR))
    assign(j, get(j, GOBPANCESTOR), envir=e)
  for(j in ls(GOCCANCESTOR))
    assign(j, get(j, GOCCANCESTOR), envir=e)
  
  stopifnot(length(ls(e))==length(ls(GOMFANCESTOR))+
            length(ls(GOBPANCESTOR))+length(ls(GOCCANCESTOR)))

  stopifnot(all(c("feature", "Name") %in% names(gff)))
  if(!"Ontology_term" %in% names(gff)) {
    stopifnot("attributes" %in% names(gff))
    gff$"Ontology_term" =getAttributeField(gff$attributes, "Ontology_term")
  }
  
  ## for each gene in 'x', get the GO classes
  stopifnot("feature" %in% names(gff))
  whg = which(gff[, "feature"]=="gene")
  mt  = match(x, gff[whg, "Name"])
  stopifnot(!any(is.na(mt)))
  gt = strsplit(gff[whg[mt], "Ontology_term"], split=",")
    
  ## extend by ancestors
  rv = sapply(gt, function(v) {
    if(any(is.na(v))) {
      stopifnot(length(v)==1)
      k = character(0)
    } else {
      k  = mget(v, e, ifnotfound=list(character(0)))
      k  = sort(unique(c(v, unlist(k))))
    }
    return(k)
  })
  names(rv) = x

  ## checks
  stopifnot(all(gff[whg[mt], "Name"] == names(rv)))
  stopifnot(!any(unlist(sapply(rv, duplicated))))

  rv
}
