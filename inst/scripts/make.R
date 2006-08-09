## --------------------------------------------------
## Create the data set davidTiling
## --------------------------------------------------
library("davidTiling")
library("affy")
options(error=recover)

celDir  = system.file("celfiles", package="davidTiling")
outDir  = "../../data"
tmpDir = tempdir()

source("~/madman/Rpacks/tilingArray/R/readCel2eSet.R")

adfData = data.frame(
    filename = I(c("09_11_04_S96_genDNA_16hrs_45C_noDMSO.cel",
      "041119_S96genDNA_re-hybe.cel",
      "041120_S96genDNA_re-hybe.cel",
      "05_04_27_2xpolyA_NAP3.cel",
      "05_04_26_2xpolyA_NAP2.cel",
      "05_04_20_2xpolyA_NAP_2to1.cel",
      "050409_totcDNA_14ug_no52.cel",
      "030505_totcDNA_15ug_affy.cel")),
    nucleicAcid = c(rep("genomic DNA", 3), rep("poly(A) RNA", 3), rep("total RNA", 2)))
rownames(adfData) = adfData$filename

pd = new("AnnotatedDataFrame",
  data=adfData,
  varMetadata = data.frame(
    labelDescription = I(c(filename="Name of the CEL file",
      nucleicAcid = "What is the sample? A factor with three levels: genomic DNA, poly(A) RNA, total RNA"))))

ed = new("MIAME", name="Lior David, Marina Granovskaia, Lars M. Steinmetz",
    lab="Stanford Genome Technology Center; European Molecular Biology Laboratory",
    contact="larsms@embl.de",
    title="A high-resolution map of transcription in the yeast genome",
    abstract="There is abundant transcription from eukaryotic genomes unaccounted for by protein coding genes. A high-resolution genome-wide survey of transcription in a well annotated genome will help relate transcriptional complexity to function. By quantifying RNA expression on both strands of the complete genome of Saccharomyces cerevisiae using a high-density oligonucleotide tiling array, this study identifies the boundary, structure, and level of coding and noncoding transcripts. A total of 85% of the genome is expressed in rich media. Apart from expected transcripts, we found operon-like transcripts, transcripts from neighboring genes not separated by intergenic regions, and genes with complex transcriptional architecture where different parts of the same gene are expressed at different levels. We mapped the positions of 3' and 5' UTRs of coding genes and identified hundreds of RNA transcripts distinct from annotated genes. These nonannotated transcripts, on average, have lower sequence conservation and lower rates of deletion phenotype than protein coding genes. Many other transcripts overlap known genes in antisense orientation, and for these pairs global correlations were discovered: UTR lengths correlated with gene function, localization, and requirements for regulation; antisense transcripts overlapped 3' UTRs more than 5' UTRs; UTRs with overlapping antisense tended to be longer; and the presence of antisense associated with gene function. These findings may suggest a regulatory role of antisense transcription in S. cerevisiae. Moreover, the data show that even this well studied genome has transcriptional complexity far beyond current annotation.",
    url="http://www.pnas.org/cgi/reprint/0601091103v1",
    pubMedIds="16569694")

## uncompress manually, the "compress" option of ReadAffy let's it crash on my computer..
for (f in pd$filename) {
  fin  = file.path(celDir, paste(f, "gz", sep="."))
  fout = file.path(tmpDir, f)
  if(!file.exists(fout)) {
    cmd = paste("gzip -dc", fin, ">", fout)
    cat(cmd, "\n")
    system(cmd)
  }
}

davidTiling = readCel2eSet(adf=pd, path=tmpDir, rotated=TRUE, experimentData=ed)

##
## save it
##
save(davidTiling, file=file.path(outDir, "davidTiling.rda"), compress=TRUE)
