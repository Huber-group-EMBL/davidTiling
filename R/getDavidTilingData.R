#' Create the data set davidTiling.
#' 
#' In version older than 1.53.0, this dataset was provided as a pre-built 
#' ExpressionSet object, available via `data(davidTiling)`. 
#' With this function, the dataset is built on-the-fly by:
#' - downloading the raw CEL files from ArrayExpress (E-TABM-14)
#' - reading the CEL files into an ExpressionSet object
#'
#' @returns An ExpressionSet object containing the davidTiling data.
#' 
#' @references 
#' David L, Huber W, Granovskaia M, et al. A high-resolution map of
#' transcription in the yeast genome. Proceedings of the National Academy of
#' Sciences of the United States of America. 2006 Apr;103(14):5320-5325.
#' \doi{10.1073/pnas.0601091103}. PMID: 16569694; PMCID: PMC1414796.
#' 
#' @seealso
#' \url{https://www.ebi.ac.uk/arrayexpress/files/E-TABM-14/}
#' 
#' @examples
#' getDavidTilingData()
#' 
#' @export
getDavidTilingData <- function() {
  ## Modified from inst/scripts/make.R

  celDir <- tools::R_user_dir("tilingArray", which = "cache")
  if (!dir.exists(celDir)) {
    dir.create(celDir, recursive = TRUE)
  }

  required_cel_files <- c(
    "09_11_04_S96_genDNA_16hrs_45C_noDMSO.cel",
    "041119_S96genDNA_re-hybe.cel",
    "041120_S96genDNA_re-hybe.cel",
    "05_04_27_2xpolyA_NAP3.cel",
    "05_04_26_2xpolyA_NAP2.cel",
    "05_04_20_2xpolyA_NAP_2to1.cel",
    "050409_totcDNA_14ug_no52.cel",
    "030505_totcDNA_15ug_affy.cel"
  )

  for (f in required_cel_files) {
    if (!file.exists(file.path(celDir, f))) {
      message("Downloading ", f, " from ArrayExpress...")
      utils::download.file(
        url = paste0(
          "https://www.ebi.ac.uk/arrayexpress/files/E-TABM-14/",
          f
        ),
        destfile = file.path(celDir, f)
      )
    } else {
      message(f, " already exists in ", celDir, ", skipping download.")
    }
  }

  adfData <- data.frame(
    filename = required_cel_files,
    nucleicAcid = rep.int(c("genomic DNA", "poly(A) RNA", "total RNA"), c(3, 3, 2)),
    row.names = required_cel_files
  )

  pd <- methods::new(
    "AnnotatedDataFrame",
    data = adfData,
    varMetadata = data.frame(
      labelDescription = c(
        filename = "Name of the CEL file",
        nucleicAcid = "What is the sample? A factor with three levels: genomic DNA, poly(A) RNA, total RNA"
      )
    )
  )

  ed <- methods::new(
    "MIAME",
    name = "Lior David, Marina Granovskaia, Lars M. Steinmetz",
    lab = "Stanford Genome Technology Center; European Molecular Biology Laboratory",
    contact = "larsms@embl.de",
    title = "A high-resolution map of transcription in the yeast genome",
    abstract = "There is abundant transcription from eukaryotic genomes unaccounted for by protein coding genes. A high-resolution genome-wide survey of transcription in a well annotated genome will help relate transcriptional complexity to function. By quantifying RNA expression on both strands of the complete genome of Saccharomyces cerevisiae using a high-density oligonucleotide tiling array, this study identifies the boundary, structure, and level of coding and noncoding transcripts. A total of 85% of the genome is expressed in rich media. Apart from expected transcripts, we found operon-like transcripts, transcripts from neighboring genes not separated by intergenic regions, and genes with complex transcriptional architecture where different parts of the same gene are expressed at different levels. We mapped the positions of 3' and 5' UTRs of coding genes and identified hundreds of RNA transcripts distinct from annotated genes. These nonannotated transcripts, on average, have lower sequence conservation and lower rates of deletion phenotype than protein coding genes. Many other transcripts overlap known genes in antisense orientation, and for these pairs global correlations were discovered: UTR lengths correlated with gene function, localization, and requirements for regulation; antisense transcripts overlapped 3' UTRs more than 5' UTRs; UTRs with overlapping antisense tended to be longer; and the presence of antisense associated with gene function. These findings may suggest a regulatory role of antisense transcription in S. cerevisiae. Moreover, the data show that even this well studied genome has transcriptional complexity far beyond current annotation.",
    url = "http://www.pnas.org/cgi/reprint/0601091103v1",
    pubMedIds = "16569694"
  )

  tilingArray::readCel2eSet(
    adf = pd,
    path = celDir,
    rotated = TRUE,
    experimentData = ed
  )
}