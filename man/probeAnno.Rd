\name{data:probeAnno}
\alias{data:probeAnno}
\alias{probeAnno}
\non_function{}
\title{An environment with probe mapping information for the Scerevisiae tiling array}
\usage{data("probeAnno")}

\description{
The environment contains probe mapping information for the Affymetrix
Scerevisiaetiling chip from 2004, which was custom-made for the Stanford
Genome Center. The chips tile the complete genome of S. cerevisiae in
steps of 8 bases, separately for each strand of each chromosome. The
two tiles for one chromosome (Watson and Crick strands) are offset by 4
bases.

In the following a brief description of the 138 elements of the
\code{probeAnno} environment.

\code{probeReverse}: a list of 8 factors, each of length 6553600,
corresponding to the rows of \code{\link{davidTiling}}. For
example, if the probe corresponding to the j-th row
in \code{\link{davidTiling}} maps to the
coding sequence of a gene, then the factor level of
\code{probeRevers$CDS[j]} is the name of the gene, and the empty string
\code{""} otherwise. This applies to samples that were hybridized to
the chip after a reverse transcription step.

\code{probeDirect}: analogous to \code{probeReverse}, but for samples
that were hybridized to the chip without a reverse transcription step.
The probes map to the opposite chromosomal strand compared to
experiments with  reverse transcription.

\code{1.+.index}: indices (from 1...6553600, corresponding to the
rows of \code{\link{davidTiling}}) of probes mapping to the Watson
strand of chromosome 1.

\code{1.+.start}, \code{1.+.end}: start and end positions in genomic
coordinates of the alignments of the probes (in the same order as in
\code{1.+.index}) to the Watson strand of chromosome 1. For 25-mers,
the values in \code{1.+.end} are those in \code{1.+.start} plus 24, but
not all probes on the array are 25-mers.

\code{1.-.unique}: specificity of the probe:
\describe{
\item{0}{has exactly one perfect match (PM) and no near-matches in the genome}
\item{1}{has exactly one PM and some near-matches}
\item{2}{has no PM but one or more near-matches}
\item{3}{has multiple PMs in the genome}}
 
\code{1.-.index}, \code{1.-.start}, \code{1.-.end}, \code{1.-.unique}:
analogous to the above, but for the Crick strand of chromosome 1.
\code{2.+.index}, \code{2.+.start}, \code{2.+.end}, \code{2.-.unique}:
analogous to the above, but for the Watson strand of chromosome 2;
and so forth. "Chromosome 17" is mitochondrial DNA.
}

\source{
Probe sequences were obtained from Affymetrix in a file called
\code{S.cerevisiae_tiling.1lq}.
The genomic sequences of the S. cerevisiae chromosomes were
downloaded  from
\code{ftp://genome-ftp.stanford.edu/pub/yeast/data_download/sequence/genomic_sequence/chromosomes/fasta}
on 7 Aug 2005, 18:16 BST in 17 files
\code{chr01.fsa}--\code{chr16.fsa}, and \code{chrmt.fsa}.
The probe sequences were matched against the chromosomal sequencences
with the program \code{MUMmer}, see the script
\code{mapProbesToGenome.sh} (in the \code{inst/scripts} directory of this
package). \code{MUMmer} results were parsed and processed into the
\code{probeAnno} environment with the script
\code{makeProbeAnno.R} (in the \code{inst/scripts} directory of this package).
}
\author{W. Huber \email{huber@ebi.ac.uk}}
\examples{
data("probeAnno")
ls(probeAnno)
str(probeAnno$"1.+.start")
}
\keyword{datasets}
