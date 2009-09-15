\name{data:davidTiling}
\alias{data:davidTiling}
\alias{davidTiling}
\non_function{}
\title{Dataset of class eSet with the raw 'CEL file' intensities}
\usage{data("davidTiling")}

\format{Intensity data for 8 arrays. The
\code{\link[Biobase]{phenoData}} slot contains the file names
and the nucleic acid type.}

\description{
The data are from an experiment that used Affymetrix
Scerevisiaetiling chips from 2004, which where custom-made for the Stanford
Genome Center. The chips tile the complete genome of S. cerevisiae in
steps of 8 bases, separately for each strand of each chromosome. The
two tiles for one chromosome (Watson and Crick strands) are offset by 4
bases.

Note that the class \code{\link[Biobase:class.eSet]{eSet}} was used instead of
\code{\link[affy:AffyBatch-class]{AffyBatch}} since the additional
overhead of 'CDF environments' in the latter is not needed here.
}

\source{
Lior David and Lars Steinmetz, both from the Stanford
Genome Center. Lars Steinmetz is also at EMBL Heidelberg.
}
\author{W. Huber \email{huber@ebi.ac.uk}}
\examples{
data("davidTiling")
dim(exprs(davidTiling))
}
\keyword{datasets}
