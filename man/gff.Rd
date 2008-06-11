\name{data:gff}
\alias{data:gff}
\alias{gff}
\non_function{}
\title{Genomic features of Saccharomyces cerevisiae}
\usage{data("gff")}
\description{A data frame with genomic features of Saccharomyces cerevisiae.}

\format{Object of class \code{data.frame}.
GFF is a file format for annotating genomes, see <insert the URL to
the documentation page for GFF at Sanger here>.
The format is essentially a rectangular table, and here
it is represented as a data frame.
}

\source{
Two GFF files were downloaded:
\code{saccharomyces_cerevisiae.gff} from
\code{ftp://genome-ftp.stanford.edu/pub/yeast/data_download/chromosomal_feature}
on 7 Aug 2005, 18:16 BST, and
\code{IGR_v24.2.p001.allowoverlap.GFF} from
\code{http://jura.wi.mit.edu/fraenkel/download/release_v24/GFF}
upon suggestion from the SGD curators on 30 Aug 2005. (Future versions of
SGD's GFF files are likely to include the latter as well).
They were parsed, combined and written into the
\code{gff} data.frame with the script
\code{makeProbeAnno.R} in the \code{inst/scripts} directory of this package.
}

\author{W. Huber \email{huber@ebi.ac.uk}}
\examples{
data("gff")
str(gff)
}
\keyword{datasets}
