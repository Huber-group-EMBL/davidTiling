\name{GOHyperG}
\alias{GOHyperG}
\title{Hypergeometric test for GO category enrichment.}
\description{
  This function is adapted from the function of the same name in
  \code{GOstats} package. Main difference is that it draws the GO
  annotations of the genes in \code{candidates} from the data.frame
  \code{gff}. It also draws a plot.
}
\usage{
GOHyperG(candidates, gff, plotmain)
}
\arguments{
  \item{candidates}{character vector}
  \item{gff}{see \code{\link{getAllGO}}}
  \item{plotmain}{character, plot title}
}

\details{The elements of \code{x} are matched against the column
  \code{gene} in \code{gff}. All are required to match.
  A list of GO terms is then extracted from the corresponding
  rows in the \code{Ontology_term} column. A gene may be annotated by
  several terms, separated by ",". Then the GO package is used to
  augment this by all ancestor terms.
}

\value{
  List of character vectors.
}

\seealso{\code{\link{getAllGO}}}
\author{W. Huber <huber@ebi.ac.uk>}

\examples{
}
\keyword{manip}
