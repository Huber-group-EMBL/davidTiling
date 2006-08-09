\name{getAllGO}
\alias{getAllGO}
\title{Get all GO categories for a list of genes.}
\description{
  The function uses the GO categories in the data.frame \code{gff} to
  obtain annotated GO categories, then the GO**ANCESTOR data in the
  \code{GO} package to add all parent terms as well.
}
\usage{
getAllGO(x, gff) 
}
\arguments{
  \item{x}{character vector.}
  \item{gff}{data.frame with columns \code{feature}, \code{Name}, and
    (\code{Ontology_term} or \code{attributes}), see details}
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

\author{W. Huber <huber@ebi.ac.uk>}

\examples{
}
\keyword{manip}
