\name{scatterWithHist}
\alias{scatterWithHist}
\title{Scatterplot with histograms of marginal distributions}
\description{
  Scatterplot with histograms of marginal distributions.
}
\usage{
scatterWithHist(x, breaks, barcols, xlab, ylab, \dots) 
}
\arguments{
  \item{x}{numeric matrix with 2 columns.}
  \item{breaks}{numeric vector with histogram breaks, see \code{\link{hist}}.}
  \item{barcols}{character vector of length 2, colors for the histogram filling.}
  \item{xlab}{characterof length 1, label for x-axis.}
  \item{ylab}{characterof length 1, label for y-axis.}
  \item{\dots}{further arguments that get passed on to \code{\link{plot}}.}
}

\details{
}

\value{
  The function is called for its side effect.
}

\author{W. Huber <huber@ebi.ac.uk>}

\examples{
  x = rexp(100)
  x = cbind(x, x+0.6*rnorm(length(x)))
  scatterWithHist(x,
     breaks=seq(min(x),max(x),length=20),
     barcols=c("mistyrose", "lightblue"),
     xlab="Daffodil", ylab="Petunia", pch=16)
}
\keyword{hplot}
