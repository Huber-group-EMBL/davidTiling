\name{showDens}
\alias{showDens}
\title{Plot function for more than one density}
\description{
  Plot function for more than one density
}
\usage{
showDens(z, breaks, xat, xtickLabels=paste(xat), col, ylab = "", ...)
}

\arguments{
  \item{z}{List: numeric vectors for computing histograms for}
  \item{breaks}{Numeric vector: breaks of the histogram}
  \item{xat}{Numeric vector: where to put the x-axis ticks}
  \item{xtickLabels}{Character vector: what to write underneath them} 
  \item{col}{Character vector: colours of the histograms}
  \item{ylab}{Character scalar: y-axis label}
  \item{\dots}{futher arguments passed on to plot}
}
\details{
  ...
}
\value{
  returns scale factor
}
\author{Wolfgang Huber \email{huber@ebi.ac.uk}}
\seealso{\code{\link{hist}}}
\examples{

showDens(list(x1=runif(100), x2=exp(runif(50))-1, x3=runif(20)),
  breaks=seq(0, 2, 0.2), xat=seq(0, 2, 0.5), col=rainbow(3), xlab="Random Numbers")

}
\keyword{hplot}
