\name{partAbund}
\alias{partAbund}

\title{
Partial abundance community matrix
}
\description{
Calculates all partial abundance community matrices following the approach proposed by Blanchet et al. (in press)
}
\usage{
partAbund(comm,countThresh=NULL)
}

\arguments{
  \item{comm}{A community \code{matrix}. If a \code{data.frame} is given instead of a \code{matrix} the function converts \code{comm} into a matrix and sends a warning message.}
  \item{countThresh}{Counting threshold. The maximum number of individuals to count by site by species. If the counting threshold is not given, community matrices with all different counting thresholds are calculated. If the value given is not an integer, this value is rounded to the lowest integer}
}

\value{
An object of class \code{partAbund}. It is essentially an array where each third dimension is community matrix with a particular counting threshold. The community matrices in this array are ordered sequentially from lowest (1) to the highest counting threshold considered.
}

\author{
F. Guillaume Blanchet
}

\examples{
data(mite)

pAbundMite<-partAbund(mite,countThresh=35)
str(pAbundMite)

}
\keyword{ multivariate }
