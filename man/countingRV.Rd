\name{countingRV}
\alias{countingRV}

\title{
Calculates RV coefficient for partial abundance matrices.
}
\description{
Calculates RV coefficients of partial abundance matrices by using as reference the complete abundance community matrix. This function is designed to be flexible with regards to the underlying distance used as \code{\link{dist}}, \code{\link{vegdist}}, \code{\link{decostand}} (through the argument \code{transfo}) and \code{\link{designdist}}. 
}
\usage{
countingRV(comm,partAbund,transfo=NULL,dist=NULL,vegdist=NULL,designdist=NULL,logbase=NULL,...)
}

\arguments{
  \item{comm}{A complete abundance community \code{matrix}. If a \code{data.frame} is given instead of a \code{matrix} the function converts \code{comm} into a matrix and sends a warning message.}
  \item{partAbund}{An object of class \code{partAbund}.}
  \item{transfo}{A character string passed to \code{\link{decostand}} defining how the data should be transformed.}
  \item{dist}{A character string passed to \code{\link{dist}} defining the distance to apply to \code{comm} and \code{partAbund} when calculating the RV coefficient.}
  \item{vegdist}{A character string passed to \code{\link{vegdist}} defining the distance to apply to \code{comm} and \code{partAbund} when calculating the RV coefficient.}
  \item{designdist}{A character string presenting an equation of the  dissimilarities passed to \code{\link{designdist}} defining the distance to apply to \code{comm} and \code{partAbund} when calculating the RV coefficient.}
  \item{logbase}{A numeric value defining the log base use for the modified Gower distance. This parameter is used only when \code{transfo} is "log" and \code{vegdist} is "altGower"}
  \item{\dots}{arguments passed to \code{\link{dist}}, \code{\link{vegdist}}, \code{\link{decostand}} (through the argument \code{transfo}) or \code{\link{designdist}} when they are used.}  
}

\details{Only a single argument between \code{dist}, \code{vegdist} and \code{designdist} an be defined. An error message is sent if multiple arguments are defined simultaneously. 

The function allow for \code{transfo} and any one of \code{dist}, \code{vegdist} and \code{designdist} to be used simultaneously. However, combining both pre-transformation and distances should only be carried out by expert users.}

\value{
A vector of RV coefficients for each partial abundance. .
}

\author{
F. Guillaume Blanchet
}

\examples{
data(mite)

### Calculate partial abundance matrix (up to 35 individuals per site per species)
pAbundMite<-partAbund(mite,countThresh=35)

### Calculate RV coefficients for the partial abundance using Hellinger pre-transformation
RVpartAbund<-countingRV(mite, pAbundMite, transfo="hellinger")

### Calculate RV coefficients for the partial abundance using the modified Gower distance with a log base of 2
RVpartAbund<-countingRV(mite, pAbundMite, transfo="log",vegdist="altGower", logbase=2)

}
\keyword{ multivariate }
