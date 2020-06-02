\name{commSimul}
\alias{commSimul}

\title{
Community simulation function
}
\description{
Simulate a community following the procedure proposed in Blanchet et al. (in press).
}
\usage{
commSimul(nquad,nsp, nind,patchrad=0.2,SAD="lognormal",sdlog=5)
}

\arguments{
  \item{nquad}{Number of quadrates in which to divide the sampling area. This should be equal to any integer squared. (See Details).}
  \item{nsp}{Integer defining the number of species to be simulated.}
  \item{nind}{Integer defining the maximum number of individuals for each species. Note that if SAD = "uniform", this value will be the same for each species. (See Details).}
  \item{patchrad}{A numeric value defining the maximum radius of a patch for a group of individuals of a single species. (See Details).}
  \item{SAD}{Character string defining the type species-abundance distribution for the community. Either "lognormal", "bstick", "uniform".}
  \item{sdlog}{This argument is active only if SAD= "lognormal". It defines the standard deviation of the lognormal distribution.}
}

\details{The number of quadrates in the sampling area (a square of unit size) is defined by dividing the area horizontally and vertically so that each quadrate is proportionally equal to the size of the sampling area. For this reason, nquad should be an integer that can be obtain from any solution of x^2. If it is not the case, a value will be calculated that approximate the desired value by rounding the square root nquad and putting it to the power of 2.

Depending on the way individuals are distributed in the sampling area (spatagg), the number of individuals per species defined by 'nind' may varies because the method use to sample individual relies on random point patterns and the value given is the average of the random point process. 

Note that when defining the patch radius through patchrad, the patch radius is sampled from a uniform distribution ranging from 0.01 to patchrad. Also, since all simulations are carried out in a square unit size, defining patchrad with a value larger than 1 is ill advised. The argument patchrad defines patch radius for individuals of each species simulated so that the patch size associated to different species within the same community can vary within the predefine range.
}

\value{
A sites (rows) by species (columns) community matrix.
}

\author{
F. Guillaume Blanchet
}

\examples{
comm1<-commSimul(10^2,20, 2000,patchrad=0.2,SAD="lognormal",sdlog=5)
}
\keyword{ multivariate }
