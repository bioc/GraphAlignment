\name{VectorToBin}
\alias{VectorToBin}
\title{Convert vector to bin representation}
\description{
  Convert a vector of real numbers to a vector of integer bin numbers 
  according to the specified lookup vector.
}
\usage{
VectorToBin(v, lookup, clamp=TRUE)
}
\arguments{
  \item{v}{input vector}
  \item{lookup}{lookup vector}
  \item{clamp}{clamp values to lookup range}
}
\value{
  The return value is a vector containing the bin numbers for the elements of the input vector.
}
\details{
  This function transforms a vector of arbitrary values into a vector of bin  numbers corresponding to the data in the input vector. Bin numbers are found using the specified lookup table (see \link{GetBinNumber}).
}
\examples{
  x <- rnorm(10)
  lookup <- 0.2*(-5:5)
  VectorToBin(x, lookup)
}
\author{Joern P. Meier, Michal Kolar, Ville Mustonen, Michael Laessig, and Johannes Berg}
\keyword{misc}
