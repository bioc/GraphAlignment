\name{MatrixToBin}
\alias{MatrixToBin}
\title{Convert matrix to bin representation}
\description{
  Convert a matrix of real numbers to a vector of integer bin numbers according to the specified lookup vector.
}
\usage{
MatrixToBin(M, lookup, clamp=TRUE)
}
\arguments{
  \item{M}{input matrix}
  \item{lookup}{lookup vector}
  \item{clamp}{clamp values to lookup range}
}
\value{
  The return value is a matrix containing the bin numbers for the elements of the input matrix.
}
\details{
  This function transforms a matrix of arbitrary values into a matrix of bin numbers corresponding to the data in the input matrix. Bin numbers are found using the specified lookup table (see \link{GetBinNumber}).
}
\examples{
  m <- matrix(rnorm(25), 5, 5)
  lookup <- 0.2*(-5:5)
  MatrixToBin(m, lookup)
}
\author{Joern P. Meier, Michal Kolar, Ville Mustonen, Michael Laessig, and Johannes Berg}
\keyword{misc}
