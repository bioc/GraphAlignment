\name{InvertPermutation}
\alias{InvertPermutation}
\title{Invert permutation}
\description{
  Invert the permutation defined by the specified vector.
}
\usage{
InvertPermutation(p)
}
\arguments{
  \item{p}{permutation vector}
}
\details{
  The return value is a permutation vector equal to the inverse of the input permutation.
}
\examples{
  p <- rank(rnorm(10))
  pInv <- InvertPermutation(p)
}
\author{Joern P. Meier, Michal Kolar, Ville Mustonen, Michael Laessig, and Johannes Berg}
\keyword{misc}
