\name{Permute}
\alias{Permute}
\title{Permute rows and columns of a matrix}
\description{
  Permute rows and columns of the input matrix.
}
\usage{
Permute(m, p, invertp=FALSE)
}
\arguments{
  \item{m}{matrix}
  \item{p}{permutation vector}
  \item{invertp}{apply inverse permutation}
}
\value{
  The return value is the row- and columnwise Permutation of the elements of the input matrix, so Permute(m, p)[i, j] is equal to m[p[i], p[j]].
}
\details{
  This function permutes rows and columns of a matrix using the specified permutation vector. The inverse of the permutation will be applied if the invertp argument is set to TRUE.
}
\examples{
  example <- GenerateExample(10, 10, 1, 0.5, TRUE, 4)
  b <- example$b
  p <- rank(rnorm(10))
  Permute(b, p, TRUE)
}
\author{Joern P. Meier, Michal Kolar, Ville Mustonen, Michael Laessig, and Johannes Berg}
\keyword{misc}
