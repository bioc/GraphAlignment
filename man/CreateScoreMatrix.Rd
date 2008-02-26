\name{CreateScoreMatrix}
\alias{CreateScoreMatrix}
\title{Create score matrix}
\description{
  Create a score matrix from two lookup vectors.
}
\usage{
CreateScoreMatrix(lookupX, lookupY)
}
\arguments{
  \item{lookupX}{x lookup vector}
  \item{lookupY}{y lookup vector}
}
\value{
  The return value is a matrix of link scores.
}
\details{
  This function creates a very simple score matrix containing the product of lookup table values for each row and column as its elements. This can be used for testing purposes.
}
\examples{
  lookupLink <- 0.1*(-10:10)
  CreateScoreMatrix(lookupLink, lookupLink)
}
\author{Joern P. Meier, Michal Kolar, Ville Mustonen, Michael Laessig, and Johannes Berg}
\keyword{misc}
