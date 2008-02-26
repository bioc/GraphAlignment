\name{Trace}
\alias{Trace}
\title{Calculate trace of a matrix}
\description{
  Calculate the trace of a square matrix.
}
\usage{
Trace(m)
}
\arguments{
  \item{m}{square matrix}
}
\details{
  The return value is the trace of the input matrix.
}
\examples{
  m <- matrix(c(1, 0, 0, 1), 2, 2)
  Trace(m)
}
\author{Joern P. Meier, Michal Kolar, Ville Mustonen, Michael Laessig, and Johannes Berg}
\keyword{misc}
