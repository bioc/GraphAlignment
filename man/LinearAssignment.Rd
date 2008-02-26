\name{LinearAssignment}
\alias{LinearAssignment}
\title{Solve linear assignment problem}
\description{
  Solve the linear assignment problem specified by the cost matrix.
}
\usage{
LinearAssignment(matrix)
}
\arguments{
  \item{matrix}{cost matrix}
}
\value{
  The return value is a permutation vector equal to the solution of the linear assignment problem specified by the cost matrix. The result is the permutation P for which MP is minimal (where M is the cost matrix).
}
\details{
  This function solves the linear assignment problem defined by the input matrix. [Jonker, Volgenant 1987]
}
\examples{
  m <- matrix(rnorm(25), 5, 5)
  px <- LinearAssignment(round(-1000 * (m / max(abs(m)))))
  m[px,]
}
\references{
  Jonker, R & Volgenant, A. (1987) Computing 38, 325--340.
}
\author{Joern P. Meier, Michal Kolar, Ville Mustonen, Michael Laessig, and Johannes Berg}
\keyword{misc}
