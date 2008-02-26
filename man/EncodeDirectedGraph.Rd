\name{EncodeDirectedGraph}
\alias{EncodeDirectedGraph}
\title{Encode directed graph}
\description{
  Encode an adjacency matrix for a directed graph into a symmetric matrix.
}
\usage{
EncodeDirectedGraph(matrix, P)
}
\arguments{
  \item{matrix}{square matrix}
  \item{P}{permutation vector}
}
\value{
  The return value is a symmetric matrix representing the encoded input matrix.
}
\details{
  This function encodes an adjacency matrix for a directed graph into a symmetric matrix. Currently only binary directed graphs are implemented. The adjacency matrix of a binary directed graph has elements 0, 1. The same graph can be represented by a symmetric adjacency matrix with elements -1, 0, 1, with the sign of the entry indicating the direction of the link. The result is obtained by setting entries (i, j) and (j, i) of the target matrix m' to 1 if entry m[i, j] = 1 and P[i] > P[j] and to -1 if m[i, j] = 1 and P[j] > P[i].
}
\examples{
  dg <- matrix(c(0, 1, 0, 0, 1, 1, 1, 0, 1), 3, 3)
  EncodeDirectedGraph(dg, 1:dim(dg)[1])
}
\author{Joern P. Meier, Michal Kolar, Ville Mustonen, Michael Laessig, and Johannes Berg}
\keyword{misc}
