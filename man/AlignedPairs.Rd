\name{AlignedPairs}
\alias{AlignedPairs}
\title{Get aligned node pairs}
\description{
  Create a matrix of pairs of aligned nodes from networks A and B using the permutation vector P, where P is in the format returned by AlignNetworks().
}
\usage{
AlignedPairs(A, B, P)
}
\arguments{
  \item{A}{adjacency matrix for network A}
  \item{B}{adjacency matrix for network B}
  \item{P}{permutation vector to be used as the alignment}
}
\value{
  The return value is a matrix with two columns. The number of rows is equal to the number of aligned node pairs. Each row in the matrix denotes a pair of aligned nodes. In each row, the first element (index 1) is the label of a node in network A, and the second element (index 2) is the label of a node in network B.
}
\details{
  This function creates a matrix containing pairs of aligned nodes from networks A and B using the permutation vector P, where P is in the format returned by \link{AlignNetworks}.
}
\examples{
  ex<-GenerateExample(dimA=22, dimB=22, filling=.5, covariance=.6,
    symmetric=TRUE, numOrths=10, correlated=seq(1,18))
  
  pinitial<-InitialAlignment(psize=34, r=ex$r, mode="reciprocal")
  
  lookupLink<-seq(-2,2,.5)
  linkParams<-ComputeLinkParameters(ex$a, ex$b, pinitial, lookupLink)
  
  lookupNode<-c(-.5,.5,1.5)
  nodeParams<-ComputeNodeParameters(dimA=22, dimB=22, ex$r,
    pinitial, lookupNode)
  
  al<-AlignNetworks(A=ex$a, B=ex$b, R=ex$r, P=pinitial,
    linkScore=linkParams$ls,
    selfLinkScore=linkParams$ls,
    nodeScore1=nodeParams$s1, nodeScore0=nodeParams$s0,
    lookupLink=lookupLink, lookupNode=lookupNode,
    bStart=.1, bEnd=30,
    maxNumSteps=50)
  
  alignedPairs<-AlignedPairs(A=ex$a, B=ex$b, al)
}
\author{Joern P. Meier, Michal Kolar, Ville Mustonen, Michael Laessig, and Johannes Berg}
\keyword{misc}
