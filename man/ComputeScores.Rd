\name{ComputeScores}
\alias{ComputeScores}
\title{Compute scores}
\description{
  Compute scores.
}
\usage{
ComputeScores(A, B, R, P, linkScore, selfLinkScore, nodeScore1,
  nodeScore0, lookupLink, lookupNode, symmetric=TRUE, clamp=TRUE)
}
\arguments{
  \item{A}{adjacency matrix for network A}
  \item{B}{adjacency matrix for network B}
  \item{R}{node similarity matrix}
  \item{P}{permutation vector to be used as the initial alignment (see \link{InitialAlignment})}
  \item{linkScore}{link score matrix (see \link{ComputeLinkParameters})}
  \item{selfLinkScore}{self link score matrix (see \link{ComputeLinkParameters})}
  \item{nodeScore1}{node score vector (s1) (see \link{ComputeNodeParameters})}
  \item{nodeScore0}{node score vector for unaligned nodes (s0) (see \link{ComputeNodeParameters})}
  \item{lookupLink}{link bin lookup table (see \link{GetBinNumber})}
  \item{lookupNode}{node bin lookup table (see \link{GetBinNumber})}
  \item{symmetric}{network symmetry flag}
  \item{clamp}{clamp values to range when performing bin lookups}
}
\value{
  The return value is a list containing the link score (sl) and the node score (sn).
}
\details{
  This function computes log-likelihood scores for an alignment using the specified scoring tables, two networks A and B and their alignment P. The total score of the alignment has two contributions, the first coming from the sequence homology (node similarity, sn) and the second from the similarity of interaction networks (sl).
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
  
  ComputeScores(A=ex$a, B=ex$b, R=ex$r, P=al,
    linkScore=linkParams$ls,
    selfLinkScore=linkParams$ls,
    nodeScore1=nodeParams$s1, nodeScore0=nodeParams$s0,
    lookupLink=lookupLink, lookupNode=lookupNode,
    symmetric=TRUE)
}
\author{Joern P. Meier, Michal Kolar, Ville Mustonen, Michael Laessig, and Johannes Berg}
\keyword{misc}
