\name{AnalyzeAlignment}
\alias{AnalyzeAlignment}
\title{Analyze an alignment}
\description{
  Analyze an alignment between two networks.
}
\usage{
AnalyzeAlignment(A, B, R, P, lookupNode, epsilon, clamp)
}
\arguments{
  \item{A}{adjacency matrix for network A}
  \item{B}{adjacency matrix for network B}
  \item{R}{node similarity matrix}
  \item{P}{permutation vector}
  \item{lookupNode}{node bin lookup vector}
  \item{epsilon}{node similarity threshold}
  \item{clamp}{clamp values to range when performing bin lookups}
}
\value{
The return value is a list containing the results. Defined values are:
  
na - number of aligned node pairs
  
nb - number of aligned node pairs where neither partner has appreciable sequence similarity with any node in the other network. Appreciable sequence similarity means that r>epsilon. Formally, this is the number of aligned node pairs (ia, ib), where no jb exists such that R[ia, jb] > epsilon and no ja such that R[ja, ib] > epsilon.

nc - number of aligned node pairs, where the partners have no appreciable sequence similarity, but one or both of them has appreciable sequence similarity with some other node in the other network. Formally, this is the number of aligned node pairs (ia, ib) with R[ia, ib] < epsilon but jb or ja exists, such that R[ia, jb] > epsilon or R[ja, ib] > epsilon.
}
\details{
  This function analyzes an alignment and returns various characteristics.
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
  
  AnalyzeAlignment(A=ex$a, B=ex$b, R=ex$r, P=al, lookupNode,
    epsilon=.5)
}
\author{Joern P. Meier, Michal Kolar, Ville Mustonen, Michael Laessig, and Johannes Berg}
\keyword{misc}
