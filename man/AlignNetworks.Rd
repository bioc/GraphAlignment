\name{AlignNetworks}
\alias{AlignNetworks}
\title{Align networks}
\description{
  Align networks A and B.
}
\usage{
AlignNetworks(A, B, R, P, linkScore, selfLinkScore, nodeScore1,
  nodeScore0, lookupLink, lookupNode, bStart, bEnd, maxNumSteps, clamp=TRUE, 
  directed=FALSE)
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
  \item{bStart}{start scaling value for simulated annealing}
  \item{bEnd}{end scaling value for simulated annealing}
  \item{maxNumSteps}{maximum number of steps}
  \item{clamp}{clamp values to range when performing bin lookups}
  \item{directed}{whether input networks should be treated as directed graphs}
}
\value{
  The return value is a permutation vector p which aligns nodes from network a with nodes from network B (including dummy nodes). The returned permutation should be read in the following way: the node i in the network A is aligned to  that node in the network B which label is at the i-th position of the permutation vector p. If the label at this position is larger than the size of the network B, the node i is not aligned.
}
\details{
  This function finds an alignment between the two input networks, specified in the form of adjacency matrices, by repeatedly calling \link{ComputeM} and \link{LinearAssignment}, up to maxNumSteps times. Simulated annealing is performed if a range is specified in the bStart and bEnd arguments. This simple procedure is described in detail in [Berg, Laessig 2006]. Different procedures can easily be implemented by the user.
  
  In each step, the matrix M is calculated from the scoring parameters and the current permutation vector P. The result is then normalized to the range [-1, 1] and, if simulated annealing is enabled, a random matrix depending on the current simulated annealing parameters is added. The linear assignment routine is used to calculate the value of P which is used to compute M in the next step.

  If the flag directed is set, directed binary networks are encoded by suitable symmetric matrices using \link{EncodeDirectedGraph}. The corresponding 3x3 matrices of the link score are computed from the 2x2 matrices given as input.

  Simulated annealing is enabled if bStart differs from bEnd. In this case, a value bStep = bEnd - bStart) / (maxNumSteps - 1) is calculated. In step n, the random matrix which is added to M is scaled by the factor 1 / [bStart + (n - 1) * bStep].
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
}
\references{
  Berg, J. & Laessig, M. (2006) Proc. Natl. Acad. Sci. USA 103, 10967-10972.
}
\author{Joern P. Meier, Michal Kolar, Ville Mustonen, Michael Laessig, and Johannes Berg}
\keyword{misc}
