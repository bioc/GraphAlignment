\name{ComputeNodeParameters}
\alias{ComputeNodeParameters}
\title{Compute node parameters}
\description{
  Compute the optimal scoring parameters (node score) for a given alignment.
}
\usage{
ComputeNodeParameters(dimA, dimB, R, P, lookupNode, clamp=TRUE)
}
\arguments{
  \item{dimA}{size of network A}
  \item{dimB}{size of network B}
  \item{R}{node similarity score matrix}
  \item{P}{permutation vector (see \link{InitialAlignment}, \link{AlignNetworks})}
  \item{lookupNode}{node bin lookup table (see \link{GetBinNumber})}
  \item{clamp}{clamp values to range when performing bin lookups}
}
\value{
  The return value is list containing the node score vectors s0 and s1.
}
\details{
  This function computes optimal node score parameters for use with \link{ComputeM} and \link{AlignNetworks}. It takes the size of the networks, a matrix of node similarities R, an initial alignment P, and the lookup table for node binning, lookupNode, as parameters.
}
\examples{
  ex<-GenerateExample(dimA=22, dimB=22, filling=.5, covariance=.6,
    symmetric=TRUE, numOrths=10, correlated=seq(1,18))
  
  pinitial<-InitialAlignment(psize=34, r=ex$r, mode="reciprocal")
  
  lookupNode<-c(-.5,.5,1.5)
  nodeParams<-ComputeNodeParameters(dimA=22, dimB=22, ex$r,
    pinitial, lookupNode)
}
\author{Joern P. Meier, Michal Kolar, Ville Mustonen, Michael Laessig, and Johannes Berg}
\keyword{misc}
