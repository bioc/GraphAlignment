\name{ComputeLinkParameters}
\alias{ComputeLinkParameters}
\title{Compute link parameters}
\description{
  Compute the optimal scoring parameters (link score) for a given alignment.
}
\usage{
ComputeLinkParameters(A, B, P, lookupLink, clamp=TRUE)
}
\arguments{
  \item{A}{adjacency matrix for network A}
  \item{B}{adjacency matrix for network B}
  \item{P}{permutation vector (see \link{InitialAlignment}, \link{AlignNetworks})}
  \item{lookupLink}{link bin lookup table (see \link{GetBinNumber})}
  \item{clamp}{clamp values to range when performing bin lookups}
}
\value{
  The return value is a list containing the self link score matrix (lsSelf) and the link score matrix (ls).
}
\details{
  This function computes optimal link score parameters for use with \link{ComputeM} and \link{AlignNetworks}. It takes two matrices as well as an initial alignment P and the lookup table for link binning, lookupLink, as parameters.
}
\examples{
  ex<-GenerateExample(dimA=22, dimB=22, filling=.5, covariance=.6,
    symmetric=TRUE, numOrths=10, correlated=seq(1,18))
  
  pinitial<-InitialAlignment(psize=34, r=ex$r, mode="reciprocal")
  
  lookupLink<-seq(-2,2,.5)
  linkParams<-ComputeLinkParameters(ex$a, ex$b, pinitial, lookupLink)
}
\author{Joern P. Meier, Michal Kolar, Ville Mustonen, Michael Laessig, and Johannes Berg}
\keyword{misc}
