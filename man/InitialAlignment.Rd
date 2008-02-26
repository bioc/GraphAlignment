\name{InitialAlignment}
\alias{InitialAlignment}
\title{Create initial alignment}
\description{
  Create a permutation for an initial alignment. Either "random" or "reciprocal" may be specified as the mode argument.
}
\usage{
InitialAlignment(psize, r=NA, mode="random")
}
\arguments{
  \item{psize}{size of the alignment}
  \item{r}{node similarity score matrix (required for mode 'reciprocal')}
  \item{mode}{type of initial alignment}
}
\value{
  The return value is a permutation vector of the specified size.
}
\details{
  To create a random initial alignment of size psize, the \link{InitialAlignment} function can be used with the mode argument set to "random". If mode is set to "reciprocal", a reciprocal best match algorithm is applied to the input matrix R to find an initial alignment. This mode requires that the psize argument is sufficiently large to allow for the addition of dummy nodes to which unaligned nodes can formally be aligned.
}
\examples{
  ex<-GenerateExample(dimA=22, dimB=22, filling=.5, covariance=.6,
    symmetric=TRUE, numOrths=10, correlated=seq(1,18))
  
  pinitial<-InitialAlignment(psize=34, r=ex$r, mode="reciprocal")
}
\author{Joern P. Meier, Michal Kolar, Ville Mustonen, Michael Laessig, and Johannes Berg}
\keyword{misc}
