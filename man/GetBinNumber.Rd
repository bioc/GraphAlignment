\name{GetBinNumber}
\alias{GetBinNumber}
\title{Get bin number}
\description{
  Get the bin number for the argument using the specified lookup vector.
}
\usage{
GetBinNumber(x, lookup, clamp=TRUE)
}
\arguments{
  \item{x}{value to be binned}
  \item{lookup}{lookup vector}
  \item{clamp}{clamp value to lookup range}
}
\value{
  This function returns the bin number corresponding to the input value. The bin number is obtained by performing a lookup in the specified lookup vector.
}
\details{
  The lookup vector defines the lower and upper boundaries for each bin. The first entry in the lookup vector is the lower boundary of the first bin, while the last value in the lookup vector is the upper boundary of the last bin. For all other entries, entry i of the lookup vector defines the upper boundary of the (i-1)-th bin and the lower boundary of the i-th bin. The number of bins is therefore n-1, where n is the length of the lookup vector. A lookup vector must have at least two elements.
  
  If clamping is enabled (clamp=TRUE), arguments which fall below the lower boundary of the first bin are treated as if they are actually in the first bin. Likewise, values which are above the upper boundary of the last bin are treated as if they are actually in the last bin. If clamping is disabled (clamp=FALSE), values outside the lookup range cause an error.
}
\examples{
  lookup <- c(-1, 0, 1)
  GetBinNumber(-0.5, lookup)
  GetBinNumber(0.5, lookup)
}
\author{Joern P. Meier, Michal Kolar, Ville Mustonen, Michael Laessig, and Johannes Berg}
\keyword{misc}
