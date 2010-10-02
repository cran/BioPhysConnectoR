\name{show.code}
\alias{show.code}
\title{Output of the Amino Acid Coding Scheme}
\description{
 Displays an encoding table for the standard amino acids.
}
\usage{
show.code(code = 0:19, offset = 0)
}
\arguments{
  \item{code}{vector of 20 integer values used for coding the standard amino acids}
  \item{offset}{offset added to the values of \code{code}}
}
\details{
This function shows the coding order of the standard amino acid
   as used in \code{aa2num()}. The offset is included.
}
\value{
No values are returned.
}
\author{Franziska Hoffgaard}
\seealso{\code{\link{aa2num}}}
\examples{
show.code()
}
\keyword{ utilities }
