\name{get.contact.list}
\alias{get.contact.list}
\title{Returns a List of Contacts for a given Contact Map }
\description{
 From a given contact map a list of contacts is computed.
}
\usage{
get.contact.list(cm, d = NULL, single = TRUE, val = 1)
}
\arguments{
  \item{cm}{binary matrix which specifies which amino acid position are in contact}
  \item{d}{vector specifying the last amino acids of each chain}
  \item{single}{bool, if TRUE only indices \eqn{i < j} are returned}
  \item{val}{numeric value, indices of matrix entries with this value are returned}
}
\details{
 All amino acid pairs with \code{cm[i,j] == val} are extracted from the \code{cm}. Each pair is listed once \eqn{[i,j] = [j,i]} (if \code{single} = TRUE) with \eqn{i < j}. If \code{d} is specified the contacts of the off-diagonal (peptide bonds) are left out, otherwise they will be returned as well. The values of \code{d} are the indices of the off-diagonal entries that do not represent a peptide bond.
}
\value{
  Return value is a list with the indices of all pairs of amino acids that fulfill the above stated condition.
}
\author{ Franziska Hoffgaard }

\seealso{\code{\link{build.contacts}}, \code{\link{simc}}}
\examples{
mat<-matrix(c(0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0,
              0, 1, 0, 1, 0, 0), ncol = 5)
get.contact.list(mat)
}
\keyword{ utilities }
