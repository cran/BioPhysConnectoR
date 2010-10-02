\name{build.contacts}
\alias{build.contacts}
\title{Determine the Contact Map and Distance Matrices}
\description{
Computation of a binary matrix specifying the contacts between each two amino acids \eqn{i} and \eqn{j} in respect to their spatial distance defined by their coordinates. Distance matrices for \eqn{\Delta x}, \eqn{\Delta y} and \eqn{\Delta z} are computed as well as the matrix containing the squared distances for all amino acid pairs.
}
\usage{
build.contacts(n, cuts, xyz)
}
\arguments{
  \item{n}{length of the amino acid sequence}
  \item{cuts}{squared cutoff}
  \item{xyz}{matrix with the x-, y- and z-coordinates of all \eqn{C_{\alpha}} atoms of the protein}
}
\details{
  If the squared distance between two \eqn{C_{\alpha}} atoms of amino acids \eqn{i} and  \eqn{j} is smaller than or equal to \code{cuts}, we assume a contact. In the contact map the value at indices \eqn{[i,j]} and \eqn{[j,i]} is set to \eqn{1}. Otherwise, the two \eqn{C_{\alpha}} atoms are not in contact, and the value is set to \eqn{0}. Per definition an atom is not in contact with itself. The contact map is a symmetric matrix. The matrices in \code{$deltas} (\code{$dx}, \code{$dy} and \code{$dz} and \code{$ds}) are symmetric matrices as well. The number of contacts between distinct amino acids is stored in \code{$cnr}.
}
\value{
Return value is a list with the following components:

  \item{$cm}{contact map}
  \item{$deltas}{list with distance matrices for x-, y- and z-direction as well as for the squared distance between all pairs of \eqn{C_{\alpha}} atoms}
  \item{$cnr}{number of contacts}
}
\author{Franziska Hoffgaard}

\seealso{\code{\link{get.contact.list}}}
\examples{
n<-10
xyz<-matrix(rep(1:10, 3), ncol = 3)
bc<-build.contacts(n, 3, xyz)
}
\keyword{ utilities }
