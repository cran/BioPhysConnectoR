\name{build.interact}
\alias{build.interact}
\title{Compute the Interaction Matrix}
\description{
Calculation of the interaction strength between all amino acids of a protein.
}
\usage{
build.interact(cseq, mj1, mj2 = mj1, d, alpha = 82)
}
\arguments{
  \item{cseq}{(coded) amino acid sequence}
  \item{mj1}{matrix for intrachain interactions between amino acids}
  \item{mj2}{matrix for interchain interactions between amino acids}
  \item{d}{vector of chain lengths}
  \item{alpha}{strength of the peptide bond}
}
\details{
 Per default only one matrix with interaction strengths between amino acids is used to compute the interaction map. But it is possible to differentiate between intrachain [Miyazawa and Jernigan (1996)] and interchain [Keskin et al. (1998)] interactions by using an additional interaction matrix. Both matrices are included in the package and can be used directly. The value of \code{d} specifies the lengths of the different chains in the sequence.

Per definition there is no interaction between an amino acid and itself, therefore the diagonal is set to zero. All entries that represent peptide bonds between two amino acides are set to \code{alpha} (only entries in the off-diagonals).
The input of \code{cseq} can be a numeric sequences \eqn{[0,19]} as well as a string vector of standard amino acid sequences in one or three-letter code.
}
\value{
 Return value is a symmetric matrix specifying the interactions between the amino acids. If \eqn{n} is the length of the sequence, the matrix dimension results in \eqn{n \times n}.
}
\references{
Miyazawa and Jernigan (1996) \emph{Journal of Molecular Biology} 256, 623--644.\cr
Hamacher (2006) \emph{Journal of Chemical Theory and Computation} 2, 873--878.\cr
Keskin, Bahar, Badretdinov, Ptitsyn and Jernigan (1998) \emph{Protein Science} 7, 2578--2586.
}
\author{Franziska Hoffgaard}
\note{
Choosing one's own matrices for inter- and intrachain interactions requires some thought:
this  funcion only deals with the standard amino acids thus it requires always  symmetric, \eqn{20 \times 20} matrices.
}
\seealso{\code{\link{aa2num}}}
\examples{
#Miyazawa/Jernigan matrix
mat<-as.matrix(read.table(system.file("mj1.txt", package = "BioPhysConnectoR")))

#Keskin matrix
mat2<-as.matrix(read.table(system.file("mj2.txt", package = "BioPhysConnectoR")))

cseq<-rep(1:5, 2)
d<-c(5, 5)
im<-build.interact(cseq = cseq, mj1 = mat, mj2 = mat2, d = d)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
