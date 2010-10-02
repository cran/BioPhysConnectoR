\name{get.bfacs}
\alias{get.bfacs}
\title{Determine B factors}
\description{
 The function calculates B factors from the inverse Hessian matrix.
}
\usage{
get.bfacs(covmat)
}
\arguments{
   \item{covmat}{inverse Hesse matrix}
}
\details{
  The B factors for each \eqn{C_{\alpha}} atom are computed from the diagonal of the covariance matrix by summing up the corresponding entries.
}
\value{
 Return value is a vector with the B factors for each \eqn{C_{\alpha}}.
}


\author{Franziska Hoffgaard}

\seealso{\code{\link{build.invhess}}}
\examples{
#extract the example pdb
pdb<-extractPDB(system.file("1KZK.pdb", package = "BioPhysConnectoR"))

#build the contact matrix for a given squared cutoff of 169
bc<-build.contacts(pdb$lca, 169, pdb$coords)

#for this example we use the fictitious value 12 as interaction strength
interaction.mat<-build.interact(cseq = rep(0, 198), mj1 = matrix(12, 20, 20),
d = pdb$chains)

#compute the inverse hessian matrix
cov.mat<-get.cov(cm = bc$cm, im = interaction.mat, deltas = bc$deltas)
bfac<-get.bfacs(cov.mat)
}
\keyword{ utilities }
