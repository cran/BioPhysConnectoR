\name{scpcp}
\Rdversion{1.1}
\alias{scpcp}
\title{
Self-Consistent Pair Contact Probability Approximation
}
\description{
The model of the Self-Consistent Pair Contact Probability (SCPCP) (Micheletti et al., 2001; Hamacher et al., 2006) computes equilibrium properties of structures with known native states. For a given contact map, extracted from a PDB file or artificially created, the fraction of native contacts, the free and internal energies are computed as well as the degree to which an amino acid is in its native state conformation. The maximum number of iteration and the preferred accuracy for the approximation can be specified.
}
\usage{
scpcp(T, R, cm, pstart = 0.5, maxiter = 2000, chains=NULL, maxtol = 1e-11,
      file = NULL, im = NULL)
}
\arguments{
  \item{T}{temperature}
  \item{R}{distance cutoff between current and native state}
  \item{cm}{contact map}
  \item{pstart}{initial probability value}
  \item{chains}{vector denoting the chain lengths}
  \item{maxiter}{maximum number of iterations}
  \item{maxtol}{tolerance}
  \item{file}{output file name}
  \item{im}{interaction matrix}
}
\value{
Returns a list with the following components
	\item{$free}{free energy}
	\item{$intern}{internal energy}
	\item{$entropy}{entropy}
	\item{$q}{fraction of native contacts}
	\item{$bfacs}{vector containing the B factors}
	\item{$pi}{probability vector}
	\item{$gmat}{resulting G matrix}
	\item{$iter}{number of iterations}
	\item{$err}{deviation of probabilities}
}
\details{
Results of each iteration as well as the final results are written into the given output file.
}
\references{
Micheletti, Banavar, Maritan (2001) \emph{Physical Review Letters} 87, 088102-1.\cr
Hamacher, Trylska, McCammon (2006)\emph{PLoS Computational Biology} 2, e10.\cr
Hamacher (2009) \emph{Eur. Biophys. J.},  in press.\cr
}
\author{
Franziska Hoffgaard
}

\examples{
p<-extractPDB(system.file("1KZK.pdb", package = "BioPhysConnectoR"))
cm<-build.contacts(n = p$lca, xyz = p$coords, cuts = 169)$cm
chains<-p$chains
im<-build.interact(cseq=p$caseq,d=chains,mj1=matrix(0.05,20,20),alpha=1)
res<-scpcp(cm = cm, maxtol = 10^(-10), T = 5, R = 0.8, im = im, 
	   maxiter = 10, pstart = 0.5, chains=chains)
}
\keyword{utilities}
