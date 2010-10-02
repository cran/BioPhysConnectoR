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
scpcp(t, r, cm, pstart = 0.5, maxiter, chains=NULL, eps = 1e-11,
      file = NULL, im = NULL)
}
\arguments{
  \item{t}{temperature}
  \item{r}{distance cutoff between current and native state}
  \item{cm}{contact map}
  \item{pstart}{initial probability value}
  \item{chains}{vector denoting the chain lengths}
  \item{maxiter}{maximum number of iterations}
  \item{eps}{tolerance}
  \item{file}{to store the results in file, the file name prefix has to be given}
  \item{im}{interaction matrix, per default each interaction strength equals 1}
}
\value{
Returns a list with the following components
	\item{$q}{fraction of native contacts}
	\item{$fe}{free energy}
	\item{$ie}{internal energy}
	\item{$p}{probability vector}
	\item{$iter}{number of iterations}
	\item{$b}{vector containing the B factors}
}
\details{
 If a file name prefix is given, the output is directed to three files. The file ending with \dQuote{_qfi.out} covers the following information: "t q fe ie i". The files \dQuote{_p.out} and \dQuote{_bfacs.out} contain "t" followed by the probabilities or the B factors for each amino acid.
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
res<-scpcp(cm = cm, eps = 10^(-10), t = 5, r = 0.8, maxiter = 100,
           pstart = 0.5, chains=chains)
}
\keyword{utilities}
