\name{simc}
\alias{simc}
\title{Computed Elastic Network Models for Switched-Off-List of Contacts}
\description{
 For each entry in the contact list the contact will be broken and the resulting covariance matrix and new B factors will be computed in the elastic network model. Furthermore the Frobenius norms between the original and the new covariance matrix can be evaluated.
 }
\usage{
simc(pdb, mj1 = NULL, mj2 = NULL, mj.avg = FALSE, cl = NULL,
     alpha = 82, cuts = 169, path = getwd(), inv2file = FALSE,
     bfacs = TRUE, frob = TRUE, loc = NULL, norm = FALSE,
     file = NULL, cluster = NULL)
}
\arguments{
  \item{pdb}{file name of the PDB}
  \item{mj1}{matrix for the intrachain interaction strengths}
  \item{mj2}{matrix for the interchain interaction strengths}
  \item{mj.avg}{logical, if TRUE only the average value of the matrices is used as value for the interaction of any two amino acids}
  \item{cl}{optional contact list to process}
  \item{alpha}{strength of the peptide bond}
  \item{cuts}{squared cutoff distance}
  \item{path}{path to the output files}
  \item{inv2file}{logical, if TRUE the inverse Hessian matrix is written to a file, otherwise it will not be stored}
  \item{bfacs}{logical, if TRUE, the B factors are written to a file, otherwise they will not be stored}
  \item{frob}{logical, if TRUE, the Frobenius norm is computed}
  \item{loc}{dimensions \verb{i1, j1, i2, j2} for a matrix subset of which the frobenius norm should be computed }
  \item{norm}{logical, if TRUE the Frobenius norm is computed for the normalized matrices}
  \item{file}{personalized file name prefix}
  \item{cluster}{snow cluster object created with \code{makeCluster()}}
}
\details{
   If no contact list is given, the full contact list is extracted from the PDB-file. Each contact (except covalent contacts) in the list is broken and the corresponding covariance matrix and B factors are computed. Those can be written into files. A user-defined contact list can be specified as well. For the computation of the Frobenius norm, different regions can be specified in \code{loc} as matrix. Each row determines a region \verb{[i1:j1,i2:j2]} to be used for the norm. The routine is parallelized for the list of contacts using \code{parLapply()} from the package \pkg{snow}. If \code{cluster} is left at its default value the computaion is carries out in serial.
   }
\value{
No values are returned.
}
\references{
Hamacher and McCammon (2005) \emph{Journal of Chemical Theory and Computation} 2, 873.\cr
Hamacher (2008) \emph{Gene} 422, 30--36.\cr
Tierney, Rossini, Li (2009) \emph{Int J Parallel Proc} 37, 78--90.\cr
}
\author{Franziska Hoffgaard, Philipp Weil}
\seealso{ \code{\link{sim}}, \code{\link{sims}} }
\examples{
\dontrun{
cl<-matrix(c(3,1,4,1,5,1,9,1,10,1,11,1,24,1,66,1,67,1,68,1),ncol=2,byrow=TRUE)
out<-simc(system.file("1KZK.pdb", package = "BioPhysConnectoR"), cuts = 169, cl=cl)

## Cluster example
makeCluster(2)->clu
out<-simc(system.file("1KZK.pdb", package = "BioPhysConnectoR"), cuts = 169, cl=cl, 
                  cluster=clu)
stopCluster(clu)
}
}
\keyword{ utilities }
