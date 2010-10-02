\name{sim}
\alias{sim}
\title{Compute the Covariance Matrices and B Factors for a List of PDBs}
\description{
 For a given list of PDB-files the respective covariance matrices and resulting B factors are computed.
}
\usage{
sim(pdbs, mj1 = NULL, mj2 = NULL, mj.avg = FALSE, alpha = 82,
    cuts = 169, path = getwd(), cluster = NULL)
}
\arguments{
  \item{pdbs}{list of PDB file names}
  \item{mj1}{matrix for the intrachain interaction strengths}
  \item{mj2}{matrix for the interchain interaction strengths}
  \item{mj.avg}{logical, if TRUE only the average value of the interaction matrix
      is used as value for the interaction of any two amino acids}
  \item{alpha}{strength of the peptide bond}
  \item{cuts}{ squared distance cutoff}
  \item{path}{ path to the output files }
  \item{cluster}{snow cluster object created with \code{makeCluster()}}
}
\details{
  Each PDB file is extracted. All features necessary for the computation of the covariance matrix
   in the elastic network model are computed. Both the covariance matrix and the B factors are computed and written to a file. The computation is parallelized using \code{parLapply()} from the package \pkg{snow}. If \code{cluster} is left at its default value the computation is carried out in serial.
}
\value{
No values are returned. The B factors and inverse Hessian matrices are written into files.
}
\references{
Hamacher and McCammon (2005) \emph{Journal of Chemical Theory and Computation} 2, 873.\cr
Tierney, Rossini, Li (2009) \emph{Int J Parallel Proc} 37, 78--90.\cr
}
\author{Franziska Hoffgaard, Philipp Weil}
\seealso{ \code{\link{sims}}, \code{\link{simc}} }
\examples{
\dontrun{
#build a list of pdb-files
pdbs<-list(system.file("1KZK.pdb", package = "BioPhysConnectoR"),
           system.file("1EBY.pdb", package = "BioPhysConnectoR"))
sim(pdbs, cuts = 169)

## Cluster example
 clu<-makeCluster(2)
 sim(pdbs, cuts = 169, cluster = clu)
 stopCluster(clu)
}
}
\keyword{ utilities }
