\name{sims}
\alias{sims}
\title{Apply a List of Different Amino Acid Sequences}
\description{
  For a set of sequences given in an alignment and a corresponding PDB file the covariance matrix for each sequence is computed based on the given molecular structure. The latter can be directed into a file. The Frobenius norms are computed upon request.
}
\usage{
sims(pdb, alignment, mj1 = NULL, mj2 = NULL, mj.avg = FALSE,
     alpha = 82, cuts = 169, path = getwd(), mimethod = "ORMI",
     gapchar = "NOGAPCHAR", inv2file = FALSE, bfacs = TRUE,
     frob = TRUE, loc = NULL, norm = FALSE, cluster = NULL)
}
\arguments{
  \item{pdb}{PDB file}
  \item{alignment}{alignment file in fasta format}
  \item{mj1}{matrix for the intrachain interaction strengths}
  \item{mj2}{matrix for the interchain interaction strengths}
  \item{mj.avg}{logical, if TRUE only the average value of each matrix is used as value for the interaction of any two amino acids}
  \item{alpha}{strength of the peptide bond}
  \item{cuts}{squared distance cutoff}
  \item{path}{path to the output files }
  \item{mimethod}{method for the computation of the mutual information}
  \item{gapchar}{character vector denoting gaps in the alignment}
  \item{inv2file}{logical, if TRUE, the inverse Hessian matrix is written to a file}
  \item{bfacs}{logical, if TRUE, the B factors are written to a file}
  \item{frob}{logical, if TRUE, the Frobenius norm is computed}
  \item{loc}{dimensions \verb{i1, j1, i2, j2} for a matrix subset of which the Frobenius norm should be computed}
  \item{norm}{logical, if TRUE the Frobenius norm is computed for the normalized matrices}
  \item{cluster}{snow cluster object created with \code{makeCluster()}}
}
\details{
  For an alignment the sequence entropy and mutual information is computed. Furthermore for each sequence in the alignment and the structure information from the PDB file the covariance matrix and B factors can be computed and the output can be written to files. The computation is parallelized using \code{parLapply()} from the package \pkg{snow}. If \code{cluster} is left at its default value the computation is carried out in serial.
}
\value{
  Return value is a list with the following components:

  \item{$entropy}{sequence entropy of the alignment}
  \item{$mi}{mutual information of the alignment}

  If the Frobenius norm is computed, the value(s) will be returned as well:

  \item{$res}{vector consisting of: row number, sequence name, Frobenius norm}

}
\references{
Hamacher and McCammon (2005) \emph{Journal of Chemical Theory and Computation} 2, 873.\cr
Grant, Rodrigues, ElSawy, McCammon, Caves, (2006) \emph{Bioinformatics} 22, 2695--2696.\cr
Tierney, Rossini, Li (2009) \emph{Int J Parallel Proc} 37, 78--90.\cr
}
\author{ Franziska Hoffgaard, Philipp Weil}
\note{Make sure your alignment contains the IUPAC standard amino acids only without any gap characters.\cr 
This function includes source code of the \pkg{bio3d} package. 
}
\seealso{\code{\link{get.entropy}}, \code{\link{get.mie}}, \code{\link{sim}}, \code{\link{simc}} }
\examples{
\dontrun{
sims(system.file("1KZK.pdb",package="BioPhysConnectoR"),system.file("align.fasta",
     package="BioPhysConnectoR"))

## Cluster example
clu<-makeCluster(2)
sims(system.file("1KZK.pdb",package="BioPhysConnectoR"),system.file("align.fasta",
     package="BioPhysConnectoR",cluster = clu))
stopCluster(clu)

}
}
\keyword{ utilities }
