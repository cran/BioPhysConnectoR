\name{extractPDB}
\alias{extractPDB}
\title{Extract Data from a PDB-File}
\description{
  This function reads a PDB file and extracts data from it.
}
\usage{
extractPDB(file.name, verbose = TRUE)
}
\arguments{
  \item{file.name}{name of the PDB file}
  \item{verbose}{logical, if TRUE prints the  output from both functions \code{read.pdb()} and \code{atom.select()}.}
}
\details{
  The input PDB file is read. All \eqn{C_{\alpha}} atoms are extracted. For each \eqn{C_{\alpha}} atom the x-, y- and z-coordinates as well as the amino acid type and the B factor are extracted. The sequence length is compared to the number of \eqn{C_{\alpha}} atoms in the PDB. For each chain of the protein the length is computed.
}
\value{
  Returns a list with the following components:\cr

  \item{pdb}{list of class \dQuote{pdb} as originally extracted by \code{read.pdb()} in the \pkg{bio3d} package}
  \item{seq}{sequence according to the \dQuote{SEQRES} entries of the PDB file}
  \item{lseq}{length of the extracted sequence \verb{seq}}
  \item{lca}{number of \eqn{C_{\alpha}} atoms}
  \item{caseq}{sequence data for the \eqn{C_{\alpha}} selection}
  \item{coords}{matrix of coordinates for each \eqn{C_{\alpha}} atom, the rows are the \eqn{C_{\alpha}}'s, the columns are x, y, z}
  \item{b}{B factor for each \eqn{C_{\alpha}}}
  \item{chains}{integer vector with the lengths of the different chains in the protein sequence}
}
\references{
Grant, Rodrigues, ElSawy, McCammon, Caves, (2006) \emph{Bioinformatics} 22, 2695--2696.
}
\author{Franziska Hoffgaard}
\note{
Alternate and insert records in the PDB file are ignored. When removing alternate records a message is printed. Please check the PDB file.

If the number of \eqn{C_{\alpha}} atoms \code{lca} and the length of the extracted sequence \code{seq} differ, a message is printed. It is not advisable to use any other parts of the BioPhysConnectoR-package until you have found an appropriate solution to deal with this problem in the pdb-file.

This function includes source code of the \pkg{bio3d} package.
}
\examples{
pdb<-extractPDB(system.file("1KZK.pdb", package = "BioPhysConnectoR"))
}
\keyword{ IO }
