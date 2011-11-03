\name{BioPhysConnectoR-package}
\Rdversion{1.1}
\alias{BioPhysConnectoR}
\alias{BioPhysConnectoR}
\docType{package}
\title{
BioPhysConnectoR
}
\description{
Functions to investigate the relation between biomolecular structures, their interactions, and the evolutionary information contained within sequence alignments of such molecules.
}
\details{
\tabular{ll}{
Package: \tab BioPhysConnectoR\cr
Version: \tab 1.6-7\cr
Date: \tab 2011-11-03\cr
Depends: \tab R (>= 2.11.0), snow, matrixcalc\cr
License: \tab GPL 2 \cr
Packaged: \tab 2010-23-09 13:50:31 UTC;weil\cr
Built: \tab R 2.12.1; x86_64-pc-linux-gnu; 2011-11-03 16:55:01 UTC; unix\cr
}

Index:
\preformatted{
aa2num                    Conversion of Amino Acids into Integer Values
BioPhysConnectoR-package  BioPhysConnectoR
build.contacts            Determine the Contact Map and Distance Matrices
build.hess                Construct the Hessian Matrix
build.interact            Compute the Interaction Matrix
build.invhess             Compute the Covariance Matrix /
                          Inverse Hessian Matrix
extractPDB                Extract Data from a PDB-File
fnorm                     Frobenius Norm of Two Matrices
get.bfacs                 Determine B factors
get.contact.list          Returns a List of Contacts for a given Contact Map
get.cov                   Compute the Covariance Matrix /
                          Inverse Hessian Matrix
get.entropy               Compute the Sequence Entropy for an Alignment
get.freqs                 Compute the Frequencies in an Alignment
get.mie                   Mutual Information
get.svd                   Singular Value Decomposition
lbpc                      List the Functions of the BioPhysConnectoR Package
mat.norm                  Normalisation of a Matrix
mat.read                  Read Matrix Data from a File
mat.sort                  Sort a Matrix According to a Specified Column
mat.write                 Writes Matrix Data to a File
read.fasta                Reads aligned or un-aligned sequences from a FASTA format file
scpcp                     Self-Consistent Pair Contact Probability
                          Approximation
show.code                 Output of the Amino Acid Coding Scheme
sim                       Compute the Correlation Matrix and B Factors
                          for a List of PDBs
simc                      Computed Elastic Network Models for
                          Switched-Off-List of Contacts
sims                      Apply a List of Different Amino Acid Sequences
}
}
\author{
Franziska Hoffgaard <hoffgaard@bio.tu-darmstadt.de>, with contributions from Philipp Weil <weil@bio.tu-darmstadt.de> and Kay Hamacher <hamacher@bio.tu-darmstadt.de>\cr

Maintainer: Franziska Hoffgaard <hoffgaard@bio.tu-darmstadt.de>
}

\references{
Hoffgaard, Weil, Hamacher (2010) \emph{BMC Bioinformatics} 11, 199. \cr
Hamacher (2006) \emph{Journal of Chemical Theory and Computation} 2, 873--878.\cr
Hamacher (2008) \emph{Gene} 422, 30--36.\cr
Hamacher (2009) \emph{Eur. Biophys. J.},  in press.\cr
Grant, Rodrigues, ElSawy, McCammon, Caves, (2006) \emph{Bioinformatics} 22, 2695--2696.\cr
Tierney, Rossini, Li (2009) \emph{Int J Parallel Proc} 37, 78--90.\cr
Novomestky (2008) \emph{matrixcalc}.\cr
Newman (2002) \emph{Physical Review Letters} 89, 208701-1 -- 208701- 4.\cr
Miyazawa, Jernigan (1996) \emph{Journal of Molecular Biology} 256, 623--644.\cr
Keskin, Bahar, Badretdinov, Ptitsyn, Jernigan (1998) \emph{Protein Science} 7, 2578--2586.\cr
Shannon (1948) \emph{The Bell System Technical Journal} 27, 379--423\cr
}
\keyword{ package }

