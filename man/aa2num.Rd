\name{aa2num}
\alias{aa2num}
\title{Conversion of Amino Acids into Integer Values }
\description{
A sequence containing standard amino acids is converted into a sequence of integer values. An offset for the code can be specified.
}
\usage{
aa2num(seq, offset = 1, code = 0:19, verbose = FALSE)
}
\arguments{
\item{seq}{string vector containing a sequence of standard amino acids}
\item{offset}{offset, added to the specified code}
\item{code}{integer vector of the values to be assigned to the amino acids}
\item{verbose}{logical, if TRUE the specific encoding is printed}
}
\details{
The default values for the code are the integer values in the range 0 to 19. It is possible to use another numerical code for the 20 amino acids. The default coding order of the amino acids is: \verb{CYS, MET, PHE, ILE, LEU, VAL, TRP, TYR, ALA, GLY, THR, SER, ASN, GLN, ASP, GLU, HIS, ARG, LYS, PRO}. Each amino acid is converted to a single element according to the \code{code}-vector. The offset \code{offset} is added. Both sequences with three-letter and one-letter code can be handled.}
\value{
Returns a sequence of integer values according to the given numerical code.
}
\references{
Grant, Rodrigues, ElSawy, McCammon, Caves, (2006) \emph{Bioinformatics} 22, 2695--2696.
}
\author{Franziska Hoffgaard}
\note{
The function assigns values only for the standard amino acids. This function includes source code of the \pkg{bio3d} package. }
\seealso{\code{\link{show.code}}}
\examples{
seq<-c("MET", "GLY", "PRO", "LYS", "ASN")
aa2num(seq)
}
\keyword{ utilities }
