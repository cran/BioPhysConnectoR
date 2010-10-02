\name{get.freqs}
\Rdversion{1.1}
\alias{freq2p}
\alias{freq1p}
\title{Compute the Frequencies in an Alignment}
\description{
One and two point frequencies of a specified column or specified columns are computed for a
given alignment.
}
\usage{
freq1p(aln, i = NULL)

freq2p(i, aln, j2 = NULL, lett = NULL, cluster = NULL)
}
\arguments{
  \item{aln}{alignment}
  \item{i}{reference column for the frequency computation}
  \item{j2}{columns to compute the two point frequency together with column \code{i}}
  \item{lett}{character vector containing the unique alphabet of the alignment}
  \item{cluster}{snow cluster object created with \code{makeCluster()}}
}
\details{
The columns \code{j2} are the specified columns (not \code{i}) to compute the frequencies and the rows represent all possible two letter pairs according to the alphabet. The computation of \code{freq2p()} is parallelized using \code{parLapply()} from the package \pkg{snow}. If \code{cluster} is left at its default value the computation is carried out in serial.
}
\value{
A matrix is returned for the \code{freq2p()}. For \code{freq1p()}, if no \code{i} is specified, a matrix containing the frequencies of the symbols in the alignment is returned. For a certain \code{i} a vector with the respective frequencies of the symbols is the result.
}

\author{
Franziska Hoffgaard, Philipp Weil
}
\references{
Tierney, Rossini, Li (2009) \emph{Int J Parallel Proc} 37, 78--90.\cr
}
\seealso{
\code{\link{get.entropy}}, \code{\link{get.mie}}
}
\examples{
seqa<-unlist(strsplit("PQITLWQRPLVTIKIGGQL",split=""))
seqb<-unlist(strsplit("PQITLWKRPLVTIRIGGQL",split=""))
seqc<-unlist(strsplit("PQITLWQRPLVTIKIGGQL",split=""))
aln<-matrix(c(seqa,seqb,seqc),nrow=3,byrow=TRUE)
f1<-freq1p(aln)
f1a<-freq1p(aln, 1)
f2<-freq2p(1, aln, 2:10)
\dontrun{
## Cluster example
clu<-makeCluster(2)
f2<-freq2p(1, aln, 2:10, cluster=clu)
stopCluster(clu)
}
}
\keyword{ utilities }
