\name{get.entropy}
\alias{get.entropy}
\alias{get.entropy2p}
\title{ Compute the Sequence Entropy for an Alignment }
\description{
  Computes the sequence entropy of an alignment. It is possible to specify which characters to omit in the computation. The joint entropy is computed using \code{get.entropy2p()}.
}
\usage{
get.entropy(aln, bool = FALSE, gapchar = "NOGAPCHAR",
            verbose = FALSE)

get.entropy2p(aln, bool = FALSE, gapchar = "NOGAPCHAR",
              verbose = FALSE)
}
\arguments{
  \item{aln}{alignment matrix}
  \item{bool}{logical, if TRUE gaps are ignored when computing the entropy of each column of the alignment}
  \item{gapchar}{character vector containing the unique set of characters
                 representing gaps in the amino acid sequence}
  \item{verbose}{logical, TRUE for getting output messages}
 }
\details{
The Shannon (1948) entropy for an alignment is computed as follows:
\deqn{H(X)=-\sum_{x\in X} p(x)\cdot\log_2(p(x))}{H(X)=-sum_x p(x)log_2(p(x))}
The joint entropy is computed for every possible column pair:
\deqn{H(X,Y)=-\sum_{x\in X}\sum_{ y\in Y} p(x,y)\cdot\log_2(p(x,y))}{H(X)=-sum_(x,y) p(x,y)log_2(p(x,y))}
where \eqn{X} and \eqn{Y} are two columns of the alignment.
}
\value{
  Return value for \code{get.entropy()} is a vector containing the entropy for each column.\cr
  Return value for \code{get.entropy2p()} is a matrix containing the joint entropies in the lower triangle.
}
\references{
Shannon (1948) \emph{The Bell System Technical Journal} 27, 379--423.
}
\author{Franziska Hoffgaard}

\seealso{\code{\link{get.mie}}}

\examples{
aln<-matrix(c("M", "H", "X", "P", "V", "-", "H", "X", "L", "V", "M", "L",
 "X", "P", "V"), 3, byrow = TRUE)
h1<-get.entropy(aln, bool = TRUE , gapchar = "-")
h2<-get.entropy(aln)

h3<-get.entropy2p(aln)
}
\keyword{ utilities }
