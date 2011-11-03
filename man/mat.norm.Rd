\name{mat.norm}
\alias{mat.norm}
\title{Normalization of a Matrix}
\description{
  A matrix is normalized by dividing each entry \eqn{[i,j]} by the square root of the product of the diagonal entries \eqn{[i,i]} and \eqn{[j,j]}. The input matrix should be a square matrix with positive diagonal entries.
}
\usage{
mat.norm(mat)
}
\arguments{
  \item{mat}{ numerical matrix to be normalized }
}
\value{
  The normalized matrix is returned.
}
\author{Franziska Hoffgaard}

\examples{
a<-matrix(runif(16, 1, 15), 4, 4)
b<-mat.norm(a)
}
\keyword{utilities}
