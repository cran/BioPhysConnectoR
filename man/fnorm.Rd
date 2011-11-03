\name{fnorm}
\alias{fnorm}
\title{Frobenius Norm of Two Matrices}
\description{
  Computation of the Frobenius norm of two matrices as the sum of the squared differences between these matrices.
}
\usage{
fnorm(mat1, mat2)
}
\arguments{
  \item{mat1}{numerical matrix}
  \item{mat2}{numerical matrix}
}
\value{
Return value is the Frobenius norm.
}
\author{Franziska Hoffgaard}
\note{ 
Both input matrices must have the same dimensions.
}
\examples{
a<-matrix(runif(15, 1, 7), 5, 3)
b<-matrix(rnorm(15, 1, 7), 5, 3)
fn<-fnorm(a, b)
}
\keyword{utilities}
