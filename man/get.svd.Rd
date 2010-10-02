\name{get.svd}
\alias{get.svd}
\title{Singular Value Decomposition}
\description{
  Computation of the singular value decomposition for a matrix. Sorts the eigenvalues in ascending order by maintaining their original order in an index array.
}
\usage{
get.svd(hessian.mat, linpack = TRUE)
}
\arguments{
  \item{hessian.mat}{input matrix is a Hessian matrix}
  \item{linpack}{logical, specifies whether LINPACK or LAPACK routines are used}
}

\value{
Return value is a list with the following components:
  \item{$v}{matrix with eigenvectors in each column}
  \item{$indx}{index array for the sorted eigenvalues}
  \item{$ev}{vector with the sorted eigenvalues}
  \item{$u}{matrix with eigenvectors in each column}
}
\author{Franziska Hoffgaard}
\note{This function uses the \code{svd()} function of \code{R}.}
\seealso{\code{\link{mat.sort}}, \code{\link{build.hess}}}
\examples{
a<-matrix(round(runif(9, 1, 5)), 3, 3)
out<-get.svd(a)
}
\keyword{ utilities }