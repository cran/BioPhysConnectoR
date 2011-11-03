\name{mat.read}
\alias{mat.read}
\title{Read Matrix Data from a File}
\description{
 A matrix is constructed from a specified input file.
}
\usage{
mat.read(file.name, ij = FALSE, sym = FALSE)
}
\arguments{
  \item{file.name}{file name}
  \item{ij}{logical, if TRUE the format in the file is: \eqn{i\quad j\quad value}, otherwise each row of the file represents one row in the matrix}
  \item{sym}{logical, if TRUE the matrix is symmetric, only important if \code{ij} = TRUE}
}
\details{
  There are two ways a matrix can be specified in the input file. If each line represents a row in the matrix (default), the matrix can simply be read. Otherwise it is possible to have an input file with the form: \eqn{i\quad j\quad value}, this means the value at position \eqn{[i,j]} of the matrix is set to \eqn{value}. If \code{ij} = TRUE and \code{sym} = TRUE only one half of the matrix needs to be provided.
}
\value{
 Return value is the matrix.
}
\author{Franziska Hoffgaard}
\keyword{ IO }