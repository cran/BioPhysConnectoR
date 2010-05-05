\name{mat.write}
\alias{mat.write}
\title{Writes Matrix Data to a File}
\description{
 Matrix data are written to a specified output file.
}
\usage{
mat.write(mat, file.name, ij = FALSE, sym = FALSE, sparse=FALSE, formatted=TRUE)
}
\arguments{
  \item{mat}{matrix which should be written}
  \item{file.name}{file name}
  \item{ij}{logical, if TRUE the format in the file will be: \eqn{i\quad j\quad value}, otherwise each row of the file represents one row in the matrix}
  \item{sym}{logical, if TRUE the matrix is symmetric, only important if \code{ij} = TRUE}
  \item{sparse}{logical, if TRUE the matrix is written in sparse format (e.g. only non-zero values are written to the file), only important if \code{ij} = TRUE}
  \item{formatted}{logical, if TRUE a blank line is inserted after all entries of one row, only important if \code{ij} = TRUE}
}
\details{
  In general there are two ways the output file can be written. First the matrix can be written to the file as \code{write.table} does and second the file can be written with the form: \eqn{i\quad j\quad value}, this means the value at position \eqn{[i,j]}. If the latter by \code{ij = TRUE} is selected, three more options can be chosen. If \code{sym = TRUE} the output file just contains the upper triangle and the diagonal of the matrix. If \code{spase = TRUE} only non-zero values are written. If \code{formatted = TRUE} a blank line seperates the entries of different rows.
}
\value{
 No values are returned.
}
\author{Philipp Weil}
\seealso{\code{\link{mat.read}}}
\keyword{IO}
