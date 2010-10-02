\name{mat.sort}
\alias{mat.sort}
\title{Sort a Matrix According to a Specified Column}
\description{
 Sorts a matrix in respect to one or more specified columns by preserving its row context. If more than one columns is given the function uses the columns in the specified order for sorting.
}
\usage{
mat.sort(mat, sort, decreasing = FALSE)
}
\arguments{
  \item{mat}{matrix to be sorted}
  \item{sort}{column indices}
   \item{decreasing}{logical, if TRUE the column is ordered decreasingly}
}
\value{
  Returns the sorted matrix.
}
\author{Franziska Hoffgaard}

\examples{
mat<-matrix(data = as.integer(runif(16, 1, 20)), nrow = 4)

#Sort the matrix in respect to the second column of mat
m<-mat.sort(mat, 2)

#Sorting mat according to more than one column
mat<-matrix(data = c(rep(3, 3), as.integer(runif(13, 1, 20))), nrow = 4)
m1<-mat.sort(mat, c(1, 2))
m2<-mat.sort(mat, c(2, 1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
