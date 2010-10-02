\name{build.hess}
\alias{build.hess}
\title{Construct the Hessian Matrix}
\description{
  The Hessian matrix is computed using the interaction matrix, the contact map, and the coordinate differences.
}
\usage{
build.hess(cm, im, deltas)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cm}{contact map}
  \item{im}{interaction map}
  \item{deltas}{list of difference matrices computed with \code{build.contacts}}
}
\details{
  For a sequence of length \eqn{n} a \eqn{3n \times 3n} matrix is computed. Only if amino acids are in contact (\code{cm} entry equals \eqn{1}) values for the Hessian matrix are computed as proposed in Hamacher (2006). If amino acids are not in contact their respective matrix entries equal \eqn{0}. Values for amino acids in contact depend on coordinate differences and interaction strengths.
}
\value{
Return value is the Hessian Matrix.
}
\references{
Hamacher (2006) \emph{Journal of Chemical Theory and Computation} 2, 873--878.
}
\author{Franziska Hoffgaard}

\seealso{\code{\link{build.interact}}, \code{\link{build.contacts}}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
