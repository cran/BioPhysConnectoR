\name{invhess}
\alias{build.invhess}
\alias{get.cov}
\title{Compute the Covariance Matrix / Inverse Hessian Matrix}
\description{
Computes the inverse Hessian matrix. The covariance matrix is computed as a pseudo-inverse derived from the eigenvalues and eigenvectors by a singular value decomposition (\code{get.svd()}) of the Hessian matrix. Otherwise, if neither the Hessian matrix nor the eigenvalues need to be stored, the inverse Hessian can directly be computed from the contact, interaction and distance matrices.
}
\usage{
build.invhess(svd_obj, singularity = 6)

get.cov(cm, im, deltas)
}
\arguments{
  \item{svd_obj}{svd object computed by \code{get.svd()} containing the eigenvector matrices, the eigenvalues and the index vector}
  \item{singularity}{number of eigenvalues equal/close to zero due to symmetries}
   \item{cm}{contact map for a protein}
  \item{im}{matrix of interaction strengths between the amino acids of the protein}
   \item{deltas}{difference matrices (x, y, z, squared) for all pairs of \eqn{C_{\alpha}} atoms as derived from \code{build.contacts()}}
}
\details{
  The calculation of the matrix omits by default the first six eigenvalues, because of translational and rotational symmetry in the model. The computation depends on the eigenvalues and -vectors. The number of eigenvalues to omit in the calculation can be specified by \code{singularity}. If the number of eigenvalues equalling zero is unknown and should be determined, the parameter \code{singularity} can be set to \verb{NULL}. The threshold for zero is set to \eqn{10^{-8}}.
}
\value{
Return value is the covariance matrix (also called inverse Hessian matrix).
}
\references{
Hamacher (2006) \emph{Journal of Chemical Theory and Computation} 2, 873--878.
}
\author{Franziska Hoffgaard}

\seealso{\code{\link{build.hess}}, \code{\link{get.svd}}}

\keyword{ utilities }
