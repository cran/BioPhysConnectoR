\name{get.mie}
\alias{get.mie}

\title{Mutual Information }
\description{
 The joint information content (mutual information) for an alignment
    is computed. Considering the gap problem there are four ways to compute it.
}
\usage{
get.mie(aln, method = "ORMI", gapchar = NULL, nullmod = NULL)
}
\arguments{
  \item{aln}{ matrix of which the mutual information will be computed }
  \item{method}{ method, that is used for the computation (see details) }
  \item{gapchar}{ symbols of the input matrix that should be handled as gaps. These symbols are omitted in the computations (see details) }
  \item{nullmod}{ integer specifying how many shuffle runs should be performed }
}
\details{
Methods: The mutual information (MI) is computed as follows:
 \deqn{I(X,Y)=\sum_{x\in X}\sum_{y\in Y} p(x,y)\cdot\log_2\left(\frac{p(x,y)}{p(x)\cdot p(y)}\right)}{I(X,Y)=sum_(x in X) sum_(y in Y) p(x,y) log_2(p(x,y)/(p(x) p(y)))}

 \deqn{I(X,Y)=H(X)+H(Y)-H(X,Y)}
\describe{
\item{ORMI}{In the ORiginal MI gaps are treated simply as any other character.}
\item{SUMI}{The SUbset MI omits for each pair of columns all the rows with at least one gap character for the computation.}
\item{DEMI}{In the Delta Entropy MI the entropies for the columns are computed separately by leaving out any gap characters. The joint entropy also considers only rows without any gap character. DEMI follows as \eqn{H(X)+H(Y)-H(X,Y)}.}
\item{ESMI}{The Enhanced Sampling MI omits the gap characters and rows with gaps for the computation of the probabilities. The information content is computed via the probabilites.}}

The null model is computed by shuffling each column content and computing the resulting MI. Over all shuffle runs the MI values in each entry of the matrix is averaged. The averaged squared MI values and the variance is computed as well.
}
\value{
 Return value is the MI matrix per default. \cr
 If a null model should be computed the returned value is a list of matrices.
	\item{\code{$mi}}{MI matrix}
	\item{\code{$nullmodel}}{MI matrix of the null model}
	\item{\code{$nullsquare}}{matrix of the averaged MI squared values for the null model}
	\item{\code{$nullvar}}{matrix of the variance of each MI value over the whole shuffle run for the null model}
}
\author{Franziska Hoffgaard}
\seealso{\code{\link{get.entropy}}, \code{\link{get.entropy2p}}, \code{\link{freq1p}}, \code{\link{freq2p}}}
\examples{
seqa<-unlist(strsplit("PQITLWQRPLVTIKIGGQL",split=""))
seqb<-unlist(strsplit("PQITLWKRPLVTIRIGGQL",split=""))
seqc<-unlist(strsplit("PQITLWQRPLVTIKIGGQL",split=""))
a<-matrix(c(seqa,seqb,seqc),nrow=3,byrow=TRUE)
mi<-get.mie(a)
mi2<-get.mie(a, method = "SUMI", gapchar = "-")
mi_null<-get.mie(a,nullmod=100)
mi2_null<-get.mie(a, method = "SUMI", gapchar = "-",nullmod=100)
}
\keyword{utilities}
