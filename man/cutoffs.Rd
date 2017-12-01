\name{cutoffs}
\alias{cutoffs}

\title{
Cut-points 
}

\description{
Calculates cut-points of ordinal variables with respect to the Gaussian copula. 
}

\usage{
cutoffs(y)
}

\arguments{
  \item{y}{
An (\eqn{n \times p}) matrix or a \code{data.frame} corresponding to the data matrix (\eqn{n} is the sample size and \eqn{p} is the number of variables). It also could be an object of class \code{"simgeno"}.
}
}

\details{ The relationship between \eqn{j}th variable and \eqn{j}th latent variable is expressed through this set of cut-points.
}

\value{
\item{cutoffs}{ A \eqn{p} by \eqn{(k + 1)} matrix representing the cut-point values under the Gaussian copula, where \eqn{k} defines the number of categories in the dataset.}
}

\author{ 
  Pariya Behrouzi and Ernst C. Wit \cr
  Maintainer: Pariya Behrouzi <pariya.behrouzi@gmail.com>
}

\references{
1. Behrouzi, P., and Wit, E. C. (2017a). Detecting Epistatic Selection with Partially Observed Genotype Data Using Copula Graphical Models. arXiv preprint, arXiv:1710.00894. \cr
2. Behrouzi, P., and Wit, E. C. (2017b). De novo construction of q-ploid linkage maps using discrete graphical models. arXiv preprint, arXiv:1710.01063. \cr
3. Behrouzi, P., and Wit, E. C. (2017c). netgwas: An R Package for Network-Based Genome-Wide Association Studies. arXiv preprint, arXiv:1710.01236.
}

\seealso{
  \code{\link{lower.upper}}, \code{\link{simgeno}} and \code{\link{netgwas-package}}.
}

\examples{
	D <- simgeno(p = 100, n = 50, k = 3)
	cutoffs(D$data)
}