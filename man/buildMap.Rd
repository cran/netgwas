\name{buildMap}
\alias{buildMap}

\title{ linkage group detection and ordering markers for class "netgwasmap"}

\description{
Implements different algorithms for detecting linkage groups and ordering markers in each linkage group.
}

\usage{
buildMap( res, opt.index, num.iso.m = NULL, use.comu = FALSE)
}

\arguments{
	\item{res}{ An object with S3 class "netgwasmap"}
	\item{opt.index}{ An index of a desired regularization parameter.}
	\item{num.iso.m}{A number of unlinked markers to be counted as isolated markers. Default is one.}
	\item{use.comu}{Using community detection algorithm to detect linkage groups. Default is FALSE.}
}

\details{
This function determines linkage groups and order markers within each linkage group for class "netgwasmap".
}

\value{
 \item{map}{ An estimated linkage map for a graph G associated with the specified regularization parameter in \code{opt.index}}
 }

\references{
1. Behrouzi, P., and Wit, E. C. (2017b). De novo construction of q-ploid linkage maps using discrete graphical models. arXiv preprint, arXiv:1710.01063. \cr
2. Behrouzi, P., and Wit, E. C. (2017c). netgwas: An R Package for Network-Based Genome-Wide Association Studies. arXiv preprint, arXiv:1710.01236.
}

\author{ 
  Pariya Behrouzi and Ernst C.Wit \cr
  Maintainer: Pariya Behrouzi \email{pariya.behrouzi@gmail.com} \cr
}

\seealso{\code{\link{netmap}}}

\examples{
\dontrun{
data(CviCol)
#Randomly change the order of markers across the genome
cvicol <- CviCol[ ,sample(1:ncol(CviCol), ncol(CviCol), replace=FALSE)]

#Constructing linkage map for Cvi x Col genotype data
out <- netmap(cvicol, cross.typ= "inbred", vis= FALSE, ncores=1); out
plot(out)
map <- out$map; map

#Visualization of other three-dimensional maps
plot(out$allres)  
#Constructing a linkage map for 5th three-dimensional map
bm <- buildMap(out, opt.index=5); bm
}
}